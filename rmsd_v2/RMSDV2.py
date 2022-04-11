import nanome
import tempfile
import time
from Bio.PDB.Structure import Structure
from Bio.PDB import Superimposer, PDBParser
from Bio import pairwise2
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import substitution_matrices
from itertools import chain
from scipy.spatial import KDTree
from nanome.util import Logs, async_callback, Matrix, ComplexUtils
from nanome.api.structure import Atom, Complex
from .menu import RMSDMenu


class RMSDV2(nanome.AsyncPluginInstance):

    def start(self):
        self.menu = RMSDMenu(self)

    @async_callback
    async def on_run(self):
        self.menu.enabled = True
        self.complexes = await self.request_complex_list()
        self.menu.render(complexes=self.complexes)

    @async_callback
    async def on_complex_list_updated(self, complexes):
        self.menu.render(complexes=complexes)

    @async_callback
    async def on_complex_added(self):
        self.complexes = await self.request_complex_list()
        await self.menu.render(complexes=self.complexes)

    @async_callback
    async def on_complex_removed(self):
        self.complexes = await self.request_complex_list()
        await self.menu.render(complexes=self.complexes)

    async def superimpose_by_entry(self, fixed_comp: Complex, moving_comps: list[Complex]):
        start_time = time.time()
        Logs.message(f"Superimposing {len(moving_comps)} structures")
        moving_comp_indices = [comp.index for comp in moving_comps]
        updated_comps = await self.request_complexes([fixed_comp.index, *moving_comp_indices])
        fixed_comp = updated_comps[0]
        moving_comps = updated_comps[1:]
        fixed_comp.locked = True
        comps_to_update = [fixed_comp]
        rmsd_results = {}
        for moving_comp in moving_comps:
            ComplexUtils.align_to(moving_comp, fixed_comp)
            parser = PDBParser(QUIET=True)
            fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
            moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
            fixed_comp.io.to_pdb(fixed_pdb.name)
            moving_comp.io.to_pdb(moving_pdb.name)
            fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
            moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)
            transform_matrix, rms = await self.superimpose(fixed_struct, moving_struct, alignment_type)
            moving_comp.set_surface_needs_redraw()
            for comp_atom in moving_comp.atoms:
                comp_atom.position = transform_matrix * comp_atom.position
            rmsd_results[moving_comp.full_name] = rms
            moving_comp.locked = True
            comps_to_update.append(moving_comp)
        await self.update_structures_deep(comps_to_update)
        end_time = time.time()
        process_time = end_time - start_time
        extra = {"process_time": process_time}
        Logs.message(
            f"Superposition completed in {round(end_time - start_time, 2)} seconds.",
            extra=extra)
        return rmsd_results

    @staticmethod
    def create_transform_matrix(superimposer: Superimposer) -> Matrix:
        """Convert rotation and transform matrix from superimposer into Nanome Matrix."""
        rot, tran = superimposer.rotran
        rot = rot.tolist()
        tran = tran.tolist()
        m = Matrix(4, 4)
        m[0][0:3] = rot[0]
        m[1][0:3] = rot[1]
        m[2][0:3] = rot[2]
        m[3][0:3] = tran
        m[3][3] = 1
        # transpose necessary because numpy and nanome matrices are opposite row/col
        m.transpose()
        return m

    async def superimpose(self, fixed_struct: Structure, moving_struct: Structure) -> tuple[Matrix, float]: 
        # Collect aligned residues
        # Align Residues based on Alpha Carbon
        mapping = self.align_structures(fixed_struct, moving_struct)
        fixed_atoms = []
        moving_atoms = []
        alpha_carbon = 'CA'
        for fixed_residue in fixed_struct.get_residues():
            fixed_id = fixed_residue.id[1]
            if fixed_id in mapping:
                fixed_atoms.append(fixed_residue[alpha_carbon])
                moving_residue_serial = mapping[fixed_id]
                moving_residue = next(
                    rez for rez in moving_struct.get_residues()
                    if rez.id[1] == moving_residue_serial)
                moving_atoms.append(moving_residue[alpha_carbon])
        assert len(moving_atoms) == len(fixed_atoms)
        Logs.message("Superimposing Structures.")
        superimposer = Superimposer()
        superimposer.set_atoms(fixed_atoms, moving_atoms)
        rms = superimposer.rms
        transform_matrix = self.create_transform_matrix(superimposer)
        Logs.message(f"RMSD: {rms}")
        return transform_matrix, rms

    async def superimpose_by_chain(self, fixed_comp: Complex, fixed_chain_name: str, moving_comp_chain_list: list):
        start_time = time.time()
        Logs.message("Superimposing by Chain.")
        moving_comp_indices = [item[0].index for item in moving_comp_chain_list]
        updated_comps = await self.request_complexes([fixed_comp.index, *moving_comp_indices])
        fixed_comp = updated_comps[0]
        moving_comps = updated_comps[1:]

        updated_moving_comps = []
        parser = PDBParser(QUIET=True)

        fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        fixed_comp.io.to_pdb(fixed_pdb.name)
        fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
        fixed_chain = next(ch for ch in fixed_struct.get_chains() if ch.id == fixed_chain_name)
        results = {}
        for i, moving_comp in enumerate(moving_comps):
            moving_chain_name = moving_comp_chain_list[i][1]
            ComplexUtils.align_to(moving_comp, fixed_comp)

            moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
            moving_comp.io.to_pdb(moving_pdb.name)
            moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)

            moving_chain = next(ch for ch in moving_struct.get_chains() if ch.id == moving_chain_name)
            transform_matrix, rms = await self.superimpose(fixed_chain, moving_chain)
            results[moving_comp.full_name] = rms
            # apply transformation to moving_comp
            moving_comp.set_surface_needs_redraw()
            for comp_atom in moving_comp.atoms:
                comp_atom.position = transform_matrix * comp_atom.position
            updated_moving_comps.append(moving_comp)

        await self.update_structures_deep(updated_moving_comps)
        end_time = time.time()
        process_time = end_time - start_time
        extra = {"process_time": process_time}
        Logs.message(
            f"Superposition completed in {round(process_time, 2)} seconds.",
            extra=extra)
        return results

    async def superimpose_by_active_site(
            self, target_reference: Complex, ligand_name: str, moving_comp_list: list[Complex], site_size=5):
        binding_site_atoms = await self.get_binding_site_atoms(target_reference, ligand_name, site_size)

    async def get_binding_site_atoms(self, target_reference: Complex, ligand_name: str, site_size=5) -> list[Atom]:
        """Identify atoms in the active site around a ligand."""
        mol = next(
            mol for i, mol in enumerate(target_reference.molecules)
            if i == target_reference.current_frame)
        target_ligands = await mol.get_ligands()
        ligand = next(ligand for ligand in target_ligands if ligand.name == ligand_name)

        # Use KDTree to find target atoms within site_size radius of ligand atoms
        ligand_atoms = chain(*[res.atoms for res in ligand.residues])
        ligand_positions = [atom.position.unpack() for atom in ligand_atoms]
        target_atoms = chain(*[ch.atoms for ch in mol.chains if not ch.name.startswith("H")])
        near_point_set = set()
        target_tree = KDTree([atom.position.unpack() for atom in target_atoms])
        target_point_indices = target_tree.query_ball_point(ligand_positions, site_size)
        for point_indices in target_point_indices:
            for point_index in point_indices:
                near_point_set.add(tuple(target_tree.data[point_index]))
        active_site_atoms = []
        for targ_atom in mol.atoms:
            if targ_atom.position.unpack() in near_point_set:
                targ_atom.selected = True
                active_site_atoms.append(targ_atom)
        Logs.message(f"{len(active_site_atoms)} atoms identified in binding site.")
        await self.update_structures_deep([target_reference])
        return active_site_atoms

    def align_structures(self, structA, structB, alignment_type='global'):
        """
        Performs a global pairwise alignment between two sequences
        using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
        as implemented in Biopython. Returns the alignment, the sequence
        identity and the residue mapping between both original sequences.

        Source: https://gist.github.com/JoaoRodrigues/e3a4f2139d10888c679eb1657a4d7080
        """

        def _get_pdb_sequence(structure):
            """
            Retrieves the AA sequence from a PDB structure.
            """
            def _aainfo(r): return (r.id[1], aa3to1.get(r.resname, "X"))
            seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
            return seq

        Logs.message(f"Running {alignment_type} alignment.")
        resseq_A = _get_pdb_sequence(structA)
        resseq_B = _get_pdb_sequence(structB)

        sequence_A = "".join([i[1] for i in resseq_A])
        sequence_B = "".join([i[1] for i in resseq_B])

        if alignment_type == 'global':
            alignment_fn = pairwise2.align.globalds
            alignment_fn_args = (
                sequence_A,
                sequence_B,
                substitution_matrices.load("BLOSUM62")
            )
            alignment_fn_kwargs = dict(
                one_alignment_only=True,
                open=-10.0,
                extend=-0.5,
                penalize_end_gaps=(False, False)
            )
        elif alignment_type == 'local':
            alignment_fn = pairwise2.align.localxx
            alignment_fn_args = (sequence_A, sequence_B)
            alignment_fn_kwargs = dict()
        alns = alignment_fn(*alignment_fn_args, **alignment_fn_kwargs)
        best_aln = alns[0]
        aligned_A, aligned_B, score, begin, end = best_aln

        # Equivalent residue numbering
        # Relative to reference
        mapping = {}
        aa_i_A, aa_i_B = 0, 0
        for aln_i, (aa_aln_A, aa_aln_B) in enumerate(zip(aligned_A, aligned_B)):
            if aa_aln_A == "-":
                if aa_aln_B != "-":
                    aa_i_B += 1
            elif aa_aln_B == "-":
                if aa_aln_A != "-":
                    aa_i_A += 1
            else:
                assert resseq_A[aa_i_A][1] == aa_aln_A
                assert resseq_B[aa_i_B][1] == aa_aln_B
                mapping[resseq_A[aa_i_A][0]] = resseq_B[aa_i_B][0]
                aa_i_A += 1
                aa_i_B += 1

        return mapping


def main():
    plugin = nanome.Plugin('RMSD V2', 'Superimpose two structures', 'alignment', False)
    plugin.set_plugin_class(RMSDV2)
    plugin.run()


if __name__ == '__main__':
    main()
