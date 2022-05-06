import nanome
import os
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
from nanome.api.structure import Complex
from .enums import AlignmentMethodEnum
from .menu import MainMenu
from .fpocket_client import FPocketClient


class SuperimposePlugin(nanome.AsyncPluginInstance):

    def start(self):
        self.menu = MainMenu(self)

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
        complexes = await self.request_complex_list()
        await self.menu.render(complexes=complexes)

    @async_callback
    async def on_complex_removed(self):
        self.complexes = await self.request_complex_list()
        await self.menu.render(complexes=self.complexes)

    async def superimpose_by_entry(self, fixed_comp_index, moving_comp_indices, alignment_method):
        start_time = time.time()
        updated_comps = await self.request_complexes([fixed_comp_index, *moving_comp_indices])
        fixed_comp = updated_comps[0]
        moving_comps = updated_comps[1:]
        Logs.message(f"Superimposing {len(moving_comps)} structures")
        fixed_comp.locked = True
        comps_to_update = [fixed_comp]
        rmsd_results = {}
        comp_count = len(moving_comps)
        for i, moving_comp in enumerate(moving_comps):
            Logs.debug(f"Superimposing Moving Complex {i}")
            ComplexUtils.align_to(moving_comp, fixed_comp)
            parser = PDBParser(QUIET=True)
            fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
            moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
            fixed_comp.io.to_pdb(fixed_pdb.name)
            moving_comp.io.to_pdb(moving_pdb.name)
            fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
            moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)

            try:
                superimposer, paired_residue_count, paired_atom_count = await self.superimpose(
                    fixed_struct, moving_struct, alignment_method)
            except Exception:
                Logs.error(f"Superimposition failed for {moving_comp.full_name}")
                continue
            rmsd_results[moving_comp.full_name] = self.format_superimposer_data(superimposer, paired_residue_count, paired_atom_count)
            # Use matrix to transform moving atoms to new position
            transform_matrix = self.create_transform_matrix(superimposer)
            for comp_atom in moving_comp.atoms:
                comp_atom.position = transform_matrix * comp_atom.position
            moving_comp.set_surface_needs_redraw()
            moving_comp.locked = True
            comps_to_update.append(moving_comp)
            self.update_loading_bar(i + 1, comp_count)

        await self.update_structures_deep(comps_to_update)
        end_time = time.time()
        process_time = end_time - start_time
        extra = {"process_time": process_time}
        Logs.message(
            f"Superposition completed in {round(end_time - start_time, 2)} seconds.",
            extra=extra)
        return rmsd_results

    async def superimpose_by_chain(self, fixed_comp_index, fixed_chain_name, moving_comp_chain_list, alignment_method):
        start_time = time.time()
        Logs.message("Superimposing by Chain.")
        moving_comp_indices = [item[0] for item in moving_comp_chain_list]
        updated_comps = await self.request_complexes([fixed_comp_index, *moving_comp_indices])
        fixed_comp = updated_comps[0]
        moving_comps = updated_comps[1:]

        fixed_comp.locked = True
        comps_to_update = [fixed_comp]
        parser = PDBParser(QUIET=True)

        fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        fixed_comp.io.to_pdb(fixed_pdb.name)
        fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
        fixed_chain = next(ch for ch in fixed_struct.get_chains() if ch.id == fixed_chain_name)
        comp_count = len(moving_comps)
        results = {}
        for i, moving_comp in enumerate(moving_comps):
            Logs.debug(f"Superimposing Moving Complex {i}")
            moving_chain_name = moving_comp_chain_list[i][1]
            ComplexUtils.align_to(moving_comp, fixed_comp)

            moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
            moving_comp.io.to_pdb(moving_pdb.name)
            moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)

            moving_chain = next(
                ch for ch in moving_struct.get_chains()
                if ch.id == moving_chain_name)

            try:
                superimposer, paired_residue_count, paired_atom_count = await self.superimpose(
                    fixed_chain, moving_chain, alignment_method)
            except Exception:
                Logs.error(f"Superimposition failed for {moving_comp.full_name}")
                continue

            transform_matrix = self.create_transform_matrix(superimposer)

            comp_data = self.format_superimposer_data(
                superimposer, paired_residue_count, paired_atom_count, moving_chain_name)
            results[moving_comp.full_name] = comp_data

            # apply transformation to moving_comp
            for comp_atom in moving_comp.atoms:
                comp_atom.position = transform_matrix * comp_atom.position
            moving_comp.locked = True
            moving_comp.set_surface_needs_redraw()
            comps_to_update.append(moving_comp)
            self.update_loading_bar(i + 1, comp_count)

        await self.update_structures_deep(comps_to_update)
        end_time = time.time()
        process_time = end_time - start_time
        extra = {"process_time": process_time}
        Logs.message(
            f"Superposition completed in {round(process_time, 2)} seconds.",
            extra=extra)
        return results

    async def superimpose_by_active_site(
            self, target_reference: int, ligand_name: str, moving_indices: list, site_size=4.5):
        # Select the binding site on the target_reference.
        updated_complexes = await self.request_complexes([target_reference, *moving_indices])
        target_reference = updated_complexes[0]
        moving_comp_list = updated_complexes[1:]
        binding_site_atoms = await self.get_binding_site_atoms(target_reference, ligand_name, site_size)
        for atom in binding_site_atoms:
            atom.selected = True
        await self.update_structures_deep([target_reference])

        # For each moving comp, select the potential binding sites.
        fpocket_client = FPocketClient()
        for moving_comp in moving_comp_list:
            with tempfile.TemporaryDirectory() as tmpdir:
                output_dir = fpocket_client.run(moving_comp, tmpdir)
                pocket_sets = fpocket_client.parse_results(moving_comp, output_dir)
            for i in range(len(pocket_sets) - 1, -1, -1):
                Logs.debug(f"Highlighting pocket {i + 1}")
                for atom in moving_comp.atoms:
                    atom.selected = atom in pocket_sets[i]
                await self.update_structures_deep([moving_comp])

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

    async def superimpose(self, fixed_struct: Structure, moving_struct: Structure, alignment_method):
        # Collect aligned residues
        # Align Residues based on Alpha Carbon
        mapping = self.align_structures(fixed_struct, moving_struct)
        fixed_atoms = []
        moving_atoms = []
        alpha_carbon = 'CA'
        skip_count = 0
        for fixed_residue in fixed_struct.get_residues():
            fixed_id = fixed_residue.id[1]
            if fixed_id not in mapping:
                continue

            new_fixed_atoms = []
            new_moving_atoms = []
            if alignment_method == AlignmentMethodEnum.ALPHA_CARBONS_ONLY:
                # Add alpha carbons.
                new_fixed_atoms.append(fixed_residue[alpha_carbon])
            else:
                # Add all heavy atoms (Non hydrogens)
                for atom in fixed_residue.get_atoms():
                    if not atom.name.startswith('H'):
                        new_fixed_atoms.append(atom)

            # Get matching atoms from moving structure.
            moving_residue_serial = mapping[fixed_id]
            moving_residue = next(
                rez for rez in moving_struct.get_residues()
                if rez.id[1] == moving_residue_serial)
            if alignment_method == AlignmentMethodEnum.ALPHA_CARBONS_ONLY:
                new_moving_atoms.append(moving_residue[alpha_carbon])
            else:
                for atom in moving_residue.get_atoms():
                    if not atom.name.startswith('H'):
                        new_moving_atoms.append(atom)

            if len(new_moving_atoms) != len(new_fixed_atoms):
                # I think we can just skip residues with differing atom counts.
                # This is an issue with Heavy atom alignment methods.
                skip_count += 1
                continue
            fixed_atoms.extend(new_fixed_atoms)
            moving_atoms.extend(new_moving_atoms)
        assert len(moving_atoms) == len(fixed_atoms), f"{len(moving_atoms)} != {len(fixed_atoms)}"
        total_residue_count = sum(1 for _ in fixed_struct.get_residues())
        paired_residue_count = len(set([atom.get_parent() for atom in fixed_atoms]))
        if skip_count > 0:
            Logs.warning(
                f"Not including {skip_count}/{total_residue_count} residue pairs "
                "in RMSD calculation due to differing atom counts.")

        Logs.message("Superimposing Structures.")
        superimposer = Superimposer()
        superimposer.set_atoms(fixed_atoms, moving_atoms)
        rms = round(superimposer.rms, 2)
        Logs.debug(f"RMSD: {rms}")
        paired_atom_count = len(fixed_atoms)
        return superimposer, paired_residue_count, paired_atom_count

    def format_superimposer_data(self, superimposer: Superimposer, paired_residue_count: int, paired_atom_count: int, chain_name=''):
        # Set up data to return to caller
        rms = round(superimposer.rms, 2)
        comp_data = {
            'rmsd': rms,
            'paired_atoms': paired_atom_count,
            'paired_residues': paired_residue_count
        }
        if chain_name:
            comp_data['chain'] = chain_name
        return comp_data

    async def get_binding_site_atoms(self, target_reference: Complex, ligand_name: str, site_size=4.5):
        """Identify atoms in the active site around a ligand."""
        mol = next(
            mol for i, mol in enumerate(target_reference.molecules)
            if i == target_reference.current_frame)
        target_ligands = await mol.get_ligands()
        ligand = next(ligand for ligand in target_ligands if ligand.name == ligand_name)
        # Use KDTree to find target atoms within site_size radius of ligand atoms
        ligand_atoms = chain(*[res.atoms for res in ligand.residues])
        binding_site_atoms = self.calculate_binding_site_atoms(target_reference, ligand_atoms)
        return binding_site_atoms

    def calculate_binding_site_atoms(self, target_reference: Complex, ligand_atoms: list, site_size=4.5):
        """Use KDTree to find target atoms within site_size radius of ligand atoms."""
        mol = next(
            mol for i, mol in enumerate(target_reference.molecules)
            if i == target_reference.current_frame)
        ligand_positions = [atom.position.unpack() for atom in ligand_atoms]
        target_atoms = chain(*[ch.atoms for ch in mol.chains if not ch.name.startswith("H")])
        target_tree = KDTree([atom.position.unpack() for atom in target_atoms])
        target_point_indices = target_tree.query_ball_point(ligand_positions, site_size)
        near_point_set = set()
        for point_indices in target_point_indices:
            for point_index in point_indices:
                near_point_set.add(tuple(target_tree.data[point_index]))
        binding_site_atoms = []

        for targ_atom in mol.atoms:
            if targ_atom.position.unpack() in near_point_set:
                binding_site_atoms.append(targ_atom)
        Logs.message(f"{len(binding_site_atoms)} atoms identified in binding site.")
        return binding_site_atoms

    def align_structures(self, structA, structB):
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

        resseq_A = _get_pdb_sequence(structA)
        resseq_B = _get_pdb_sequence(structB)

        sequence_A = "".join([i[1] for i in resseq_A])
        sequence_B = "".join([i[1] for i in resseq_B])
        alns = pairwise2.align.globalds(
            sequence_A,
            sequence_B,
            substitution_matrices.load("BLOSUM62"),
            one_alignment_only=True,
            open=-10.0,
            extend=-0.5,
            penalize_end_gaps=(False, False))
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

    def update_loading_bar(self, current, total):
        self.menu.update_loading_bar(current, total)


def main():
    default_description = 'Superimpose two or more structures'
    description = os.environ.get("PLUGIN_DESCRIPTION", "") or default_description
    plugin = nanome.Plugin('Superimpose', description, 'alignment', False)
    plugin.set_plugin_class(SuperimposePlugin)
    plugin.run()


if __name__ == '__main__':
    main()
