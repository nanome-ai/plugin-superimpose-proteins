import nanome
import tempfile
from Bio.PDB import Superimposer, PDBParser
from Bio import pairwise2
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import substitution_matrices
from nanome.util import Logs, async_callback, Matrix, ComplexUtils
from nanome.util.enums import NotificationTypes

from .menu import RMSDMenu


class RMSDV2(nanome.AsyncPluginInstance):

    def start(self):
        self.menu = RMSDMenu(self)

    @async_callback
    async def on_run(self):
        self.menu.enabled = True
        complexes = await self.request_complex_list()
        self.menu.render(complexes=complexes)

    @async_callback
    async def on_complex_list_updated(self, complexes):
        self.menu.render(complexes=complexes)

    @async_callback
    async def on_complex_added(self):
        complexes = await self.request_complex_list()
        default_values = len(complexes) == 2
        await self.menu.render(complexes=complexes, default_values=default_values)

    @async_callback
    async def on_complex_removed(self):
        complexes = await self.request_complex_list()
        await self.menu.render(complexes=complexes)

    async def msa_superimpose(self, fixed_comp, moving_comps, alignment_type='global'):
        moving_comp_indices = [comp.index for comp in moving_comps]
        updated_comps = await self.request_complexes([fixed_comp.index, *moving_comp_indices])
        fixed_comp = updated_comps[0]
        moving_comps = updated_comps[1:]
        fixed_comp.locked = True
        comps_to_update = [fixed_comp]
        rmsd_results = {}
        for moving_comp in moving_comps:
            updated_moving_comp, rms = await self.superimpose(fixed_comp, moving_comp, alignment_type)
            rmsd_results[moving_comp.full_name] = rms
            updated_moving_comp.locked = True
            comps_to_update.append(updated_moving_comp)
        await self.update_structures_deep(comps_to_update)
        return rmsd_results

    @staticmethod
    def create_transform_matrix(superimposer):
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

    async def superimpose(self, fixed_comp, moving_comp, alignment_type='global'):
        Logs.message(f"Superimposing {moving_comp.full_name} onto {fixed_comp.full_name}.")
        ComplexUtils.align_to(moving_comp, fixed_comp)
        parser = PDBParser(QUIET=True)
        fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        fixed_comp.io.to_pdb(fixed_pdb.name)
        moving_comp.io.to_pdb(moving_pdb.name)
        fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
        moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)

        Logs.message("Aligning Structures.")
        mapping = self.align_sequences(fixed_struct, moving_struct, alignment_type)

        # Collect aligned residues
        # Align Residues based on Alpha Carbon
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
        Logs.message(f"RMSD: {rms}")
        transform_matrix = self.create_transform_matrix(superimposer)
        # apply transformation to moving_comp
        moving_comp.set_surface_needs_redraw()
        for comp_atom in moving_comp.atoms:
            comp_atom.position = transform_matrix * comp_atom.position
        return moving_comp, rms

    async def superimpose_by_chain(self, fixed_comp, fixed_chain_name, moving_comp, moving_chain_name):
        fixed_comp, moving_comp = await self.request_complexes([fixed_comp.index, moving_comp.index])
        Logs.message("Parsing Complexes to BioPython Structures.")
        parser = PDBParser(QUIET=True)
        ComplexUtils.align_to(moving_comp, fixed_comp)

        fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        fixed_comp.io.to_pdb(fixed_pdb.name)
        moving_comp.io.to_pdb(moving_pdb.name)
        fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
        moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)

        Logs.message("Aligning Structures.")
        fixed_chain = next(ch for ch in fixed_struct.get_chains() if ch.id == fixed_chain_name)
        moving_chain = next(ch for ch in moving_struct.get_chains() if ch.id == moving_chain_name)
        mapping = self.align_sequences(fixed_chain, moving_chain)

        # Collect aligned residues
        # Align Residues based on Alpha Carbon
        fixed_atoms = []
        moving_atoms = []
        alpha_carbon = 'CA'
        for fixed_residue in fixed_chain.get_residues():
            fixed_id = fixed_residue.id[1]
            if fixed_id in mapping:
                fixed_atoms.append(fixed_residue[alpha_carbon])
                moving_residue_serial = mapping[fixed_id]
                moving_residue = next(
                    rez for rez in moving_chain.get_residues()
                    if rez.id[1] == moving_residue_serial)
                moving_atoms.append(moving_residue[alpha_carbon])
        assert len(moving_atoms) == len(fixed_atoms)
        Logs.message("Superimposing Structures.")
        superimposer = Superimposer()
        superimposer.set_atoms(fixed_atoms, moving_atoms)
        # superimposer.apply(moving_struct.get_atoms())
        rms = superimposer.rms
        Logs.message(f'RMSD: {rms}')

        m = self.create_transform_matrix(superimposer)
        # apply transformation to moving_comp
        moving_comp.set_surface_needs_redraw()
        for comp_atom in moving_comp.atoms:
            comp_atom.position = m * comp_atom.position
        await self.update_structures_deep([moving_comp])
        self.send_notification(NotificationTypes.message, f'RMSD: {rms}')
        return {moving_comp.full_name: rms}

    def align_sequences(self, structA, structB, alignment_type='global'):
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

        if alignment_type == 'global':
            Logs.message("Using Global Alignment")
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
            Logs.message("Using Local Alignment")
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
