import tempfile
from Bio.PDB import Superimposer, PDBParser
from Bio import SeqIO, pairwise2
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import substitution_matrices

import nanome
from nanome.util import Logs, async_callback


class RMSDV2(nanome.AsyncPluginInstance):

    def start(self):
        pass

    @async_callback
    async def on_run(self):
        Logs.message("Retrieving complexes")
        complexes = await self.request_complex_list()
        complexes = await self.request_complexes([comp.index for comp in complexes])

        Logs.message("Parsing Complexes to BioPython Structures.")
        parser = PDBParser(QUIET=True)
        fixed_comp = complexes[0]
        moving_comp = complexes[1]
        fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        fixed_comp.io.to_pdb(fixed_pdb.name)
        moving_comp.io.to_pdb(moving_pdb.name)
        fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
        moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)
        Logs.message("Aligning Structures.")
        mapping = self.align_sequences(fixed_struct, moving_struct)

        fixed_atoms = []
        moving_atoms = []
        for fixed_atom in fixed_struct.get_atoms():
            if fixed_atom.serial_number in mapping:
                fixed_atoms.append(fixed_atom)
                moving_atom_serial = mapping[fixed_atom.serial_number] 
                moving_atom = next(
                    atom for atom in moving_struct.get_atoms()
                    if atom.serial_number == moving_atom_serial)
                moving_atoms.append(moving_atom)

        # Logs.message("Selecting Chains")
        # fixed_chain = next(fixed_struct.get_chains())
        # moving_chain = next(moving_struct.get_chains())
        # fixed_chain_name = 'A'
        # moving_chain_name = 'A'
        Logs.message("Superimposing Structures.")
        superimposer = Superimposer()
        superimposer.set_atoms(fixed_atoms, moving_atoms)
        Logs.message(superimposer.rms)

    def align_sequences(self, structA, structB):
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

            _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, "X"))
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
            penalize_end_gaps=(False, False),
        )

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
