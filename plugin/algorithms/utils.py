import time
from Bio import pairwise2
from Bio.PDB import Superimposer
from Bio.PDB.Structure import Structure
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.PDB.Polypeptide import is_aa

from nanome.util import Matrix, Logs

def format_superimposer_data(superimposer: Superimposer, paired_residue_count: int, paired_atom_count: int, chain_name=''):
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

def align_structures(structA, structB):
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

    start_time = time.time()
    Logs.message("Calculating alignment")
    resseq_A = _get_pdb_sequence(structA)
    resseq_B = _get_pdb_sequence(structB)

    sequence_A = "".join([i[1] for i in resseq_A])
    sequence_B = "".join([i[1] for i in resseq_B])
    alns = pairwise2.align.globalms(sequence_A, sequence_B, 2, -1, -10, -.5)
    best_aln = alns[0]
    aligned_A, aligned_B, score, begin, end = best_aln

    # Equivalent residue numbering
    # Relative to reference
    mapping = {}
    aa_i_A, aa_i_B = 0, 0
    for aa_aln_A, aa_aln_B in zip(aligned_A, aligned_B):
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
    end_time = time.time()
    Logs.message(f"Alignment completed in {round(end_time - start_time, 2)} seconds.")
    return mapping


def superimpose(fixed_struct: Structure, moving_struct: Structure, overlay_method):
    """Align residues from each Structure, and calculate RMSD."""
    paired_res_id_mapping = align_structures(fixed_struct, moving_struct)
    fixed_atoms = []
    moving_atoms = []
    alpha_carbon_name = 'CA'
    skipped = 0
    for fixed_id, moving_id in paired_res_id_mapping.items():
        fixed_residue = next(
            res for res in fixed_struct.get_residues()
            if res.id[1] == fixed_id)
        moving_residue = next(
            res for res in moving_struct.get_residues()
            if res.id[1] == moving_id)

        new_fixed_atoms = []
        new_moving_atoms = []
        if overlay_method.name == 'ALPHA_CARBONS_ONLY':
            # Add alpha carbons.
            try:
                fixed_alpha_carbon = fixed_residue[alpha_carbon_name]
                moving_alpha_carbon = moving_residue[alpha_carbon_name]
            except KeyError:
                Logs.debug(f"Skipping Residue {fixed_id}, missing alpha carbon.")
                skipped += 1
                continue
            else:
                new_fixed_atoms.append(fixed_alpha_carbon)
                new_moving_atoms.append(moving_alpha_carbon)
        else:
            # Add all heavy atoms (Non hydrogens)
            for atom in fixed_residue.get_atoms():
                if not atom.name.startswith('H'):
                    new_fixed_atoms.append(atom)
            for atom in moving_residue.get_atoms():
                if not atom.name.startswith('H'):
                    new_moving_atoms.append(atom)

        if len(new_moving_atoms) != len(new_fixed_atoms):
            # We can skip residues with differing atom counts.
            # This is an issue with Heavy atom alignment methods.
            skipped += 1
            continue
        fixed_atoms.extend(new_fixed_atoms)
        moving_atoms.extend(new_moving_atoms)

    assert len(moving_atoms) == len(fixed_atoms), f"{len(moving_atoms)} != {len(fixed_atoms)}"
    paired_residue_count = len(paired_res_id_mapping) - skipped
    superimposer = Superimposer()
    superimposer.set_atoms(fixed_atoms, moving_atoms)
    rms = round(superimposer.rms, 2)
    Logs.debug(f"RMSD: {rms}")
    paired_atom_count = len(fixed_atoms)
    return superimposer, paired_residue_count, paired_atom_count