
import tempfile
from Bio.PDB import PDBParser
import tempfile
from nanome.util import Logs
from nanome.api.structure import Complex

from . import utils

__all__ = ["superimpose_by_chain"]


def superimpose_by_chain(fixed_comp: Complex, fixed_chain_name: str, moving_comp: Complex, moving_chain_name: str, overlay_method):
    parser = PDBParser(QUIET=True)
    fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
    fixed_comp.io.to_pdb(fixed_pdb.name)
    fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
    fixed_chain = next(ch for ch in fixed_struct.get_chains() if ch.id == fixed_chain_name)

    moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
    moving_comp.io.to_pdb(moving_pdb.name)
    moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)
    moving_chain = next(
        ch for ch in moving_struct.get_chains()
        if ch.id == moving_chain_name)

    try:
        superimposer, paired_residue_count, paired_atom_count = utils.superimpose_structures(
            fixed_chain, moving_chain, overlay_method)
    except Exception:
        Logs.error(f"Superimposition failed for {moving_comp.full_name}")
        return None, None

    transform_matrix = utils.create_transform_matrix(superimposer)
    comp_data = utils.format_superimposer_data(
        superimposer, paired_residue_count, paired_atom_count, moving_chain_name)
    return transform_matrix, comp_data
