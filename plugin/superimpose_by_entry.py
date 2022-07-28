from Bio.PDB import PDBParser
import tempfile
from nanome.util import Logs
from nanome.api.structure import Complex

from . import utils

__all__ = ["superimpose_by_entry"]


def superimpose_by_entry(fixed_comp: Complex, moving_comp: Complex, overlay_method):
    parser = PDBParser(QUIET=True)
    fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
    moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
    fixed_comp.io.to_pdb(fixed_pdb.name)
    moving_comp.io.to_pdb(moving_pdb.name)
    fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
    moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)

    try:
        superimposer, paired_residue_count, paired_atom_count = utils.superimpose_structures(
            fixed_struct, moving_struct, overlay_method)
    except Exception:
        Logs.error(f"Superimposition failed for {moving_comp.full_name}")
        return None, None

    rmsd_results = utils.format_superimposer_data(superimposer, paired_residue_count, paired_atom_count)
    transform_matrix = utils.create_transform_matrix(superimposer)
    return transform_matrix, rmsd_results