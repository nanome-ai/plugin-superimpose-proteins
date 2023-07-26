from Bio.PDB import PDBParser
import tempfile
from nanome.util import Logs
from nanome.api.structure import Complex

from . import utils

__all__ = ["superimpose_by_selection"]


def superimpose_by_selection(fixed_comp: Complex, moving_comp: Complex, overlay_method):
    parser = PDBParser(QUIET=True)
    fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
    moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
    fixed_comp.io.to_pdb(fixed_pdb.name)
    moving_comp.io.to_pdb(moving_pdb.name)

    # Reduce structures to only selected atoms
    fixed_selected_residues = [res for res in fixed_comp.residues if any([atom.selected for atom in res.atoms])]
    moving_selected_residues = [res for res in moving_comp.residues if any([atom.selected for atom in res.atoms])]
    fixed_selected_comp = utils.extract_binding_site(fixed_comp, fixed_selected_residues, comp_name=fixed_comp.name)
    moving_selected_comp = utils.extract_binding_site(moving_comp, moving_selected_residues, comp_name=moving_comp.name)

    fixed_struct = parser.get_structure(fixed_selected_comp.full_name, fixed_pdb.name)
    moving_struct = parser.get_structure(moving_selected_comp.full_name, moving_pdb.name)

    try:
        superimposer, paired_residue_count, paired_atom_count = utils.superimpose_structures(
            fixed_struct, moving_struct, overlay_method)
    except Exception:
        Logs.warning("Failed to superimpose entries")
        return None, None

    rmsd_results = utils.format_superimposer_data(superimposer, paired_residue_count, paired_atom_count)
    transform_matrix = utils.create_transform_matrix(superimposer)
    return transform_matrix, rmsd_results
