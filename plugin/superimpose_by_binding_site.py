import os
import tempfile
from Bio.PDB import Superimposer
from nanome.util import Logs, ComplexUtils
from nanome.api.structure import Complex
from nanome.api import structure

from .fpocket_client import FPocketClient
from .site_motif_client import SiteMotifClient
from . import utils


async def superimpose_by_binding_site(fixed_comp, moving_comps, fixed_binding_site_comp, plugin_instance):
    fpocket_client = FPocketClient()
    sitemotif_client = SiteMotifClient()
    temp_dir = tempfile.TemporaryDirectory()

    fixed_binding_site_pdb = tempfile.NamedTemporaryFile(dir=temp_dir.name, suffix='.pdb')
    fixed_binding_site_comp.io.to_pdb(path=fixed_binding_site_pdb.name)
    fixed_pdb = fixed_binding_site_pdb.name

    pocket_residue_pdbs = []
    plugin_instance.update_submit_btn_text('Finding Pockets...')
    for moving_comp in moving_comps:
        fpocket_results = fpocket_client.run(moving_comp, temp_dir.name)
        pocket_pdbs = fpocket_client.get_pocket_pdb_files(fpocket_results)
        comp_residue_pdbs = utils.clean_fpocket_pdbs(pocket_pdbs, moving_comp)
        pocket_residue_pdbs.extend(comp_residue_pdbs)

    align_output_file = os.path.join(temp_dir.name, 'align_output.txt')
    plugin_instance.update_submit_btn_text('Aligning Pockets...')
    sitemotif_client.run(fixed_pdb, pocket_residue_pdbs, align_output_file)
    
    fixed_binding_site_pdb.close()
    output_data = {}
    binding_site_comps = []
    for moving_comp in moving_comps:
        pdb1, _, alignment = sitemotif_client.find_match(moving_comp.index, align_output_file)
        if os.path.basename(fixed_pdb) == pdb1:
            comp1 = fixed_comp
            comp2 = moving_comp
        else:
            comp1 = moving_comp
            comp2 = fixed_comp

        comp1_atoms, comp2_atoms = sitemotif_client.parse_residue_pairs(comp1, comp2, alignment)
        comp1_bp_atoms = utils.convert_atoms_to_biopython(comp1_atoms)
        comp2_bp_atoms = utils.convert_atoms_to_biopython(comp2_atoms)
        superimposer = Superimposer()
        if comp1 == fixed_comp:
            superimposer.set_atoms(comp1_bp_atoms, comp2_bp_atoms)
        else:
            superimposer.set_atoms(comp2_bp_atoms, comp1_bp_atoms)

        rms = round(superimposer.rms, 2)
        Logs.debug(f"RMSD: {rms}")
        paired_atom_count = len(comp1_atoms)
        paired_residue_count = paired_atom_count
        rmsd_results = utils.format_superimposer_data(superimposer, paired_residue_count, paired_atom_count)
        transform_matrix = utils.create_transform_matrix(superimposer)
        output_data[moving_comp.index] = (transform_matrix, rmsd_results)

        # Create a separate complex for the binding site.
        if moving_comp == comp1:
            comp_atoms = comp1_atoms
        else:
            comp_atoms = comp2_atoms
        binding_site_residues = set(atom.residue for atom in comp_atoms)
        binding_site_comp = utils.extract_binding_site(moving_comp, binding_site_residues)
        binding_site_comp.position = moving_comp.position
        binding_site_comp.rotation = moving_comp.rotation
        # Make sure the binding aligns with its original position.
        # ComplexUtils.align_to(binding_site_comp, moving_comp)

        binding_site_comps.append(binding_site_comp)
        
    temp_dir.cleanup()
    created_binding_site_comps = await plugin_instance.add_to_workspace(binding_site_comps)
    for old_bsc, new_bsc, moving_comp in zip(binding_site_comps, created_binding_site_comps, moving_comps):
        new_bsc.position = old_bsc.position
        new_bsc.rotation = old_bsc.rotation

        transform_matrix = output_data[moving_comp.index][0]
        for comp_atom in new_bsc.atoms:
            new_position = transform_matrix * comp_atom.position
            comp_atom.position = new_position
    await plugin_instance.update_structures_deep(created_binding_site_comps)
    return output_data
