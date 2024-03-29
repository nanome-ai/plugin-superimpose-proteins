import itertools
import os
import tempfile
from Bio.PDB import Superimposer
from nanome.util import Logs

from .fpocket_client import FPocketClient
from .site_motif_client import SiteMotifClient
from . import utils


async def superimpose_by_binding_site(fixed_comp, moving_comps, fixed_binding_site_comp, extract_binding_sites, overlay_method: str, plugin_instance):
    fpocket_client = FPocketClient()
    sitemotif_client = SiteMotifClient()
    temp_dir = tempfile.TemporaryDirectory()

    fixed_binding_site_pdb = tempfile.NamedTemporaryFile(dir=temp_dir.name, suffix='.pdb')
    fixed_binding_site_comp.io.to_pdb(path=fixed_binding_site_pdb.name)
    fixed_binding_site_comp.position = fixed_comp.position
    fixed_binding_site_comp.rotation = fixed_comp.rotation
    fixed_binding_site_comp.locked = True
    fixed_pdb = fixed_binding_site_pdb.name

    pocket_residue_pdbs = []
    plugin_instance.update_submit_btn_text('Finding Pockets...')
    for moving_comp in moving_comps:
        fpocket_results = fpocket_client.run(moving_comp, temp_dir.name)
        pocket_pdbs = fpocket_client.get_pocket_pdb_files(fpocket_results)

        comp_residue_pdbs = []
        if pocket_pdbs:
            comp_residue_pdbs = utils.clean_fpocket_pdbs(pocket_pdbs, moving_comp)
        else:
            # If no pockets were found, use entire complex as pocket
            Logs.debug("Manually adding pdb to Fpocket results")
            filepath = os.path.join(fpocket_results, 'pockets', f'{moving_comp.index}_pocket1_atm.pdb')
            moving_comp.io.to_pdb(path=filepath)
            comp_residue_pdbs = [filepath]
            moving_comp.whole_structure_alignment = True
        pocket_residue_pdbs.extend(comp_residue_pdbs)

    align_output_file = os.path.join(temp_dir.name, 'align_output.txt')
    plugin_instance.update_submit_btn_text('Aligning Pockets...')
    sitemotif_client.run(fixed_pdb, pocket_residue_pdbs, align_output_file)

    fixed_binding_site_pdb.close()
    output_data = {}
    binding_site_comps = []
    if extract_binding_sites:
        binding_site_comps.append(fixed_binding_site_comp)

    for moving_comp in moving_comps:
        pdb1, _, alignment = sitemotif_client.find_match(moving_comp.index, align_output_file)
        if os.path.basename(fixed_pdb) == pdb1:
            comp1 = fixed_comp
            comp2 = moving_comp
        else:
            comp1 = moving_comp
            comp2 = fixed_comp

        comp1_atoms, comp2_atoms = sitemotif_client.parse_residue_pairs(comp1, comp2, alignment, overlay_method)
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
        paired_residue_count = len(set([atm.residue for atm in comp1_atoms]))
        rmsd_results = utils.format_superimposer_data(superimposer, paired_residue_count, paired_atom_count)
        transform_matrix = utils.create_transform_matrix(superimposer)
        output_data[moving_comp.index] = (transform_matrix, rmsd_results)

        # If we're using a pocket within a larger struct, Create a separate complex for the binding site.
        whole_structure_alignment = any([
            getattr(moving_comp, 'whole_structure_alignment', False),
            'binding site' in moving_comp.full_name.lower()
        ])
        if not whole_structure_alignment and extract_binding_sites:
            if moving_comp == comp1:
                comp_atoms = comp1_atoms
            else:
                comp_atoms = comp2_atoms
            binding_site_residues = set(atom.residue for atom in comp_atoms)
            comp_name = f'{moving_comp.name} Run {plugin_instance.run_index} binding site'
            binding_site_comp = utils.extract_binding_site(
                moving_comp, binding_site_residues, comp_name=comp_name)
            binding_site_comp.position = moving_comp.position
            binding_site_comp.rotation = moving_comp.rotation
            binding_site_comp.locked = True
            binding_site_comps.append(binding_site_comp)

    temp_dir.cleanup()
    created_binding_site_comps = []
    if binding_site_comps:
        created_binding_site_comps = await plugin_instance.add_to_workspace(binding_site_comps)

    original_comps = [fixed_comp] + moving_comps
    for old_bsc, new_bsc, original_comp in zip(binding_site_comps, created_binding_site_comps, original_comps):
        new_bsc.position = old_bsc.position
        new_bsc.rotation = old_bsc.rotation
        try:
            transform_matrix = output_data[original_comp.index][0]
        except KeyError:
            continue
        for comp_atom in new_bsc.atoms:
            new_position = transform_matrix * comp_atom.position
            comp_atom.position = new_position
            new_bsc.boxed = False
        new_bsc.set_surface_needs_redraw()
    plugin_instance.update_structures_shallow(created_binding_site_comps)
    plugin_instance.update_structures_shallow(itertools.chain(*[comp.atoms for comp in created_binding_site_comps]))
    return output_data
