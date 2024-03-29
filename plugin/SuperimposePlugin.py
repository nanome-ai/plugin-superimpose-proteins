import nanome
import os
import tempfile
import itertools
from nanome.util import Logs, async_callback, ComplexUtils
from nanome.api.structure import Complex
from nanome.util import enums

from . import __version__
from . import utils
from .superimpose_by_chain import superimpose_by_chain
from .superimpose_by_entry import superimpose_by_entry
from .superimpose_by_binding_site import superimpose_by_binding_site
from .superimpose_by_selection import superimpose_by_selection
from .menu import MainMenu

PDBOPTIONS = Complex.io.PDBSaveOptions()
PDBOPTIONS.write_bonds = True


class SuperimposePlugin(nanome.AsyncPluginInstance):

    def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.menu = MainMenu(self)
        self.complexes = []
        self.run_index = 0

    def on_stop(self):
        self.temp_dir.cleanup()

    @async_callback
    async def on_run(self):
        self.menu.open_menu()
        if not self.complexes:
            self.set_plugin_list_button(enums.PluginListButtonType.run, text='Loading...', usable=False)
            workspace = await self.request_workspace()
            self.complexes = workspace.complexes
        self.menu._menu.title = 'Superimpose Proteins'
        self.menu.render(force_enable=True)
        self.set_plugin_list_button(enums.PluginListButtonType.run, text='Run', usable=True)

    @async_callback
    async def on_complex_list_changed(self):
        workspace = await self.request_workspace()
        self.complexes = workspace.complexes
        await self.menu.render()

    async def superimpose_by_entry(self, fixed_comp_index, moving_comp_indices, overlay_method):
        updated_comps = await self.request_complexes([fixed_comp_index, *moving_comp_indices])
        self.run_index += 1
        fixed_comp = updated_comps[0]
        moving_comps = updated_comps[1:]
        fixed_comp.locked = True
        fixed_comp.boxed = False
        comps_to_update = [fixed_comp]
        rmsd_results = {}
        comp_count = len(moving_comps)
        for i, moving_comp in enumerate(moving_comps):
            Logs.debug(f"Starting Structure {i + 1}")
            ComplexUtils.align_to(moving_comp, fixed_comp)
            transform_matrix, rmsd_data = superimpose_by_entry(fixed_comp, moving_comp, overlay_method)
            rmsd_results[moving_comp.full_name] = rmsd_data
            # Use matrix to transform moving atoms to new position
            for comp_atom in moving_comp.atoms:
                comp_atom.position = transform_matrix * comp_atom.position
            moving_comp.set_surface_needs_redraw()
            moving_comp.locked = True
            moving_comp.boxed = False
            comps_to_update.append(moving_comp)
            self.update_loading_bar(i + 1, comp_count)

        await self.update_structures_deep(comps_to_update)
        # Due to a bug in nanome-core, if a complex is unlocked, we need to
        # make a separate call to remove box from around complexes.
        self.update_structures_shallow(comps_to_update)

        # Update comps in stored complex list
        for i in range(len(self.complexes)):
            comp_index = self.complexes[i].index
            updated_comp = next(
                (updated_comp for updated_comp in comps_to_update
                 if updated_comp.index == comp_index), None)
            if updated_comp:
                self.complexes[i] = updated_comp
        return rmsd_results

    async def superimpose_by_chain(self, fixed_comp_index, fixed_chain_name, moving_comp_chain_list, overlay_method):
        self.run_index += 1
        moving_comp_indices = [item[0] for item in moving_comp_chain_list]
        updated_comps = await self.request_complexes([fixed_comp_index, *moving_comp_indices])
        fixed_comp = updated_comps[0]
        moving_comps = updated_comps[1:]
        fixed_comp.locked = True
        fixed_comp.boxed = False
        comps_to_update = [fixed_comp]
        comp_count = len(moving_comps)
        rmsd_results = {}
        for i, moving_comp in enumerate(moving_comps):
            Logs.debug(f"Superimposing Moving Complex {i + 1}")
            moving_chain_name = moving_comp_chain_list[i][1]
            ComplexUtils.align_to(moving_comp, fixed_comp)

            transform_matrix, rmsd_data = superimpose_by_chain(
                fixed_comp, fixed_chain_name, moving_comp, moving_chain_name, overlay_method)
            rmsd_results[moving_comp.full_name] = rmsd_data

            # apply transformation to moving_comp
            for comp_atom in moving_comp.atoms:
                comp_atom.position = transform_matrix * comp_atom.position

            moving_comp.locked = True
            moving_comp.boxed = False
            moving_comp.set_surface_needs_redraw()
            comps_to_update.append(moving_comp)
            self.update_loading_bar(i + 1, comp_count)

        await self.update_structures_deep(comps_to_update)
        # Due to a bug in nanome-core, if a complex is unlocked, we need to
        # make a separate call to remove box from around complexes.
        self.update_structures_shallow(comps_to_update)

        # Update comps in stored complex list
        for i in range(len(self.complexes)):
            comp_index = self.complexes[i].index
            updated_comp = next(
                (updated_comp for updated_comp in comps_to_update
                 if updated_comp.index == comp_index), None)
            if updated_comp:
                self.complexes[i] = updated_comp
        return rmsd_results

    async def superimpose_by_binding_site(
            self, fixed_index: int, ligand_index: int, moving_indices: list,
            overlay_method, extract_binding_sites=False):
        self.run_index += 1
        updated_complexes = await self.request_complexes([fixed_index, *moving_indices])
        fixed_comp = updated_complexes[0]
        moving_comp_list = updated_complexes[1:]
        # If -1 passed as ligand index, use the whole complex as the binding site
        if ligand_index == -1:
            fixed_binding_site_residues = list(fixed_comp.residues)
        else:
            site_size = 5  # Number of angstroms around the ligand to consider as the binding site
            fixed_binding_site_residues = await self.get_binding_site_residues(fixed_comp, ligand_index, site_size)

        fixed_bs_name = f'{fixed_comp.name} Run {self.run_index} binding site'
        fixed_binding_site_comp = utils.extract_binding_site(fixed_comp, fixed_binding_site_residues, fixed_bs_name)
        fixed_binding_site_pdb = tempfile.NamedTemporaryFile(dir=self.temp_dir.name, suffix='.pdb')
        fixed_binding_site_comp.io.to_pdb(path=fixed_binding_site_pdb.name)

        fixed_comp.locked = True
        fixed_comp.boxed = False
        comps_to_update = [fixed_comp]
        rmsd_results = {}
        for moving_comp in moving_comp_list:
            ComplexUtils.align_to(moving_comp, fixed_comp)

        superimpose_data = await superimpose_by_binding_site(
            fixed_comp, moving_comp_list, fixed_binding_site_comp,
            extract_binding_sites, overlay_method.name, self)
        fixed_binding_site_pdb.close()
        self.update_submit_btn_text('Updating Workspace...')
        for comp_index, data in superimpose_data.items():
            moving_comp = next(comp for comp in moving_comp_list if comp.index == comp_index)
            transform_matrix, rmsd_data = data
            rmsd_results[moving_comp.full_name] = rmsd_data
            moving_comp = next(comp for comp in moving_comp_list if comp.index == comp_index)
            for comp_atom in moving_comp.atoms:
                new_position = transform_matrix * comp_atom.position
                comp_atom.position = new_position
            moving_comp.locked = True
            moving_comp.boxed = False
            moving_comp.set_surface_needs_redraw()
            comps_to_update.append(moving_comp)

        self.update_structures_shallow(comps_to_update)
        self.update_structures_shallow(itertools.chain(*[comp.atoms for comp in comps_to_update]))
        # Due to a bug in nanome-core, if a complex is unlocked, we need to
        # make a second call to remove box from around complexes.
        self.update_structures_shallow(comps_to_update)
        # Update comps in stored complex list
        for i in range(len(self.complexes)):
            comp_index = self.complexes[i].index
            updated_comp = next(
                (updated_comp for updated_comp in comps_to_update
                 if updated_comp.index == comp_index), None)
            if updated_comp:
                self.complexes[i] = updated_comp
        return rmsd_results

    async def get_binding_site_residues(self, target_reference: Complex, ligand_index: int, site_size=5):
        """Identify atoms in the active site around a ligand."""
        mol = next(
            mol for i, mol in enumerate(target_reference.molecules)
            if i == target_reference.current_frame)
        target_ligands = await mol.get_ligands()
        ligand = target_ligands[ligand_index]
        ligand_atoms = itertools.chain(*[res.atoms for res in ligand.residues])
        binding_site_atoms = utils.calculate_binding_site_atoms(target_reference, ligand_atoms, site_size)
        residue_set = set()
        for atom in binding_site_atoms:
            residue_set.add(atom.residue)
        return residue_set

    def update_loading_bar(self, current, total):
        self.menu.update_loading_bar(current, total)

    def update_submit_btn_text(self, new_text):
        if getattr(self, 'menu', None):
            self.menu.update_submit_btn_text(new_text)

    async def superimpose_by_selection(self, fixed_comp_index, moving_comp_indices, overlay_method):
        updated_comps = await self.request_complexes([fixed_comp_index, *moving_comp_indices])
        # Validate that each complex has at least one selected atom
        fixed_comp, *moving_comps = updated_comps
        if not any([atom.selected for atom in fixed_comp.atoms]):
            msg = f"Fixed Structure {fixed_comp.full_name} has no selected atoms"
            self.send_notification(enums.NotificationTypes.error, msg)
            Logs.warning(msg)
            raise ValueError(f"Make sure all structures have selected atoms.")
        self.run_index += 1
        fixed_comp.locked = True
        fixed_comp.boxed = False
        comps_to_update = [fixed_comp]
        rmsd_results = {}
        comp_count = len(moving_comps)
        for i, moving_comp in enumerate(moving_comps):
            Logs.debug(f"Starting Structure {i + 1}")
            ComplexUtils.align_to(moving_comp, fixed_comp)
            transform_matrix, rmsd_data = superimpose_by_selection(fixed_comp, moving_comp, overlay_method)
            rmsd_results[moving_comp.full_name] = rmsd_data
            # Use matrix to transform moving atoms to new position
            for comp_atom in moving_comp.atoms:
                comp_atom.position = transform_matrix * comp_atom.position
            moving_comp.set_surface_needs_redraw()
            moving_comp.locked = True
            moving_comp.boxed = False
            comps_to_update.append(moving_comp)
            self.update_loading_bar(i + 1, comp_count)

        await self.update_structures_deep(comps_to_update)
        # Due to a bug in nanome-core, if a complex is unlocked, we need to
        # make a separate call to remove box from around complexes.
        self.update_structures_shallow(comps_to_update)

        # Update comps in stored complex list
        for i in range(len(self.complexes)):
            comp_index = self.complexes[i].index
            updated_comp = next(
                (updated_comp for updated_comp in comps_to_update
                 if updated_comp.index == comp_index), None)
            if updated_comp:
                self.complexes[i] = updated_comp
        return rmsd_results


def main():
    plugin_title = 'Superimpose Proteins'
    default_description = f'Overlay multiple proteins in 3D space for visual comparison and calculate Root Mean Square Deviation (RMSD) values for a range of similarity from identical (RMSD =0) to very different (RMSD>10).\n\nVersion {__version__}'
    description = os.environ.get("PLUGIN_DESCRIPTION", "") or default_description
    plugin = nanome.Plugin(plugin_title, description, 'Alignment', False, version=__version__)
    plugin.set_plugin_class(SuperimposePlugin)
    plugin.run()


if __name__ == '__main__':
    main()
