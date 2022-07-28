import nanome
import os
import sys
import tempfile
import time
from Bio.PDB import Superimposer, PDBParser
from itertools import chain
from scipy.spatial import KDTree
from nanome.util import Logs, async_callback, Matrix, ComplexUtils
from nanome.api.structure import Complex
from nanome.util import enums

from . import __version__
from . import algorithms
from . import utils
from .menu import MainMenu
from .fpocket_client import FPocketClient
from .site_motif_client import SiteMotifClient

PDBOPTIONS = Complex.io.PDBSaveOptions()
PDBOPTIONS.write_bonds = True


class SuperimposePlugin(nanome.AsyncPluginInstance):

    def start(self):
        self.menu = MainMenu(self)
        self.temp_dir = tempfile.TemporaryDirectory()
        self.complexes = []

    def on_stop(self):
        self.temp_dir.cleanup()

    @async_callback
    async def on_run(self):
        self.menu.enabled = True
        if not self.complexes:
            self.set_plugin_list_button(enums.PluginListButtonType.run, text='Loading...', usable=False)
            workspace = await self.request_workspace()
            self.complexes = workspace.complexes
        self.menu.render(force_enable=True)
        self.set_plugin_list_button(enums.PluginListButtonType.run, text='Run', usable=True)

    @async_callback
    async def on_complex_list_updated(self, complexes):
        self.complexes = complexes
        await self.menu.render()

    @async_callback
    async def on_complex_added(self):
        workspace = await self.request_workspace()
        self.complexes = workspace.complexes
        await self.menu.render()

    @async_callback
    async def on_complex_removed(self):
        workspace = await self.request_workspace()
        self.complexes = workspace.complexes
        await self.menu.render()

    async def superimpose_by_entry(self, fixed_comp_index, moving_comp_indices, overlay_method):
        updated_comps = await self.request_complexes([fixed_comp_index, *moving_comp_indices])
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
            transform_matrix, rmsd_data = algorithms.superimpose_by_entry(fixed_comp, moving_comp, overlay_method)
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
        moving_comp_indices = [item[0] for item in moving_comp_chain_list]
        updated_comps = await self.request_complexes([fixed_comp_index, *moving_comp_indices])
        fixed_comp = updated_comps[0]
        moving_comps = updated_comps[1:]

        fixed_comp.locked = True
        fixed_comp.boxed = False
        comps_to_update = [fixed_comp]
        parser = PDBParser(QUIET=True)

        fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        fixed_comp.io.to_pdb(fixed_pdb.name)
        fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
        fixed_chain = next(ch for ch in fixed_struct.get_chains() if ch.id == fixed_chain_name)
        comp_count = len(moving_comps)
        results = {}
        for i, moving_comp in enumerate(moving_comps):
            Logs.debug(f"Superimposing Moving Complex {i + 1}")
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
                    fixed_chain, moving_chain, overlay_method)
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
        return results

    async def superimpose_by_binding_site(
            self, fixed_index: int, ligand_name: str, moving_indices: list, site_size=5):
        # Select the binding site on the fixed_index.
        updated_complexes = await self.request_complexes([fixed_index, *moving_indices])
        fixed_comp = updated_complexes[0]
        moving_comp_list = updated_complexes[1:]
        fixed_binding_site_residues = await self.get_binding_site_residues(fixed_comp, ligand_name, site_size)
        fixed_binding_site_comp = utils.extract_binding_site(fixed_comp, fixed_binding_site_residues)

        fixed_binding_site_pdb = tempfile.NamedTemporaryFile(dir=self.temp_dir.name, suffix='.pdb')
        fixed_binding_site_comp.io.to_pdb(path=fixed_binding_site_pdb.name)

        fpocket_client = FPocketClient()
        sitemotif_client = SiteMotifClient()

        fixed_comp.locked = True
        comps_to_update = [fixed_comp]
        comp_count = len(moving_comp_list)
        rmsd_results = {}
        for i, moving_comp in enumerate(moving_comp_list):
            ComplexUtils.align_to(moving_comp, fixed_comp)
            Logs.debug(f"Identifying binding sites for moving comp {i + 1}")
            fpocket_results = fpocket_client.run(moving_comp, self.temp_dir.name)
            pocket_pdbs = fpocket_client.get_pocket_pdb_files(fpocket_results)
            pocket_residue_pdbs = utils.clean_fpocket_pdbs(pocket_pdbs, moving_comp)

            fixed_pdb = fixed_binding_site_pdb.name
            pdb1, pdb2, alignment = sitemotif_client.find_match(fixed_pdb, pocket_residue_pdbs)
            if fixed_pdb == pdb1:
                comp1 = fixed_comp
                comp2 = moving_comp
            else:
                comp1 = moving_comp
                comp2 = fixed_comp

            # Get biopython representation of alpha carbon atoms, and pass to superimposer
            comp1_atoms, comp2_atoms = sitemotif_client.parse_residue_pairs(comp1, comp2, alignment)
            comp1_bp_atoms = utils.convert_atoms_to_biopython(comp1_atoms)
            comp2_bp_atoms = utils.convert_atoms_to_biopython(comp2_atoms)
            superimposer = Superimposer()
            if comp1 == fixed_comp:
                superimposer.set_atoms(comp1_bp_atoms, comp2_bp_atoms)
            else:
                superimposer.set_atoms(comp2_bp_atoms, comp1_bp_atoms)
            # Select all alpha carbons used in the superimpose
            comp1_atoms_selected = 0
            comp2_atoms_selected = 0
            for atom in comp1.atoms:
                atom.selected = atom in comp1_atoms
                if atom.selected:
                    comp1_atoms_selected += 1
                    atom.atom_mode = enums.AtomRenderingMode.BallStick
                    atom.residue.ribboned = False
            for atom in comp2.atoms:
                atom.selected = atom in comp2_atoms
                if atom.selected:
                    comp2_atoms_selected += 1
                    atom.atom_mode = enums.AtomRenderingMode.BallStick
                    atom.residue.ribboned = False

            rms = round(superimposer.rms, 2)
            Logs.debug(f"Comp1 atoms selected: {comp1_atoms_selected}")
            Logs.debug(f"Comp2 atoms selected: {comp2_atoms_selected}")
            Logs.debug(f"RMSD: {rms}")
            paired_atom_count = len(comp1_atoms)
            paired_residue_count = paired_atom_count
            rmsd_results[moving_comp.full_name] = self.format_superimposer_data(superimposer, paired_residue_count, paired_atom_count)
            transform_matrix = self.create_transform_matrix(superimposer)
            for comp_atom in moving_comp.atoms:
                new_position = transform_matrix * comp_atom.position
                comp_atom.position = new_position
            moving_comp.locked = True
            moving_comp.boxed = False
            moving_comp.set_surface_needs_redraw()
            comps_to_update.append(moving_comp)
            self.update_loading_bar(i + 1, comp_count)

        await self.update_structures_deep(comps_to_update)
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

    async def get_binding_site_residues(self, target_reference: Complex, ligand_name: str, site_size=4.5):
        """Identify atoms in the active site around a ligand."""
        mol = next(
            mol for i, mol in enumerate(target_reference.molecules)
            if i == target_reference.current_frame)
        target_ligands = await mol.get_ligands()
        ligand = next(ligand for ligand in target_ligands if ligand.name == ligand_name)
        # Use KDTree to find target atoms within site_size radius of ligand atoms
        ligand_atoms = chain(*[res.atoms for res in ligand.residues])
        binding_site_atoms = self.calculate_binding_site_atoms(target_reference, ligand_atoms)
        residue_set = set()
        for atom in binding_site_atoms:
            residue_set.add(atom.residue)
        return residue_set

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
        return binding_site_atoms

    def update_loading_bar(self, current, total):
        self.menu.update_loading_bar(current, total)


def main():
    plugin_title = 'Superimpose Proteins'
    default_description = f'Overlay multiple proteins in 3D space for visual comparison and calculate Root Mean Square Deviation (RMSD) values for a range of similarity from identical (RMSD =0) to very different (RMSD>10).\n\nVersion {__version__}'
    description = os.environ.get("PLUGIN_DESCRIPTION", "") or default_description
    plugin = nanome.Plugin(plugin_title, description, 'Alignment', False)
    plugin.set_plugin_class(SuperimposePlugin)
    plugin.run()


if __name__ == '__main__':
    main()
