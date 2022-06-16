import copy
import nanome
import os
import tempfile
import time

from Bio.PDB.Structure import Structure
from Bio.PDB import Superimposer, PDBParser
from Bio import pairwise2
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import substitution_matrices
from itertools import chain
from scipy.spatial import KDTree
from nanome.util import Logs, async_callback, Matrix, ComplexUtils
from nanome.api.structure import Complex
from .enums import OverlayMethodEnum
from .menu import MainMenu
from .fpocket_client import FPocketClient
from .site_motif_client import SiteMotifClient
import sys


def extract_binding_site(comp, binding_site_atoms):
    # Copy comp, and remove all residues that are not part of the binding site
    orig_recursion_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(100000)
    new_comp = copy.deepcopy(comp)
    sys.setrecursionlimit(orig_recursion_limit)
    new_comp.name = f'{comp.name} binding site'
    new_comp.index = -1
    new_comp.set_surface_needs_redraw()

    binding_site_atom_indices = iter(a.index for a in binding_site_atoms)
    comp_residues = list(new_comp.residues)
    for i in range(len(comp_residues) - 1, -1, -1):
        res = comp_residues[i]
        binding_atoms_in_res = [a for a in res.atoms if a.index in binding_site_atom_indices]
        if not binding_atoms_in_res:
            res.chain.remove_residue(res)
        else:
            for atom in binding_atoms_in_res:
                atom.index = -1
    return new_comp


class SuperimposePlugin(nanome.AsyncPluginInstance):

    def start(self):
        self.menu = MainMenu(self)
        self.temp_dir = tempfile.TemporaryDirectory()
        self.complexes = []

    @async_callback
    async def on_run(self):
        self.menu.enabled = True
        self.complexes = await self.request_complex_list()
        self.menu.render(complexes=self.complexes, force_enable=True)

    @async_callback
    async def on_complex_list_updated(self, complexes):
        self.complexes = complexes
        await self.menu.render(complexes=self.complexes)

    @async_callback
    async def on_complex_added(self):
        self.complexes = await self.request_complex_list()
        await self.menu.render(complexes=self.complexes)

    @async_callback
    async def on_complex_removed(self):
        self.complexes = await self.request_complex_list()
        await self.menu.render(complexes=self.complexes)

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
            parser = PDBParser(QUIET=True)
            fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
            moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
            fixed_comp.io.to_pdb(fixed_pdb.name)
            moving_comp.io.to_pdb(moving_pdb.name)
            fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
            moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)

            try:
                superimposer, paired_residue_count, paired_atom_count = await self.superimpose(
                    fixed_struct, moving_struct, overlay_method)
            except Exception:
                Logs.error(f"Superimposition failed for {moving_comp.full_name}")
                continue
            rmsd_results[moving_comp.full_name] = self.format_superimposer_data(superimposer, paired_residue_count, paired_atom_count)
            # Use matrix to transform moving atoms to new position
            transform_matrix = self.create_transform_matrix(superimposer)
            for comp_atom in moving_comp.atoms:
                comp_atom.position = transform_matrix * comp_atom.position
            moving_comp.set_surface_needs_redraw()
            moving_comp.locked = True
            moving_comp.boxed = False
            comps_to_update.append(moving_comp)
            self.update_loading_bar(i + 1, comp_count)

        await self.update_structures_deep(comps_to_update)
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
            self, fixed_index: int, ligand_name: str, moving_indices: list, site_size=4.5):
        # Select the binding site on the fixed_index.
        updated_complexes = await self.request_complexes([fixed_index, *moving_indices])
        fixed_comp = updated_complexes[0]
        moving_comp_list = updated_complexes[1:]
        fixed_binding_site_atoms = await self.get_binding_site_atoms(fixed_comp, ligand_name, site_size)
        fixed_binding_site_comp = extract_binding_site(fixed_comp, fixed_binding_site_atoms)

        fixed_binding_site_pdb = tempfile.NamedTemporaryFile(dir=self.temp_dir.name, suffix='.pdb')
        fixed_binding_site_comp.io.to_pdb(path=fixed_binding_site_pdb.name)

        fpocket_client = FPocketClient()
        sitemotif_client = SiteMotifClient()
        for moving_comp in moving_comp_list:
            print(f"Calculating binding site for {moving_comp.full_name}")
            fpocket_results = fpocket_client.run(moving_comp, self.temp_dir.name)
            pocket_pdbs = fpocket_client.get_pocket_pdb_files(fpocket_results)
            matching_pocket = sitemotif_client.find_match(fixed_binding_site_pdb.name, pocket_pdbs)
        return {}

    @staticmethod
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

    async def superimpose(self, fixed_struct: Structure, moving_struct: Structure, overlay_method):
        """Align residues from each Structure, and calculate RMS"""
        paired_res_id_mapping = self.align_structures(fixed_struct, moving_struct)
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
            if overlay_method == OverlayMethodEnum.ALPHA_CARBONS_ONLY:
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

    def format_superimposer_data(self, superimposer: Superimposer, paired_residue_count: int, paired_atom_count: int, chain_name=''):
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

    async def get_binding_site_atoms(self, target_reference: Complex, ligand_name: str, site_size=4.5):
        """Identify atoms in the active site around a ligand."""
        mol = next(
            mol for i, mol in enumerate(target_reference.molecules)
            if i == target_reference.current_frame)
        target_ligands = await mol.get_ligands()
        ligand = next(ligand for ligand in target_ligands if ligand.name == ligand_name)
        # Use KDTree to find target atoms within site_size radius of ligand atoms
        ligand_atoms = chain(*[res.atoms for res in ligand.residues])
        binding_site_atoms = self.calculate_binding_site_atoms(target_reference, ligand_atoms)
        return binding_site_atoms

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
        Logs.message(f"{len(binding_site_atoms)} atoms identified in binding site.")
        return binding_site_atoms

    def align_structures(self, structA, structB):
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
        alns = pairwise2.align.globalds(
            sequence_A,
            sequence_B,
            substitution_matrices.load("BLOSUM62"),
            one_alignment_only=True,
            open=-10.0,
            extend=-0.5,
            penalize_end_gaps=(False, False))
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

    def update_loading_bar(self, current, total):
        self.menu.update_loading_bar(current, total)


def main():
    default_description = 'Superimpose two or more structures'
    description = os.environ.get("PLUGIN_DESCRIPTION", "") or default_description
    plugin = nanome.Plugin('Superimpose', description, 'alignment', False)
    plugin.set_plugin_class(SuperimposePlugin)
    plugin.run()


if __name__ == '__main__':
    main()
