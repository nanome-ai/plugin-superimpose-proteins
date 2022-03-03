import tempfile
from os import path
from tkinter import Menu
from Bio.PDB import Superimposer, PDBParser
from Bio import pairwise2
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import substitution_matrices

import nanome
from nanome.api import ui
from nanome.api.structure import Complex, Chain
from nanome.util import Logs, async_callback, Vector3, ComplexUtils

BASE_PATH = path.dirname(f'{path.realpath(__file__)}')
MENU_PATH = path.join(BASE_PATH, 'menu.json')


class RMSDMenu:

    def __init__(self, plugin_instance):
        super().__init__()
        self._menu = ui.Menu.io.from_json(MENU_PATH)
        self.plugin = plugin_instance

    @property
    def ln_struct1_complex(self):
        return self._menu.root.find_node('ln_struct1_complex')
    
    @property
    def ln_struct1_chain(self):
        return self._menu.root.find_node('ln_struct1_chain')

    @property
    def ln_struct2_complex(self):
        return self._menu.root.find_node('ln_struct2_complex')
    
    @property
    def ln_struct2_chain(self):
        return self._menu.root.find_node('ln_struct2_chain')

    @async_callback
    async def render(self, complexes=None, default_values=False):
        complexes = complexes or []
        self.complexes = complexes

        # for comp in self.complexes:
        #     comp.register_complex_updated_callback(self.on_complex_updated)

        # Find the first complex with selected atoms, and make that the default.
        # I guess that works for now.
        default_complex = None
        if complexes:
            default_complex = complexes[0]

        self.display_structures(complexes, self.ln_struct1_complex, default_structure=default_complex)
        if default_complex:
            await self.display_chains(default_complex, self.ln_struct1_chain)

        self.display_structures(complexes, self.ln_struct2_complex, default_structure=default_complex)
        if default_complex:
            await self.display_chains(default_complex, self.ln_struct2_chain)

        # self.display_structures(complexes, self.ln_ligands)
        # Determine whether we should currently be showing the ligand dropdown.
        # enable_ligands_node = self.btn_show_all_interactions.selected
        # self.toggle_ln_ligands_visibility(enable_ligands_node)
        # self.dd_complexes = self.ln_complexes.get_content()
        # self.dd_ligands = self.ln_ligands.get_content()
        # self.dd_ligands.register_item_clicked_callback(self.update_dropdown)
        # self.dd_complexes.register_item_clicked_callback(self.toggle_complex)
        self.plugin.update_menu(self._menu)

    def display_structures(self, complexes, layoutnode, default_structure=False):
        """Create dropdown of complexes, and add to provided layoutnode."""
        dropdown_items = self.create_structure_dropdown_items(complexes)
        dropdown = ui.Dropdown()
        dropdown.max_displayed_items = 12
        dropdown.items = dropdown_items

        # set default item selected.
        if default_structure:
            for ddi in dropdown.items:
                select_ddi = False
                if isinstance(default_structure, Complex):
                    select_ddi = ddi.complex.index == default_structure.index

                if select_ddi:
                    ddi.selected = True
                    break

        layoutnode.set_content(dropdown)
        self.plugin.update_node(layoutnode)
    
    async def display_chains(self, complex, layoutnode):
        """Create dropdown of chains, and add to provided layoutnode."""
        dropdown = ui.Dropdown()
        dropdown_items = []
        if sum(1 for ch in complex.chains) == 0:
            # get deep complex
            complex = (await self.plugin.request_complexes([complex.index]))[0]
        chain_names = [ch.name for ch in complex.chains]
        for chain_name in chain_names:
            ddi = ui.DropdownItem(chain_name)
            dropdown_items.append(ddi)
        dropdown.items = dropdown_items
        layoutnode.set_content(dropdown)
        self.plugin.update_node(layoutnode)

    def create_structure_dropdown_items(self, complexes):
        """Generate list of buttons corresponding to provided complexes."""
        complex_ddis = []
        ddi_labels = set()

        for struct in complexes:
            struct_name = struct.full_name

            # Make sure we have a unique name for every structure
            ddi_label = struct_name
            if ddi_label in ddi_labels:
                num = 1
                while ddi_label in ddi_labels:
                    ddi_label = f'{struct_name} {{{num}}}'
                    num += 1

            ddi_labels.add(ddi_label)
            ddi = ui.DropdownItem(ddi_label)
            ddi.complex = struct
            complex_ddis.append(ddi)

        return complex_ddis


class RMSDV2(nanome.AsyncPluginInstance):

    def start(self):
        self.menu = RMSDMenu(self)

    @async_callback
    async def on_run(self):
        self.menu.enabled = True
        complexes = await self.request_complex_list()
        Logs.message('RMSDV2 Run.')
        self.menu.render(complexes=complexes)
        return
        Logs.message("Retrieving complexes")
        complexes = await self.request_complex_list()
        complexes = await self.request_complexes([comp.index for comp in complexes])

        Logs.message("Parsing Complexes to BioPython Structures.")
        parser = PDBParser(QUIET=True)
        fixed_comp = complexes[0]
        moving_comp = complexes[1]
        ComplexUtils.align_to(moving_comp, fixed_comp)

        fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        fixed_comp.io.to_pdb(fixed_pdb.name)
        moving_comp.io.to_pdb(moving_pdb.name)
        fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
        moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)
        # Logs.message("Selecting Chains")
        # fixed_chain = next(fixed_struct.get_chains())
        # moving_chain = next(moving_struct.get_chains())
        # fixed_chain_name = 'A'
        # moving_chain_name = 'A'

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

        Logs.message("Superimposing Structures.")
        superimposer = Superimposer()
        superimposer.set_atoms(fixed_atoms, moving_atoms)
        superimposer.apply([*fixed_atoms, *moving_atoms])
        Logs.message(f'RMSD: {superimposer.rms}')
        Logs.message("Updating Workspace")
        # Apply changes to Workspace
        for comp_atom in moving_comp.atoms:
            struc_atom = next((a for a in moving_atoms if comp_atom.serial == a.serial_number), None)
            if not struc_atom:
                continue
            comp_atom.position = Vector3(*struc_atom.coord)

        for comp_atom in fixed_comp.atoms:
            struc_atom = next((a for a in moving_atoms if comp_atom.serial == a.serial_number), None)
            if not struc_atom:
                continue
            comp_atom.position = Vector3(*struc_atom.coord)
        await self.update_structures_deep([moving_struct])

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
