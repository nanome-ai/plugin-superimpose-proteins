import tempfile
from os import path
from tkinter import Menu
from Bio.PDB import Superimposer, PDBParser
from Bio import pairwise2
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import substitution_matrices
from functools import partial
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

    @property
    def btn_submit(self):
        return self._menu.root.find_node('ln_submit').get_content()

    @async_callback
    async def render(self, complexes=None, default_values=False):
        complexes = complexes or []
        self.complexes = complexes

        self.display_structures(complexes, self.ln_struct1_complex)
        dd_struct1_complex = self.ln_struct1_complex.get_content()
        dd_struct1_complex.register_item_clicked_callback(
            partial(self.update_chain_dropdown, self.ln_struct1_chain))
        
        self.display_structures(complexes, self.ln_struct2_complex)
        dd_struct2_complex = self.ln_struct2_complex.get_content()
        dd_struct2_complex.register_item_clicked_callback(
            partial(self.update_chain_dropdown, self.ln_struct2_chain))
        self.btn_submit.register_pressed_callback(self.submit)
        self.plugin.update_menu(self._menu)

    @async_callback
    async def update_chain_dropdown(self, ln_chain_dropdown, complex_dropdown, complex_dd_item):
        """Update chain dropdown to reflect changes in complex."""
        comp = complex_dd_item.complex
        Logs.message("Updating Chain dropdown")
        await self.display_chains(comp, ln_chain_dropdown)

    def display_structures(self, complexes, layoutnode, default_structure=False):
        """Create dropdown of complexes, and add to provided layoutnode."""
        dropdown_items = self.create_structure_dropdown_items(complexes)
        dropdown = ui.Dropdown()
        dropdown.max_displayed_items = len(dropdown_items)
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
        dropdown.max_displayed_items = len(dropdown_items)
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

    @async_callback
    async def submit(self, btn):
        dd_fixed_comp = self.ln_struct1_complex.get_content()
        dd_fixed_chain = self.ln_struct1_chain.get_content()
        dd_moving_comp = self.ln_struct2_complex.get_content()
        dd_moving_chain = self.ln_struct2_chain.get_content()

        fixed_comp = next(
            item.complex for item in dd_fixed_comp.items if item.selected
        )
        fixed_chain = next(
            item.name for item in dd_fixed_chain.items
            if item.selected
        )
        moving_comp = next(
            item.complex for item in dd_moving_comp.items if item.selected
        )
        moving_chain = next(
            item.name for item in dd_moving_chain.items
            if item.selected
        )
        await self.plugin.superimpose(fixed_comp, fixed_chain, moving_comp, moving_chain)


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

    async def superimpose(self, fixed_comp, fixed_chain_name, moving_comp, moving_chain_name):
        fixed_comp, moving_comp = await self.request_complexes([fixed_comp.index, moving_comp.index])
        Logs.message("Parsing Complexes to BioPython Structures.")
        parser = PDBParser(QUIET=True)
        ComplexUtils.align_to(moving_comp, fixed_comp)

        fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        fixed_comp.io.to_pdb(fixed_pdb.name)
        moving_comp.io.to_pdb(moving_pdb.name)
        fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
        moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)

        Logs.message("Aligning Structures.")
        fixed_chain = next(ch for ch in fixed_struct.get_chains() if ch.id == fixed_chain_name)
        moving_chain = next(ch for ch in moving_struct.get_chains() if ch.id == moving_chain_name)
        mapping = self.align_sequences(fixed_chain, moving_chain)

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
