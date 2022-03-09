import tempfile
from os import path
from Bio.PDB import Superimposer, PDBParser
from Bio import pairwise2
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import substitution_matrices
from functools import partial
import nanome
from nanome.api import ui
from nanome.api.structure import Complex
from nanome.util import Logs, async_callback, Matrix, ComplexUtils
from nanome.util.enums import NotificationTypes

BASE_PATH = path.dirname(f'{path.realpath(__file__)}')
MENU_PATH = path.join(BASE_PATH, 'menu2.json')


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
    def ln_main_list(self):
        return self._menu.root.find_node('ln_main_list')

    @property
    def ln_comparator_list(self):
        return self._menu.root.find_node('ln_comparator_list')

    @property
    def btn_submit(self):
        return self._menu.root.find_node('ln_submit').get_content()

    @property
    def lbl_rmsd_value(self):
        return self._menu.root.find_node('ln_rmsd_value').get_content()

    @async_callback
    async def render(self, complexes=None, default_values=False):
        complexes = complexes or []
        self.complexes = complexes

        self.display_structures(complexes, self.ln_main_list)
        # dd_struct1_complex = self.ln_main_list.get_content()
        # dd_struct1_complex.register_item_clicked_callback(
        #     partial(self.update_chain_dropdown, self.ln_struct1_chain))
        
        self.display_structures(complexes, self.ln_comparator_list)
        # dd_struct2_complex = self.ln_comparator_list.get_content()
        # dd_struct2_complex.register_item_clicked_callback(
        #     partial(self.update_chain_dropdown, self.ln_struct2_chain))
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
        self.populate_complex_list(complexes, layoutnode)
        # dropdown_items = self.create_structure_dropdown_items(complexes)
        # ui_list = ui.UIList()
        # dropdown.max_displayed_items = len(dropdown_items)
        # dropdown.items = dropdown_items

        # set default item selected.
        # if default_structure:
        #     for ddi in dropdown.items:
        #         select_ddi = False
        #         if isinstance(default_structure, Complex):
        #             select_ddi = ddi.complex.index == default_structure.index

        #         if select_ddi:
        #             ddi.selected = True
        #             break

        # layoutnode.set_content(ui_list)
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
        Logs.message("Submit button Pressed.")
        self.btn_submit.unusable = True
        self.plugin.update_content(self.btn_submit)
        dd_fixed_comp = self.ln_struct1_complex.get_content()
        dd_fixed_chain = self.ln_struct1_chain.get_content()
        dd_moving_comp = self.ln_struct2_complex.get_content()
        dd_moving_chain = self.ln_struct2_chain.get_content()

        fixed_comp = next(
            (item.complex for item in dd_fixed_comp.items if item.selected), None
        )
        fixed_chain = next(
            (item.name for item in dd_fixed_chain.items if item.selected), None
        )
        moving_comp = next(
            (item.complex for item in dd_moving_comp.items if item.selected), None
        )
        moving_chain = next(
            (item.name for item in dd_moving_chain.items if item.selected), None
        )
        if not any([fixed_comp, fixed_chain, moving_comp, moving_chain]):
            message = "Please select a structure and chain."
            Logs.error(message)
            self.plugin.send_notification(NotificationTypes.error, message)
            self.btn_submit.unusable = False
            self.plugin.update_content(self.btn_submit)
            return
        rmsd_result = await self.plugin.superimpose(fixed_comp, fixed_chain, moving_comp, moving_chain)

        self.btn_submit.unusable = False
        self.lbl_rmsd_value.text_value = rmsd_result
        self.plugin.update_content(self.lbl_rmsd_value, self.btn_submit)
        Logs.message("Superposition completed.")

    def populate_complex_list(self, complex_list, layoutnode):
        ui_list = ui.UIList()
        ui_list.display_rows = min(len(complex_list), 5)
        complex_ln_list = []
        for complex in complex_list:
            # Create button representing each complex in the workspace.
            item = ui.LayoutNode()
            btn_ln = item.create_child_node()
            btn = btn_ln.add_new_button()
            btn.text.active = True
            btn.text.value.set_all(complex.full_name)
            btn.complex = complex
            complex_ln_list.append(item)
            # btn.register_pressed_callback(self.select_complex)
        ui_list.items = complex_ln_list
        layoutnode.set_content(ui_list)
        self.plugin.update_node(layoutnode)

class RMSDV2(nanome.AsyncPluginInstance):

    def start(self):
        self.menu = RMSDMenu(self)

    @async_callback
    async def on_run(self):
        self.menu.enabled = True
        complexes = await self.request_complex_list()
        Logs.message('RMSDV2 Run.')
        self.menu.render(complexes=complexes)
    
    @async_callback
    async def on_complex_list_updated(self, complexes):
        self.menu.render(complexes=complexes)

    @async_callback
    async def on_complex_added(self):
        complexes = await self.request_complex_list()
        await self.menu.render(complexes=complexes, default_values=True)

    @async_callback
    async def on_complex_removed(self):
        complexes = await self.request_complex_list()
        await self.menu.render(complexes=complexes)

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

        # Collect aligned residues
        # Align Residues based on Alpha Carbon
        fixed_atoms = []
        moving_atoms = []
        alpha_carbon = 'CA'
        for fixed_residue in fixed_chain.get_residues():
            fixed_id = fixed_residue.id[1]
            if fixed_id in mapping:
                fixed_atoms.append(fixed_residue[alpha_carbon])
                moving_residue_serial = mapping[fixed_id] 
                moving_residue = next(
                    rez for rez in moving_chain.get_residues()
                    if rez.id[1] == moving_residue_serial)
                moving_atoms.append(moving_residue[alpha_carbon])
        assert len(moving_atoms) == len(fixed_atoms)
        Logs.message("Superimposing Structures.")
        superimposer = Superimposer()
        superimposer.set_atoms(fixed_atoms, moving_atoms)        
        # superimposer.apply(moving_struct.get_atoms())
        rms = superimposer.rms
        Logs.message(f'RMSD: {rms}')
        
        # convert numpy rot + tran matrices to 4x4 nanome matrix
        rot, tran = superimposer.rotran
        rot = rot.tolist()
        tran = tran.tolist()
        m = Matrix(4, 4)
        m[0][0:3] = rot[0]
        m[1][0:3] = rot[1]
        m[2][0:3] = rot[2]
        m[3][0:3] = tran
        m[3][3] = 1
        Logs.debug(f"Matrix m = {str(m)}")
        # transpose necessary because numpy and nanome matrices are opposite row/col
        m.transpose()
        # apply transformation to moving_comp
        moving_comp.set_surface_needs_redraw()
        for comp_atom in moving_comp.atoms:
            comp_atom.position = m * comp_atom.position
        await self.update_structures_deep([moving_comp])
        self.send_notification(NotificationTypes.message, f'RMSD: {rms}')
        return rms

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
