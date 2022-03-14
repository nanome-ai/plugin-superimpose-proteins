from itertools import chain
import tempfile
from os import path
from Bio.PDB import Superimposer, PDBParser
from Bio import pairwise2
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import substitution_matrices
import nanome
from nanome.api import ui
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
    def ln_fixed_struct(self):
        return self._menu.root.find_node('ln_fixed_struct')

    @property
    def ln_moving_structs(self):
        return self._menu.root.find_node('ln_moving_structs')

    @property
    def ln_fixed_selection(self):
        return self._menu.root.find_node('ln_fixed_selection')

    @property
    def ln_moving_selections(self):
        return self._menu.root.find_node('ln_moving_selection')

    @property
    def ln_main_list(self):
        return self._menu.root.find_node('ln_main_list')

    @property
    def btn_global(self):
        return self._menu.root.find_node('ln_btn_global').get_content()
    
    @property
    def btn_local(self):
        return self._menu.root.find_node('ln_btn_local').get_content()

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

        self.btn_global.register_pressed_callback(self.select_alignment_type)
        self.btn_local.register_pressed_callback(self.select_alignment_type)

        self.set_complex_dropown(complexes, self.ln_fixed_struct)
        dd_fixed = self.ln_fixed_struct.get_content()
        dd_fixed.register_item_clicked_callback(self.handle_fixed_structure_selected)


        self.set_complex_dropown(complexes, self.ln_moving_structs, multi_select=True)
        
        dd_moving = self.ln_moving_structs.get_content()
        dd_moving.register_item_clicked_callback(self.handle_moving_structures_selected)
        self.btn_submit.register_pressed_callback(self.submit)
        self.plugin.update_menu(self._menu)

    def select_alignment_type(self, btn):
        btn_group = [self.btn_global, self.btn_local]
        btns_to_update = []
        if not btn.selected:
            btn.selected = True
            btns_to_update.append(btn)
            # Turn off other btn if it was selected
            for grp_btn in btn_group:
                if grp_btn != btn and grp_btn.selected:
                    grp_btn.selected = False
                    btns_to_update.append(grp_btn)
        if btns_to_update:
            self.plugin.update_content(*btns_to_update)

    @async_callback
    async def handle_fixed_structure_selected(self, dropdown, ddi):
        """Callback for when a complex is selected in a dropdown."""
        ln_selection = self.ln_fixed_selection
        comp = ddi.complex
        self.create_list_of_chain_dropdowns([comp], ln_selection)

    @async_callback
    async def handle_moving_structures_selected(self, dropdown, ddi):
        """Callback for when a complex is selected in a dropdown."""
        ln_selection = self.ln_moving_selections

        if not hasattr(dropdown, '_selected_items'):
            dropdown._selected_items = []

        if ddi in dropdown._selected_items:
            # Deselect item
            ddi.selected = False
            dropdown._selected_items.remove(ddi)

        if ddi.selected and ddi not in dropdown._selected_items:
            dropdown._selected_items.append(ddi)
        # Reselect selected items
        for ddi in dropdown._selected_items:
            ddi.selected = True
        self.plugin.update_content(dropdown)
        selected_complexes = [ddi.complex for ddi in dropdown._selected_items]
        self.create_list_of_chain_dropdowns(selected_complexes, ln_selection)



    @async_callback
    async def create_list_of_chain_dropdowns(self, comp_list, layoutnode):
        ui_list = ui.UIList()
        ui_list.max_displayed_items = 5
        ui_list.display_rows = 2
        list_items = []
        for comp in comp_list:
            ln_listitem = ui.LayoutNode()
            ln_listitem.layout_orientation = 1
            
            ln_label = ln_listitem.create_child_node()
            ln_label.set_size_ratio(0.25)
            lbl = ui.Label(comp.full_name)
            lbl.text_auto_size = False
            lbl.text_size = 0.3
            ln_label.set_content(lbl)
            
            ln_dd = ln_listitem.create_child_node()
            ln_dd.forward_dist = 0.004
            # ln_dd.padding = (0.03, 0.03, 0.03, 0.03)

            chain_dd = await self.create_chain_dropdown(comp)
            ln_dd.set_content(chain_dd)
            list_items.append(ln_listitem)
        ui_list.items = list_items
        layoutnode.set_content(ui_list)
        self.plugin.update_node(layoutnode)

    @async_callback
    async def update_chain_dropdown(self, ln_chain_dropdown, complex_dropdown, complex_dd_item):
        """Update chain dropdown to reflect changes in complex."""
        comp = complex_dd_item.complex
        Logs.message("Updating Chain dropdown")
        await self.display_chains(comp, ln_chain_dropdown)

    async def create_chain_dropdown(self, complex):
        """Create dropdown of chains, and add to provided layoutnode."""
        dropdown = ui.Dropdown()
        dropdown_items = []
        if sum(1 for ch in complex.chains) == 0:
            # get deep complex
            complex = (await self.plugin.request_complexes([complex.index]))[0]

        # Filter out hetatm chains which have an H appended to beginning of name
        chain_names = [
            ch.name for ch in complex.chains
            if not ch.name.startswith('H')
            or len(ch.name) == 1]
        for chain_name in chain_names:
            ddi = ui.DropdownItem(chain_name)
            dropdown_items.append(ddi)
        dropdown.max_displayed_items = len(dropdown_items)
        dropdown.items = dropdown_items
        return dropdown

    def create_structure_dropdown_items(self, complexes, multi_select=False):
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
            ddi.close_on_selected = not multi_select
            ddi.complex = struct
            complex_ddis.append(ddi)

        return complex_ddis

    @async_callback
    async def submit(self, btn):
        Logs.message("Submit button Pressed.")
        self.btn_submit.unusable = True
        self.plugin.update_content(self.btn_submit)

        reference_list = self.ln_main_list.get_content()
        comparator_list = self.ln_comparator_list.get_content()

        reference_item = next(item for item in reference_list.items if item.get_content().selected)
        comparator_items = [item for item in comparator_list.items if item.get_content().selected]

        alignment_type = 'global' if self.btn_global.selected else 'local'
        reference_comp = reference_item.get_content().complex
        comparator_comps = [item.get_content().complex for item in comparator_items]
        rmsd_result = await self.plugin.msa_superimpose(reference_comp, comparator_comps, alignment_type)

        self.btn_submit.unusable = False
        self.lbl_rmsd_value.text_value = rmsd_result
        self.plugin.update_content(self.lbl_rmsd_value, self.btn_submit)
        Logs.message("Superposition completed.")

    @property
    def prefab_complex_dropdown_btn(self):
        if not hasattr(self, '_prefab_complex_dropdown_btn'):
            ln = ui.LayoutNode()
            btn = ln.add_new_button()
            btn.selected = False
            btn.toggle_on_press = True
            # self._prefab_complex_dropdown_btn = btn
            self._prefab_complex_dropdown_btn = ln
        return self._prefab_complex_dropdown_btn

    def set_complex_dropown(self, complexes, layoutnode, multi_select=False):
        """Create dropdown of complexes, and add to provided layoutnode."""
        dropdown_items = self.create_structure_dropdown_items(complexes, multi_select=multi_select)
        dropdown = ui.Dropdown()
        dropdown.max_displayed_items = len(dropdown_items)
        dropdown.items = dropdown_items

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
    def complex_btn_selected(self, btn, single_selection=False, ui_list=None):
        btn.selected = not btn.selected
        btns_to_update = [btn]
        if btn.selected and single_selection and ui_list:
            for item in ui_list.items:
                item_btn = item.get_content()
                if item_btn != btn and item_btn.selected:
                    item_btn.selected = False
                    btns_to_update.append(item_btn)
        self.plugin.update_content(*btns_to_update)


class RMSDV2(nanome.AsyncPluginInstance):

    def start(self):
        self.menu = RMSDMenu(self)

    @async_callback
    async def on_run(self):
        self.menu.enabled = True
        complexes = await self.request_complex_list()
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

    async def msa_superimpose(self, reference_comp, comparator_comps, alignment_type='global'):
        moving_comp_indices = [comp.index for comp in comparator_comps]
        if alignment_type not in ['global', 'local']:
            raise ValueError(f'Alignment type must be either "global" or "local"')

        updated_comps = await self.request_complexes([reference_comp.index, *moving_comp_indices])
        fixed_comp = updated_comps[0]
        moving_comps = updated_comps[1:]
        fixed_comp.locked = True
        comps_to_update = [fixed_comp]
        for moving_comp in moving_comps:
            updated_moving_comp = await self.superimpose(fixed_comp, moving_comp, alignment_type)
            updated_moving_comp.locked = True
            comps_to_update.append(updated_moving_comp)
        await self.update_structures_deep(comps_to_update)

    @staticmethod
    def create_transform_matrix(superimposer):
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
            
    async def superimpose(self, fixed_comp, moving_comp, alignment_type='global'):
        Logs.message(f"Superimposing {moving_comp.full_name} onto {fixed_comp.full_name}.")
        ComplexUtils.align_to(moving_comp, fixed_comp)
        parser = PDBParser(QUIET=True)
        fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        fixed_comp.io.to_pdb(fixed_pdb.name)
        moving_comp.io.to_pdb(moving_pdb.name)
        fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
        moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)

        Logs.message("Aligning Structures.")
        mapping = self.align_sequences(fixed_struct, moving_struct, alignment_type)

        # Collect aligned residues
        # Align Residues based on Alpha Carbon
        fixed_atoms = []
        moving_atoms = []
        alpha_carbon = 'CA'
        for fixed_residue in fixed_struct.get_residues():
            fixed_id = fixed_residue.id[1]
            if fixed_id in mapping:
                fixed_atoms.append(fixed_residue[alpha_carbon])
                moving_residue_serial = mapping[fixed_id] 
                moving_residue = next(
                    rez for rez in moving_struct.get_residues()
                    if rez.id[1] == moving_residue_serial)
                moving_atoms.append(moving_residue[alpha_carbon])
        assert len(moving_atoms) == len(fixed_atoms)
        Logs.message("Superimposing Structures.")
        superimposer = Superimposer()
        superimposer.set_atoms(fixed_atoms, moving_atoms)        
        rms = superimposer.rms
        Logs.message(f"RMSD: {rms}")
        transform_matrix = self.create_transform_matrix(superimposer)
        # apply transformation to moving_comp
        moving_comp.set_surface_needs_redraw()
        for comp_atom in moving_comp.atoms:
            comp_atom.position = transform_matrix * comp_atom.position
        return moving_comp

    async def superimpose_by_chain(self, fixed_comp, fixed_chain_name, moving_comp, moving_chain_name):
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
        
        m = self.create_transform_matrix(superimposer)
        # apply transformation to moving_comp
        moving_comp.set_surface_needs_redraw()
        for comp_atom in moving_comp.atoms:
            comp_atom.position = m * comp_atom.position
        await self.update_structures_deep([moving_comp])
        self.send_notification(NotificationTypes.message, f'RMSD: {rms}')
        return rms

    def align_sequences(self, structA, structB, alignment_type='global'):
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

        if alignment_type == 'global':
            Logs.message("Using Global Alignment")
            alignment_fn = pairwise2.align.globalds
            alignment_fn_args = (
                sequence_A,
                sequence_B,
                substitution_matrices.load("BLOSUM62")
            )
            alignment_fn_kwargs = dict(
                one_alignment_only=True,
                open=-10.0,
                extend=-0.5,
                penalize_end_gaps=(False, False)
            )
        elif alignment_type == 'local':
            Logs.message("Using Local Alignment")
            alignment_fn = pairwise2.align.localxx
            alignment_fn_args = (sequence_A, sequence_B)
            alignment_fn_kwargs = dict()
        alns = alignment_fn(*alignment_fn_args, **alignment_fn_kwargs)
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
