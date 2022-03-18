from os import path

from nanome.api import ui
from nanome.util import Logs, async_callback


BASE_PATH = path.dirname(f'{path.realpath(__file__)}')
MENU_PATH = path.join(BASE_PATH, 'menu2.json')

class SelectionModeController:

    def __init__(self, plugin, menu):
        self.plugin = plugin
        self._menu = menu
        for btn in self.mode_selection_btn_group:
            # btn.toggle_on_press = True
            btn.register_pressed_callback(self.on_mode_selected)

    @property
    def mode_selection_btn_group(self):
        return [self.btn_global_align, self.btn_align_by_chain, self.btn_align_by_binding_site]

    @property
    def btn_global_align(self):
        return self._menu.root.find_node('ln_btn_global_align').get_content()

    @property
    def btn_align_by_chain(self):
        return self._menu.root.find_node('ln_btn_align_by_chain').get_content()

    @property
    def btn_align_by_binding_site(self):
        return self._menu.root.find_node('ln_btn_align_by_binding_site').get_content()

    @property
    def global_align_panel(self):
        return self._menu.root.find_node('Global Panel')
    
    @property
    def chain_align_panel(self):
        return self._menu.root.find_node('Chain Panel')

    def on_mode_selected(self, btn):
        btn.selected = True
        btns_to_update = [btn]
        for group_item in self.mode_selection_btn_group:
            if btn._content_id != group_item._content_id:
                group_item.selected = False
                btns_to_update.append(group_item)
        if btn.name == 'btn_global_align':
            self.global_align_panel.enabled = True
            self.chain_align_panel.enabled = False
        elif btn.name == 'btn_align_by_chain':
            self.chain_align_panel.enabled = True
            self.global_align_panel.enabled = False
        self.plugin.update_node(self._menu.root)
        # self.plugin.update_content(btns_to_update)



class RMSDMenu:

    def __init__(self, plugin_instance):
        super().__init__()
        self._menu = ui.Menu.io.from_json(MENU_PATH)
        self.selection_mode_controller = SelectionModeController(plugin_instance, self._menu)
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
    def ln_rmsd_value(self):
        return self._menu.root.find_node('ln_rmsd_value')

    @async_callback
    async def render(self, complexes=None, default_values=False):
        complexes = complexes or []
        self.complexes = complexes
        self.btn_global.register_pressed_callback(self.select_alignment_type)
        self.btn_local.register_pressed_callback(self.select_alignment_type)

        self.set_complex_dropdown(complexes, self.ln_fixed_struct)
        dd_fixed = self.ln_fixed_struct.get_content()
        dd_fixed.register_item_clicked_callback(self.handle_fixed_structure_selected)

        self.set_complex_dropdown(complexes, self.ln_moving_structs, multi_select=True)
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
            ln_dd.name = f'{comp.full_name} chain'
            ln_dd.forward_dist = 0.004

            chain_dd = await self.create_chain_dropdown(comp)
            ln_dd.set_content(chain_dd)
            list_items.append(ln_listitem)
        ui_list.items = list_items
        layoutnode.set_content(ui_list)
        self.plugin.update_node(layoutnode)

    async def create_chain_dropdown(self, complex):
        """Create dropdown of chains, and add to provided layoutnode."""
        dropdown = ui.Dropdown()
        dropdown.complex = complex
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

        fixed_comp = next(ddi.complex for ddi in self.ln_fixed_struct.get_content().items if ddi.selected)
        moving_comps = [ddi.complex for ddi in self.ln_moving_structs.get_content().items if ddi.selected]

        rmsd_results = await self.plugin.msa_superimpose(fixed_comp, moving_comps)
        self.btn_submit.unusable = False

        results_list = self.render_rmsd_results(rmsd_results)
        self.ln_rmsd_value.set_content(results_list)
        self.plugin.update_node(self.ln_rmsd_value)
        self.plugin.update_content(self.btn_submit)
        Logs.message("Superposition completed.")

    def render_rmsd_results(self, rmsd_results):
        """Render rmsd results in a list of labels."""
        results_list = ui.UIList()
        results_list.max_displayed_items = 10
        results_list.display_rows = 10
        results_list.items = []

        for comp_name, rms in rmsd_results.items():
            ln_rmsd_result = ui.LayoutNode()
            ln_rmsd_result.layout_orientation = 1
            ln_rmsd_result.set_size_ratio(0.25)
            ln_rmsd_result.set_content(ui.Label(f'{comp_name}: {rms:.2f}'))
            results_list.items.append(ln_rmsd_result)

        return results_list

    def set_complex_dropdown(self, complexes, layoutnode, multi_select=False):
        """Create dropdown of complexes, and add to provided layoutnode."""
        dropdown_items = self.create_structure_dropdown_items(complexes, multi_select=multi_select)
        dropdown = ui.Dropdown()
        dropdown.max_displayed_items = len(dropdown_items)
        dropdown.items = dropdown_items

        layoutnode.set_content(dropdown)
        self.plugin.update_node(layoutnode)
