import asyncio
from os import path
from nanome.api import ui
from nanome.util import Logs, async_callback, Color
from nanome.util.enums import NotificationTypes, SizingTypes, PaddingTypes


BASE_PATH = path.dirname(f'{path.realpath(__file__)}')
MENU_PATH = path.join(BASE_PATH, 'menu.json')
INFO_ICON_PATH = path.join(BASE_PATH, 'info_icon.png')


class SelectionModeController:

    def __init__(self, plugin, menu):
        self.plugin = plugin
        self._menu = menu
        self.current_mode = 'global'
        for btn in self.mode_selection_btn_group:
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
    def entry_align_panel(self):
        return self._menu.root.find_node('Entry Panel')

    @property
    def chain_align_panel(self):
        return self._menu.root.find_node('Chain Panel')

    def set_default_state(self):
        self.btn_global_align.selected = True
        self.btn_align_by_chain.selected = False
        self.btn_align_by_binding_site.selected = False
        self.entry_align_panel.enabled = True
        self.chain_align_panel.enabled = False

    @async_callback
    async def on_mode_selected(self, btn):
        btn.selected = True
        btns_to_update = [btn]
        for group_item in self.mode_selection_btn_group:
            if btn._content_id != group_item._content_id:
                group_item.selected = False
                btns_to_update.append(group_item)
        if btn.name == 'btn_global_align':
            Logs.message("Switched to entry mode")
            self.current_mode = 'global'
            self.entry_align_panel.enabled = True
            self.chain_align_panel.enabled = False
        elif btn.name == 'btn_align_by_chain':
            Logs.message("Switched to chain mode.")
            self.current_mode = 'chain'
            self.chain_align_panel.enabled = True
            self.entry_align_panel.enabled = False
            # Get deep complexes if necessary
            for comp in self.plugin.complexes:
                if sum(1 for _ in comp.chains) == 0:
                    self.btn_align_by_chain.unusable = True
                    self.plugin.update_content(self.btn_align_by_chain)
                    comp_indices = [cmp.index for cmp in self.plugin.complexes]
                    # Use event loop because on_mode_selected is called in __init__
                    self.plugin.complexes = await self.plugin.request_complexes(comp_indices)
                    # This is kinda iffy, but it works
                    await self.plugin.menu.render(complexes=self.plugin.complexes)
                    self.btn_align_by_chain.unusable = False
                    break
        self.plugin.update_menu(self._menu)


class EntryAlignController:

    def __init__(self, plugin, menu):
        self.plugin = plugin
        self._menu = menu

    def render(self, complexes=None):
        complexes = complexes or []
        default_fixed = None
        default_moving = None
        if len(complexes) == 2:
            default_fixed = complexes[0]
            default_moving = complexes[1]
        self.set_complex_dropdown(complexes, self.ln_fixed_struct, default_comp=default_fixed)
        self.populate_moving_comp_list(complexes, default_comp=default_moving)

    @property
    def root(self):
        return self._menu.root.find_node('Entry Panel')

    @property
    def ln_moving_comp_list(self):
        return self.root.find_node('ln_moving_comp_list')

    @property
    def ln_fixed_struct(self):
        return self.root.find_node('ln_fixed_struct')

    @property
    def ln_moving_structs(self):
        return self.root.find_node('ln_moving_structs')

    @property
    def ln_fixed_selection(self):
        return self.root.find_node('ln_fixed_selection')

    @property
    def ln_moving_selections(self):
        return self.root.find_node('ln_moving_selection')

    def get_fixed_complex(self):
        return next((ddi.complex for ddi in self.ln_fixed_struct.get_content().items if ddi.selected), None)

    def get_moving_complexes(self):
        comps = []
        for item in self._menu.root.find_node('ln_moving_comp_list').get_content().items:
            btn = item.get_children()[0].get_content()
            if btn.selected:
                comps.append(btn.comp)
        return comps

    def set_complex_dropdown(self, complexes, layoutnode, multi_select=False, default_comp=None):
        """Create dropdown of complexes, and add to provided layoutnode."""
        dropdown_items = self.create_structure_dropdown_items(complexes, multi_select=multi_select, default_comp=default_comp)
        dropdown = ui.Dropdown()
        dropdown.max_displayed_items = len(dropdown_items)
        dropdown.items = dropdown_items
        layoutnode.set_content(dropdown)
        self.plugin.update_node(layoutnode)
        if multi_select:
            dropdown.register_item_clicked_callback(self.multi_select_item_selected)

    def create_structure_dropdown_items(self, complexes, multi_select=False, default_comp=None):
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
            if ddi.complex == default_comp:
                ddi.selected = True
            complex_ddis.append(ddi)

        return complex_ddis

    @async_callback
    async def multi_select_item_selected(self, dropdown, ddi):
        """Callback for when a complex is selected in a dropdown."""

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

    def populate_moving_comp_list(self, complexes, default_comp=None):
        green = Color(36, 184, 177)
        comp_list = self.ln_moving_comp_list.get_content()
        comp_list.items = []
        for comp in complexes:
            ln = ui.LayoutNode()
            ln.padding_type = PaddingTypes.fixed

            ln.forward_dist = 0.05
            ln.layout_orientation = 1
            btn_ln = ln.create_child_node(ui.LayoutNode())
            lbl_ln = ln.create_child_node(ui.LayoutNode())
            lbl_ln.padding = (0.02, 0.00, 0.0, 0.0)

            btn_ln.sizing_type = SizingTypes.ratio
            btn_ln.sizing_value = 0.1
            btn = btn_ln.add_new_button()
            btn.text.value.set_all("")
            btn.mesh.active = True
            btn.mesh.enabled.set_all(False)
            btn.mesh.enabled.selected = True
            btn.mesh.enabled.hover = True
            btn.mesh.enabled.highlighted = True
            btn.mesh.enabled.selected_highlighted = True
            btn.mesh.color.selected = green
            btn.mesh.color.highlighted = green
            btn.mesh.color.selected_highlighted = green

            btn.comp = comp
            if comp == default_comp:
                btn.selected = True

            btn.toggle_on_press = True
            lbl_ln.add_new_label(comp.full_name)
            comp_list.items.append(ln)

        self.plugin.update_node(self.ln_moving_comp_list)


class ChainAlignController:

    def __init__(self, plugin, menu):
        self.plugin = plugin
        self._menu = menu

    @property
    def root(self):
        return self._menu.root.find_node('Chain Panel')

    @property
    def ln_moving_comp_list(self):
        return self.root.find_node('ln_moving_comp_list')

    @property
    def ln_fixed_struct(self):
        return self.root.find_node('ln_fixed_struct')

    @property
    def ln_moving_struct(self):
        return self.root.find_node('ln_moving_struct')

    @property
    def ln_fixed_chain(self):
        return self.root.find_node('ln_fixed_chain')

    @property
    def ln_moving_chain(self):
        return self.root.find_node('ln_moving_chain')

    @async_callback
    async def render(self, complexes=None):
        complexes = complexes or []

        default_fixed = None
        default_moving = None
        if len(complexes) == 2:
            default_fixed = complexes[0]
            default_moving = complexes[1]

        fixed_dropdown = self.ln_fixed_struct.get_content()
        dropdown_items = self.create_structure_dropdown_items(complexes, default_comp=default_fixed)
        fixed_dropdown.items = dropdown_items
        fixed_dropdown.max_displayed_items = len(dropdown_items)

        if default_fixed:
            chain_dropdown = self.ln_fixed_chain.get_content()
            chain_dropdown.items = self.create_chain_dropdown_items(default_fixed)

        if len(self.plugin.complexes) == 2 and chain_dropdown.items:
            chain_dropdown.items[0].selected = True

        fixed_dropdown.register_item_clicked_callback(self.update_fixed_chain_dropdown)
        self.populate_moving_comp_list(complexes, default_comp=default_moving)

    def get_fixed_complex(self):
        return next((ddi.complex for ddi in self.ln_fixed_struct.get_content().items if ddi.selected), None)

    def get_moving_complexes_and_chains(self):
        comp_chain_list = []
        lst = self.ln_moving_comp_list.get_content()
        for item in lst.items:
            btn = item.get_children()[0].get_content()
            chain_dd = item.get_children()[1].get_children()[1].get_content()
            if not btn.selected:
                continue

            selected_comp = btn.comp
            selected_chain = None
            for chain_ddi in chain_dd.items:
                if chain_ddi.selected:
                    selected_chain = chain_ddi.name
                    break
                print('here')
            comp_chain_list.append((selected_comp, selected_chain))
        return comp_chain_list

    def get_fixed_chain(self):
        return next((ddi.name for ddi in self.ln_fixed_chain.get_content().items if ddi.selected), None)

    def get_moving_chains(self):
        return next((ddi.name for ddi in self.ln_moving_chain.get_content().items if ddi.selected), None)

    def update_fixed_chain_dropdown(self, dropdown, ddi):
        comp = ddi.complex
        chain_dropdown = self.ln_fixed_chain.get_content()
        chain_dropdown.items = self.create_chain_dropdown_items(comp)
        self.plugin.update_content(chain_dropdown)

    @staticmethod
    def create_chain_dropdown_items(comp, set_default=False):
        """Update chain dropdown to reflect changes in complex."""
        dropdown_items = []
        # Filter out hetatom chains (HA, HB, etc)
        chain_names = [
            ch.name for ch in comp.chains
            if not ch.name.startswith('H') or len(ch.name) < 2
        ]
        for chain_name in chain_names:
            ddi = ui.DropdownItem(chain_name)
            dropdown_items.append(ddi)
        if set_default and dropdown_items:
            dropdown_items[0].selected = True
        return dropdown_items

    def create_structure_dropdown_items(self, complexes, default_comp=None):
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
            if ddi.complex == default_comp:
                ddi.selected = True
            complex_ddis.append(ddi)

        return complex_ddis

    def populate_moving_comp_list(self, complexes, default_comp=None):
        green = Color(36, 184, 177)
        comp_list = self.ln_moving_comp_list.get_content()
        comp_list.items = []
        for comp in complexes:
            ln = ui.LayoutNode()
            ln.forward_dist = 0.05
            ln.layout_orientation = 1
            btn_ln = ln.create_child_node(ui.LayoutNode())
            info_ln = ln.create_child_node(ui.LayoutNode())
            info_ln.padding = (0.02, 0.00, 0.0, 0.0)

            btn_ln.sizing_type = SizingTypes.ratio
            btn_ln.sizing_value = 0.1
            btn = btn_ln.add_new_button()
            btn.text.value.set_all("")
            btn.mesh.active = True
            btn.mesh.enabled.set_all(False)
            btn.mesh.enabled.selected = True
            btn.mesh.enabled.hover = True
            btn.mesh.enabled.highlighted = True
            btn.mesh.enabled.selected_highlighted = True
            btn.mesh.color.selected = green
            btn.mesh.color.highlighted = green
            btn.mesh.color.selected_highlighted = green
            btn.comp = comp
            if btn.comp == default_comp:
                btn.selected = True
            btn.toggle_on_press = True

            info_ln.layout_orientation = 1
            lbl_ln = info_ln.create_child_node(ui.LayoutNode())
            lbl_ln.add_new_label(comp.full_name)
            btn_list_ln = info_ln.create_child_node(ui.LayoutNode())

            comp_dd = ui.Dropdown()
            set_default = bool(default_comp)
            comp_dd.items = self.create_chain_dropdown_items(comp, set_default=set_default)
            comp_dd = btn_list_ln.set_content(comp_dd)
            comp_list.items.append(ln)

        self.plugin.update_node(self.ln_moving_comp_list)


class RMSDMenu:

    def __init__(self, plugin_instance):
        super().__init__()
        self._menu = ui.Menu.io.from_json(MENU_PATH)
        self.plugin = plugin_instance
        self.selection_mode_controller = SelectionModeController(plugin_instance, self._menu)
        self.global_align_controller = EntryAlignController(plugin_instance, self._menu)
        self.chain_align_controller = ChainAlignController(plugin_instance, self._menu)
        self.selection_mode_controller.set_default_state()
        self.btn_color_override.toggle_on_press = True
        self.btn_color_override.switch.active = True
        self.ln_info_img.add_new_image(INFO_ICON_PATH)

    @property
    def btn_submit(self):
        return self._menu.root.find_node('ln_submit').get_content()

    @property
    def ln_rmsd_value(self):
        return self._menu.root.find_node('ln_rmsd_value')

    @property
    def lst_rmsd_results(self):
        return self._menu.root.find_node('ln_rmsd_results').get_content()

    @property
    def btn_color_override(self):
        return self._menu.root.find_node('ln_color_override').get_content()

    @property
    def ln_info_img(self):
        return self._menu.root.find_node('ln_info_img')

    @property
    def ln_moving_comp_list(self):
        return self._menu.root.find_node('ln_moving_comp_list')

    @async_callback
    async def render(self, complexes=None):
        complexes = complexes or []
        self.global_align_controller.render(complexes)
        await self.chain_align_controller.render(complexes)
        self.btn_submit.register_pressed_callback(self.submit)
        self.plugin.update_menu(self._menu)

    @async_callback
    async def submit(self, btn):
        Logs.message("Submit button Pressed.")
        self.btn_submit.unusable = True
        self.plugin.update_content(self.btn_submit)
        current_mode = self.selection_mode_controller.current_mode
        rmsd_results = None
        if current_mode == 'global':
            fixed_comp = self.global_align_controller.get_fixed_complex()
            moving_comps = self.global_align_controller.get_moving_complexes()
            if not all([fixed_comp, moving_comps]):
                msg = "Please select all complexes."
                Logs.warning(msg)
                self.plugin.send_notification(NotificationTypes.warning, msg)
            else:
                rmsd_results = await self.plugin.msa_superimpose(fixed_comp, moving_comps)
        if current_mode == 'chain':
            fixed_comp = self.chain_align_controller.get_fixed_complex()
            fixed_chain = self.chain_align_controller.get_fixed_chain()
            moving_comp_chain_list = self.chain_align_controller.get_moving_complexes_and_chains()
            # moving_chain = self.chain_align_controller.get_moving_chains()
            if not all([fixed_comp, fixed_chain, moving_comp_chain_list]):
                msg = "Please select all complexes and chains."
                Logs.warning(msg)
                self.plugin.send_notification(NotificationTypes.error, msg)
            else:
                rmsd_results = await self.plugin.superimpose_by_chain(fixed_comp, fixed_chain, moving_comp_chain_list)
        if rmsd_results:
            self.render_rmsd_results(rmsd_results)
        self.btn_submit.unusable = False
        self.plugin.update_content(self.btn_submit)

    def render_rmsd_results(self, rmsd_results):
        """Render rmsd results in a list of labels."""
        new_menu = ui.Menu()
        new_menu.index = 100
        new_menu.title = "RMSD Values"
        ln = ui.LayoutNode()
        results_list = ui.UIList()
        for name, rms_val in rmsd_results.items():
            item = ui.LayoutNode()
            item.add_new_label(f"{name}: {rms_val:.2f}")
            results_list.items.append(item)
        ln.set_content(results_list)
        new_menu.root.add_child(ln)
        new_menu.enabled = True
        self.plugin.update_menu(new_menu)
