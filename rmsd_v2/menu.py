import asyncio
from os import path
from nanome.api import ui
from nanome.util import Logs, async_callback
from nanome.util.enums import NotificationTypes


BASE_PATH = path.dirname(f'{path.realpath(__file__)}')
MENU_PATH = path.join(BASE_PATH, 'menu_json', 'newMenu.json')
MENU_ITEM_PATH = path.join(BASE_PATH, 'menu_json', 'menu_item.json')
INFO_ICON_PATH = path.join(BASE_PATH, 'assets', 'info_icon.png')


class SelectionModeController:

    def __init__(self, plugin, menu):
        self.plugin = plugin
        self._menu = menu
        self.current_mode = 'global'
        for btn in self.mode_selection_btn_group:
            btn.register_pressed_callback(self.on_mode_selected)

    @property
    def mode_selection_btn_group(self):
        return [self.btn_global_align, self.btn_align_by_chain, self.btn_align_by_active_site]

    @property
    def btn_global_align(self):
        return self._menu.root.find_node('ln_btn_global_align').get_content()

    @property
    def btn_align_by_chain(self):
        return self._menu.root.find_node('ln_btn_align_by_chain').get_content()

    @property
    def btn_align_by_active_site(self):
        return self._menu.root.find_node('ln_btn_align_by_active_site').get_content()

    @property
    def entry_align_panel(self):
        return self._menu.root.find_node('Entry Panel')

    @property
    def chain_align_panel(self):
        return self._menu.root.find_node('Chain Panel')

    # @async_callback
    def on_mode_selected(self, btn, update=True, log=True):
        btn.selected = True
        btns_to_update = [btn]
        for group_item in self.mode_selection_btn_group:
            if btn._content_id != group_item._content_id:
                group_item.selected = False
                btns_to_update.append(group_item)
        if btn.name == 'btn_global_align':
            if log:
                Logs.message("Switched to entry mode")
            self.current_mode = 'global'
            self.entry_align_panel.enabled = True
            self.chain_align_panel.enabled = False
            self.active_site_panel.enabled = False
        elif btn.name == 'btn_align_by_chain':
            if log:
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
                    loop = asyncio.get_event_loop()
                    self.plugin.complexes = loop.run_until_complete(
                        self.plugin.request_complexes(comp_indices))
                    # This is kinda iffy, but it works
                    loop.run_until_complete(self.plugin.menu.render(complexes=self.plugin.complexes))
                    self.btn_align_by_chain.unusable = False
                    break
        if update:
            self.plugin.update_menu(self._menu)


class EntryAlignController:

    def __init__(self, plugin, menu):
        self.plugin = plugin
        self._menu = menu

    def render(self, complexes=None):
        complexes = complexes or []
        default_moving = None
        if len(complexes) == 2:
            default_moving = complexes[1]
        # self.set_complex_dropdown(complexes, self.ln_fixed_struct, default_comp=default_fixed)
        self.populate_moving_comp_list(complexes, default_comp=default_moving)

    @property
    def root(self):
        return self._menu.root.find_node('Entry Panel')

    @property
    def ln_moving_comp_list(self):
        return self.root.find_node('ln_moving_comp_list')

    @property
    def ln_fixed_struct(self):
        return self.panel_root.find_node('ln_fixed_struct')

    @property
    def ln_moving_structs(self):
        return self.panel_root.find_node('ln_moving_structs')

    @property
    def ln_fixed_selection(self):
        return self.panel_root.find_node('ln_fixed_selection')

    @property
    def ln_moving_selections(self):
        return self.panel_root.find_node('ln_moving_selection')

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
        comp_list = self.ln_moving_comp_list.get_content()
        for comp in complexes:
            ln = ui.LayoutNode.io.from_json(MENU_ITEM_PATH)
            btn_fixed = ln.find_node('btn_fixed').get_content()
            btn_fixed.register_pressed_callback(self.btn_fixed_clicked)
            btn_moving = ln.find_node('btn_moving').get_content()
            struct_name = ln.find_node('lbl_struct_name').get_content()
            struct_name.text_value = comp.full_name
            btn_fixed.toggle_on_press = True
            btn_moving.toggle_on_press = True
            comp_list.items.append(ln)
            continue
        self.plugin.update_node(self.ln_moving_comp_list)

    def btn_fixed_clicked(self, btn):
        """Only one fixed strcuture can be selected at a time."""
        btns_to_update = [btn]
        for menu_item in self.ln_moving_comp_list.get_content().items:
            btn_fixed = menu_item.find_node('btn_fixed').get_content()
            btn_moving = menu_item.find_node('btn_moving').get_content()
            if btn_fixed == btn:
                # Fixed structure cannot also be a moving structure.
                if not btn.selected:
                    btn_moving.unusable = False
                else:
                    btn_moving.selected = False
                    btn_moving.unusable = True
            else:
                btn_fixed.selected = False
                btn_moving.unusable = False
            btns_to_update.append(btn_fixed)
            btns_to_update.append(btn_moving)
        self.plugin.update_content(*btns_to_update)


class RMSDMenu:

    def __init__(self, plugin_instance):
        super().__init__()
        self._menu = ui.Menu.io.from_json(MENU_PATH)
        self.plugin = plugin_instance
        self.selection_mode_controller = SelectionModeController(plugin_instance, self._menu)
        self.global_align_controller = EntryAlignController(plugin_instance, self._menu)
        # Make sure Global Align Panel is always default
        default_mode = self.selection_mode_controller.btn_global_align
        self.selection_mode_controller.on_mode_selected(default_mode, update=False, log=False)
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
                rmsd_results = await self.plugin.superimpose_by_entry(fixed_comp, moving_comps)
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
