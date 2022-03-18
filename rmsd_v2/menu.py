from os import path
from functools import partial
from nanome.api import ui
from nanome.util import Logs, async_callback


BASE_PATH = path.dirname(f'{path.realpath(__file__)}')
MENU_PATH = path.join(BASE_PATH, 'menu.json')


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
    def global_align_panel(self):
        return self._menu.root.find_node('Global Panel')

    @property
    def chain_align_panel(self):
        return self._menu.root.find_node('Chain Panel')

    def on_mode_selected(self, btn, update=True):
        btn.selected = True
        btns_to_update = [btn]
        for group_item in self.mode_selection_btn_group:
            if btn._content_id != group_item._content_id:
                group_item.selected = False
                btns_to_update.append(group_item)
        if btn.name == 'btn_global_align':
            self.current_mode = 'global'
            self.global_align_panel.enabled = True
            self.chain_align_panel.enabled = False
        elif btn.name == 'btn_align_by_chain':
            self.current_mode = 'chain'
            self.chain_align_panel.enabled = True
            self.global_align_panel.enabled = False
        if update:
            self.plugin.update_menu(self._menu)


class GlobalAlignController:

    def __init__(self, plugin, menu):
        self.plugin = plugin
        self._menu = menu

    def render(self, complexes=None):
        complexes = complexes or []
        self.set_complex_dropdown(complexes, self.ln_fixed_struct)
        self.set_complex_dropdown(complexes, self.ln_moving_structs, multi_select=True)

    @property
    def root(self):
        return self._menu.root.find_node('Global Panel')

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
        return next(ddi.complex for ddi in self.ln_fixed_struct.get_content().items if ddi.selected)

    def get_moving_complexes(self):
        return [ddi.complex for ddi in self.ln_moving_structs.get_content().items if ddi.selected]

    def set_complex_dropdown(self, complexes, layoutnode, multi_select=False):
        """Create dropdown of complexes, and add to provided layoutnode."""
        dropdown_items = self.create_structure_dropdown_items(complexes, multi_select=multi_select)
        dropdown = ui.Dropdown()
        dropdown.max_displayed_items = len(dropdown_items)
        dropdown.items = dropdown_items

        layoutnode.set_content(dropdown)
        self.plugin.update_node(layoutnode)

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


class ChainAlignController:

    def __init__(self, plugin, menu):
        self.plugin = plugin
        self._menu = menu

    @property
    def root(self):
        return self._menu.root.find_node('Chain Panel')

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
        fixed_dropdown = self.ln_fixed_struct.get_content()
        moving_dropdown = self.ln_moving_struct.get_content()

        dropdown_items = self.create_structure_dropdown_items(complexes)
        fixed_dropdown.items = dropdown_items
        fixed_dropdown.max_displayed_items = len(dropdown_items)

        dropdown_items = self.create_structure_dropdown_items(complexes)
        moving_dropdown.items = dropdown_items
        moving_dropdown.max_displayed_items = len(dropdown_items)

        fixed_dropdown.register_item_clicked_callback(
            partial(self.update_chain_dropdown, self.ln_fixed_chain))
        moving_dropdown.register_item_clicked_callback(
            partial(self.update_chain_dropdown, self.ln_moving_chain))

        # self.ln_fixed_struct.set_content(fixed_dropdown)
        # self.ln_moving_struct.set_content(moving_dropdown)
        # self.plugin.update_node(self.ln_fixed_struct, self.ln_moving_struct)

    def get_fixed_complex(self):
        return next(ddi.complex for ddi in self.ln_fixed_struct.get_content().items if ddi.selected)

    def get_moving_complex(self):
        return next(ddi.complex for ddi in self.ln_moving_struct.get_content().items if ddi.selected)

    def get_fixed_chain(self):
        return next(ddi.name for ddi in self.ln_fixed_chain.get_content().items if ddi.selected)

    def get_moving_chain(self):
        return next(ddi.name for ddi in self.ln_moving_chain.get_content().items if ddi.selected)

    @async_callback
    async def update_chain_dropdown(self, ln_chain_dropdown, complex_dropdown, complex_dd_item):
        """Update chain dropdown to reflect changes in complex."""
        comp = complex_dd_item.complex
        Logs.message("Updating Chain dropdown")
        dropdown = ln_chain_dropdown.get_content()
        dropdown_items = []
        if sum(1 for ch in comp.chains) == 0:
            # get deep complex
            comp = (await self.plugin.request_complexes([comp.index]))[0]
        # Filter out hetatom chains (HA, HB, etc)
        chain_names = [
            ch.name for ch in comp.chains
            if not ch.name.startswith('H') or len(ch.name) < 2
        ]
        for chain_name in chain_names:
            ddi = ui.DropdownItem(chain_name)
            dropdown_items.append(ddi)
        dropdown.max_displayed_items = len(dropdown_items)
        dropdown.items = dropdown_items
        self.plugin.update_node(ln_chain_dropdown)

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


class RMSDMenu:

    def __init__(self, plugin_instance):
        super().__init__()
        self._menu = ui.Menu.io.from_json(MENU_PATH)
        self.plugin = plugin_instance
        self.selection_mode_controller = SelectionModeController(plugin_instance, self._menu)
        self.global_align_controller = GlobalAlignController(plugin_instance, self._menu)
        self.chain_align_controller = ChainAlignController(plugin_instance, self._menu)
        # Make sure Global Align Panel is always open
        default_mode = self.selection_mode_controller.btn_global_align
        self.selection_mode_controller.on_mode_selected(default_mode, update=False)

    @property
    def btn_submit(self):
        return self._menu.root.find_node('ln_submit').get_content()

    @property
    def ln_rmsd_value(self):
        return self._menu.root.find_node('ln_rmsd_value')

    @async_callback
    async def render(self, complexes=None):
        complexes = complexes or []
        self.complexes = complexes
        self.global_align_controller.render(complexes)
        self.chain_align_controller.render(complexes)
        self.btn_submit.register_pressed_callback(self.submit)
        self.plugin.update_menu(self._menu)

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

    @async_callback
    async def submit(self, btn):
        Logs.message("Submit button Pressed.")
        self.btn_submit.unusable = True
        self.plugin.update_content(self.btn_submit)
        current_mode = self.selection_mode_controller.current_mode
        if current_mode == 'global':
            fixed_comp = self.global_align_controller.get_fixed_complex()
            moving_comps = self.global_align_controller.get_moving_complexes()
            rmsd_results = await self.plugin.msa_superimpose(fixed_comp, moving_comps)
        if current_mode == 'chain':
            fixed_comp = self.chain_align_controller.get_fixed_complex()
            fixed_chain = self.chain_align_controller.get_fixed_chain()
            moving_comp = self.chain_align_controller.get_moving_complex()
            moving_chain = self.chain_align_controller.get_moving_chain()
            rmsd_results = await self.plugin.superimpose_by_chain(fixed_comp, fixed_chain, moving_comp, moving_chain)
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
