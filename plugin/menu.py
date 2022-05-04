from os import path
from nanome.api import ui
from nanome.util import Logs, async_callback, Color
from nanome.util.enums import NotificationTypes


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


BASE_PATH = path.dirname(f'{path.realpath(__file__)}')
MENU_PATH = path.join(BASE_PATH, 'menu_json', 'menu.json')
COMP_LIST_ITEM_PATH = path.join(BASE_PATH, 'menu_json', 'comp_list_item.json')
RMSD_MENU_PATH = path.join(BASE_PATH, 'menu_json', 'rmsd_menu.json')
RMSD_TABLE_ENTRY = path.join(BASE_PATH, 'menu_json', 'rmsd_list_entry.json')

GEAR_ICON_PATH = path.join(BASE_PATH, 'assets', 'gear.png')


class MainMenu:

    def __init__(self, plugin_instance):
        super().__init__()
        self._menu = ui.Menu.io.from_json(MENU_PATH)
        self.plugin = plugin_instance
        self.btn_advanced.icon.value.set_all(GEAR_ICON_PATH)

        self.current_mode = 'global'
        for btn in self.mode_selection_btn_group:
            btn.register_pressed_callback(self.on_mode_selected)
        self.btn_rmsd_table.register_pressed_callback(self.open_rmsd_table)
        self.btn_docs.register_pressed_callback(self.open_docs_page)

    @property
    def btn_submit(self):
        return self._menu.root.find_node('ln_submit').get_content()

    @property
    def ln_btn_rmsd_table(self):
        return self._menu.root.find_node('btn_rmsd_table')

    @property
    def btn_rmsd_table(self):
        return self.ln_btn_rmsd_table.get_content()

    @property
    def ln_moving_comp_list(self):
        return self._menu.root.find_node('ln_moving_comp_list')
    
    @property
    def lbl_moving_structures(self):
        return self._menu.root.find_node('lbl_moving_structures').get_content()

    @async_callback
    async def render(self, complexes=None):
        complexes = complexes or []
        current_mode = self.current_mode
        self.populate_comp_list(complexes, current_mode)
        self.btn_submit.register_pressed_callback(self.submit)
        self.plugin.update_menu(self._menu)

    @async_callback
    async def submit(self, btn):
        Logs.message("Submit button Pressed.")
        self.btn_submit.unusable = True
        self.plugin.update_content(self.btn_submit)
        current_mode = self.current_mode
        rmsd_results = None
        fixed_comp_index = self.get_fixed_comp_index() or 0
        if current_mode == 'global':
            moving_comp_indices = self.get_moving_comp_indices()
            if not all([fixed_comp_index, moving_comp_indices]):
                msg = "Please select all complexes."
                Logs.warning(msg)
                self.plugin.send_notification(NotificationTypes.error, msg)
            else:
                rmsd_results = await self.plugin.superimpose_by_entry(fixed_comp_index, moving_comp_indices)
        if current_mode == 'chain':
            fixed_chain = self.get_fixed_chain()
            moving_comp_chain_list = self.get_moving_comp_indices_and_chains()
            if not all([fixed_comp_index, fixed_chain, moving_comp_chain_list]):
                msg = "Please select all complexes and chains."
                Logs.warning(msg)
                self.plugin.send_notification(NotificationTypes.error, msg)
            else:
                rmsd_results = await self.plugin.superimpose_by_chain(fixed_comp_index, fixed_chain, moving_comp_chain_list)
        if rmsd_results:
            fixed_name = next(comp.full_name for comp in self.plugin.complexes if comp.index == fixed_comp_index)
            if current_mode == 'chain':
                fixed_name = f'{fixed_name} Chain {fixed_chain}'
            self.render_rmsd_results(rmsd_results, fixed_name)
            self.ln_btn_rmsd_table.enabled = True
        self.btn_submit.unusable = False
        self.plugin.update_node(self.ln_btn_rmsd_table)
        self.plugin.update_content(self.btn_submit)

    def render_rmsd_results(self, rmsd_results, fixed_comp_name):
        """Render rmsd results in a list."""
        new_menu = ui.Menu.io.from_json(RMSD_MENU_PATH)
        new_menu.index = 200

        comp_header_lbl = new_menu.root.find_node('comp_name_header').get_content()
        comp_header_lbl.text_value = comp_header_lbl.text_value.replace('<fixed>', fixed_comp_name)
        results_list = new_menu.root.find_node('results_list').get_content()
        list_items = []
        row_color1 = Color(21, 26, 37)
        row_color2 = Color(42, 52, 63)
        for i, comp_name in enumerate(rmsd_results, 1):
            results_data = rmsd_results[comp_name]
            rmsd_val = results_data['rmsd']
            paired_residue_count = results_data['paired_residues']
            if 'chain' in results_data:
                comp_name = f'{comp_name} Chain {results_data["chain"]}'

            item = ui.LayoutNode().io.from_json(RMSD_TABLE_ENTRY)
            item_mesh = item.add_new_mesh()
            item_mesh.mesh_color = row_color1 if i % 2 == 0 else row_color2

            item.get_children()[0].get_content().text_value = i
            item.get_children()[1].get_content().text_value = comp_name
            item.get_children()[2].get_content().text_value = rmsd_val
            item.get_children()[3].get_content().text_value = paired_residue_count
            list_items.append(item)

        results_list.items = list_items
        new_menu.enabled = False
        self.rmsd_menu = new_menu

    @property
    def root(self):
        return self._menu.root.find_node('Entry Panel')

    def get_fixed_comp_index(self):
        for item in self._menu.root.find_node('ln_moving_comp_list').get_content().items:
            btn_fixed = item.find_node('btn_fixed').get_content()
            struct_name = item.find_node('lbl_struct_name').get_content().text_value
            if btn_fixed.selected:
                comp = next(comp for comp in self.plugin.complexes if comp.full_name == struct_name)
                return comp.index

    def get_moving_comp_indices(self):
        comps = []
        for item in self._menu.root.find_node('ln_moving_comp_list').get_content().items:
            btn_moving = item.find_node('btn_moving').get_content()
            struct_name = item.find_node('lbl_struct_name').get_content().text_value
            if btn_moving.selected:
                comp = next(comp for comp in self.plugin.complexes if comp.full_name == struct_name)
                comps.append(comp.index)
        return comps

    def get_fixed_chain(self):
        for item in self._menu.root.find_node('ln_moving_comp_list').get_content().items:
            btn_fixed = item.find_node('btn_fixed').get_content()
            if btn_fixed.selected:
                dd_chain = item.find_node('dd_chain').get_content()
                selected_ddi = next((ddi for ddi in dd_chain.items if ddi.selected), None)
                return getattr(selected_ddi, 'name', '')

    def populate_comp_list(self, complexes, mode='global', default_comp=None):
        comp_list = self.ln_moving_comp_list.get_content()
        comp_list.items = []
        set_default_values = len(complexes) == 2
        for i, comp in enumerate(complexes):
            ln = ui.LayoutNode.io.from_json(COMP_LIST_ITEM_PATH)
            btn_fixed = ln.find_node('btn_fixed').get_content()
            btn_moving = ln.find_node('btn_moving').get_content()
            lbl_struct_name = ln.find_node('lbl_struct_name').get_content()
            ln_lbl_chain_count = ln.find_node('lbl_chain_count')
            
            lbl_chain_count = ln_lbl_chain_count.get_content()
            btn_fixed.register_pressed_callback(self.btn_fixed_clicked)
            btn_moving.register_pressed_callback(self.btn_moving_clicked)

            lbl_struct_name.text_value = comp.full_name
            btn_fixed.toggle_on_press = True
            btn_moving.toggle_on_press = True
            chain_count = sum(
                1 for ch in comp.chains
                if not ch.name.startswith('H')
                or len(ch.name) < 2)
            lbl_chain_count.text_value = f'{chain_count} Chain(s)'

            ln_dd_chain = ln.find_node('dd_chain')
            ln_dd_chain.enabled = mode == 'chain'
            # Set up chain dropdown if in chain mode
            if mode == 'chain':
                ln_lbl_chain_count.enabled = True
                dd_chain = ln_dd_chain.get_content()
                comp = next(
                    cmp for cmp in self.plugin.complexes
                    if cmp.full_name == lbl_struct_name.text_value)
                dd_chain.items = create_chain_dropdown_items(comp)
            # Set default selections if required.
            if set_default_values and i == 0:
                btn_fixed.selected = True
                btn_moving.selected = False
                if ln_dd_chain.enabled:
                    dd_chain.items[0].selected = True
            elif set_default_values and i == 1:
                btn_fixed.selected = False
                btn_moving.selected = True
                if ln_dd_chain.enabled:
                    dd_chain.items[0].selected = True
            comp_list.items.append(ln)

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
    
    def btn_moving_clicked(self, btn):
        """Only one fixed strcuture can be selected at a time."""
        btns_to_update = [btn]
        selected_count = 0
        for menu_item in self.ln_moving_comp_list.get_content().items:
            btn_moving = menu_item.find_node('btn_moving').get_content()
            if btn_moving.selected:
                selected_count += 1
        label_text = 'Moving Structures'
        if selected_count > 0:
            label_text += f' ({selected_count})'
        self.lbl_moving_structures.text_value = label_text
        self.plugin.update_content(self.lbl_moving_structures, *btns_to_update)

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
    def btn_docs(self):
        return self._menu.root.find_node('btn_docs').get_content()

    @property
    def btn_advanced(self):
        return self._menu.root.find_node('btn_advanced').get_content()

    @async_callback
    async def on_mode_selected(self, btn, update=True, log=True):
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
            self.render(complexes=self.plugin.complexes)
        elif btn.name == 'btn_align_by_chain':
            if log:
                Logs.message("Switched to chain mode.")
            self.current_mode = 'chain'
            # Get deep complexes if necessary
            for comp in self.plugin.complexes:
                if sum(1 for _ in comp.chains) == 0:
                    self.btn_align_by_chain.unusable = True
                    self.plugin.update_content(self.btn_align_by_chain)
                    comp_indices = [cmp.index for cmp in self.plugin.complexes]
                    self.plugin.complexes = await self.plugin.request_complexes(comp_indices)
                    self.btn_align_by_chain.unusable = False
                    break
        await self.plugin.menu.render(complexes=self.plugin.complexes)
        if update:
            self.plugin.update_menu(self._menu)

    def get_moving_comp_indices_and_chains(self):
        comp_chain_list = []
        lst = self.ln_moving_comp_list.get_content()
        for item in lst.items:
            btn = item.find_node('btn_moving').get_content()
            chain_dd = item.find_node('dd_chain').get_content()
            comp_name = item.find_node('lbl_struct_name').get_content().text_value
            if not btn.selected:
                continue

            selected_comp = next((
                comp for comp in self.plugin.complexes
                if comp.full_name == comp_name
            ), None)
            selected_chain = None
            for chain_ddi in chain_dd.items:
                if chain_ddi.selected:
                    selected_chain = chain_ddi.name
                    break
                print('here')
            comp_chain_list.append((selected_comp.index, selected_chain))
        return comp_chain_list

    def open_rmsd_table(self, btn):
        self.rmsd_menu.enabled = True
        self.plugin.update_menu(self.rmsd_menu)

    def open_docs_page(self, btn):
        docs_url = 'https://docs.nanome.ai/plugins/cheminteractions.html'
        self.plugin.open_url(docs_url)
