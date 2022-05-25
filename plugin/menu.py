import functools
import time
from os import path
from nanome.api import ui
from nanome.util import Logs, async_callback, Color
from nanome.util.enums import NotificationTypes, VertAlignOptions
from .enums import AlignmentModeEnum, AlignmentMethodEnum

BASE_PATH = path.dirname(f'{path.realpath(__file__)}')
MENU_PATH = path.join(BASE_PATH, 'menu_json', 'menu.json')
COMP_LIST_ITEM_PATH = path.join(BASE_PATH, 'menu_json', 'comp_list_item.json')
RMSD_MENU_PATH = path.join(BASE_PATH, 'menu_json', 'rmsd_menu.json')
RMSD_TABLE_ENTRY = path.join(BASE_PATH, 'menu_json', 'rmsd_list_entry.json')

GEAR_ICON_PATH = path.join(BASE_PATH, 'assets', 'gear.png')
GOLD_PIN_ICON_PATH = path.join(BASE_PATH, 'assets', 'TargetReferenceIcon.png')
DASHED_PIN_ICON_PATH = path.join(BASE_PATH, 'assets', 'TargetReferenceHoverIcon.png')
TRANSPARENCY_PATH = path.join(BASE_PATH, 'assets', 'transparent.png')
LOAD_ICON_PATH = path.join(BASE_PATH, 'assets', 'LoadIcon.png')
EXPORT_ICON_PATH = path.join(BASE_PATH, 'assets', 'Export.png')

DOCS_URL = 'https://docs.nanome.ai/plugins/cheminteractions.html'


def create_chain_buttons(comp, set_default=False):
    """Update chain dropdown to reflect changes in complex."""
    list_items = []
    # Filter out hetatom chains (HA, HB, etc)
    chain_names = [
        ch.name for ch in comp.chains
        if not ch.name.startswith('H') or len(ch.name) < 2
    ]
    for chain_name in chain_names:
        btn_ln = ui.LayoutNode()
        new_btn = btn_ln.add_new_button(chain_name)
        new_btn.toggle_on_press = True
        list_items.append(btn_ln)
    if set_default and list_items:
        list_items[0].selected = True
    return list_items


class MainMenu:

    def __init__(self, plugin_instance):
        super().__init__()
        self._menu = ui.Menu.io.from_json(MENU_PATH)
        self.plugin = plugin_instance
        self.btn_advanced.icon.value.set_all(GEAR_ICON_PATH)
        self.rmsd_menus = []

        self.load_icon.file_path = LOAD_ICON_PATH
        self.current_mode = AlignmentModeEnum.ENTRY
        for btn in self.mode_selection_btn_group:
            btn.register_pressed_callback(self.on_mode_selected)
        self.btn_rmsd_table.register_pressed_callback(self.open_rmsd_menu)
        self.btn_docs.register_pressed_callback(self.open_docs_page)
        self.btn_select_all.register_pressed_callback(
            functools.partial(self.toggle_all_moving_complexes, True))
        self.btn_deselect_all.register_pressed_callback(
            functools.partial(self.toggle_all_moving_complexes, False))

    @property
    def btn_submit(self):
        return self._menu.root.find_node('ln_submit').get_content()

    @property
    def ln_btn_align_by_binding_site(self):
        return self._menu.root.find_node('ln_btn_align_by_binding_site')

    @property
    def ln_binding_site_mode(self):
        return self._menu.root.find_node('ln_binding_site_mode')

    @property
    def btn_align_by_binding_site(self):
        return self.ln_btn_align_by_binding_site.get_content()

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
    def ln_empty_list(self):
        return self._menu.root.find_node('ln_empty_list')

    @property
    def load_icon(self):
        return self._menu.root.find_node('ln_load_icon').get_content()

    @property
    def lbl_moving_structures(self):
        return self._menu.root.find_node('lbl_moving_structures').get_content()

    @property
    def dd_align_using(self):
        return self._menu.root.find_node('ln_align_using').get_content()

    @property
    def ln_loading_bar(self):
        return self._menu.root.find_node('ln_loading_bar')

    @property
    def loading_bar(self):
        return self.ln_loading_bar.get_content()

    @property
    def btn_select_all(self):
        return self._menu.root.find_node('ln_btn_select_all').get_content()

    @property
    def btn_deselect_all(self):
        return self._menu.root.find_node('ln_btn_deselect_all').get_content()

    @async_callback
    async def render(self, complexes=None):
        complexes = complexes or []
        self.ln_binding_site_mode.enabled = False
        self.populate_comp_list(complexes, self.current_mode)
        self.check_if_ready_to_submit()
        self.btn_submit.register_pressed_callback(self.submit)
        self.plugin.update_menu(self._menu)

    @async_callback
    async def submit(self, btn):
        current_mode = self.current_mode
        self.btn_submit.unusable = True
        self.btn_submit.text.value.unusable = "Calculating..."
        self.plugin.update_content(self.btn_submit)
        fixed_comp_index = self.get_fixed_comp_index() or 0

        # Get alignment method based on dropdown selection
        heavy_atoms_method = 'Heavy atoms only'
        alignment_method = None
        selected_alignment_method = next((
            item.name for item in self.dd_align_using.items
            if item.selected), None)

        if selected_alignment_method == heavy_atoms_method:
            alignment_method = AlignmentMethodEnum.HEAVY_ATOMS_ONLY
        else:
            alignment_method = AlignmentMethodEnum.ALPHA_CARBONS_ONLY
        Logs.message("Submit button Pressed.")

        self.ln_loading_bar.enabled = True
        self.loading_bar.percentage = 0
        self.plugin.update_node(self.ln_loading_bar)

        start_time = time.time()
        rmsd_results = None
        moving_comp_count = 0
        try:
            if current_mode == AlignmentModeEnum.ENTRY:
                moving_comp_indices = self.get_moving_comp_indices()
                moving_comp_count = len(moving_comp_indices)
                Logs.message(f"Superimposing {moving_comp_count} structures by {current_mode.name.lower()}, using {alignment_method.name.lower()}")
                rmsd_results = await self.plugin.superimpose_by_entry(fixed_comp_index, moving_comp_indices, alignment_method)
            elif current_mode == AlignmentModeEnum.CHAIN:
                fixed_chain = self.get_fixed_chain()
                moving_comp_chain_list = self.get_moving_comp_indices_and_chains()
                moving_comp_count = len(moving_comp_chain_list)
                Logs.message(f"Superimposing {moving_comp_count} structures by {current_mode.name.lower()}, using {alignment_method.name.lower()}")
                rmsd_results = await self.plugin.superimpose_by_chain(fixed_comp_index, fixed_chain, moving_comp_chain_list, alignment_method)
            elif current_mode == AlignmentModeEnum.BINDING_SITE:
                ligand_name = self.get_binding_site_ligand()
                moving_comp_indices = self.get_moving_comp_indices()
                moving_comp_count = len(moving_comp_indices)
                if not all([fixed_comp_index, ligand_name, moving_comp_indices]):
                    msg = "Please select all complexes and chains."
                    Logs.warning(msg)
                    self.plugin.send_notification(NotificationTypes.error, msg)
                else:
                    Logs.message(f"Superimposing {moving_comp_count} structures by {current_mode.name.lower()}, using {alignment_method.name.lower()}")
                    rmsd_results = await self.plugin.superimpose_by_binding_site(
                        fixed_comp_index, ligand_name, moving_comp_indices)
        except Exception as e:
            rmsd_results = {}
            Logs.error("Error calculating Superposition.")

        if rmsd_results:
            fixed_name = next(comp.full_name for comp in self.plugin.complexes if comp.index == fixed_comp_index)
            if current_mode == AlignmentModeEnum.CHAIN:
                fixed_name = f'{fixed_name} Chain {fixed_chain}'
            self.render_rmsd_results(rmsd_results, fixed_name)
            self.ln_btn_rmsd_table.enabled = True
        self.btn_submit.unusable = False
        self.btn_submit.text.value.unusuable = "Superimpose"
        self.ln_loading_bar.enabled = False
        self.plugin.update_node(self.ln_btn_rmsd_table, self.ln_loading_bar)
        self.plugin.update_content(self.btn_submit)
        end_time = time.time()
        # Log data about run
        elapsed_time = round(end_time - start_time, 2)
        log_extra = {
            'alignment_method': selected_alignment_method,
            'alignment_mode': current_mode.name,
            'moving_complexes': moving_comp_count,
            'elapsed_time': elapsed_time
        }
        msg = f"Superimpose completed in {elapsed_time} seconds."
        Logs.message(msg, extra=log_extra)
        self.plugin.send_notification(NotificationTypes.success, msg)

    def render_rmsd_results(self, rmsd_results, fixed_comp_name):
        """Render rmsd results in a list."""
        rmsd_menu = RMSDMenu(self.plugin)
        self.rmsd_menus.append(rmsd_menu)
        rmsd_menu.index = 255 - len(self.rmsd_menus)
        rmsd_menu.render(rmsd_results, fixed_comp_name, run_number=len(self.rmsd_menus))

    def get_fixed_comp_index(self):
        for item in self._menu.root.find_node('ln_moving_comp_list').get_content().items:
            ln_btn_fixed = item.find_node('ln_btn_fixed')
            if not ln_btn_fixed:
                # This is usually the row with the label for hidden complexes.
                # Just skip.
                continue
            btn_fixed = ln_btn_fixed.get_content()
            if btn_fixed.selected:
                comp = next((
                    comp for comp in self.plugin.complexes
                    if comp.index == item.comp_index), None)
                return getattr(comp, 'index', None)

    def get_moving_comp_indices(self):
        comps = []
        for item in self._menu.root.find_node('ln_moving_comp_list').get_content().items:
            if not item.find_node('ln_btn_moving'):
                continue
            btn_moving = item.find_node('ln_btn_moving').get_content()
            if btn_moving.selected:
                comp = next((
                    comp for comp in self.plugin.complexes
                    if comp.index == item.comp_index), None)
                if comp:
                    comps.append(comp.index)
        return comps

    def get_binding_site_ligand(self):
        for item in self._menu.root.find_node('ln_moving_comp_list').get_content().items:
            btn_fixed = item.find_node('ln_btn_fixed').get_content()
            if btn_fixed.selected:
                dd_chain = item.find_node('dd_chain').get_content()
                selected_ddi = next((ddi for ddi in dd_chain.items if ddi.selected), None)
                return getattr(selected_ddi, 'name', '')

    def get_fixed_chain(self):
        for item in self.ln_moving_comp_list.get_content().items:
            btn_fixed = item.find_node('ln_btn_fixed').get_content()
            if btn_fixed.selected:
                dd_chain = item.find_node('dd_chain').get_content()
                selected_ddi = next((ddi for ddi in dd_chain.items if ddi.selected), None)
                return getattr(selected_ddi, 'name', '')

    def populate_comp_list(self, complexes, mode=AlignmentModeEnum.ENTRY):
        comp_list = self.ln_moving_comp_list.get_content()
        set_default_values = len(complexes) == 2
        visible_items = []
        hidden_items = []
        if len(complexes) == 0:
            self.ln_moving_comp_list.enabled = False
            self.ln_empty_list.enabled = True
            return
        else:
            self.ln_moving_comp_list.enabled = True
            self.ln_empty_list.enabled = False

        for i, comp in enumerate(complexes):
            ln = ui.LayoutNode.io.from_json(COMP_LIST_ITEM_PATH)
            ln.comp_index = comp.index
            btn_fixed = ln.find_node('ln_btn_fixed').get_content()
            btn_fixed.selected = False
            btn_fixed.icon.value.set_each(
                selected=GOLD_PIN_ICON_PATH,
                highlighted=DASHED_PIN_ICON_PATH,
                idle=TRANSPARENCY_PATH,
                selected_highlighted=GOLD_PIN_ICON_PATH)
            btn_moving = ln.find_node('ln_btn_moving').get_content()
            lbl_struct_name = ln.find_node('lbl_struct_name').get_content()
            ln_chain_list = ln.find_node('ln_chain_list')
            list_chain_btns = ln_chain_list.get_content()
            ln_chain_list.remove_content()

            ln_dd_chain = ln.find_node('dd_chain')
            ln_dd_chain.enabled = False
            dd_chain = ln_dd_chain.get_content()

            btn_fixed.register_pressed_callback(self.btn_fixed_clicked)
            btn_moving.register_pressed_callback(functools.partial(self.btn_moving_clicked, dd_chain))

            lbl_struct_name.text_value = comp.full_name

            overflow_size = 15
            if len(lbl_struct_name.text_value) > overflow_size:
                letters_to_keep = overflow_size - 3
                lbl_struct_name.text_value = lbl_struct_name.text_value[:letters_to_keep] + '...'

            btn_fixed.toggle_on_press = True
            btn_moving.toggle_on_press = True

            # Set up chain dropdown if in chain mode
            if mode == AlignmentModeEnum.CHAIN:
                ln_chain_list.enabled = True
                dd_chain.register_item_clicked_callback(
                    functools.partial(self.chain_selected_callback, btn_fixed, btn_moving))
                comp = next(
                    cmp for cmp in self.plugin.complexes
                    if cmp.index == ln.comp_index)
                ln_btns = create_chain_buttons(comp)
                for ln_btn in ln_btns:
                    ln_chain_list.add_child(ln_btn)
                
            # Set default selections if required.
            if set_default_values and i == 0:
                btn_fixed.selected = True
                btn_moving.selected = False
                if ln_dd_chain.enabled and dd_chain.items:
                    dd_chain.items[0].selected = True
            elif set_default_values and i == 1:
                btn_fixed.selected = False
                btn_moving.selected = True
                if ln_dd_chain.enabled and dd_chain.items:
                    dd_chain.items[0].selected = True

            if comp.visible:
                visible_items.append(ln)
            else:
                hidden_items.append(ln)
                continue

        comp_list.items = visible_items
        if hidden_items:
            hidden_item_header = ui.LayoutNode()
            hidden_item_header.set_padding(left=0.02)
            label = hidden_item_header.add_new_label(f"Hidden Items ({len(hidden_items)})")
            label.text_auto_size = False
            label.text_size = .3
            label.text_vertical_align = VertAlignOptions.Middle
            comp_list.items.append(hidden_item_header)
            comp_list.items.extend(hidden_items)

        comp_list.display_rows = min(max(len(comp_list.items), 4), 6)
        self.plugin.update_node(self.ln_moving_comp_list)

    def chain_selected_callback(self, btn_fixed, btn_moving, dd, ddi):
        content_to_update = [dd]
        if not btn_fixed.selected and not btn_moving.selected:
            btn_moving.selected = True
            content_to_update.append(btn_moving)
        self.plugin.update_content(*content_to_update)
        self.check_if_ready_to_submit()

    @async_callback
    async def btn_fixed_clicked(self, btn):
        """Only one fixed strcuture can be selected at a time."""
        btns_to_update = [btn]
        for menu_item in self.ln_moving_comp_list.get_content().items:
            if not menu_item.find_node('ln_btn_fixed'):
                continue
            btn_fixed = menu_item.find_node('ln_btn_fixed').get_content()
            btn_moving = menu_item.find_node('ln_btn_moving').get_content()
            ln_dd_chain = menu_item.find_node('dd_chain')
            if btn_fixed == btn:
                if btn_fixed.selected:
                    # Make sure the other button is not selected
                    btn_moving.selected = False
                    btn_moving.unusable = True
                    if self.current_mode == AlignmentModeEnum.BINDING_SITE:
                        ln_dd_chain.enabled = True
                        dd_ligand = ln_dd_chain.get_content()
                        dd_ligand.register_item_clicked_callback(self.check_if_ready_to_submit)
                        comp_name = menu_item.find_node('lbl_struct_name').get_content().text_value
                        comp = next(cmp for cmp in self.plugin.complexes if cmp.full_name == comp_name)
                        dd_ligand.items = await self.create_ligand_dropdown_items(comp)
                else:
                    btn_moving.unusable = False
            else:
                btn_fixed.selected = False
                btn_moving.unusable = False
                if self.current_mode == AlignmentModeEnum.BINDING_SITE:
                    ln_dd_chain.enabled = False

            btns_to_update.append(btn_fixed)
            btns_to_update.append(btn_moving)
        self.update_selection_counter()
        self.check_if_ready_to_submit()
        self.plugin.update_node(self.ln_moving_comp_list)

    async def create_ligand_dropdown_items(self, comp):
        # Get ligands for binding site dropdown
        mol = next(mo for mo in comp.molecules)
        ligands = await mol.get_ligands()
        dropdown_items = []
        for lig in ligands:
            dropdown_items.append(ui.DropdownItem(lig.name))
        return dropdown_items

    def btn_moving_clicked(self, dd_chain, btn):
        btns_to_update = [btn]
        selected_count = 0
        if not btn.selected and any(ddi.selected for ddi in dd_chain.items):
            for ddi in dd_chain.items:
                ddi.selected = False
            btns_to_update.append(dd_chain)

        for menu_item in self.ln_moving_comp_list.get_content().items:
            ln = menu_item.find_node('ln_btn_moving')
            if not ln:
                continue
            btn_moving = ln.get_content()
            if btn_moving.selected:
                selected_count += 1
        self.update_selection_counter()
        self.check_if_ready_to_submit()
        self.plugin.update_content(self.lbl_moving_structures, self.btn_submit, *btns_to_update)

    def update_selection_counter(self):
        counter = 0
        for menu_item in self.ln_moving_comp_list.get_content().items:
            if not menu_item.find_node('ln_btn_fixed'):
                continue
            btn_fixed = menu_item.find_node('ln_btn_fixed').get_content()
            btn_moving = menu_item.find_node('ln_btn_moving').get_content()
            if btn_fixed.selected:
                counter += 1
            if btn_moving.selected:
                counter += 1
        if counter == 0:
            new_text = f'Superimpose'
        else:
            new_text = f'Superimpose ({counter})'

        self.btn_submit.text.value.idle = new_text
        self.btn_submit.text.value.highlighted = new_text

    @property
    def mode_selection_btn_group(self):
        return [self.btn_entry_align, self.btn_align_by_chain, self.btn_align_by_binding_site]

    @property
    def btn_entry_align(self):
        return self._menu.root.find_node('ln_btn_entry_align').get_content()

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
        if btn.name == 'btn_entry_align':
            if log:
                Logs.message("Switched to entry mode")
            self.current_mode = AlignmentModeEnum.ENTRY
            self.render(complexes=self.plugin.complexes)
        elif btn.name == 'btn_align_by_chain':
            Logs.message("Switched to chain mode.")
            self.current_mode = AlignmentModeEnum.CHAIN
            # Get deep complexes if necessary
        elif btn.name == 'btn_align_by_binding_site':
            self.current_mode = AlignmentModeEnum.BINDING_SITE

        if self.current_mode in [AlignmentModeEnum.CHAIN, AlignmentModeEnum.BINDING_SITE]:
            for comp in self.plugin.complexes:
                if sum(1 for _ in comp.chains) == 0:
                    btn.unusable = True
                    self.plugin.update_content(btn)
                    comp_indices = [cmp.index for cmp in self.plugin.complexes]
                    self.plugin.complexes = await self.plugin.request_complexes(comp_indices)
                    btn.unusable = False
                    break

        await self.plugin.menu.render(complexes=self.plugin.complexes)
        if update:
            self.plugin.update_menu(self._menu)

    def get_moving_comp_indices_and_chains(self):
        comp_chain_list = []
        lst = self.ln_moving_comp_list.get_content()
        for item in lst.items:
            ln_btn = item.find_node('ln_btn_moving')
            if not ln_btn:
                continue
            btn = ln_btn.get_content()
            chain_dd = item.find_node('dd_chain').get_content()
            comp_index = item.comp_index
            if not btn.selected:
                continue

            selected_comp = next((
                comp for comp in self.plugin.complexes
                if comp.index == comp_index
            ), None)
            selected_chain = None
            for chain_ddi in chain_dd.items:
                if chain_ddi.selected:
                    selected_chain = chain_ddi.name
                    break
            comp_chain_list.append((selected_comp.index, selected_chain))
        return comp_chain_list

    def open_rmsd_menu(self, btn):
        if self.rmsd_menus:
            self.rmsd_menu = self.rmsd_menus[-1]
            self.rmsd_menu.enabled = True
            self.rmsd_menu.update()

    def open_docs_page(self, btn):
        self.plugin.open_url(DOCS_URL)

    def update_loading_bar(self, current, total):
        self.loading_bar.percentage = current / total
        self.plugin.update_content(self.loading_bar)

    def check_if_ready_to_submit(self):
        """Enable or disable submit button based on if required fields are selected."""
        fixed_comp_index = self.get_fixed_comp_index()
        if not fixed_comp_index:
            return False
        ready_to_submit = False
        if self.current_mode == AlignmentModeEnum.ENTRY:
            moving_comp_indices = self.get_moving_comp_indices()
            ready_to_submit = all([fixed_comp_index, moving_comp_indices])
        elif self.current_mode == AlignmentModeEnum.CHAIN:
            fixed_chain = self.get_fixed_chain()
            moving_comp_chain_list = self.get_moving_comp_indices_and_chains()
            moving_comps_selected = moving_comp_chain_list and all([
                bool(comp_index and chain) for comp_index, chain in moving_comp_chain_list
            ])
            if fixed_comp_index and fixed_chain and moving_comps_selected:
                ready_to_submit = True
        elif self.current_mode == AlignmentModeEnum.BINDING_SITE:
            fixed_ligand = self.get_fixed_chain()
            moving_comps_selected = bool(self.get_moving_comp_indices())
            if fixed_comp_index and fixed_ligand and moving_comps_selected:
                ready_to_submit = True

        self.btn_submit.unusable = not ready_to_submit
        self.plugin.update_content(self.btn_submit)

    def toggle_all_moving_complexes(self, value: bool, btn: ui.Button):
        """Select or deselect all complexes as moving complexes"""
        for item in self.ln_moving_comp_list.get_content().items:
            ln_btn_moving = item.find_node('ln_btn_moving')
            if not ln_btn_moving:
                continue
            btn_moving = ln_btn_moving.get_content()
            if not btn_moving.unusable:
                btn_moving.selected = value
            dd_chain = item.find_node('dd_chain').get_content()
            if dd_chain.items:
                dd_chain.items[0].selected = value
        self.plugin.update_node(self.ln_moving_comp_list)
        self.check_if_ready_to_submit()


class RMSDMenu(ui.Menu):

    def __init__(self, plugin_instance):
        super().__init__()
        self._menu = ui.Menu.io.from_json(RMSD_MENU_PATH)
        self.plugin = plugin_instance
        self._menu.enabled = False
        self.img_export.file_path = EXPORT_ICON_PATH
        self.btn_docs.register_pressed_callback(self.open_docs_page)
        self.btn_export.register_pressed_callback(self.export_as_csv)

    @property
    def btn_docs(self):
        return self._menu.root.find_node('btn_docs').get_content()

    @property
    def btn_export(self):
        return self._menu.root.find_node('ln_btn_export').get_content()

    @property
    def img_export(self):
        return self._menu.root.find_node('ln_img_export').get_content()

    def render(self, rmsd_results, fixed_comp_name, run_number=0):
        self.rmsd_results = rmsd_results
        results_list = self._menu.root.find_node('results_list').get_content()
        list_items = []
        row_color_dark = Color(21, 26, 37)
        row_color_light = Color(42, 52, 63)

        # Fixed comp is the first row.
        item = ui.LayoutNode().io.from_json(RMSD_TABLE_ENTRY)
        item_mesh = item.add_new_mesh()
        item_mesh.mesh_color = row_color_light
        item.get_children()[0].add_new_image(GOLD_PIN_ICON_PATH)
        item.get_children()[1].get_content().text_value = fixed_comp_name
        item.get_children()[2].get_content().text_value = '--'
        item.get_children()[3].get_content().text_value = '--'
        item.get_children()[4].get_content().text_value = '--'
        list_items.append(item)

        # Add moving comps and results to the table.
        for i, comp_name in enumerate(rmsd_results, 1):
            results_data = rmsd_results[comp_name]
            rmsd_val = results_data['rmsd']
            paired_atom_count = results_data['paired_atoms']
            paired_residue_count = results_data['paired_residues']
            if 'chain' in results_data:
                comp_name = f'{comp_name} Chain {results_data["chain"]}'

            item = ui.LayoutNode().io.from_json(RMSD_TABLE_ENTRY)
            item_mesh = item.add_new_mesh()
            item_mesh.mesh_color = row_color_light if i % 2 == 0 else row_color_dark

            item.get_children()[0].get_content().text_value = i
            item.get_children()[1].get_content().text_value = comp_name
            item.get_children()[2].get_content().text_value = rmsd_val
            item.get_children()[3].get_content().text_value = paired_residue_count
            item.get_children()[4].get_content().text_value = paired_atom_count
            list_items.append(item)
        results_list.items = list_items
        self._menu.title = f"RMSD Run {run_number}"

    def open_docs_page(self, btn):
        self.plugin.open_url(DOCS_URL)

    def update(self):
        self.plugin.update_menu(self._menu)

    @property
    def enabled(self):
        return self._menu._enabled

    @enabled.setter
    def enabled(self, value):
        self._menu._enabled = value

    @property
    def index(self):
        return self._menu.index

    @index.setter
    def index(self, value):
        self._menu.index = value

    def export_as_csv(self, btn):
        Logs.message("Exporting RMSD results to CSV...")
