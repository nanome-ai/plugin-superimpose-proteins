import asyncio
import functools
import time
import os
from nanome.api import ui
from nanome.util import Logs, async_callback, Color
from nanome.util.enums import NotificationTypes, ScalingOptions
from .enums import AlignmentModeEnum, OverlayMethodEnum

FEATURE_FLAG_BINDING_SITE = os.environ.get('FEATURE_FLAG_BINDING_SITE', "").lower() in ("yes", "true", "t", "1")

BASE_PATH = os.path.dirname(f'{os.path.realpath(__file__)}')
MENU_PATH = os.path.join(BASE_PATH, 'menu_json', 'menu.json')
COMP_LIST_ITEM_PATH = os.path.join(BASE_PATH, 'menu_json', 'comp_list_item.json')
RMSD_MENU_PATH = os.path.join(BASE_PATH, 'menu_json', 'rmsd_menu.json')
RMSD_TABLE_ENTRY = os.path.join(BASE_PATH, 'menu_json', 'rmsd_list_entry.json')

GEAR_ICON_PATH = os.path.join(BASE_PATH, 'assets', 'gear.png')
GOLD_PIN_ICON_PATH = os.path.join(BASE_PATH, 'assets', 'gold_pin.png')
GREY_PIN_ICON_PATH = os.path.join(BASE_PATH, 'assets', 'grey_pin.png')
LOAD_ICON_PATH = os.path.join(BASE_PATH, 'assets', 'LoadIcon.png')
EXPORT_ICON_PATH = os.path.join(BASE_PATH, 'assets', 'Export.png')


DOCS_URL = 'https://docs.nanome.ai/plugins/superimpose.html'


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
        self.dd_run_history.register_item_clicked_callback(self.open_rmsd_menu)
        self.btn_docs.register_pressed_callback(self.open_docs_page)
        self.btn_select_all.register_pressed_callback(
            functools.partial(self.toggle_all_moving_complexes, True))
        self.btn_deselect_all.register_pressed_callback(
            functools.partial(self.toggle_all_moving_complexes, False))

        overlay_button_group = [self.btn_alpha_carbons, self.btn_heavy_atoms]
        overlay_method_selected_callback = functools.partial(
            self.overlay_method_selected, overlay_button_group
        )
        self.btn_alpha_carbons.register_pressed_callback(overlay_method_selected_callback)
        self.btn_heavy_atoms.register_pressed_callback(overlay_method_selected_callback)
        self.btn_submit.register_pressed_callback(self.submit)

    @property
    def btn_submit(self):
        return self._menu.root.find_node('ln_submit').get_content()

    @property
    def ln_binding_site_mode(self):
        return self._menu.root.find_node('ln_binding_site_mode')

    @property
    def ln_binding_site_mode(self):
        return self._menu.root.find_node('ln_binding_site_mode')

    @property
    def btn_align_by_binding_site(self):
        return self.ln_btn_align_by_binding_site.get_content()

    @property
    def ln_run_history(self):
        return self._menu.root.find_node('ln_run_history')

    @property
    def dd_run_history(self):
        return self.ln_run_history.get_content()

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
    def btn_alpha_carbons(self):
        return self._menu.root.find_node('ln_btn_alpha_carbons').get_content()

    @property
    def btn_heavy_atoms(self):
        return self._menu.root.find_node('ln_btn_heavy_atoms').get_content()

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
    async def render(self, force_enable=False):

        self.ln_binding_site_mode.enabled = FEATURE_FLAG_BINDING_SITE
        self._menu.root.find_node('binding_site_spacer').enabled = not FEATURE_FLAG_BINDING_SITE
        await self.populate_comp_list()
        comp_list = self.ln_moving_comp_list.get_content()

        if len(comp_list.items) >= 2:
            first_item = comp_list.items[0]
            first_fixed_btn = first_item.find_node('ln_btn_fixed').get_content()
            first_fixed_btn.selected = True
            self.btn_fixed_pressed(first_fixed_btn)
        if len(comp_list.items) == 2:
            second_item = comp_list.items[1]
            second_moving_btn = second_item.find_node('ln_btn_moving').get_content()
            ln_chain_list = second_item.find_node('ln_chain_list')
            second_moving_btn.selected = True
            self.btn_moving_pressed(ln_chain_list, second_moving_btn)

        self.check_if_ready_to_submit()
        if force_enable:
            self._menu.enabled = True
        self.plugin.update_menu(self._menu)

    @async_callback
    async def submit(self, btn):
        current_mode = self.current_mode
        self.btn_submit.unusable = True
        original_unusable_text = self.btn_submit.text.value.unusable
        self.btn_submit.text.value.unusable = "Calculating..."
        self.plugin.update_content(self.btn_submit)
        fixed_comp_index = self.get_fixed_comp_index() or 0

        # Get alignment method based on dropdown selection
        heavy_atoms_method = 'btn_heavy_atoms'
        overlay_method = None

        selected_overlay_method = next(
            btn.name for btn in
            [self.btn_alpha_carbons, self.btn_heavy_atoms]
            if btn.selected)

        if selected_overlay_method.lower() == heavy_atoms_method:
            overlay_method = OverlayMethodEnum.HEAVY_ATOMS_ONLY
        else:
            overlay_method = OverlayMethodEnum.ALPHA_CARBONS_ONLY
        Logs.message("Submit button Pressed.")

        self.ln_loading_bar.enabled = True
        self.loading_bar.percentage = 0
        self.plugin.update_node(self.ln_loading_bar)

        start_time = time.time()
        rmsd_results = None
        moving_comp_count = 0
        run_successful = False
        try:
            if current_mode == AlignmentModeEnum.ENTRY:
                moving_comp_indices = self.get_moving_comp_indices()
                moving_comp_count = len(moving_comp_indices)
                Logs.message(f"Superimposing {moving_comp_count + 1} structures by {current_mode.name.lower()}, using {overlay_method.name.lower()}")
                rmsd_results = await self.plugin.superimpose_by_entry(fixed_comp_index, moving_comp_indices, overlay_method)
            elif current_mode == AlignmentModeEnum.CHAIN:
                fixed_chain = self.get_fixed_chain()
                moving_comp_chain_list = self.get_moving_comp_indices_and_chains()
                moving_comp_count = len(moving_comp_chain_list)
                Logs.message(f"Superimposing {moving_comp_count + 1} structures by {current_mode.name.lower()}, using {overlay_method.name.lower()}")
                rmsd_results = await self.plugin.superimpose_by_chain(fixed_comp_index, fixed_chain, moving_comp_chain_list, overlay_method)
            elif current_mode == AlignmentModeEnum.BINDING_SITE:
                ligand_name = self.get_binding_site_ligand()
                moving_comp_indices = self.get_moving_comp_indices()
                moving_comp_count = len(moving_comp_indices)
                if not all([fixed_comp_index, ligand_name, moving_comp_indices]):
                    msg = "Please select all complexes and ligand."
                    Logs.warning(msg)
                    self.plugin.send_notification(NotificationTypes.error, msg)
                else:
                    Logs.message(f"Superimposing {moving_comp_count + 1} structures by {current_mode.name.lower()}, using {overlay_method.name.lower()}")
                    rmsd_results = await self.plugin.superimpose_by_binding_site(
                        fixed_comp_index, ligand_name, moving_comp_indices)
            run_successful = True
        except Exception as e:
            rmsd_results = {}
            Logs.error("Error calculating Superposition.")

        if rmsd_results:
            fixed_name = next(comp.full_name for comp in self.plugin.complexes if comp.index == fixed_comp_index)
            if current_mode == AlignmentModeEnum.CHAIN:
                fixed_name = f'{fixed_name} Chain {fixed_chain}'
            self.render_rmsd_results(rmsd_results, fixed_name)
        self.btn_submit.unusable = False
        self.btn_submit.text.value.unusable = original_unusable_text
        self.ln_loading_bar.enabled = False
        self.plugin.update_node(self.ln_loading_bar)
        self.plugin.update_content(self.btn_submit)
        end_time = time.time()
        # Log data about run
        elapsed_time = round(end_time - start_time, 2)
        log_extra = {
            'overlay_method': overlay_method.name,
            'alignment_mode': current_mode.name,
            'moving_complexes': moving_comp_count,
            'elapsed_time': elapsed_time
        }
        if run_successful:
            msg = f"Superimpose completed in {elapsed_time} seconds."
            Logs.message(msg, extra=log_extra)
            self.plugin.send_notification(NotificationTypes.success, msg)
        else:
            msg = f"Superimpose failed after {elapsed_time} seconds."
            Logs.error(msg, extra=log_extra)
            self.plugin.send_notification(NotificationTypes.error, msg)

    def render_rmsd_results(self, rmsd_results, fixed_comp_name):
        """Render rmsd results in a list."""
        rmsd_menu = RMSDMenu(self.plugin)
        self.rmsd_menus.append(rmsd_menu)
        rmsd_menu.index = 255 - len(self.rmsd_menus)
        rmsd_menu.render(rmsd_results, fixed_comp_name, run_number=len(self.rmsd_menus))
        run_number = len(self.rmsd_menus)
        ddi = ui.DropdownItem(f"Run {run_number}")
        ddi.run_number = run_number
        self.dd_run_history.items.insert(0, ddi)
        self.dd_run_history.permanent_title = f"RMSD Tables ({run_number})"
        self.plugin.update_content(self.dd_run_history)

    def get_fixed_comp_index(self):
        for item in self.ln_moving_comp_list.get_content().items:
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
        for item in self.ln_moving_comp_list.get_content().items:
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
        for item in self.ln_moving_comp_list.get_content().items:
            ln_btn_fixed = item.find_node('ln_btn_fixed')
            if not ln_btn_fixed:
                continue
            btn_fixed = ln_btn_fixed.get_content()
            if btn_fixed.selected:
                dd_ligands = item.find_node('dd_ligands').get_content()
                selected_ddi = next((ddi for ddi in dd_ligands.items if ddi.selected), None)
                return getattr(selected_ddi, 'name', '')

    def get_fixed_chain(self):
        for item in self.ln_moving_comp_list.get_content().items:
            ln_btn_fixed = item.find_node('ln_btn_fixed')
            if not ln_btn_fixed:
                continue
            btn_fixed = ln_btn_fixed.get_content()
            if btn_fixed.selected:
                # Check if a chain has been selected
                ln_chain_btns = item.find_node('ln_chain_list').get_children()
                chain_btns = [ln.get_content() for ln in ln_chain_btns if ln.get_content()]
                selected_chain = next((
                    chain_btn.text.value.idle
                    for chain_btn in chain_btns
                    if chain_btn.selected), '')
                return selected_chain

    async def populate_comp_list(self):
        """Create list items representing each protein complex in the workspace."""
        complexes = self.plugin.complexes
        comp_list = self.ln_moving_comp_list.get_content()
        comp_list.display_rows = 4
        list_items = []

        template_list_item = self.create_template_list_item()
        # Create list items for each complex based on template
        for comp in complexes:
            # Skip any complexes that are only hetatoms.
            if not any(not atm.is_het for atm in comp.atoms):
                continue

            menu_item = template_list_item.clone()
            comp_index = comp.index
            menu_item.comp_index = comp_index
            await self.configure_comp_list_item(menu_item)
            list_items.append(menu_item)

        # If less than 2 proteins available, show empty list message on menu instead of list.
        if len(list_items) < 2:
            self.ln_moving_comp_list.enabled = False
            self.ln_empty_list.enabled = True
            return
        else:
            self.ln_moving_comp_list.enabled = True
            self.ln_empty_list.enabled = False

        comp_list.items = list_items
        self.plugin.update_node(self.ln_moving_comp_list)

    async def refresh_comp_list(self):
        comp_list = self.ln_moving_comp_list.get_content()

        if self.current_mode == AlignmentModeEnum.BINDING_SITE:
            self.btn_align_by_binding_site.unusable = True
            self.plugin.update_content(self.btn_align_by_binding_site)

        for menu_item in comp_list.items:
            await self.configure_comp_list_item(menu_item)

        if self.current_mode == AlignmentModeEnum.BINDING_SITE:
            self.btn_align_by_binding_site.unusable = False
            self.plugin.update_content(self.btn_align_by_binding_site)
        self.plugin.update_content(comp_list)

    async def configure_comp_list_item(self, menu_item):
        """Configure a list item for a complex based on current mode."""
        btn_fixed = menu_item.find_node('ln_btn_fixed').get_content()
        btn_moving = menu_item.find_node('ln_btn_moving').get_content()
        ln_chain_list = menu_item.find_node('ln_chain_list')
        chain_label = menu_item.find_node('select chains').get_content()
        ln_chain_selection = menu_item.find_node('chain_selection')
        ln_ligand_selection = menu_item.find_node('ligand_selection')
        ln_dd_ligands = ln_ligand_selection.find_node('dd_ligands')

        lbl_struct_name = menu_item.find_node('lbl_struct_name').get_content()
        comp = next((cmp for cmp in self.plugin.complexes if cmp.index == menu_item.comp_index), None)
        if not comp:
            Logs.warning(f"Could not find complex '{menu_item.comp_index}' in workspace.")
            return
        # Set up labels and text overflows.
        if self.current_mode == AlignmentModeEnum.CHAIN:
            chain_label.text_value = 'Select Chain'
        else:
            chain_label.text_value = 'Chains'
        lbl_struct_name.text_value = comp.full_name
        overflow_size = 25
        if len(lbl_struct_name.text_value) > overflow_size:
            letters_to_keep = overflow_size - 3
            lbl_struct_name.text_value = lbl_struct_name.text_value[:letters_to_keep] + '...'

        # Set up chain buttons
        if not ln_chain_list.children:
            ln_chain_list.children = self.create_chain_buttons(comp)
        # Chain buttons should only be clickable if in Chain mode.
        chain_btns = [ln_btn.get_content() for ln_btn in ln_chain_list.get_children() if ln_btn.get_content()]
        for btn in chain_btns:
            btn.unusable = self.current_mode != AlignmentModeEnum.CHAIN

        # Show or hide chain section
        ln_chain_selection.enabled = self.current_mode != AlignmentModeEnum.BINDING_SITE
        ln_ligand_selection.enabled = all([
            self.current_mode == AlignmentModeEnum.BINDING_SITE
            and btn_fixed.selected
        ])

        btn_fixed.toggle_on_press = True
        btn_moving.toggle_on_press = True
        btn_moving.register_pressed_callback(functools.partial(self.btn_moving_pressed, ln_chain_list))

        if self.current_mode == AlignmentModeEnum.BINDING_SITE:
            # Set up ligand dropdown if we're in binding site mode
            dd_ligands = ui.Dropdown()
            dd_ligands.max_displayed_items = 5
            dd_ligands.items = await self.create_ligand_dropdown_items(comp)
            dd_ligands.permanent_title = True
            dd_ligands.permanent_title = 'Select Ligand' if dd_ligands.items else 'No Ligands'
            ln_dd_ligands.set_content(dd_ligands)

    def create_template_list_item(self):
        """Create a template list item for the complex list."""
        list_item = ui.LayoutNode.io.from_json(COMP_LIST_ITEM_PATH)
        # Fixed button
        btn_fixed = list_item.find_node('ln_btn_fixed').get_content()
        btn_fixed.selected = False
        btn_fixed.toggle_on_press = True
        btn_fixed.icon.value.set_each(
            selected=GOLD_PIN_ICON_PATH,
            highlighted=GREY_PIN_ICON_PATH,
            idle=GREY_PIN_ICON_PATH,
            selected_highlighted=GOLD_PIN_ICON_PATH)
        btn_fixed.register_pressed_callback(self.btn_fixed_pressed)
        btn_fixed.toggle_on_press = True
        # Moving button
        btn_moving = list_item.find_node('ln_btn_moving').get_content()
        btn_moving.toggle_on_press = True
        return list_item

    def chain_selected_callback(self, comp_index, btn_group, pressed_btn):
        # One item in button group selected at a time.
        Logs.debug(f"Chain selected: {pressed_btn.text.value.idle}")
        comp = next(comp for comp in self.plugin.complexes if comp.index == comp_index)
        btns_to_update = [pressed_btn]
        for grp_btn in btn_group:
            if grp_btn is not pressed_btn and grp_btn.selected:
                grp_btn.selected = False
                btns_to_update.append(grp_btn)

        # If the corresponding moving or fixed button hasn't been selected, select it.
        for menu_item in self.ln_moving_comp_list.get_content().items:
            ln_chain_list = menu_item.find_node('ln_chain_list')
            chain_btns = [ln.get_content() for ln in ln_chain_list.get_children() if ln.get_content()]
            if pressed_btn in chain_btns:
                btn_fixed = menu_item.find_node('ln_btn_fixed').get_content()
                btn_moving = menu_item.find_node('ln_btn_moving').get_content()
                if not btn_fixed.selected and not btn_moving.unusable:
                    btn_moving.selected = pressed_btn.selected
                    btns_to_update.append(btn_moving)
                elif btn_fixed.selected and not pressed_btn.selected:
                    # Deselect fixed btn
                    btn_fixed.selected = False
                    btn_moving.unusable = False
                    btns_to_update.append(btn_fixed)
                    btns_to_update.append(btn_moving)

        self.plugin.update_content(btns_to_update)
        self.check_if_ready_to_submit()
        if self.current_mode == AlignmentModeEnum.CHAIN:
            chain_name = pressed_btn.text.value.idle
            self.toggle_chain_atoms_selected(comp, chain_name, pressed_btn.selected)

    @async_callback
    async def btn_fixed_pressed(self, pressed_btn):
        """Handle selections for fixed button."""
        Logs.message("Fixed protein selected")
        content_to_update = [pressed_btn]
        selected_comp_chain_btn = None
        deselected_comp_chain_btn = None
        comp_list = self.ln_moving_comp_list.get_content()
        for menu_item in comp_list.items:
            if not menu_item.find_node('ln_btn_fixed'):
                continue
            btn_fixed = menu_item.find_node('ln_btn_fixed').get_content()
            btn_moving = menu_item.find_node('ln_btn_moving').get_content()
            ln_chain_list = menu_item.find_node('ln_chain_list')
            chain_btns = [ln.get_content() for ln in ln_chain_list.get_children() if ln.get_content()]
            # Handle row who's fixed selection is being toggled.
            if btn_fixed is pressed_btn:
                if btn_fixed.selected:
                    # If Button being selected, disable moving button, and select first chain.
                    btn_moving.selected = False
                    btn_moving.unusable = True
                    default_selection_index = 0
                    chain_button_selections = [ch_btn.selected for ch_btn in chain_btns]
                    selected_chain_index = next(
                        (i for i, btn_sel in enumerate(chain_button_selections) if btn_sel),
                        default_selection_index)
                    for i, ch_btn in enumerate(chain_btns):
                        if i == selected_chain_index:
                            selected_comp_chain_btn = ch_btn
                        selected_val = i == selected_chain_index
                        if ch_btn.selected != selected_val:
                            ch_btn.selected = selected_val
                            content_to_update.append(ch_btn)
                else:
                    # If Button being deselected, enable moving button, and deselect all chains.
                    btn_moving.unusable = False
                    if chain_btns:
                        selected_comp_chain_btn = next((btn for btn in chain_btns if btn.selected), None)
                        if selected_comp_chain_btn:
                            selected_comp_chain_btn.selected = False
                            content_to_update.append(selected_comp_chain_btn)
            # Handle row other than the one that fixed selection is being toggled.
            else:
                btn_fixed.selected = False
                btn_moving.unusable = False
                # Unless moving structure is selected, deselect all chains
                if not btn_moving.selected:
                    for ch_btn in [btn for btn in chain_btns if btn.selected]:
                        ch_btn.selected = False
                        content_to_update.append(ch_btn)
                        deselected_comp_chain_btn = ch_btn

            content_to_update.append(btn_fixed)
            content_to_update.append(btn_moving)
            if self.current_mode == AlignmentModeEnum.BINDING_SITE:
                ln_chain_selection = menu_item.find_node('chain_selection')
                ln_chain_selection.enabled = False
                ln_lig_selection = menu_item.find_node('ligand_selection')
                ln_lig_selection.enabled = btn_fixed.selected
        self.update_selection_counter()
        self.check_if_ready_to_submit()
        self.plugin.update_content(comp_list)
        # Update chain selections. Save to the end for performance.
        if selected_comp_chain_btn:
            self.toggle_chain_button(selected_comp_chain_btn)
        if deselected_comp_chain_btn:
            self.toggle_chain_button(deselected_comp_chain_btn)

    async def create_ligand_dropdown_items(self, comp):
        # Get ligands for binding site dropdown
        mol = next(mo for mo in comp.molecules)
        try:
            ligands = await asyncio.wait_for(mol.get_ligands(), 10)
        except asyncio.TimeoutError:
            Logs.warning("get_ligands timeout error")
            ligands = []
        dropdown_items = []
        for lig in ligands:
            dropdown_items.append(ui.DropdownItem(lig.name))
        return dropdown_items

    def btn_moving_pressed(self, ln_chain_list, btn_moving):
        Logs.message("Moving protein selected")
        btns_to_update = [btn_moving]
        chain_btns = [
            ln.get_content() for ln in ln_chain_list.get_children()
            if ln.get_content()
        ]
        # If moving struct being selected, and no chains are already selected, select the first one
        if btn_moving.selected and chain_btns and not any(ch_btn.selected for ch_btn in chain_btns):
            first_chain_btn = chain_btns[0]
            first_chain_btn.selected = True
            self.toggle_chain_button(first_chain_btn)
            btns_to_update.append(first_chain_btn)
        # If moving struct being deselected, and any chains are already selected, deselect all
        elif not btn_moving.selected and any(ch_btn.selected for ch_btn in chain_btns):
            for ch_btn in [btn for btn in chain_btns if btn.selected]:
                ch_btn.selected = False
                self.toggle_chain_button(ch_btn)
                btns_to_update.append(ch_btn)

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
    async def on_mode_selected(self, mode_btn):
        mode_btn.selected = True
        btns_to_update = [mode_btn]
        for group_item in self.mode_selection_btn_group:
            if mode_btn._content_id != group_item._content_id:
                group_item.selected = False
                btns_to_update.append(group_item)
        if mode_btn.name == 'btn_entry_align':
            Logs.message("Switched to entry mode")
            self.current_mode = AlignmentModeEnum.ENTRY
        elif mode_btn.name == 'btn_align_by_chain':
            Logs.message("Switched to chain mode.")
            self.current_mode = AlignmentModeEnum.CHAIN
        elif mode_btn.name == 'btn_align_by_binding_site':
            Logs.message("Switched to binding site mode.")
            self.current_mode = AlignmentModeEnum.BINDING_SITE
        await self.refresh_comp_list()
        self.plugin.update_menu(self._menu)

    def get_moving_comp_indices_and_chains(self):
        comp_chain_list = []
        lst = self.ln_moving_comp_list.get_content()
        for item in lst.items:
            ln_btn = item.find_node('ln_btn_moving')
            if not ln_btn:
                continue
            btn_moving = ln_btn.get_content()
            if not btn_moving.selected:
                continue
            # Check if a chain has been selected
            ln_chain_btns = item.find_node('ln_chain_list').get_children()
            chain_btns = [ln.get_content() for ln in ln_chain_btns if ln.get_content()]
            selected_chain_btns = [ch_btn for ch_btn in chain_btns if ch_btn.selected]
            selected_chain = next((
                chain_btn.text.value.idle
                for chain_btn in selected_chain_btns), '')
            comp_index = item.comp_index
            comp_chain_list.append((comp_index, selected_chain))
        return comp_chain_list

    def open_rmsd_menu(self, dd, ddi):
        ddi.selected = False
        self.plugin.update_content(dd)
        run_number = ddi.run_number
        Logs.message(f"Opening results for run {run_number}")
        if self.rmsd_menus:
            rmsd_menu = next((m for m in self.rmsd_menus if m.run_number == run_number), None)
            rmsd_menu._menu.enabled = True
            rmsd_menu.update()

    def open_docs_page(self, btn):
        self.plugin.open_url(DOCS_URL)

    def update_loading_bar(self, current, total):
        self.loading_bar.percentage = current / total
        self.plugin.update_content(self.loading_bar)

    def check_if_ready_to_submit(self):
        """Enable or disable submit button based on if required fields are selected."""
        ready_to_submit = False
        fixed_comp_index = self.get_fixed_comp_index()
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
        content_to_update = []
        chain_btns_to_toggle = []

        for item in self.ln_moving_comp_list.get_content().items:
            ln_btn_moving = item.find_node('ln_btn_moving')
            if not ln_btn_moving:
                continue

            btn_moving = ln_btn_moving.get_content()
            ln_chain_btns = item.find_node('ln_chain_list').get_children()
            chain_btns = [ln.get_content() for ln in ln_chain_btns if ln.get_content()]
            btn_fixed = item.find_node('ln_btn_fixed').get_content()

            # If deselecting all, deselect fixed selection.
            if not value and btn_fixed.selected:
                btn_fixed.selected = False
                btn_moving.selected = False
                btn_moving.unusable = False
                content_to_update.append(btn_fixed)

            btn_moving.selected = value
            content_to_update.append(btn_moving)
            # Select first chain, or deselect all chains
            if chain_btns and not any([btn.selected for btn in chain_btns]):
                ch_btn = chain_btns[0]
                ch_btn.selected = value
                chain_btns_to_toggle.append(ch_btn)
                content_to_update.append(ch_btn)
            elif not value and any([btn.selected for btn in chain_btns]):
                for ch_btn in chain_btns:
                    if ch_btn.selected:
                        ch_btn.selected = False
                        chain_btns_to_toggle.append(ch_btn)
                        content_to_update.append(ch_btn)
        self.plugin.update_content(*content_to_update)
        self.check_if_ready_to_submit()
        for ch_btn in chain_btns_to_toggle:
            self.toggle_chain_button(ch_btn)

    def toggle_chain_button(self, chain_btn: ui.Button):
        if self.current_mode != AlignmentModeEnum.CHAIN:
            return
        comp_index = chain_btn.comp_index
        chain_name = chain_btn.text.value.idle
        comp = next((
            comp for comp in self.plugin.complexes
            if comp.index == comp_index), None)
        self.toggle_chain_atoms_selected(comp, chain_name, chain_btn.selected)

    def create_chain_buttons(self, comp, set_default=False):
        """Create list of LayoutNodes containing a button for each chain in comp."""
        list_items = []
        # Filter out hetatom chains (HA, HB, etc)
        chain_names = sorted([
            ch.name for ch in comp.chains
            if not ch.name.startswith('H') or len(ch.name) < 2
        ])
        for chain_name in chain_names:
            btn_ln = ui.LayoutNode()
            btn_ln.set_padding(top=0.01, down=0.01)
            new_btn = btn_ln.add_new_button(chain_name)
            new_btn.text.min_size = 0.1
            new_btn.text.max_size = 0.25
            new_btn.toggle_on_press = True
            new_btn.unusable = not self.current_mode == AlignmentModeEnum.CHAIN
            list_items.append(btn_ln)
        if set_default and list_items:
            list_items[0].selected = True
        # Set up button group callback
        btn_list = [ln.get_content() for ln in list_items]
        for ln_btn in list_items:
            btn = ln_btn.get_content()
            btn.comp_index = comp.index
            selected_callback_fn = functools.partial(
                self.chain_selected_callback, comp.index, btn_list
            )
            btn.register_pressed_callback(selected_callback_fn)

        # Add padding when less than 4 chains
        while len(list_items) < 4:
            list_items.append(ui.LayoutNode())
        return list_items

    def toggle_chain_atoms_selected(self, comp, chain_name, value: bool):
        """Select or deselect all atoms in a chain."""
        chain = next((ch for ch in comp.chains if ch.name == chain_name), None)
        if not chain:
            Logs.warning(f"Chain {chain_name} not found in comp. Skipping atom.")
            return

        atoms_to_update = []
        for chain in comp.chains:
            # Select or deselect all atoms in the provided chain
            if chain.name == chain_name:
                chain_atms_to_update = (atm for atm in chain.atoms if atm.selected != value)
                for atm in chain_atms_to_update:
                    atm.selected = value
                    atoms_to_update.append(atm)
            else:
                # Deselect all atoms not on selected chain
                selected_nonchain_atms = (atm for atm in chain.atoms if atm.selected)
                for atm in selected_nonchain_atms:
                    atm.selected = False
                    atoms_to_update.append(atm)
        if atoms_to_update:
            Logs.debug(f"Updating {len(atoms_to_update)} atom selections")
            self.plugin.update_structures_shallow(atoms_to_update)

    def overlay_method_selected(self, btn_group, selected_btn):
        """Callback for when an overlay method is selected."""
        btns_to_update = []
        for btn in btn_group:
            if btn == selected_btn and not btn.selected:
                btn.selected = True
                btns_to_update.append(btn)
            elif btn != selected_btn and btn.selected:
                btn.selected = False
                btns_to_update.append(btn)
        if btns_to_update:
            self.plugin.update_content(*btns_to_update)


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
        self.run_number = run_number
        results_list = self._menu.root.find_node('results_list').get_content()
        list_items = []
        row_color_dark = Color(21, 26, 37)
        row_color_light = Color(42, 52, 63)

        # Fixed comp is the first row.
        item = ui.LayoutNode().io.from_json(RMSD_TABLE_ENTRY)
        item_mesh = item.add_new_mesh()
        item_mesh.mesh_color = row_color_light
        image_ln = item.get_children()[0]
        image_ln.set_padding(top=0.025, down=0.025, left=0.025, right=0.025)
        image = image_ln.add_new_image(GOLD_PIN_ICON_PATH)
        image.scaling_options = ScalingOptions.fill
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
    def index(self):
        return self._menu.index

    @index.setter
    def index(self, value):
        self._menu.index = value

    def export_as_csv(self, btn):
        Logs.message("Exporting RMSD results to CSV...")
