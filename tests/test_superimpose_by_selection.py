import asyncio
import os
from nanome.api.structure import Complex, Substructure
from nanome.util.enums import OverlayMethodEnum

from plugin.superimpose_by_selection import superimpose_by_selection


fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures')

async def test_superimpose_by_selection(self):
    # Load kinases
    kinase_folder = os.path.join(fixtures_dir, 'kinases')
    complex_list = []
    for i, pdb_file in enumerate(os.listdir(kinase_folder)):
        new_comp = Complex.io.from_pdb(path=os.path.join(kinase_folder, pdb_file))
        new_comp.full_name = pdb_file.split('.')[0]
        new_comp.index = i + 1
        complex_list.append(new_comp)

    fixed_comp_name = '2OIB'
    fixed_comp = next(comp for comp in complex_list if comp.full_name == fixed_comp_name)
    moving_comp_indices = [cmp.index for cmp in complex_list if cmp.index != fixed_comp.index]

    # For now just testing one moving comp

    result = await superimpose_by_selection(
        fixed_comp,
        moving_comp,
        OverlayMethodEnum.HEAVY_ATOMS_ONLY,
        self.plugin_instance
    )
    self.assertEqual(result, {})