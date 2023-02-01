import asyncio
import os
import unittest
from unittest.mock import patch
from random import randint
from unittest.mock import MagicMock

import nanome
from nanome.api.structure import Complex, Substructure
from plugin.SuperimposePlugin import SuperimposePlugin
from plugin.enums import OverlayMethodEnum

fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures')


def run_awaitable(awaitable, *args, **kwargs):
    loop = asyncio.get_event_loop()
    if loop.is_running:
        loop = asyncio.new_event_loop()
    result = loop.run_until_complete(awaitable(*args, **kwargs))
    loop.close()
    return result


class PluginFunctionTestCase(unittest.TestCase):

    def setUp(self):
        file_4hhb = f'{fixtures_dir}/4hhb.pdb'
        file_1mbo = f'{fixtures_dir}/1mbo.pdb'
        self.complex_4hhb = Complex.io.from_pdb(path=file_4hhb)
        self.complex_1mbo = Complex.io.from_pdb(path=file_1mbo)
        for atom in self.complex_4hhb.atoms:
            atom.index = randint(1000000000, 9999999999)
        for atom in self.complex_1mbo.atoms:
            atom.index = randint(1000000000, 9999999999)
        self.plugin_instance = SuperimposePlugin()
        self.plugin_instance.start()
        self.plugin_instance._network = MagicMock()
        nanome._internal.network.PluginNetwork._instance = MagicMock()

    def tearDown(self):
        self.plugin_instance.on_stop()

    def test_superimpose_by_chain(self):
        # Make sure clean_complex function returns valid pdb can be parsed into a Complex structure.
        chain_name_4hhb = 'A'
        chain_name_1mbo = 'A'

        request_complexes_mock = MagicMock()
        fut = asyncio.Future()
        fut.set_result([self.complex_4hhb, self.complex_1mbo])
        request_complexes_mock.return_value = fut
        self.plugin_instance.request_complexes = request_complexes_mock

        update_structures_mock = MagicMock()
        update_fut = asyncio.Future()
        update_fut.set_result([self.complex_1mbo])
        update_structures_mock.return_value = update_fut
        self.plugin_instance.update_structures_deep = update_structures_mock

        alignment_method = OverlayMethodEnum.ALPHA_CARBONS_ONLY
        moving_comp_chain_list = [(self.complex_1mbo.index, chain_name_1mbo)]
        result = run_awaitable(
            self.plugin_instance.superimpose_by_chain,
            self.complex_4hhb.index, chain_name_4hhb,
            moving_comp_chain_list,
            alignment_method
        )
        expected_output = {
            self.complex_1mbo.full_name: {'rmsd': 1.97, 'paired_atoms': 141, 'paired_residues': 141, 'chain': 'A'}
        }
        self.assertEqual(result, expected_output)

    def test_superimpose_by_entry(self):
        request_complexes_mock = MagicMock()
        fut = asyncio.Future()
        fut.set_result([self.complex_4hhb, self.complex_1mbo])
        request_complexes_mock.return_value = fut
        self.plugin_instance.request_complexes = request_complexes_mock
        
        update_structures_mock = MagicMock()
        update_fut = asyncio.Future()
        update_fut.set_result([self.complex_1mbo])
        update_structures_mock.return_value = update_fut
        self.plugin_instance.update_structures_deep = update_structures_mock

        alignment_method = OverlayMethodEnum.ALPHA_CARBONS_ONLY
        result = run_awaitable(
            self.plugin_instance.superimpose_by_entry,
            self.complex_4hhb.index,
            [self.complex_1mbo.index],
            alignment_method
        )
        expected_result = {self.complex_1mbo.full_name: {'paired_atoms': 141, 'paired_residues': 141, 'rmsd': 1.78}}
        self.assertEqual(result, expected_result)

    def test_superimpose_by_entry_heavy_atoms(self):
        request_complexes_mock = MagicMock()
        fut = asyncio.Future()
        fut.set_result([self.complex_4hhb, self.complex_1mbo])
        request_complexes_mock.return_value = fut
        self.plugin_instance.request_complexes = request_complexes_mock

        update_structures_mock = MagicMock()
        update_fut = asyncio.Future()
        update_fut.set_result([self.complex_1mbo])
        update_structures_mock.return_value = update_fut
        self.plugin_instance.update_structures_deep = update_structures_mock

        alignment_method = OverlayMethodEnum.HEAVY_ATOMS_ONLY
        result = run_awaitable(
            self.plugin_instance.superimpose_by_entry,
            self.complex_4hhb.index,
            [self.complex_1mbo.index],
            alignment_method
        )
        expected_result = {self.complex_1mbo.full_name: {'paired_atoms': 366, 'paired_residues': 46, 'rmsd': 1.93}}
        self.assertEqual(result, expected_result)

    def test_superimpose_by_chain_heavy_atoms(self):
        # Make sure clean_complex function returns valid pdb can be parsed into a Complex structure.
        chain_name_4hhb = 'A'
        chain_name_1mbo = 'A'

        request_complexes_mock = MagicMock()
        fut = asyncio.Future()
        fut.set_result([self.complex_4hhb, self.complex_1mbo])
        request_complexes_mock.return_value = fut
        self.plugin_instance.request_complexes = request_complexes_mock

        update_structures_mock = MagicMock()
        update_fut = asyncio.Future()
        update_fut.set_result([self.complex_1mbo])
        update_structures_mock.return_value = update_fut
        self.plugin_instance.update_structures_deep = update_structures_mock

        alignment_method = OverlayMethodEnum.HEAVY_ATOMS_ONLY
        moving_comp_chain_list = [(self.complex_1mbo.index, chain_name_1mbo)]
        result = run_awaitable(
            self.plugin_instance.superimpose_by_chain,
            self.complex_4hhb.index, chain_name_4hhb,
            moving_comp_chain_list,
            alignment_method
        )
        expected_output = {self.complex_1mbo.full_name: {'rmsd': 1.91, 'paired_atoms': 378, 'paired_residues': 47, 'chain': 'A'}}
        self.assertEqual(result, expected_output)

    @unittest.skip("Sitemotif not working from tests")
    @patch('nanome.api.structure.Molecule.get_ligands')
    @patch('nanome.api.plugin_instance.PluginInstance.request_complexes')
    def test_superimpose_by_binding_site(self, request_complexes_mock, get_ligands_mock, *mocks):
        # Load kinases
        kinase_folder = os.path.join(fixtures_dir, 'kinases')
        complex_list = []
        for i, pdb_file in enumerate(os.listdir(kinase_folder)):
            new_comp = Complex.io.from_pdb(path=os.path.join(kinase_folder, pdb_file))
            new_comp.full_name = pdb_file.split('.')[0]
            new_comp.index = i + 1
            complex_list.append(new_comp)

        fut = asyncio.Future()
        fut.set_result(complex_list)
        request_complexes_mock.return_value = fut
        fixed_comp_name = '2OIB'
        fixed_comp = next(comp for comp in complex_list if comp.full_name == fixed_comp_name)
        moving_comp_indices = [cmp.index for cmp in complex_list if cmp.index != fixed_comp.index]

        # Build mock substructure.
        ligand_index = 0
        substruct = Substructure()
        substruct._name = 'GLN#341 : TPO#345'
        substruct._residues = [
            res for res in next(fixed_comp.molecules).residues
            if res.serial >= 341 and res.serial <= 345
        ]
        fut = asyncio.Future()
        fut.set_result([substruct])
        get_ligands_mock.return_value = fut

        result = run_awaitable(
            self.plugin_instance.superimpose_by_binding_site,
            fixed_comp.index,
            ligand_index,
            moving_comp_indices,
            OverlayMethodEnum.HEAVY_ATOMS_ONLY,
            self.plugin_instance
        )
        self.assertEqual(result, {})
