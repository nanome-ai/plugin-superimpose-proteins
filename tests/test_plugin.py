import asyncio
import os
import unittest
from unittest.mock import patch
from random import randint

from unittest.mock import MagicMock
from nanome.api.structure import Complex
from plugin.Superimpose import SuperimposePlugin


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

    @patch('nanome._internal._network.PluginNetwork._instance')
    @patch('nanome.api.plugin_instance.PluginInstance.update_structures_deep')
    @patch('nanome.api.plugin_instance.PluginInstance.request_complexes')
    def test_superimpose_by_chain(self, request_complexes_mock, update_structures_mock, *mocks):
        # Make sure clean_complex function returns valid pdb can be parsed into a Complex structure.
        chain_name_4hhb = 'A'
        chain_name_1mbo = 'A'

        fut = asyncio.Future()
        fut.set_result([self.complex_4hhb, self.complex_1mbo])
        request_complexes_mock.return_value = fut

        update_fut = asyncio.Future()
        update_fut.set_result([self.complex_1mbo])
        update_structures_mock.return_value = update_fut

        moving_comp_chain_list = [(self.complex_1mbo.index, chain_name_1mbo)]
        result = run_awaitable(
            self.plugin_instance.superimpose_by_chain,
            self.complex_4hhb.index, chain_name_4hhb,
            moving_comp_chain_list
        )
        expected_result = 1.95456
        expected_output = {self.complex_1mbo.full_name: expected_result}
        self.assertEqual(result, expected_output)

    @patch('nanome._internal._network.PluginNetwork._instance')
    @patch('nanome.api.plugin_instance.PluginInstance.update_structures_deep')
    @patch('nanome.api.plugin_instance.PluginInstance.request_complexes')
    def test_superimpose_by_entry(self, request_complexes_mock, update_structures_mock, *mocks):
        fut = asyncio.Future()
        fut.set_result([self.complex_4hhb, self.complex_1mbo])
        request_complexes_mock.return_value = fut
        update_fut = asyncio.Future()
        update_fut.set_result([self.complex_1mbo])
        update_structures_mock.return_value = update_fut

        result = run_awaitable(
            self.plugin_instance.superimpose_by_entry,
            self.complex_4hhb.index,
            [self.complex_1mbo.index]
        )
        expected_result = {self.complex_1mbo.full_name: 27.69646}
        self.assertEqual(result, expected_result)
