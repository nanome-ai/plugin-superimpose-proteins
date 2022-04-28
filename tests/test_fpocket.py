import os
import tempfile
import unittest
from random import randint


from nanome.api.structure import Complex
from rmsd_v2.fpocket_client import FPocketClient


fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures')


class FPocketClientTestCase(unittest.TestCase):

    def setUp(self):
        file_4hhb = f'{fixtures_dir}/4hhb.pdb'
        file_1mbo = f'{fixtures_dir}/1mbo.pdb'
        self.complex_4hhb = Complex.io.from_pdb(path=file_4hhb)
        self.complex_1mbo = Complex.io.from_pdb(path=file_1mbo)
        for atom in self.complex_4hhb.atoms:
            atom.index = randint(1000000000, 9999999999)
        for atom in self.complex_1mbo.atoms:
            atom.index = randint(1000000000, 9999999999)
        self.client = FPocketClient()

    def test_fpocket(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = self.client.run(self.complex_4hhb, tmpdir)
            pocket_list = self.client.parse_results(self.complex_4hhb, output_dir)
        self.assertEqual(len(pocket_list), 24)

