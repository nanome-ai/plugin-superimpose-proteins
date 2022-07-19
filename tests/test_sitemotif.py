import os
import tempfile
import unittest
from random import randint


from nanome.api.structure import Complex
from plugin.site_motif_client import SiteMotifClient


fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures')


@unittest.skip("Sitemotif not working from tests")
class SiteMotifClientTestCase(unittest.TestCase):

    def setUp(self):
        self.pdb1 = f'{fixtures_dir}/kinases/2OIB.pdb'
        self.pdb2 = f'{fixtures_dir}/kinases/1JAM.pdb'
        self.comp1 = Complex.io.from_pdb(path=self.pdb1)
        self.comp2 = Complex.io.from_pdb(path=self.pdb2)
        self.client = SiteMotifClient()

    def test_find_match(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb1, pdb2, alignment = self.client.find_match(self.pdb1, [self.pdb2])
            self.assertTrue(len(alignment) > 1)
