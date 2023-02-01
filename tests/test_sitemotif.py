import os
import tempfile
import unittest
from random import randint


from nanome.api.structure import Complex
from plugin.site_motif_client import SiteMotifClient


fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures')


# @unittest.skip("Sitemotif not working from tests")
class SiteMotifClientTestCase(unittest.TestCase):

    def setUp(self):
        self.pdb1 = f'{fixtures_dir}/kinases/2OIB.pdb'
        self.pdb2 = f'{fixtures_dir}/kinases/1JAM.pdb'
        self.comp1 = Complex.io.from_pdb(path=self.pdb1)
        self.comp2 = Complex.io.from_pdb(path=self.pdb2)
        self.client = SiteMotifClient()

    @unittest.skip("Sitemotif not working from tests")
    def test_sitemotif_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self.client.run(self.pdb1, [self.pdb2], tmpdir)
            align_output_file = f'{tmpdir}/align_output.txt'
            with open(align_output_file, 'r') as f:
                lines = f.readlines()
                self.assertTrue(len(lines) > 1)
