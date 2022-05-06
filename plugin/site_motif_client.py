import os
import shutil
import subprocess
import tempfile


class SiteMotifClient:

    def find_match(self, binding_site_pdb:str, pocket_pdbs:list):
        """Finds the pocket that best matches provided binding site."""
        with tempfile.TemporaryDirectory() as motif_dir:
            shutil.copy(binding_site_pdb, motif_dir)
            for pdb in pocket_pdbs:
                shutil.copy(pdb, motif_dir)
                pass

            # Run pairs.py
            ouput = subprocess.run(
                ["Pairs", motif_dir],
                cwd=motif_dir)
        

            motif_path = '/home/mike/workspace/rmsd-2/site-motif/bin/site-motif'
            output = subprocess.run(
                [motif_path, motif_dir],
                shell=True,
                stdout=subprocess.PIPE)
            print('done?')
            
