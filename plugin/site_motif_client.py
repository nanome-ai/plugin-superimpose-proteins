import shutil
import subprocess


class SiteMotifClient:

    def find_match(self, binding_site_pdb: str, pocket_pdbs: list):
        """Finds the pocket that best matches provided binding site."""
        site_folder = '/tmp/site_folder'
        # with tempfile.TemporaryDirectory() as site_folder:
        shutil.copy(binding_site_pdb, site_folder)
        for pdb in pocket_pdbs:
            shutil.copy(pdb, site_folder)

        site_motif_path = 'site-motif/bin/site-motif'
        output = subprocess.run(
            [site_motif_path, site_folder],
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        if output.returncode != 0:
            raise Exception(output.stderr.decode('utf-8'))