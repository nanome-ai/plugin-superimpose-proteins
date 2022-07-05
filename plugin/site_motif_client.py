import os
import shutil
import site_motif
import subprocess
import tempfile


class SiteMotifClient:

    def find_match(self, binding_site_pdb: str, pocket_pdbs: list):
        """Finds the pocket that best matches provided binding site."""
        with tempfile.TemporaryDirectory() as output_dir:
            sites_dir = tempfile.TemporaryDirectory(dir=output_dir)
            shutil.copy(binding_site_pdb, sites_dir.name)
            for pdb in pocket_pdbs:
                shutil.copy(pdb, sites_dir.name)
            pairs_filepath = site_motif.write_pairs(sites_dir.name, output_dir)
            # Remove pairs not involving the fixed binding site
            with open(pairs_filepath, "r") as f:
                lines = f.readlines()
            with open(pairs_filepath, "w") as f:
                for line in lines:
                    binding_site_filename = os.path.basename(binding_site_pdb)
                    if binding_site_filename in line:
                        f.write(line)
            pdb_size_filepath = site_motif.write_pdb_size(sites_dir.name, output_dir)
            script_path = f'{os.path.dirname(site_motif.__file__)}/pocket_matrix_mpi7.py'
            command = f"mpiexec -n 4 python {script_path} {sites_dir.name} {pairs_filepath} {pdb_size_filepath} {output_dir}"
            output = subprocess.run(command.split())
            align_output = f'{output_dir}/align_output.txt'
            sites_dir.cleanup()
            if not os.path.exists(align_output):
                raise Exception('No align output found')
            with open(align_output, 'r') as f:
                lines = f.readlines()
                if not lines:
                    raise Exception("No lines found in align output")
                for line in lines:
                    if 'Best match' in line:
                        return line.split()[0]
                raise Exception('No best match found')
