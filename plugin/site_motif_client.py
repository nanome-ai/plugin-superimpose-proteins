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
            for pdb in pocket_pdbs:
                shutil.copy(pdb, sites_dir.name)
            pairs_filepath = site_motif.write_pairs(sites_dir.name, output_dir)
            pdb_size_filepath = site_motif.write_pdb_size(sites_dir.name, output_dir)
            script_path = f'{os.path.dirname(site_motif.__file__)}/pocket_matrix_mpi7.py'
            command = f"mpiexec -n 4 python {script_path} {sites_dir.name} {pairs_filepath} {pdb_size_filepath} {output_dir}"
            output = subprocess.run(
                command.split(),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            align_output = f'{output_dir}/align_output.txt'
            sites_dir.cleanup()
            if output.returncode != 0:
                raise Exception(output.stderr.decode('utf-8'))
            elif not os.path.exists(align_output):
                raise Exception('No align output found')
            with open(align_output, 'r') as f:
                lines = f.readlines()
                if not lines:
                    raise Exception("No lines found in align output")
                for line in lines:
                    if 'Best match' in line:
                        return line.split()[0]
                raise Exception('No best match found')
