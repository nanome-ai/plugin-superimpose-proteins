import os
import shutil
import site_motif
import subprocess
import tempfile
from nanome.util import Logs


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
            
            try:
                subprocess.run(command.split(), timeout=60)
            except subprocess.TimeoutExpired:
                Logs.warning("SiteMotif: Timeout. SiteMotif may not have finished.")
                pass
            align_output = f'{output_dir}/align_output.txt'
            sites_dir.cleanup()
            with open(align_output, "r") as f:
                with open('align_output.txt', 'w') as w:
                    w.write(f.read())
            if not os.path.exists(align_output):
                raise Exception('No align output found')
            with open(align_output, 'r') as f:
                lines = f.readlines()
                if not lines:
                    raise Exception("No lines found in align output")
            
            max_paired_res = -1
            residue_alignment = ''
            pdb_1 = ''
            pdb_2 = ''
            for line in lines:
                line_split = line.split('\t')
                if line_split[0] == line_split[1]:
                    continue
                aligned_residues = line_split[-1]
                paired_residue_count = len(aligned_residues.split(' '))
                if paired_residue_count > max_paired_res:
                    max_paired_res = paired_residue_count
                    pdb_1 = line_split[0]
                    pdb_2 = line_split[1]
                    residue_alignment = aligned_residues

            # Replace pdb names with full paths
            for pdb_file in [binding_site_pdb, *pocket_pdbs]:
                if pdb_1 in pdb_file:
                    pdb_1 = pdb_file
                if pdb_2 in pdb_file:
                    pdb_2 = pdb_file
            return pdb_1, pdb_2, residue_alignment
    
    def parse_atom_pairs(self, pdb_1, pdb_2, alignment):
        pass
                
