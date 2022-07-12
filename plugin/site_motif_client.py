import os
import shutil
import site_motif
import subprocess
import tempfile
from nanome.util import Logs
from nanome.api.structure import Complex
from Bio.PDB import PDBParser

class SiteMotifClient:

    @staticmethod
    def find_match(binding_site_pdb: str, pocket_pdbs: list):
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
            # Uncomment for debugging
            # with open(align_output, "r") as f:
            #     with open('align_output.txt', 'w') as w:
            #         w.write(f.read())
            if not os.path.exists(align_output):
                raise Exception('No align output found')
            with open(align_output, 'r') as f:
                lines = f.readlines()
                if not lines:
                    raise Exception("No lines found in align output")

            # Find the alignment that contains the longest fragment match
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
            Logs.message("Longest fragment match: " + str(max_paired_res))
            # Replace pdb names with full paths
            for pdb_file in [binding_site_pdb, *pocket_pdbs]:
                if pdb_1 in pdb_file:
                    pdb_1 = pdb_file
                if pdb_2 in pdb_file:
                    pdb_2 = pdb_file
                if pdb_1 and pdb_2:
                    break
            return pdb_1, pdb_2, residue_alignment
    
    @staticmethod
    def parse_residue_pairs(comp1: Complex, comp2: Complex, alignment: str):
        """Parses aligned residues from complexes.
        
        alignment format ex. ARG-A-334 GLN-A-155_VAL-A-343 LYS-A-198_ILE-A-364 ALA-A-151
        where comp1 is the first residue and comp2 is the second residue
        """
        alignment_pairings = alignment.strip().split(' ')
        parser = PDBParser(QUIET=True)
        comp1_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        comp2_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        comp1.io.to_pdb(comp1_pdb.name)
        comp2.io.to_pdb(comp2_pdb.name)
        
        struct1 = parser.get_structure(comp1.full_name, comp1_pdb.name)
        struct2 = parser.get_structure(comp2.full_name, comp2_pdb.name)
        comp1_atom_list = []
        comp2_atom_list = []
        for residue_pair in alignment_pairings:
            residues = residue_pair.split('_')
            if len(residues) != 2:
                # For now skip single residues returned, or anything weird.
                continue
            res1, res2 = residues
            res1_name, res1_chain, res1_serial = res1.split('-')
            res2_name, res2_chain, res2_serial = res2.split('-')
            comp1_res = None
            for rez in struct1.get_residues():
                rez_chain_name = rez.get_parent().get_id()
                rez_serial = rez.get_id()[1]
                rez_name = rez.get_resname()
                if rez_chain_name == res1_chain and rez_serial == int(res1_serial) and rez_name == res1_name:
                    comp1_res = rez
                    break

            comp2_res = None
            for rez in struct2.get_residues():
                rez_chain_name = rez.get_parent().get_id()
                rez_serial = rez.get_id()[1]
                rez_name = rez.get_resname()
                if rez_chain_name == res2_chain and rez_serial == int(res2_serial) and rez_name == res2_name:
                    comp2_res = rez
                    break

            # Get alpha carbon positions for each paired residue
            ca1 = next(atom for atom in comp1_res.get_atoms() if atom.name == 'CA')
            ca2 = next(atom for atom in comp2_res.get_atoms() if atom.name == 'CA')
            comp1_atom_list.append(ca1)
            comp2_atom_list.append(ca2)
        return comp1_atom_list, comp2_atom_list
  
                
