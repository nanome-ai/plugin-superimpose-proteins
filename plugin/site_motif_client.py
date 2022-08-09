import os
import shutil
import site_motif
import subprocess
import tempfile
from nanome.util import Logs
from nanome.api.structure import Complex


class SiteMotifClient:

    @staticmethod
    def run(binding_site_pdb: str, pocket_pdbs: list, align_output_file: str):
        """Runs SiteMotif and writes results to align_output_file."""
        with tempfile.TemporaryDirectory() as output_dir:
            sites_dir = tempfile.TemporaryDirectory(dir=output_dir)
            shutil.copy(binding_site_pdb, sites_dir.name)
            for pdb in pocket_pdbs:
                shutil.copy(pdb, sites_dir.name)
            pairs_filepath = site_motif.write_pairs(sites_dir.name, output_dir, reference_pdb=binding_site_pdb)
            pdb_size_filepath = site_motif.write_pdb_size(sites_dir.name, output_dir)
            script_path = f'{os.path.dirname(site_motif.__file__)}/pocket_matrix_mpi7.py'
            command = f"mpiexec -n 4 python {script_path} {sites_dir.name} {pairs_filepath} {pdb_size_filepath} {output_dir}"
            try:
                subprocess.run(command.split(), timeout=60)
            except subprocess.TimeoutExpired:
                Logs.warning("SiteMotif: Timeout. SiteMotif may not have finished.")

            align_output = f'{output_dir}/align_output.txt'
            sites_dir.cleanup()
            shutil.copy(align_output, align_output_file)
            # Uncomment for debugging
            # shutil.copy(align_output, 'align_output.txt')
            if not os.path.exists(align_output_file):
                raise Exception('No align output found')

    @staticmethod
    def find_match(comp_index: int, align_output_file: str):
        """Finds the pocket that best matches provided binding site."""
        # Find the alignment that contains the longest fragment match
        max_paired_res = -1
        residue_alignment = ''
        pdb_1 = ''
        pdb_2 = ''
        lines = []
        with open(align_output_file, 'r') as f:
            lines = f.readlines()
            if not lines:
                raise Exception("No lines found in align output")

        for line in lines:
            line_split = line.split('\t')
            if str(comp_index) not in line_split[0] and str(comp_index) not in line_split[1]:
                continue
            if line_split[0] == line_split[1]:
                continue
            aligned_residues = line_split[-1].strip()
            paired_residue_count = len(aligned_residues.split(' '))
            if paired_residue_count > max_paired_res:
                max_paired_res = paired_residue_count
                pdb_1 = line_split[0]
                pdb_2 = line_split[1]
                residue_alignment = aligned_residues
        Logs.message(f"Longest fragment match: {max_paired_res}")
        return pdb_1, pdb_2, residue_alignment

    @staticmethod
    def parse_residue_pairs(comp1: Complex, comp2: Complex, alignment: str):
        """Parses aligned residues from complexes.

        alignment format ex. ARG-A-334 GLN-A-155_VAL-A-343 LYS-A-198_ILE-A-364 ALA-A-151
        where comp1 is the first residue and comp2 is the second residue
        """
        alignment_pairings = alignment.strip().split(' ')
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
            comp2_res = None
            for rez in comp1.residues:
                rez_chain_name = rez.chain.name
                rez_serial = rez.serial
                rez_name = rez.name
                if rez_chain_name == res1_chain and rez_serial == int(res1_serial) and rez_name == res1_name:
                    comp1_res = rez
                    break

            for rez in comp2.residues:
                rez_chain_name = rez.chain.name
                rez_serial = rez.serial
                rez_name = rez.name
                if rez_chain_name == res2_chain and rez_serial == int(res2_serial) and rez_name == res2_name:
                    comp2_res = rez
                    break

            if not comp1_res:
                Logs.warning(f"Could not find {res1_name} {res1_chain} {res1_serial} on {comp1.full_name}")
                continue
            if not comp2_res:
                Logs.warning(f"Could not find {res2_name} {res2_chain} {res2_serial} on {comp2.full_name}")
                continue
            if comp1_res.name == comp2_res.name and len(list(comp1_res.atoms)) == len(list(comp2_res.atoms)):
                # If whole residues are the same, add all atoms
                comp1_res_atoms = sorted(comp1_res.atoms, key=lambda x: x.name)
                comp2_res_atoms = sorted(comp2_res.atoms, key=lambda x: x.name)
                comp1_atom_list.extend(comp1_res_atoms)
                comp2_atom_list.extend(comp2_res_atoms)
            else:
                # Find alpha carbons from residues and add to list
                ca1 = next(atom for atom in comp1_res.atoms if atom.name == 'CA')
                ca2 = next(atom for atom in comp2_res.atoms if atom.name == 'CA')
                comp1_atom_list.append(ca1)
                comp2_atom_list.append(ca2)
        return comp1_atom_list, comp2_atom_list
