import tempfile
import os
import subprocess
from nanome.util import Logs

class FPocketClient:
    

    def get_pockets(self, comp):
        with tempfile.TemporaryDirectory() as tmpdir:
            comp_pdb = tempfile.NamedTemporaryFile(dir=tmpdir, suffix=".pdb")
            comp_pdb_path = comp_pdb.name
            comp_filename = comp_pdb.name.split('/')[-1].split('.pdb')[0]
            output_folder = os.path.join(tmpdir, f"{comp_filename}_out")
            comp.io.to_pdb(path=comp_pdb_path)
            completed_process = subprocess.run(
                ["fpocket", "-f",  comp_pdb_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            comp_pdb.close()
            Logs.message(f"fpocket output: {completed_process.returncode}")
            # Convert atoms
            pocket_folder = os.path.join(output_folder, 'pockets')
            # Every pocket returns a .pdb of atoms and a .pqr of vertices
            # Get atoms that are part of pocket
            output_pdb_files = [fi for fi in os.listdir(pocket_folder) if fi.endswith('.pdb')]
            list_of_pocket_serials = []
            for i in range(1, len(output_pdb_files)):
                pocket_file = os.path.join(pocket_folder, f'pocket{i}_atm.pdb')
                pocket_serials = []
                with open(pocket_file) as fd:
                    for line in fd:
                        if not line.startswith("ATOM"):
                            continue
                        # Collect serials for pocket atoms
                        row = line.split()
                        pocket_atom_serial = int(row[1])
                        pocket_serials.append(pocket_atom_serial)
                list_of_pocket_serials.append(pocket_serials)

            # Return list of atoms for each pocket
            pocket_atom_lists = []
            for pocket_serial_list in list_of_pocket_serials:
                pocket_atoms = []
                for pocket_serial in pocket_serial_list:
                    atom = next(a for a in comp.atoms if a.serial == pocket_serial)
                    pocket_atoms.append(atom)
                pocket_atom_lists.append(pocket_atoms)
            return pocket_atom_lists
