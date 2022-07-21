import tempfile
from Bio.PDB import PDBParser

from nanome.api.structure import Complex, Molecule, Chain
from nanome.util import Logs


PDBOPTIONS = Complex.io.PDBSaveOptions()
PDBOPTIONS.write_bonds = True


def extract_binding_site(comp, binding_site_residues):
    """Copy comp, and remove all residues that are not part of the binding site."""
    new_comp = Complex()
    new_mol = Molecule()
    new_comp.add_molecule(new_mol)
    new_comp.name = f'{comp.name} binding site'
    new_comp.index = -1

    binding_site_residue_indices = [r.index for r in binding_site_residues]
    # Logs.debug(f'Binding site residues: {len(binding_site_residues)}')
    for ch in comp.chains:
        reses_on_chain = [res for res in ch.residues if res.index in binding_site_residue_indices]
        if reses_on_chain:
            new_ch = Chain()
            new_ch.name = ch.name
            new_ch.residues = reses_on_chain
            new_mol.add_chain(new_ch)
    # Logs.debug(f'New comp residues: {len(list(new_comp.residues))}')
    return new_comp


def clean_fpocket_pdbs(fpocket_pdbs, comp: Complex):
    """Add full residue data to pdb files."""
    Logs.debug(f"Cleaning {len(fpocket_pdbs)} fpocket pdbs")
    for i, pocket_pdb in enumerate(fpocket_pdbs):
        pocket_residues = set()
        with open(pocket_pdb) as f:
            for line in f:
                if line.startswith('ATOM'):
                    chain_name = line[21]
                    res_serial = int(line[22:26])
                    chain = next(ch for ch in comp.chains if ch.name == chain_name)
                    residue = next(rez for rez in chain.residues if rez.serial == res_serial)
                    pocket_residues.add(residue)
        pocket_comp = extract_binding_site(comp, pocket_residues)
        pocket_comp.io.to_pdb(path=pocket_pdb, options=PDBOPTIONS)
    return fpocket_pdbs


def convert_atoms_to_biopython(atom_list: list):
    """Converts atoms to biopython format."""
    parser = PDBParser(QUIET=True)
    comp = atom_list[0].complex
    comp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
    comp.io.to_pdb(comp_pdb.name)
    struct1 = parser.get_structure(comp.full_name, comp_pdb.name)
    # struct2 = parser.get_structure(comp2.full_name, comp2_pdb.name)
    bp_atom_list = []
    for atom in atom_list:
        res_serial = atom.residue.serial
        chain_name = atom.chain.name
        atom_name = atom.name  # should always be CA
        bp_atom = next(atm for atm in struct1.get_atoms() if all([
            atm.name == atom_name,
            atm.get_parent().id[1] == res_serial,
            atm.get_parent().get_parent().id == chain_name
        ]))
        bp_atom_list.append(bp_atom)
    return bp_atom_list

