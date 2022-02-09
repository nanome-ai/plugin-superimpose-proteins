import nanome
from nanome.api.ui import Menu
from nanome.util import Logs, async_callback
from Bio.PDB import Superimposer, PDBParser
from Bio import SeqIO, pairwise2
import tempfile

class RMSDV2(nanome.AsyncPluginInstance):

    def start(self):
        pass

    @async_callback
    async def on_run(self):
        Logs.message("Retrieving complexes")
        complexes = await self.request_complex_list()
        complexes = await self.request_complexes([comp.index for comp in complexes])

        Logs.message("Parsing Complexes to BioPython Structures.")
        parser = PDBParser(QUIET=True)
        fixed_comp = complexes[0]
        moving_comp = complexes[1]
        fixed_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        moving_pdb = tempfile.NamedTemporaryFile(suffix=".pdb")
        fixed_comp.io.to_pdb(fixed_pdb.name)
        moving_comp.io.to_pdb(moving_pdb.name)
        # fixed_struct = parser.get_structure(fixed_comp.full_name, fixed_pdb.name)
        # moving_struct = parser.get_structure(moving_comp.full_name, moving_pdb.name)
        # fixed_chain = next(fixed_struct.get_chains())
        # moving_chain = next(moving_struct.get_chains())
        Logs.message("Selecting Chains")
        fixed_chain_name = 'A'
        moving_chain_name = 'A'
        for record in SeqIO.parse(fixed_pdb.name, 'pdb-atom'):
            if record.annotations.get('chain') == fixed_chain_name:
                fixed_seq = record.seq
                break
        
        for record in SeqIO.parse(moving_pdb.name, 'pdb-atom'):
            if record.annotations.get('chain') == moving_chain_name:
                moving_seq = record.seq
                break
        
        Logs.message("Aligning chains")
        alignment = pairwise2.align.globalxx(fixed_seq, moving_seq)
        # self.align_chains(fixed_chain, moving_chain)
        Logs.message("SuperImposing Structures.")
        superimposer = Superimposer()
        # superimposer.set_atoms(list(fixed_chain.get_atoms()), list(moving_chain.get_atoms()))
        Logs.message(superimposer.rms)
        self.update_menu(self.menu)
    
    def align_chains(self, fixed_chain, moving_chain):
        pass


def main():
    plugin = nanome.Plugin('RMSD V2', 'Superimpose two structures', 'alignment', False)
    plugin.set_plugin_class(RMSDV2)
    plugin.run()


if __name__ == '__main__':
    main()
