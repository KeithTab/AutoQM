from Bio import PDB
from pathlib import Path
from functools import lru_cache, wraps
from multiprocessing import Pool, cpu_count
from typing import Dict, List, Generator, Tuple
import time
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from utils.standardize import standardize_pdb_file


def timing_decorator(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start
        if kwargs.get('verbose', False):
            print(f"  ⏱ {func.__name__} completed in {elapsed:.2f}s")
        return result
    return wrapper


def error_handler(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            return {
                'status': 'error',
                'error_message': str(e),
                'function': func.__name__
            }
    return wrapper


class PDBProcessor:
    __slots__ = ('parser', 'io')
    
    def __init__(self):
        self.parser = PDB.PDBParser(QUIET=True)
        self.io = PDB.PDBIO()
    
    class ProteinSelect(PDB.Select):
        __slots__ = ()
        
        STANDARD_RESIDUES = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        def accept_residue(self, residue):
            return residue.resname in self.STANDARD_RESIDUES
        
        def accept_atom(self, atom):
            return atom.element != 'H'
    
    class LigandSelect(PDB.Select):
        __slots__ = ()
        
        STANDARD_RESIDUES = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        def accept_residue(self, residue):
            return residue.resname not in self.STANDARD_RESIDUES and residue.resname != 'HOH'
    
    @error_handler
    @timing_decorator
    def add_hydrogens_to_ligand(self, ligand_pdb_path: Path, verbose: bool = False) -> bool:
        ligand_pdb_path = Path(ligand_pdb_path)
        
        with open(ligand_pdb_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM') or line.startswith('ATOM'):
                    residue_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    break
            else:
                residue_name = 'LIG'
                chain_id = 'A'
        
        mol = Chem.MolFromPDBFile(str(ligand_pdb_path), removeHs=False, sanitize=False)
        if mol is None:
            if verbose:
                print(f"  ⚠ RDKit failed to read {ligand_pdb_path.name}")
            return False
        
        rdDetermineBonds.DetermineBonds(mol, charge=0)
        
        mol_with_h = Chem.AddHs(mol, addCoords=True)
        
        conf = mol_with_h.GetConformer()
        heavy_atom_indices = []
        for atom in mol_with_h.GetAtoms():
            if atom.GetAtomicNum() != 1:
                heavy_atom_indices.append(atom.GetIdx())
        
        if AllChem.MMFFHasAllMoleculeParams(mol_with_h):
            props = AllChem.MMFFGetMoleculeProperties(mol_with_h)
            ff = AllChem.MMFFGetMoleculeForceField(mol_with_h, props)
            
            for idx in heavy_atom_indices:
                pos = conf.GetAtomPosition(idx)
                ff.AddFixedPoint(idx)
            
            ff.Minimize(maxIts=200)
        else:
            ff = AllChem.UFFGetMoleculeForceField(mol_with_h)
            
            for idx in heavy_atom_indices:
                ff.AddFixedPoint(idx)
            
            ff.Minimize(maxIts=200)
        
        temp_path = ligand_pdb_path.with_suffix('.tmp')
        writer = Chem.PDBWriter(str(temp_path))
        writer.write(mol_with_h)
        writer.close()
        
        fixed_lines = []
        with open(temp_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM') or line.startswith('ATOM'):
                    line = line[:17] + f'{residue_name:>3}' + ' ' + chain_id + line[22:]
                fixed_lines.append(line)
        
        with open(ligand_pdb_path, 'w') as f:
            f.writelines(fixed_lines)
        
        temp_path.unlink(missing_ok=True)
        
        if verbose:
            print(f"  ✓ Added hydrogens to {ligand_pdb_path.name} (residue: {residue_name}, chain: {chain_id})")
        return True
    
    def _count_atoms_generator(self, structure) -> Generator[Tuple[str, int], None, None]:
        total = protein = protein_no_h = ligand = 0
        
        STANDARD_RESIDUES = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    is_protein = residue.resname in STANDARD_RESIDUES
                    is_ligand = residue.resname not in STANDARD_RESIDUES and residue.resname != 'HOH'
                    
                    for atom in residue:
                        total += 1
                        if is_protein:
                            protein += 1
                            if atom.element != 'H':
                                protein_no_h += 1
                        elif is_ligand:
                            ligand += 1
        
        yield ('total_atoms', total)
        yield ('protein_atoms', protein)
        yield ('protein_atoms_no_h', protein_no_h)
        yield ('ligand_atoms', ligand)
    
    @error_handler
    @timing_decorator
    def process_pdb(self, input_pdb_path: Path, output_protein_path: Path, 
                    output_ligand_path: Path, verbose: bool = False) -> Dict:
        input_pdb_path = Path(input_pdb_path)
        output_protein_path = Path(output_protein_path)
        output_ligand_path = Path(output_ligand_path)
        
        output_protein_path.parent.mkdir(parents=True, exist_ok=True)
        output_ligand_path.parent.mkdir(parents=True, exist_ok=True)
        
        STANDARD_RESIDUES = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        protein_lines = []
        ligand_lines = []
        total_atoms = protein_atoms = protein_no_h = ligand_atoms = 0
        
        with open(input_pdb_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    resname = line[17:20].strip()
                    element = line[76:78].strip() if len(line) > 77 else line[12:14].strip()[0]
                    
                    total_atoms += 1
                    
                    if resname in STANDARD_RESIDUES:
                        if element != 'H':
                            protein_lines.append(line)
                            protein_no_h += 1
                        protein_atoms += 1
                    elif resname != 'HOH':
                        ligand_lines.append(line)
                        ligand_atoms += 1
        
        with open(output_protein_path, 'w') as f:
            f.writelines(protein_lines)
            if protein_lines:
                f.write('END\n')
        
        standardize_pdb_file(output_protein_path, verbose=verbose)
        
        with open(output_ligand_path, 'w') as f:
            f.writelines(ligand_lines)
            if ligand_lines:
                f.write('END\n')
        
        h_added = self.add_hydrogens_to_ligand(output_ligand_path, verbose=verbose)
        
        stats = {
            'input_file': str(input_pdb_path),
            'total_atoms': total_atoms,
            'protein_atoms': protein_atoms,
            'protein_atoms_no_h': protein_no_h,
            'ligand_atoms': ligand_atoms,
            'hydrogens_added': h_added,
            'output_protein_file': str(output_protein_path),
            'output_ligand_file': str(output_ligand_path),
            'status': 'success'
        }
        
        return stats
    
    @staticmethod
    def _process_single_file(args: Tuple) -> Dict:
        input_path, output_protein, output_ligand, verbose = args
        processor = PDBProcessor()
        
        try:
            stats = processor.process_pdb(input_path, output_protein, output_ligand, verbose)
            
            if verbose:
                base_name = input_path.stem
                print(f"Processed: {base_name}")
                print(f"  Protein atoms (no H): {stats['protein_atoms_no_h']}")
                print(f"  Ligand atoms: {stats['ligand_atoms']}")
                print(f"  Hydrogens added: {'Yes' if stats.get('hydrogens_added') else 'No'}")
            
            return stats
            
        except Exception as e:
            return {
                'input_file': str(input_path),
                'status': 'error',
                'error_message': str(e)
            }
    
    def _prepare_batch_args(self, input_files: List[Path], output_dir: Path, 
                           verbose: bool) -> Generator[Tuple, None, None]:
        for input_file in input_files:
            input_path = Path(input_file)
            base_name = input_path.stem
            
            file_output_dir = output_dir / base_name
            file_output_dir.mkdir(parents=True, exist_ok=True)
            
            output_protein = file_output_dir / f"{base_name}_protein.pdb"
            output_ligand = file_output_dir / f"{base_name}_ligand.pdb"
            
            yield (input_path, output_protein, output_ligand, verbose)
    
    @timing_decorator
    def process_batch(self, input_files: List[Path], output_dir: Path, 
                     verbose: bool = False, parallel: bool = True, 
                     n_processes: int = None) -> List[Dict]:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        args_list = list(self._prepare_batch_args(input_files, output_dir, verbose))
        
        if parallel and len(args_list) > 1:
            if n_processes is None:
                n_processes = min(cpu_count(), len(args_list))
            
            if verbose:
                print(f"Processing {len(args_list)} files using {n_processes} processes...")
            
            with Pool(processes=n_processes) as pool:
                results = pool.map(self._process_single_file, args_list)
        else:
            results = list(map(self._process_single_file, args_list))
        
        return results


def main():
    processor = PDBProcessor()


if __name__ == '__main__':
    main()
