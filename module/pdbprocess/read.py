"""
Optimized PDB Processing Module
"""

from Bio import PDB
from pathlib import Path
import subprocess
import tempfile
from functools import lru_cache, wraps
from multiprocessing import Pool, cpu_count
from typing import Dict, List, Generator, Tuple
import time


def timing_decorator(func):
    """Decorator to measure function execution time"""
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
    """Decorator for consistent error handling"""
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
    """Optimized PDB processor with parallel processing"""
    
    __slots__ = ('parser', 'io', '_obabel_available')
    
    def __init__(self):
        self.parser = PDB.PDBParser(QUIET=True)
        self.io = PDB.PDBIO()
        self._obabel_available = None
    
    @property
    @lru_cache(maxsize=1)
    def obabel_available(self) -> bool:
        """Check if obabel is available (cached)"""
        if self._obabel_available is None:
            try:
                result = subprocess.run(
                    ['obabel', '--version'],
                    capture_output=True,
                    timeout=5
                )
                self._obabel_available = result.returncode == 0
            except (FileNotFoundError, subprocess.TimeoutExpired):
                self._obabel_available = False
        return self._obabel_available
    
    class ProteinSelect(PDB.Select):
        """Optimized protein selector"""
        __slots__ = ()
        
        def accept_residue(self, residue):
            return residue.id[0] == ' '
        
        def accept_atom(self, atom):
            return atom.element != 'H'
    
    class LigandSelect(PDB.Select):
        """Optimized ligand selector"""
        __slots__ = ()
        
        def accept_residue(self, residue):
            return residue.id[0] != ' ' and residue.resname != 'HOH'
    
    @error_handler
    @timing_decorator
    def add_hydrogens_to_ligand(self, ligand_pdb_path: Path, verbose: bool = False) -> bool:
        """
        Add hydrogens to ligand using obabel with optimized I/O
        
        Args:
            ligand_pdb_path: Path to ligand PDB file
            verbose: Enable verbose output
            
        Returns:
            True if successful, False otherwise
        """
        if not self.obabel_available:
            return False
        
        ligand_pdb_path = Path(ligand_pdb_path)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_file:
            tmp_path = Path(tmp_file.name)
        
        try:
            result = subprocess.run(
                ['obabel', str(ligand_pdb_path), '-O', str(tmp_path), '-h'],
                capture_output=True,
                text=True,
                timeout=60
            )
            
            if result.returncode == 0 and tmp_path.exists():
                tmp_path.replace(ligand_pdb_path)
                return True
            
            tmp_path.unlink(missing_ok=True)
            return False
                
        except (subprocess.TimeoutExpired, Exception):
            tmp_path.unlink(missing_ok=True)
            return False
    
    def _count_atoms_generator(self, structure) -> Generator[Tuple[str, int], None, None]:
        """Generator to count atoms efficiently"""
        total = protein = protein_no_h = ligand = 0
        
        for model in structure:
            for chain in model:
                for residue in chain:
                    is_protein = residue.id[0] == ' '
                    is_ligand = residue.id[0] != ' ' and residue.resname != 'HOH'
                    
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
        """Process PDB file with optimized I/O operations"""
        input_pdb_path = Path(input_pdb_path)
        output_protein_path = Path(output_protein_path)
        output_ligand_path = Path(output_ligand_path)
        
        # Create parent directories efficiently
        output_protein_path.parent.mkdir(parents=True, exist_ok=True)
        output_ligand_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Parse structure once
        structure = self.parser.get_structure('complex', str(input_pdb_path))
        
        # Initialize stats using generator
        stats = {'input_file': str(input_pdb_path)}
        stats.update(dict(self._count_atoms_generator(structure)))
        
        # Save protein and ligand
        self.io.set_structure(structure)
        self.io.save(str(output_protein_path), self.ProteinSelect())
        
        self.io.set_structure(structure)
        self.io.save(str(output_ligand_path), self.LigandSelect())
        
        # Add hydrogens if available
        h_added = self.add_hydrogens_to_ligand(output_ligand_path, verbose=verbose)
        
        stats.update({
            'hydrogens_added': h_added,
            'output_protein_file': str(output_protein_path),
            'output_ligand_file': str(output_ligand_path),
            'status': 'success'
        })
        
        return stats
    
    @staticmethod
    def _process_single_file(args: Tuple) -> Dict:
        """Static method for parallel processing"""
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
        """Generator to prepare arguments for batch processing"""
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
        """Process multiple PDB files with optional parallel processing"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Prepare arguments using generator
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
