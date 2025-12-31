from pathlib import Path
import subprocess
import shutil
import os
import tempfile
from multiprocessing import Pool, cpu_count
from functools import lru_cache, wraps
from typing import Dict, List, Optional, Generator
import time
from utils.charge import detect_ligand_charge

def timing_decorator(func):

    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start
        if kwargs.get('verbose', False):
            print(f'  ⏱ {func.__name__} completed in {elapsed:.2f}s')
        return result
    return wrapper

def error_handler(func):

    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            return {'status': 'error', 'error_message': str(e), 'function': func.__name__}
    return wrapper

class TopologyGenerator:
    __slots__ = ('_acpype_path', '_leap_template', '_leap_sys2_template', '_obabel_path')

    def __init__(self):
        self._acpype_path = shutil.which('acpype')
        if not self._acpype_path:
            raise RuntimeError('acpype not found in PATH')
        self._leap_template = None
        self._leap_sys2_template = None
        self._obabel_path = shutil.which('obabel')

    @property
    @lru_cache(maxsize=1)
    def leap_template(self) -> str:
        if self._leap_template is None:
            self._leap_template = """source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams frcmod.ionsjc_tip3p
set default PBradii mbondi3

rec=loadpdb {protein_pdb}
loadamberparams {frcmod}
lig=loadmol2 {mol2}

gascomplex= combine {{rec lig}}

savepdb gascomplex gas.complex.pdb
saveamberparm gascomplex gas.complex.parm7 gas.complex.rst7
saveamberparm rec gas.receptor.parm7 gas.receptor.rst7
saveamberparm lig gas.ligand.parm7 gas.ligand.rst7

solvcomplex= combine {{rec lig}}
solvatebox solvcomplex TIP3PBOX 12.0
addions solvcomplex Cl- 0
addions solvcomplex Na+ 0

savepdb solvcomplex wet.complex.pdb
charge solvcomplex
check solvcomplex
saveamberparm solvcomplex wet.complex.prmtop wet.complex.rst7
quit
"""
        return self._leap_template

    @property
    @lru_cache(maxsize=1)
    def leap_sys2_template(self) -> str:
        if self._leap_sys2_template is None:
            self._leap_sys2_template = """source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams frcmod.ionsjc_tip3p
set default PBradii mbondi3

rec=loadpdb {protein_pdb}
loadamberparams {frcmod}

savepdb rec gas.sys2.pdb
saveamberparm rec gas.sys2.parm7 gas.sys2.rst7

solvcomplex= combine {{rec}}
solvatebox solvcomplex TIP3PBOX 12.0
addions solvcomplex Cl- 0
addions solvcomplex Na+ 0

savepdb solvcomplex wet.sys2.pdb
charge solvcomplex
check solvcomplex
saveamberparm solvcomplex wet.sys2.prmtop wet.sys2.rst7
quit
"""
        return self._leap_sys2_template

    def _collect_output_files(self, working_dir: Path) -> Dict[str, List[str]]:
        return {'topology_files': [str(f) for f in working_dir.rglob('*.itp')], 'gro_files': [str(f) for f in working_dir.rglob('*.gro')], 'top_files': [str(f) for f in working_dir.rglob('*.top')], 'frcmod_files': [str(f) for f in working_dir.rglob('*.frcmod')], 'mol2_files': [str(f) for f in working_dir.rglob('*.mol2') if 'gaff2' in f.name]}

    def _add_hydrogens_to_ligand(self, ligand_pdb_path: Path, verbose: bool=False) -> bool:
        if not self._obabel_path:
            if verbose:
                print('  ⚠ obabel not found; skipping ligand hydrogen addition')
            return False
        ligand_pdb_path = Path(ligand_pdb_path)
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp_file:
            tmp_path = Path(tmp_file.name)
        try:
            result = subprocess.run([self._obabel_path, str(ligand_pdb_path), '-O', str(tmp_path), '-h'], capture_output=True, text=True, timeout=120)
            if result.returncode == 0 and tmp_path.exists():
                tmp_path.replace(ligand_pdb_path)
                if verbose:
                    print(f'  ✓ Added hydrogens to {ligand_pdb_path.name} via obabel')
                return True
            if verbose:
                print(f'  ⚠ obabel failed to add hydrogens for {ligand_pdb_path.name}')
            tmp_path.unlink(missing_ok=True)
            return False
        except subprocess.TimeoutExpired:
            if verbose:
                print(f'  ⚠ obabel hydrogenation timed out for {ligand_pdb_path.name}')
            tmp_path.unlink(missing_ok=True)
            return False
        except Exception as exc:
            if verbose:
                print(f'  ⚠ obabel error for {ligand_pdb_path.name}: {exc}')
            tmp_path.unlink(missing_ok=True)
            return False

    @error_handler
    @timing_decorator
    def generate_topology(self, ligand_pdb_path: Path, output_dir: Optional[Path]=None, verbose: bool=False) -> Dict:
        ligand_pdb_path = Path(ligand_pdb_path)
        if not ligand_pdb_path.exists():
            raise FileNotFoundError(f'Ligand PDB file not found: {ligand_pdb_path}')
        output_dir = Path(output_dir) if output_dir else ligand_pdb_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        ligand_name = ligand_pdb_path.stem
        working_dir = output_dir / f'{ligand_name}_acpype'
        hydrogens_added = self._add_hydrogens_to_ligand(ligand_pdb_path, verbose)
        charge = detect_ligand_charge(ligand_pdb_path, verbose=verbose)
        if working_dir.exists():
            frcmod_files = list(working_dir.rglob('*.frcmod'))
            mol2_files = [f for f in working_dir.rglob('*.mol2') if 'gaff2' in f.name]
            if frcmod_files and mol2_files:
                if verbose:
                    print(f'Topology exists for {ligand_name}, regenerating leap.in...')
                self._generate_leap_input(output_dir, working_dir, frcmod_files, mol2_files, verbose)
                file_dict = self._collect_output_files(working_dir)
                return {'status': 'success', 'ligand_name': ligand_name, 'input_file': str(ligand_pdb_path), 'output_dir': str(working_dir), 'hydrogens_added': hydrogens_added, **file_dict, 'skipped': True}
        working_dir.mkdir(parents=True, exist_ok=True)
        cmd = ['acpype', '-i', str(ligand_pdb_path), '-n', str(charge), '-o', 'gmx']
        try:
            if verbose:
                print(f'Generating topology for {ligand_name}...')
            result = subprocess.run(cmd, cwd=str(working_dir), capture_output=True, text=True, check=True)
            if verbose and result.stdout:
                print(result.stdout)
            file_dict = self._collect_output_files(working_dir)
            self._generate_leap_input(ligand_pdb_path.parent, working_dir, [Path(f) for f in file_dict['frcmod_files']], [Path(f) for f in file_dict['mol2_files']], verbose)
            return {'status': 'success', 'ligand_name': ligand_name, 'input_file': str(ligand_pdb_path), 'output_dir': str(working_dir), 'hydrogens_added': hydrogens_added, **file_dict}
        except subprocess.CalledProcessError as e:
            return {'status': 'error', 'ligand_name': ligand_name, 'input_file': str(ligand_pdb_path), 'hydrogens_added': hydrogens_added, 'error_message': e.stderr if e.stderr else str(e)}
        except Exception as e:
            return {'status': 'error', 'ligand_name': ligand_name, 'input_file': str(ligand_pdb_path), 'hydrogens_added': hydrogens_added, 'error_message': str(e)}

    @timing_decorator
    def _generate_leap_input(self, parent_dir: Path, acpype_dir: Path, frcmod_files: List[Path], mol2_files: List[Path], verbose: bool=False) -> bool:
        parent_dir = Path(parent_dir)
        acpype_dir = Path(acpype_dir)
        protein_pdb = next(parent_dir.glob('*_protein.pdb'), None)
        if not protein_pdb:
            if verbose:
                print(f'Warning: No protein PDB file found in {parent_dir}')
            return False
        if not frcmod_files or not mol2_files:
            if verbose:
                print(f'Warning: Missing frcmod or mol2 files in {acpype_dir}')
            return False
        protein_rel = os.path.relpath(protein_pdb, parent_dir)
        frcmod_rel = os.path.relpath(frcmod_files[0], parent_dir)
        mol2_rel = os.path.relpath(mol2_files[0], parent_dir)
        leap_content = self.leap_template.format(protein_pdb=protein_rel, frcmod=frcmod_rel, mol2=mol2_rel)
        leap_file = parent_dir / 'leap.in'
        leap_file.write_text(leap_content)
        leap_sys2_content = self.leap_sys2_template.format(protein_pdb=protein_rel, frcmod=frcmod_rel)
        leap_sys2_file = parent_dir / 'leap_sys2.in'
        leap_sys2_file.write_text(leap_sys2_content)
        if verbose:
            print(f'Generated leap.in at {leap_file}')
            print(f'Generated leap_sys2.in at {leap_sys2_file}')
        return True

    def process_batch(self, ligand_pdb_files, output_dir=None, verbose=False):
        results = []
        for ligand_file in ligand_pdb_files:
            ligand_path = Path(ligand_file)
            if verbose:
                print(f'Processing {ligand_path.name}...')
            result = self.generate_topology(ligand_path, output_dir, verbose)
            results.append(result)
            if verbose:
                if result['status'] == 'success':
                    print(f"  Success: topology files in {result['output_dir']}")
                else:
                    print(f"  Error: {result['error_message']}")
        return results

    def _generate_topology_wrapper(self, args):
        ligand_pdb, verbose = args
        ligand_name = ligand_pdb.stem.replace('_ligand', '')
        output_subdir = ligand_pdb.parent
        return self.generate_topology(ligand_pdb, output_subdir, verbose=verbose)

    def _ligand_files_generator(self, output_dir: Path) -> Generator[Path, None, None]:
        yield from output_dir.rglob('*_ligand.pdb')

    @timing_decorator
    def process_from_directory(self, output_dir: Path, verbose: bool=False, parallel: bool=False, n_processes: Optional[int]=None) -> List[Dict]:
        output_dir = Path(output_dir)
        ligand_files = list(self._ligand_files_generator(output_dir))
        if not ligand_files:
            return []
        if parallel and len(ligand_files) > 1:
            if n_processes is None:
                n_processes = min(cpu_count(), len(ligand_files))
            if verbose:
                print(f'Processing {len(ligand_files)} ligands using {n_processes} processes...')
            args_list = [(ligand_pdb, verbose) for ligand_pdb in ligand_files]
            with Pool(processes=n_processes) as pool:
                results = pool.map(self._generate_topology_wrapper, args_list)
            return results
        else:
            return list(map(lambda ligand_pdb: self.generate_topology(ligand_pdb, ligand_pdb.parent, verbose=verbose), ligand_files))

def main():
    generator = TopologyGenerator()
if __name__ == '__main__':
    main()