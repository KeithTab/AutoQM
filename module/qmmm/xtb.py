from pathlib import Path
import subprocess
import shutil
from collections import Counter
from multiprocessing import Pool, cpu_count
from functools import lru_cache, wraps
from typing import Dict, List, Tuple, Optional, Set, Generator
import time
from utils.charge import *

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

class XTBProcessor:
    __slots__ = ('amino_acid_charges', '_standard_amino_acids', '_obabel_path')

    def __init__(self):
        self.amino_acid_charges = {'ALA': 0, 'ARG': +1, 'ASN': 0, 'ASP': -1, 'CYS': 0, 'GLN': 0, 'GLU': -1, 'GLY': 0, 'HIS': 0, 'ILE': 0, 'LEU': 0, 'LYS': +1, 'MET': 0, 'PHE': 0, 'PRO': 0, 'SER': 0, 'THR': 0, 'TRP': 0, 'TYR': 0, 'VAL': 0, 'HIE': 0, 'HID': 0, 'HIP': +1, 'CYX': 0}
        self._standard_amino_acids = frozenset(self.amino_acid_charges.keys())
        self._obabel_path = shutil.which('obabel')

    def run_tleap(self, leap_input_path):
        leap_input_path = Path(leap_input_path)
        work_dir = leap_input_path.parent

        def _run_single_leap(leap_file):
            try:
                return subprocess.run(['tleap', '-f', leap_file.name], cwd=str(work_dir), capture_output=True, text=True, timeout=300)
            except subprocess.TimeoutExpired:
                raise RuntimeError(f'tleap command timed out for {leap_file.name}')
            except Exception as err:
                raise RuntimeError(f'tleap failed for {leap_file.name}: {err}')
        try:
            _run_single_leap(leap_input_path)
            gas_complex_pdb = work_dir / 'gas.complex.pdb'
            if not gas_complex_pdb.exists():
                return {'status': 'error', 'error_message': 'gas.complex.pdb not generated', 'output_dir': str(work_dir)}
            leap_sys2_path = leap_input_path.with_name('leap_sys2.in')
            leap_sys2_status = 'skipped'
            if leap_sys2_path.exists():
                try:
                    _run_single_leap(leap_sys2_path)
                    leap_sys2_status = 'success'
                except RuntimeError as err:
                    leap_sys2_status = f'error: {err}'
            return {'status': 'success', 'gas_complex_pdb': str(gas_complex_pdb), 'output_dir': str(work_dir), 'leap_sys2_status': leap_sys2_status}
        except RuntimeError as err:
            return {'status': 'error', 'error_message': str(err), 'output_dir': str(work_dir)}

    def _parse_pdb_atoms(self, pdb_path: Path) -> Tuple[List[Tuple], List[Dict], Dict]:
        ligand_coords = []
        protein_atoms = []
        residue_info = {}
        with open(pdb_path, 'r') as f:
            for line in f:
                if not line.startswith('ATOM'):
                    continue
                resname = line[17:20].strip()
                resnum = int(line[22:26].strip())
                coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                if resname not in self._standard_amino_acids:
                    ligand_coords.append(coords)
                else:
                    protein_atoms.append({'resnum': resnum, 'resname': resname, 'coords': coords})
                    residue_info.setdefault(resnum, resname)
        return (ligand_coords, protein_atoms, residue_info)

    @timing_decorator
    def get_nearby_residues(self, pdb_path: Path, ligand_resname: str='LIG', distance_cutoff: float=2.5, verbose: bool=False) -> List[Tuple[int, str]]:
        pdb_path = Path(pdb_path)
        try:
            ligand_coords, protein_atoms, residue_info = self._parse_pdb_atoms(pdb_path)
            if not ligand_coords or not protein_atoms:
                return []
            nearby_residues: Set[int] = set()
            cutoff_sq = distance_cutoff ** 2
            for lx, ly, lz in ligand_coords:
                for atom in protein_atoms:
                    px, py, pz = atom['coords']
                    dist_sq = (lx - px) ** 2 + (ly - py) ** 2 + (lz - pz) ** 2
                    if dist_sq <= cutoff_sq:
                        nearby_residues.add(atom['resnum'])
            return [(resnum, residue_info[resnum]) for resnum in sorted(nearby_residues)]
        except Exception as e:
            return []

    def get_ligand_resnum(self, pdb_path):
        pdb_path = Path(pdb_path)
        standard_amino_acids = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HIE', 'HID', 'HIP', 'CYX'}
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        resname = line[17:20].strip()
                        if resname not in standard_amino_acids:
                            resnum = int(line[22:26].strip())
                            return resnum
            return None
        except Exception as e:
            return None

    def find_common_residues(self, residue_lists):
        if not residue_lists:
            return []
        residue_counter = Counter()
        residue_names = {}
        for residues_with_names in residue_lists:
            for resnum, resname in residues_with_names:
                residue_counter[resnum] += 1
                residue_names[resnum] = resname
        num_lists = len(residue_lists)
        common_residues = [(res, residue_names[res]) for res, count in residue_counter.items() if count == num_lists]
        common_residues.sort(key=lambda x: x[0])
        return common_residues

    def calculate_qm_charge(self, residues_with_names, work_dir=None):
        ligand_charge = 0
        if work_dir:
            work_dir = Path(work_dir)
            ligand_pdb = next(work_dir.glob('*_ligand.pdb'), None)
            if ligand_pdb:
                ligand_charge = detect_ligand_charge(ligand_pdb, verbose=False)
        total_charge = ligand_charge
        for resnum, resname in residues_with_names:
            charge = self.amino_acid_charges.get(resname, 0)
            total_charge += charge
        return total_charge

    def format_qmmask(self, residue_numbers_with_names: List[Tuple[int, str]]) -> str:
        if not residue_numbers_with_names:
            return ''
        residue_numbers = [resnum for resnum, _ in residue_numbers_with_names]

        def generate_ranges():
            start = end = residue_numbers[0]
            for num in residue_numbers[1:]:
                if num == end + 1:
                    end = num
                else:
                    yield (f'{start}-{end}' if start != end else str(start))
                    start = end = num
            yield (f'{start}-{end}' if start != end else str(start))
        return ','.join(generate_ranges())

    def run_sander_min(self, work_dir):
        work_dir = Path(work_dir)
        required_files = ['min.mdin', 'wet.complex.prmtop', 'wet.complex.rst7']
        for req_file in required_files:
            if not (work_dir / req_file).exists():
                return {'status': 'error', 'error_message': f'Required file {req_file} not found', 'output_dir': str(work_dir)}
        mdin_file = work_dir / 'min.mdin'
        spin_attempts = self._spin_attempts(mdin_file)
        last_error = 'Sander minimization failed'
        last_stdout = None
        last_stderr = None
        for idx, spin_value in enumerate(spin_attempts):
            if idx > 0:
                self._set_spin_value(mdin_file, spin_value)
            try:
                result = subprocess.run(['sander', '-O', '-i', 'min.mdin', '-o', '01.min.mdout', '-p', 'wet.complex.prmtop', '-c', 'wet.complex.rst7', '-ref', 'wet.complex.rst7', '-x', '01.min.trj', '-inf', '01.min.mdinfo', '-r', '01.min.rst7'], cwd=str(work_dir), capture_output=True, text=True, timeout=600)
            except subprocess.TimeoutExpired:
                last_error = 'Sander minimization timed out'
                if idx == len(spin_attempts) - 1:
                    return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir)}
                continue
            except Exception as e:
                last_error = str(e)
                if idx == len(spin_attempts) - 1:
                    return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir)}
                continue
            mdout_file = work_dir / '01.min.mdout'
            rst7_file = work_dir / '01.min.rst7'
            if result.returncode == 0 and mdout_file.exists() and rst7_file.exists():
                return {'status': 'success', 'mdout': str(mdout_file), 'rst7': str(rst7_file), 'output_dir': str(work_dir), 'spin_used': spin_value}
            last_error = 'Sander minimization output files not generated'
            last_stdout = result.stdout
            last_stderr = result.stderr
        return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir), 'stdout': last_stdout, 'stderr': last_stderr}

    def run_sander_sp(self, work_dir):
        work_dir = Path(work_dir)
        required_files = ['sp.mdin', 'wet.complex.prmtop', '01.min.rst7']
        for req_file in required_files:
            if not (work_dir / req_file).exists():
                return {'status': 'error', 'error_message': f'Required file {req_file} not found', 'output_dir': str(work_dir)}
        mdin_file = work_dir / 'sp.mdin'
        spin_attempts = self._spin_attempts(mdin_file)
        last_error = 'Sander SP calculation failed'
        last_stdout = None
        last_stderr = None
        for idx, spin_value in enumerate(spin_attempts):
            if idx > 0:
                self._set_spin_value(mdin_file, spin_value)
            try:
                result = subprocess.run(['sander', '-O', '-i', 'sp.mdin', '-o', 'sp.mdout', '-p', 'wet.complex.prmtop', '-c', '01.min.rst7', '-ref', '01.min.rst7', '-y', '01.min.trj', '-inf', 'sp.mdinfo', '-r', 'sp.rst7', '-x', 'sp.trj'], cwd=str(work_dir), capture_output=True, text=True, timeout=600)
            except subprocess.TimeoutExpired:
                last_error = 'Sander SP calculation timed out'
                if idx == len(spin_attempts) - 1:
                    return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir)}
                continue
            except Exception as e:
                last_error = str(e)
                if idx == len(spin_attempts) - 1:
                    return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir)}
                continue
            mdout_file = work_dir / 'sp.mdout'
            rst7_file = work_dir / 'sp.rst7'
            if result.returncode == 0 and mdout_file.exists() and rst7_file.exists():
                return {'status': 'success', 'mdout': str(mdout_file), 'rst7': str(rst7_file), 'output_dir': str(work_dir), 'spin_used': spin_value}
            last_error = 'Sander SP output files not generated'
            last_stdout = result.stdout
            last_stderr = result.stderr
        return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir), 'stdout': last_stdout, 'stderr': last_stderr}

    def run_sander_min_sys2(self, work_dir):
        work_dir = Path(work_dir)
        required_files = ['min_sys2.mdin', 'wet.sys2.prmtop', 'wet.sys2.rst7']
        for req_file in required_files:
            if not (work_dir / req_file).exists():
                return {'status': 'error', 'error_message': f'Required file {req_file} not found', 'output_dir': str(work_dir)}
        mdin_file = work_dir / 'min_sys2.mdin'
        spin_attempts = self._spin_attempts(mdin_file)
        last_error = 'Sander SYS2 minimization failed'
        last_stdout = None
        last_stderr = None
        for idx, spin_value in enumerate(spin_attempts):
            if idx > 0:
                self._set_spin_value(mdin_file, spin_value)
            try:
                result = subprocess.run(['sander', '-O', '-i', 'min_sys2.mdin', '-o', '02.min.mdout', '-p', 'wet.sys2.prmtop', '-c', 'wet.sys2.rst7', '-ref', 'wet.sys2.rst7', '-inf', '02.min.mdinfo', '-r', '02.min.rst7', '-x', '02.min.trj'], cwd=str(work_dir), capture_output=True, text=True, timeout=600)
            except subprocess.TimeoutExpired:
                last_error = 'Sander SYS2 minimization timed out'
                if idx == len(spin_attempts) - 1:
                    return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir)}
                continue
            except Exception as e:
                last_error = str(e)
                if idx == len(spin_attempts) - 1:
                    return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir)}
                continue
            mdout_file = work_dir / '02.min.mdout'
            rst7_file = work_dir / '02.min.rst7'
            if result.returncode == 0 and mdout_file.exists() and rst7_file.exists():
                return {'status': 'success', 'mdout': str(mdout_file), 'rst7': str(rst7_file), 'output_dir': str(work_dir), 'spin_used': spin_value}
            last_error = 'Sander SYS2 minimization output files not generated'
            last_stdout = result.stdout
            last_stderr = result.stderr
        return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir), 'stdout': last_stdout, 'stderr': last_stderr}

    def run_sander_sp_sys2(self, work_dir):
        work_dir = Path(work_dir)
        required_files = ['sp_sys2.mdin', 'wet.sys2.prmtop', 'wet.sys2.rst7', '02.min.rst7', '02.min.trj']
        for req_file in required_files:
            if not (work_dir / req_file).exists():
                return {'status': 'error', 'error_message': f'Required file {req_file} not found', 'output_dir': str(work_dir)}
        mdin_file = work_dir / 'sp_sys2.mdin'
        spin_attempts = self._spin_attempts(mdin_file)
        last_error = 'Sander SYS2 SP calculation failed'
        last_stdout = None
        last_stderr = None
        for idx, spin_value in enumerate(spin_attempts):
            if idx > 0:
                self._set_spin_value(mdin_file, spin_value)
            try:
                result = subprocess.run(['sander', '-O', '-i', 'sp_sys2.mdin', '-o', 'sp_sys2.mdout', '-p', 'wet.sys2.prmtop', '-c', 'wet.sys2.rst7', '-ref', '02.min.rst7', '-inf', 'sp_sys2.mdinfo', '-r', 'sp_sys2.rst7', '-x', 'sp_sys2.trj', '-y', '02.min.trj'], cwd=str(work_dir), capture_output=True, text=True, timeout=600)
            except subprocess.TimeoutExpired:
                last_error = 'Sander SYS2 SP calculation timed out'
                if idx == len(spin_attempts) - 1:
                    return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir)}
                continue
            except Exception as e:
                last_error = str(e)
                if idx == len(spin_attempts) - 1:
                    return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir)}
                continue
            mdout_file = work_dir / 'sp_sys2.mdout'
            rst7_file = work_dir / 'sp_sys2.rst7'
            if result.returncode == 0 and mdout_file.exists() and rst7_file.exists():
                return {'status': 'success', 'mdout': str(mdout_file), 'rst7': str(rst7_file), 'output_dir': str(work_dir), 'spin_used': spin_value}
            last_error = 'Sander SYS2 SP output files not generated'
            last_stdout = result.stdout
            last_stderr = result.stderr
        return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir), 'stdout': last_stdout, 'stderr': last_stderr}

    def pdb_to_xyz(self, pdb_path, xyz_path):
        pdb_path = Path(pdb_path)
        xyz_path = Path(xyz_path)
        atoms = []
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        element = line[76:78].strip()
                        if not element:
                            atom_name = line[12:16].strip()
                            element = ''.join([c for c in atom_name if c.isalpha()])[:2]
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        atoms.append((element, x, y, z))
            with open(xyz_path, 'w') as f:
                f.write(f'{len(atoms)}\n')
                f.write('Converted from PDB\n')
                for element, x, y, z in atoms:
                    f.write(f'{element:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n')
            return True
        except Exception as e:
            return False

    def _convert_ligand_pdb_to_xyz(self, ligand_pdb: Path, xyz_path: Path, verbose: bool=False) -> Dict:
        ligand_pdb = Path(ligand_pdb)
        xyz_path = Path(xyz_path)
        if not self._obabel_path:
            return {'status': 'error', 'error_message': 'obabel executable not found in PATH'}
        try:
            result = subprocess.run([self._obabel_path, str(ligand_pdb), '-O', str(xyz_path)], capture_output=True, text=True, timeout=120)
        except subprocess.TimeoutExpired:
            return {'status': 'error', 'error_message': 'obabel conversion timed out'}
        except Exception as exc:
            return {'status': 'error', 'error_message': f'obabel conversion failed: {exc}'}
        if result.returncode == 0 and xyz_path.exists():
            if verbose:
                print(f'  ✓ Generated {xyz_path.name} from {ligand_pdb.name} via obabel')
            return {'status': 'success', 'xyz_file': str(xyz_path)}
        return {'status': 'error', 'error_message': 'obabel conversion returned non-zero exit status', 'stdout': result.stdout, 'stderr': result.stderr}

    def run_xtb_ligand(self, work_dir, verbose: bool=False):
        work_dir = Path(work_dir)
        ligand_pdb = next(iter(sorted(work_dir.glob('*_ligand.pdb'))), None)
        if not ligand_pdb or not ligand_pdb.exists():
            return {'status': 'error', 'error_message': 'Ligand PDB file not found for XTB ligand calculation'}
        qmmm_xyz = work_dir / 'qmmm.xyz'
        conversion = self._convert_ligand_pdb_to_xyz(ligand_pdb, qmmm_xyz, verbose=verbose)
        if conversion['status'] != 'success':
            return conversion
        ligand_charge = detect_ligand_charge(ligand_pdb, verbose=False)
        uhf_attempts = [2, 1] if ligand_charge >= 0 else [1, 2]
        last_error = 'XTB ligand calculation failed'
        last_stdout = None
        last_stderr = None
        for uhf in uhf_attempts:
            try:
                result = subprocess.run(['xtb', 'qmmm.xyz', '-chrg', str(ligand_charge), '-uhf', str(uhf)], cwd=str(work_dir), capture_output=True, text=True, timeout=600)
            except subprocess.TimeoutExpired:
                last_error = 'XTB ligand calculation timed out'
                last_stdout = None
                last_stderr = None
                continue
            except Exception as exc:
                last_error = str(exc)
                last_stdout = None
                last_stderr = None
                continue
            ligand_out = work_dir / 'qm_ligand.out'
            ligand_out.write_text(result.stdout)
            if result.returncode == 0:
                if verbose:
                    print(f'  ✓ XTB ligand calculation succeeded with uhf={uhf}')
                return {'status': 'success', 'output_file': str(ligand_out), 'uhf_used': uhf}
            last_error = 'XTB ligand command returned non-zero exit status'
            last_stdout = result.stdout
            last_stderr = result.stderr
        return {'status': 'error', 'error_message': last_error, 'stdout': last_stdout, 'stderr': last_stderr}

    def run_xtb_qm(self, work_dir, qm_charge=0):
        work_dir = Path(work_dir)
        qmmm_pdb = work_dir / 'qmmm_region.pdb'
        if not qmmm_pdb.exists():
            return {'status': 'error', 'error_message': 'qmmm_region.pdb not found', 'output_dir': str(work_dir)}
        qmmm_xyz = work_dir / 'qmmm_region.xyz'
        if not self.pdb_to_xyz(qmmm_pdb, qmmm_xyz):
            return {'status': 'error', 'error_message': 'Failed to convert PDB to XYZ', 'output_dir': str(work_dir)}
        for uhf_value in [1, 2]:
            try:
                result = subprocess.run(['xtb', 'qmmm_region.xyz', '-chrg', str(qm_charge), '-uhf', str(uhf_value)], cwd=str(work_dir), capture_output=True, text=True, timeout=600)
                if result.returncode == 0:
                    qm_out = work_dir / 'qm.out'
                    with open(qm_out, 'w') as f:
                        f.write(result.stdout)
                    if qm_out.exists():
                        return {'status': 'success', 'qm_out': str(qm_out), 'output_dir': str(work_dir), 'uhf_used': uhf_value}
                elif uhf_value == 2:
                    qm_out = work_dir / 'qm.out'
                    with open(qm_out, 'w') as f:
                        f.write(result.stdout)
                    return {'status': 'error', 'error_message': f'XTB failed with both uhf=1 and uhf=2', 'output_dir': str(work_dir), 'stdout': result.stdout, 'stderr': result.stderr}
            except subprocess.TimeoutExpired:
                if uhf_value == 2:
                    return {'status': 'error', 'error_message': 'XTB calculation timed out', 'output_dir': str(work_dir)}
            except Exception as e:
                if uhf_value == 2:
                    return {'status': 'error', 'error_message': str(e), 'output_dir': str(work_dir)}
        return {'status': 'error', 'error_message': 'XTB calculation failed unexpectedly', 'output_dir': str(work_dir)}

    def run_sander_vdwm(self, work_dir):
        work_dir = Path(work_dir)
        required_files = ['vdwm.mdin', 'wet.complex.prmtop', '01.min.rst7']
        for req_file in required_files:
            if not (work_dir / req_file).exists():
                return {'status': 'error', 'error_message': f'Required file {req_file} not found', 'output_dir': str(work_dir)}
        mdin_file = work_dir / 'vdwm.mdin'
        spin_attempts = self._spin_attempts(mdin_file)
        last_error = 'Sander VDWM calculation failed'
        last_stdout = None
        last_stderr = None
        for idx, spin_value in enumerate(spin_attempts):
            if idx > 0:
                self._set_spin_value(mdin_file, spin_value)
            try:
                result = subprocess.run(['sander', '-O', '-i', 'vdwm.mdin', '-o', 'vdwm.mdout', '-p', 'wet.complex.prmtop', '-c', '01.min.rst7', '-ref', '01.min.rst7', '-y', '01.min.trj', '-inf', 'vdwm.mdinfo', '-r', 'vdwm.rst7', '-x', 'vdwm.trj'], cwd=str(work_dir), capture_output=True, text=True, timeout=600)
            except subprocess.TimeoutExpired:
                last_error = 'Sander VDWM calculation timed out'
                if idx == len(spin_attempts) - 1:
                    return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir)}
                continue
            except Exception as e:
                last_error = str(e)
                if idx == len(spin_attempts) - 1:
                    return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir)}
                continue
            mdout_file = work_dir / 'vdwm.mdout'
            rst7_file = work_dir / 'vdwm.rst7'
            if result.returncode == 0 and mdout_file.exists() and rst7_file.exists():
                return {'status': 'success', 'mdout': str(mdout_file), 'rst7': str(rst7_file), 'output_dir': str(work_dir), 'spin_used': spin_value}
            last_error = 'Sander VDWM output files not generated'
            last_stdout = result.stdout
            last_stderr = result.stderr
        return {'status': 'error', 'error_message': last_error, 'output_dir': str(work_dir), 'stdout': last_stdout, 'stderr': last_stderr}

    def _run_sander_wrapper(self, args):
        work_dir, qm_charge, verbose = args
        if verbose:
            print(f'Running sander minimization in {work_dir.name}...')
        min_result = self.run_sander_min(work_dir)
        if min_result['status'] != 'success':
            return {'status': 'error', 'min_result': min_result, 'output_dir': str(work_dir), 'ligand_xtb_result': None}
        if verbose:
            print(f'Running sander SP calculation in {work_dir.name}...')
        sp_result = self.run_sander_sp(work_dir)
        if sp_result['status'] != 'success':
            return {'status': 'error', 'min_result': min_result, 'sp_result': sp_result, 'output_dir': str(work_dir), 'ligand_xtb_result': None}
        if verbose:
            print(f'Running sander VDWM calculation in {work_dir.name}...')
        vdwm_result = self.run_sander_vdwm(work_dir)
        if vdwm_result['status'] != 'success':
            return {'status': 'error', 'min_result': min_result, 'sp_result': sp_result, 'vdwm_result': vdwm_result, 'output_dir': str(work_dir), 'ligand_xtb_result': None}
        if verbose:
            print(f'Running SYS2 minimization in {work_dir.name}...')
        sys2_min_result = self.run_sander_min_sys2(work_dir)
        if sys2_min_result['status'] != 'success':
            return {'status': 'error', 'min_result': min_result, 'sp_result': sp_result, 'vdwm_result': vdwm_result, 'sys2_min_result': sys2_min_result, 'output_dir': str(work_dir), 'ligand_xtb_result': None}
        if verbose:
            print(f'Running SYS2 SP calculation in {work_dir.name}...')
        sys2_sp_result = self.run_sander_sp_sys2(work_dir)
        if sys2_sp_result['status'] != 'success':
            return {'status': 'error', 'min_result': min_result, 'sp_result': sp_result, 'vdwm_result': vdwm_result, 'sys2_min_result': sys2_min_result, 'sys2_sp_result': sys2_sp_result, 'output_dir': str(work_dir), 'ligand_xtb_result': None}
        if verbose:
            print(f'Running ligand-only XTB calculation in {work_dir.name}...')
        ligand_xtb_result = self.run_xtb_ligand(work_dir, verbose=verbose)
        if ligand_xtb_result['status'] != 'success':
            return {'status': 'error', 'min_result': min_result, 'sp_result': sp_result, 'vdwm_result': vdwm_result, 'sys2_min_result': sys2_min_result, 'sys2_sp_result': sys2_sp_result, 'ligand_xtb_result': ligand_xtb_result, 'output_dir': str(work_dir)}
        if verbose:
            print(f'Running XTB QM calculation in {work_dir.name}...')
        xtb_result = self.run_xtb_qm(work_dir, qm_charge=qm_charge)
        return {'status': 'success' if xtb_result['status'] == 'success' else 'error', 'min_result': min_result, 'sp_result': sp_result, 'vdwm_result': vdwm_result, 'sys2_min_result': sys2_min_result, 'sys2_sp_result': sys2_sp_result, 'ligand_xtb_result': ligand_xtb_result, 'xtb_result': xtb_result, 'output_dir': str(work_dir)}

    def run_sander_batch(self, output_dir, qm_charge=0, verbose=False, parallel=True, n_processes=None):
        output_dir = Path(output_dir)
        subdirs = [d for d in output_dir.rglob('*') if d.is_dir() and (d / 'min.mdin').exists() and (d / 'wet.complex.prmtop').exists()]
        if not subdirs:
            return {'status': 'error', 'error_message': 'No subdirectories with required files found'}
        if parallel:
            if n_processes is None:
                n_processes = min(cpu_count(), len(subdirs))
            if verbose:
                print(f'Running sander in parallel using {n_processes} processes...')
            args_list = [(subdir, qm_charge, False) for subdir in subdirs]
            with Pool(processes=n_processes) as pool:
                results = pool.map(self._run_sander_wrapper, args_list)
        else:
            results = []
            for subdir in subdirs:
                result = self._run_sander_wrapper((subdir, qm_charge, verbose))
                results.append(result)
        success_count = sum((1 for r in results if r['status'] == 'success'))
        error_count = sum((1 for r in results if r['status'] == 'error'))
        if verbose:
            for result in results:
                if result['status'] == 'success':
                    print(f"Completed: {Path(result['output_dir']).name} (minimization + SP + VDWM + XTB)")
                else:
                    error_msg = result.get('error_message', 'Unknown error')
                    if 'min_result' in result and result['min_result']['status'] == 'error':
                        error_msg = f"Minimization failed: {result['min_result']['error_message']}"
                    elif 'sp_result' in result and result['sp_result']['status'] == 'error':
                        error_msg = f"SP calculation failed: {result['sp_result']['error_message']}"
                    elif 'vdwm_result' in result and result['vdwm_result']['status'] == 'error':
                        error_msg = f"VDWM calculation failed: {result['vdwm_result']['error_message']}"
                    elif 'sys2_min_result' in result and result['sys2_min_result']['status'] == 'error':
                        error_msg = f"SYS2 minimization failed: {result['sys2_min_result']['error_message']}"
                    elif 'sys2_sp_result' in result and result['sys2_sp_result']['status'] == 'error':
                        error_msg = f"SYS2 SP calculation failed: {result['sys2_sp_result']['error_message']}"
                    elif 'ligand_xtb_result' in result and result['ligand_xtb_result'] and (result['ligand_xtb_result']['status'] == 'error'):
                        error_msg = f"Ligand XTB failed: {result['ligand_xtb_result']['error_message']}"
                    elif 'xtb_result' in result and result['xtb_result']['status'] == 'error':
                        error_msg = f"XTB calculation failed: {result['xtb_result']['error_message']}"
                    print(f"Failed: {Path(result['output_dir']).name} - {error_msg}")
        return {'status': 'success' if success_count > 0 else 'error', 'success_count': success_count, 'error_count': error_count, 'total_count': len(subdirs), 'results': results}

    def update_mdin_with_qmmask(self, mdin_path: Path, qmmask: str, qmcharge: int=0) -> bool:
        mdin_path = Path(mdin_path)
        try:
            content = mdin_path.read_text()
            lines = content.splitlines(keepends=True)
            updated_lines = [f" qmmask = ':{qmmask}',\n" if 'qmmask' in line and '=' in line else f' qmcharge = {qmcharge},\n' if 'qmcharge' in line and '=' in line else line for line in lines]
            mdin_path.write_text(''.join(updated_lines))
            return True
        except Exception as e:
            return False

    def _get_spin_value(self, mdin_path: Path) -> Optional[int]:
        mdin_path = Path(mdin_path)
        try:
            for line in mdin_path.read_text().splitlines():
                stripped = line.strip().lower()
                if stripped.startswith('spin') and '=' in line:
                    try:
                        value_part = line.split('=', 1)[1]
                        value = value_part.replace(',', ' ').split()[0]
                        return int(value)
                    except (ValueError, IndexError):
                        continue
        except Exception:
            return None
        return None

    def _set_spin_value(self, mdin_path: Path, spin_value: int) -> bool:
        mdin_path = Path(mdin_path)
        try:
            lines = mdin_path.read_text().splitlines(keepends=True)
        except Exception:
            return False
        updated_lines = []
        updated = False
        for line in lines:
            stripped = line.strip().lower()
            if stripped.startswith('spin') and '=' in line:
                updated_lines.append(f' spin = {spin_value},\n')
                updated = True
            else:
                updated_lines.append(line)
        if not updated:
            return False
        try:
            mdin_path.write_text(''.join(updated_lines))
            return True
        except Exception:
            return False

    def _spin_attempts(self, mdin_path: Path) -> List[int]:
        current_spin = self._get_spin_value(mdin_path) or 1
        alternate_spin = 2 if current_spin == 1 else 1
        if alternate_spin == current_spin:
            return [current_spin]
        return [current_spin, alternate_spin]

    def process_directory(self, output_dir, distance_cutoff=None, verbose=False, parallel=True, n_processes=None):
        output_dir = Path(output_dir)
        subdirs = [d for d in output_dir.rglob('*') if d.is_dir() and (d / 'leap.in').exists()]
        if not subdirs:
            return {'status': 'error', 'error_message': 'No subdirectories with leap.in found'}
        if distance_cutoff is None:
            if parallel:
                if n_processes is None:
                    n_processes = min(cpu_count(), len(subdirs))
                if verbose:
                    print(f'Running tleap in parallel using {n_processes} processes...')
                with Pool(processes=n_processes) as pool:
                    results = pool.starmap(self._run_tleap_only, [(subdir,) for subdir in subdirs])
            else:
                results = []
                for subdir in subdirs:
                    result = self._run_tleap_only(subdir)
                    results.append(result)
            tleap_success = sum((1 for r in results if r['status'] == 'success'))
            return {'status': 'success', 'total_count': len(subdirs), 'success_count': tleap_success, 'results': results}
        if parallel:
            if n_processes is None:
                n_processes = min(cpu_count(), len(subdirs))
            if verbose:
                print(f'Running tleap in parallel using {n_processes} processes...')
            args_list = [(subdir, distance_cutoff, False) for subdir in subdirs]
            with Pool(processes=n_processes) as pool:
                results = pool.map(self._run_tleap_and_analyze, args_list)
        else:
            results = []
            for subdir in subdirs:
                result = self._run_tleap_and_analyze((subdir, distance_cutoff, verbose))
                results.append(result)
        per_directory_masks = []
        for result in results:
            if result['status'] != 'success':
                continue
            residues = sorted(result['residues'], key=lambda x: x[0])
            residues_with_ligand = list(residues)
            ligand_resnum = result.get('ligand_resnum')
            if ligand_resnum is not None:
                residues_with_ligand.append((ligand_resnum, 'LIG'))
                residues_with_ligand.sort(key=lambda x: x[0])
            qmmask_local = self.format_qmmask(residues_with_ligand) if residues_with_ligand else ''
            sys2_residues = [(resnum, resname) for resnum, resname in residues_with_ligand if resname != 'LIG']
            qmmask_sys2_local = self.format_qmmask(sys2_residues) if sys2_residues else ''
            qm_charge_local = self.calculate_qm_charge(residues_with_ligand, work_dir=result.get('work_dir')) if residues_with_ligand else 0
            if verbose:
                print(f"{result['subdir']}: Found {len(residues_with_ligand)} residues for QM mask")
            per_directory_masks.append({'subdir': result['subdir'], 'residues': residues_with_ligand, 'qmmask': qmmask_local, 'qmmask_sys2': qmmask_sys2_local, 'qm_charge': qm_charge_local})
        if not per_directory_masks:
            return {'status': 'error', 'error_message': 'No successful processing', 'results': results}
        success_count = 0
        for subdir in subdirs:
            subdir_data = next((entry for entry in per_directory_masks if entry['subdir'] == subdir.name), None)
            if not subdir_data:
                continue
            mdin_files = [('min.mdin', True, subdir_data['qmmask']), ('min_sys2.mdin', False, subdir_data['qmmask_sys2']), ('sp.mdin', False, subdir_data['qmmask']), ('sp_sys2.mdin', False, subdir_data['qmmask_sys2'])]
            for filename, count_success, mask_value in mdin_files:
                mdin_path = subdir / filename
                if not mdin_path.exists():
                    continue
                if not mask_value:
                    continue
                success = self.update_mdin_with_qmmask(mdin_path, mask_value, qmcharge=subdir_data['qm_charge'])
                if success and verbose:
                    print(f'Updated {subdir.name}/{filename}')
                if success and count_success:
                    success_count += 1
        return {'status': 'success', 'success_count': success_count, 'total_count': len(subdirs), 'results': results, 'per_directory_masks': per_directory_masks}

    def _run_tleap_only(self, subdir):
        leap_in = subdir / 'leap.in'
        result = self.run_tleap(leap_in)
        if result['status'] == 'success':
            return {'status': 'success', 'subdir': subdir.name}
        return {'status': 'error', 'subdir': subdir.name, 'error_message': result.get('error_message', 'Unknown error')}

    def _run_tleap_and_analyze(self, args):
        subdir, distance_cutoff, verbose = args
        leap_in = subdir / 'leap.in'
        if verbose:
            print(f'Processing {subdir.name}...')
        result = self.run_tleap(leap_in)
        if result['status'] == 'success':
            gas_pdb = Path(result['gas_complex_pdb'])
            nearby_residues = self.get_nearby_residues(gas_pdb, distance_cutoff=distance_cutoff)
            ligand_resnum = self.get_ligand_resnum(gas_pdb)
            return {'status': 'success', 'subdir': subdir.name, 'work_dir': subdir, 'residues': nearby_residues, 'ligand_resnum': ligand_resnum}
        return {'status': 'error', 'subdir': subdir.name, 'error_message': result.get('error_message', 'Unknown error')}

def main():
    processor = XTBProcessor()
if __name__ == '__main__':
    main()