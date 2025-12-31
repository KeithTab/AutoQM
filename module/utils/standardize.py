from pathlib import Path
from typing import Optional
from pdbfixer import PDBFixer
from openmm.app import PDBFile


def standardize_pdb_file(pdb_path: Path, output_path: Path = None, 
                         verbose: bool = False) -> bool:
    pdb_path = Path(pdb_path)
    if output_path is None:
        output_path = pdb_path
    else:
        output_path = Path(output_path)
    
    if not pdb_path.exists():
        if verbose:
            print(f"  ⚠ PDB file not found: {pdb_path}")
        return False
    
    try:
        fixer = PDBFixer(filename=str(pdb_path))
        
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        
        with open(output_path, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
        
        _fix_terminal_atoms(output_path)
        _fix_disulfide_bonds(output_path, verbose=verbose)
        
        if verbose:
            print(f"  ✓ Standardized {pdb_path.name} using PDBFixer")
        
        return True
        
    except Exception as e:
        if verbose:
            print(f"  ⚠ PDBFixer failed for {pdb_path.name}: {e}")
        return False


def _fix_terminal_atoms(pdb_path: Path):
    lines = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                if ' OXT ' in line:
                    line = line.replace(' OXT ', ' O   ')
                elif ' OC1 ' in line:
                    line = line.replace(' OC1 ', ' O   ')
                elif ' OC2 ' in line:
                    line = line.replace(' OC2 ', ' OXT ')
            lines.append(line)
    
    with open(pdb_path, 'w') as f:
        f.writelines(lines)


def _fix_disulfide_bonds(pdb_path: Path, distance_cutoff: float = 2.5, verbose: bool = False):
    lines = []
    cys_residues = {}
    
    with open(pdb_path, 'r') as f:
        for line in f:
            lines.append(line)
            if line.startswith(('ATOM', 'HETATM')):
                resname = line[17:20].strip()
                if resname == 'CYS':
                    atom_name = line[12:16].strip()
                    if atom_name == 'SG':
                        chain = line[21].strip()
                        resnum = int(line[22:26].strip())
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        key = f"{chain}_{resnum}"
                        cys_residues[key] = {'x': x, 'y': y, 'z': z, 'chain': chain, 'resnum': resnum}
    
    bonded_cys = set()
    cutoff_sq = distance_cutoff ** 2
    
    for key1, cys1 in cys_residues.items():
        for key2, cys2 in cys_residues.items():
            if key1 >= key2:
                continue
            dx = cys1['x'] - cys2['x']
            dy = cys1['y'] - cys2['y']
            dz = cys1['z'] - cys2['z']
            dist_sq = dx*dx + dy*dy + dz*dz
            
            if dist_sq <= cutoff_sq:
                bonded_cys.add(key1)
                bonded_cys.add(key2)
                if verbose:
                    print(f"  ✓ Detected disulfide bond: {key1} - {key2} ({dist_sq**0.5:.2f} Å)")
    
    if bonded_cys:
        new_lines = []
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                resname = line[17:20].strip()
                if resname == 'CYS':
                    chain = line[21].strip()
                    resnum = int(line[22:26].strip())
                    key = f"{chain}_{resnum}"
                    if key in bonded_cys:
                        line = line[:17] + 'CYX' + line[20:]
            new_lines.append(line)
        
        with open(pdb_path, 'w') as f:
            f.writelines(new_lines)
        
        if verbose:
            print(f"  ✓ Renamed {len(bonded_cys)} CYS to CYX (disulfide bonded)")


def standardize_pdb_content(pdb_content: str) -> Optional[str]:
    import tempfile
    
    try:
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            tmp.write(pdb_content)
            tmp_path = tmp.name
        
        try:
            fixer = PDBFixer(filename=tmp_path)
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            
            output = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
            output_path = output.name
            output.close()
            
            with open(output_path, 'w') as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
            
            with open(output_path, 'r') as f:
                result = f.read()
            
            Path(tmp_path).unlink(missing_ok=True)
            Path(output_path).unlink(missing_ok=True)
            
            return result
            
        except Exception:
            Path(tmp_path).unlink(missing_ok=True)
            return pdb_content
            
    except Exception:
        return pdb_content
