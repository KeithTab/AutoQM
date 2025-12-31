from pathlib import Path
from typing import Optional, Set, Tuple
from Bio import PDB

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    PDBFIXER_AVAILABLE = True
except ImportError:
    PDBFIXER_AVAILABLE = False
    print("Warning: PDBFixer not available. Install with: conda install -c conda-forge pdbfixer")


def standardize_pdb_file(pdb_path: Path, output_path: Path = None, 
                         verbose: bool = False) -> bool:
    if not PDBFIXER_AVAILABLE:
        if verbose:
            print("  ⚠ PDBFixer not available, skipping standardization")
        return False
    
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


def _detect_disulfide_bonds_biopython(pdb_path: Path, verbose: bool = False) -> Set[Tuple[str, int]]:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', str(pdb_path))
    
    bonded_cys = set()
    ssbond_records = set()
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('SSBOND'):
                try:
                    chain1 = line[15].strip()
                    res1 = int(line[17:21].strip())
                    chain2 = line[29].strip()
                    res2 = int(line[31:35].strip())
                    ssbond_records.add((chain1, res1))
                    ssbond_records.add((chain2, res2))
                    if verbose:
                        print(f"  ✓ Found SSBOND record: {chain1}_{res1} - {chain2}_{res2}")
                except:
                    pass
    
    if ssbond_records:
        bonded_cys = ssbond_records
    else:
        cys_sg_atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() == 'CYS':
                        if 'SG' in residue:
                            sg_atom = residue['SG']
                            cys_sg_atoms.append({
                                'atom': sg_atom,
                                'chain': chain.id,
                                'resnum': residue.id[1],
                                'key': (chain.id, residue.id[1])
                            })
        
        ns = PDB.NeighborSearch([item['atom'] for item in cys_sg_atoms])
        
        for i, cys1 in enumerate(cys_sg_atoms):
            for cys2 in cys_sg_atoms[i+1:]:
                distance = cys1['atom'] - cys2['atom']
                if 1.8 <= distance <= 2.5:
                    bonded_cys.add(cys1['key'])
                    bonded_cys.add(cys2['key'])
                    if verbose:
                        print(f"  ✓ Detected disulfide bond: {cys1['chain']}_{cys1['resnum']} - {cys2['chain']}_{cys2['resnum']} ({distance:.2f} Å)")
    
    return bonded_cys


def _fix_disulfide_bonds(pdb_path: Path, verbose: bool = False):
    try:
        bonded_cys = _detect_disulfide_bonds_biopython(pdb_path, verbose=verbose)
    except Exception as e:
        if verbose:
            print(f"  ⚠ BioPython disulfide detection failed: {e}, using fallback")
        bonded_cys = set()
    
    if not bonded_cys:
        return
    
    lines = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                resname = line[17:20].strip()
                if resname == 'CYS':
                    chain = line[21].strip()
                    resnum = int(line[22:26].strip())
                    if (chain, resnum) in bonded_cys:
                        line = line[:17] + 'CYX' + line[20:]
            lines.append(line)
    
    with open(pdb_path, 'w') as f:
        f.writelines(lines)
    
    if verbose:
        print(f"  ✓ Renamed {len(bonded_cys)} CYS to CYX (disulfide bonded)")


def standardize_pdb_content(pdb_content: str) -> Optional[str]:
    if not PDBFIXER_AVAILABLE:
        return pdb_content
    
    try:
        import tempfile
        
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
