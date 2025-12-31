from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops


def detect_ligand_charge(ligand_pdb_path: Path, verbose: bool = False) -> int:
    ligand_pdb_path = Path(ligand_pdb_path)
    
    if not ligand_pdb_path.exists():
        if verbose:
            print(f"  ⚠ Ligand file not found: {ligand_pdb_path}")
        return 0
    
    mol = Chem.MolFromPDBFile(str(ligand_pdb_path), removeHs=False, sanitize=True)
    
    if mol is None:
        if verbose:
            print(f"  ⚠ Failed to read ligand: {ligand_pdb_path.name}")
        return 0
    
    total_charge = rdmolops.GetFormalCharge(mol)
    
    if verbose:
        print(f"  ✓ Detected charge for {ligand_pdb_path.name}: {total_charge:+d}")
    
    return total_charge


def detect_batch_charges(ligand_pdb_files, verbose: bool = False):
    charges = {}
    
    for ligand_path in ligand_pdb_files:
        ligand_path = Path(ligand_path)
        charge = detect_ligand_charge(ligand_path, verbose=verbose)
        charges[str(ligand_path)] = charge
    
    return charges
