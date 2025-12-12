"""
Optimized XTB Analysis Module
Analyzes output files from QMMM calculations and exports results to Excel.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Generator
import pandas as pd
from functools import lru_cache


class XTBAnalyzer:
    """Optimized analyzer for XTB QMMM calculation results"""
    
    __slots__ = ('output_dir', '_regex_patterns', '_eh_to_kcal')
    
    # Class-level compiled regex patterns for better performance
    QM_ENERGY_PATTERN = re.compile(r'total energy\s+([-+]?\d+\.\d+)\s+Eh', re.IGNORECASE)
    VDWAALS_PATTERN = re.compile(r'VDWAALS\s*=\s*([-+]?\d+\.\d+)', re.IGNORECASE)
    XTBESCF_PATTERN = re.compile(r'XTBESCF\s*=\s*([-+]?\d+\.\d+)', re.IGNORECASE)
    EPTOT_PATTERN = re.compile(r'EPtot\s*=\s*([-+]?\d+\.\d+)', re.IGNORECASE)
    
    def __init__(self, output_dir: Path):
        """
        Initialize the analyzer with pre-compiled regex patterns
        
        Args:
            output_dir: Path to the output directory containing calculation results
        """
        self.output_dir = Path(output_dir)
        self._eh_to_kcal = 627.509474
        
    def extract_qm_energy(self, qm_out_file: Path) -> Optional[float]:
        """
        Extract total energy from qm.out file using compiled regex
        
        Args:
            qm_out_file: Path to qm.out file
            
        Returns:
            Total energy in Eh (Hartree), or None if not found
        """
        if not qm_out_file.exists():
            return None
            
        try:
            # Read file once
            content = qm_out_file.read_text()
            match = self.QM_ENERGY_PATTERN.search(content)
            
            return float(match.group(1)) if match else None
                
        except Exception as e:
            print(f"Error reading {qm_out_file}: {e}")
            return None
    
    def extract_vdwm_energy(self, vdwm_mdout_file: Path) -> Optional[float]:
        """
        Extract VDWAALS energy from vdwm.mdout file using compiled regex
        
        Args:
            vdwm_mdout_file: Path to vdwm.mdout file
            
        Returns:
            VDWAALS energy in kcal/mol, or None if not found
        """
        if not vdwm_mdout_file.exists():
            return None
            
        try:
            content = vdwm_mdout_file.read_text()
            match = self.VDWAALS_PATTERN.search(content)
            
            return float(match.group(1)) if match else None
                        
        except Exception as e:
            print(f"Error reading {vdwm_mdout_file}: {e}")
            return None
    
    def extract_sp_energies(self, sp_mdout_file: Path) -> Tuple[Optional[float], Optional[float]]:
        """
        Extract VDWAALS and XTBESCF energies from sp.mdout file using compiled regex
        
        Args:
            sp_mdout_file: Path to sp.mdout file
            
        Returns:
            Tuple of (VDWAALS, XTBESCF) energies in kcal/mol, or (None, None) if not found
        """
        if not sp_mdout_file.exists():
            return None, None
            
        try:
            content = sp_mdout_file.read_text()
            
            # Use compiled regex patterns for efficiency
            vdw_match = self.VDWAALS_PATTERN.search(content)
            xtb_match = self.XTBESCF_PATTERN.search(content)
            
            vdwaals = float(vdw_match.group(1)) if vdw_match else None
            xtbescf = float(xtb_match.group(1)) if xtb_match else None
            
            return vdwaals, xtbescf
                    
        except Exception as e:
            print(f"Error reading {sp_mdout_file}: {e}")
            return None, None

    def extract_eptot(self, mdout_file: Path) -> Optional[float]:
        """Extract EPtot from an mdout file."""
        if not mdout_file.exists():
            return None
        try:
            content = mdout_file.read_text()
            match = self.EPTOT_PATTERN.search(content)
            return float(match.group(1)) if match else None
        except Exception as e:
            print(f"Error reading {mdout_file}: {e}")
            return None
    
    def analyze_directory(self, dir_path: Path) -> Dict[str, Optional[float]]:
        """
        Analyze all output files in a single directory.
        
        Args:
            dir_path: Path to the directory containing output files
            
        Returns:
            Dictionary with extracted energy values
        """
        result = {
            'directory': dir_path.name,
            'qm_energy_Eh': None,
            'qm_ligand_Eh': None,
            'vdwm_vdwaals_kcal': None,
            'sp_vdwaals_kcal': None,
            'sp_xtbescf_kcal': None,
            'sp_eptot_kcal': None,
            'sp_sys2_eptot_kcal': None
        }
        
        # Extract QM energy from qm.out
        qm_out = dir_path / 'qm.out'
        result['qm_energy_Eh'] = self.extract_qm_energy(qm_out)
        qm_lig_out = dir_path / 'qm_ligand.out'
        result['qm_ligand_Eh'] = self.extract_qm_energy(qm_lig_out)
        
        # Extract VDWAALS from vdwm.mdout
        vdwm_mdout = dir_path / 'vdwm.mdout'
        result['vdwm_vdwaals_kcal'] = self.extract_vdwm_energy(vdwm_mdout)
        
        # Extract VDWAALS and XTBESCF from sp.mdout
        sp_mdout = dir_path / 'sp.mdout'
        sp_vdwaals, sp_xtbescf = self.extract_sp_energies(sp_mdout)
        result['sp_vdwaals_kcal'] = sp_vdwaals
        result['sp_xtbescf_kcal'] = sp_xtbescf
        result['sp_eptot_kcal'] = self.extract_eptot(sp_mdout)
        
        sp_sys2_mdout = dir_path / 'sp_sys2.mdout'
        result['sp_sys2_eptot_kcal'] = self.extract_eptot(sp_sys2_mdout)
        
        return result
    
    def _find_output_dirs(self) -> Generator[Path, None, None]:
        """Generator to find directories with output files"""
        for subdir in sorted(self.output_dir.rglob('*')):
            if subdir.is_dir() and any([
                (subdir / 'qm.out').exists(),
                (subdir / 'qm_ligand.out').exists(),
                (subdir / 'sp.mdout').exists(),
                (subdir / 'sp_sys2.mdout').exists(),
                (subdir / 'vdwm.mdout').exists()
            ]):
                yield subdir
    
    def analyze_all(self, verbose: bool = True) -> pd.DataFrame:
        """
        Analyze all subdirectories using optimized DataFrame operations
        
        Args:
            verbose: Whether to print progress messages
            
        Returns:
            pandas DataFrame with all results
        """
        # Use generator to collect results
        results = []
        for subdir in self._find_output_dirs():
            if verbose:
                print(f"Analyzing: {subdir.relative_to(self.output_dir)}")
            results.append(self.analyze_directory(subdir))
        
        # Create DataFrame once
        df = pd.DataFrame(results)
        
        if df.empty:
            return df
        
        # Vectorized operations for better performance
        if 'qm_energy_Eh' in df.columns:
            df['qm_energy_kcal'] = df['qm_energy_Eh'] * self._eh_to_kcal
        if 'qm_ligand_Eh' in df.columns:
            df['qm_ligand_kcal'] = df['qm_ligand_Eh'] * self._eh_to_kcal
        
        # Calculate all energies in one go using vectorized operations
        if all(col in df.columns for col in ['sp_xtbescf_kcal', 'qm_energy_kcal']):
            df['E(QM/MM)_elec'] = df['sp_xtbescf_kcal'] - df['qm_energy_kcal']
        
        if all(col in df.columns for col in ['sp_vdwaals_kcal', 'vdwm_vdwaals_kcal']):
            df['E(QM/MM)_LJ'] = df['sp_vdwaals_kcal'] - df['vdwm_vdwaals_kcal']
        
        if all(col in df.columns for col in ['E(QM/MM)_elec', 'E(QM/MM)_LJ']):
            df['E(QM/MM)'] = df['E(QM/MM)_elec'] + df['E(QM/MM)_LJ']

        # New composite energy: EPtot(sp) - EPtot(sp_sys2) - QM ligand (kcal)
        if all(col in df.columns for col in ['sp_eptot_kcal', 'sp_sys2_eptot_kcal', 'qm_ligand_kcal']):
            df['E(QM/MM)_new'] = df['sp_eptot_kcal'] - df['sp_sys2_eptot_kcal'] - df['qm_ligand_kcal']
        
        # Reorder columns efficiently
        column_order = [
            'directory',
            'qm_energy_Eh', 'qm_energy_kcal',
            'qm_ligand_Eh', 'qm_ligand_kcal',
            'sp_eptot_kcal', 'sp_sys2_eptot_kcal',
            'sp_vdwaals_kcal', 'sp_xtbescf_kcal', 'vdwm_vdwaals_kcal',
            'E(QM/MM)_elec', 'E(QM/MM)_LJ', 'E(QM/MM)', 'E(QM/MM)_new'
        ]
        
        existing_cols = [col for col in column_order if col in df.columns]
        
        return df[existing_cols]
    
    def save_to_excel(self, df: pd.DataFrame, output_file: Optional[Path] = None) -> Path:
        """
        Save analysis results to Excel file.
        
        Args:
            df: DataFrame containing analysis results
            output_file: Path to output Excel file (default: output_dir/analysis_results.xlsx)
            
        Returns:
            Path to the saved Excel file
        """
        if output_file is None:
            output_file = self.output_dir / 'analysis_results.xlsx'
        else:
            output_file = Path(output_file)
        
        # Create Excel writer with formatting
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Energy Analysis', index=False)
            
            # Get workbook and worksheet
            workbook = writer.book
            worksheet = writer.sheets['Energy Analysis']
            
            # Auto-adjust column widths
            for column in worksheet.columns:
                max_length = 0
                column_letter = column[0].column_letter
                
                for cell in column:
                    try:
                        if cell.value:
                            max_length = max(max_length, len(str(cell.value)))
                    except:
                        pass
                
                adjusted_width = min(max_length + 2, 50)
                worksheet.column_dimensions[column_letter].width = adjusted_width
        
        return output_file
    
    def run_analysis(self, verbose: bool = True) -> Path:
        """
        Run complete analysis and save results to Excel.
        
        Args:
            verbose: Whether to print progress messages
            
        Returns:
            Path to the saved Excel file
        """
        if verbose:
            print(f"Starting analysis of: {self.output_dir}")
            print("-" * 60)
        
        # Analyze all directories
        df = self.analyze_all(verbose=verbose)
        
        if df.empty:
            print("No output files found!")
            return None
        
        # Save to Excel
        excel_file = self.save_to_excel(df)
        
        if verbose:
            print("-" * 60)
            print(f"Analysis complete!")
            print(f"Total directories analyzed: {len(df)}")
            print(f"Results saved to: {excel_file}")
            print("\nSummary statistics:")
            print(df.describe())
        
        return excel_file


def main():
    """Command-line interface for the analyzer."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Analyze XTB QMMM calculation results'
    )
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Path to output directory containing calculation results'
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        help='Path to output Excel file (default: analysis_results.xlsx in input directory)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print verbose output'
    )
    
    args = parser.parse_args()
    
    # Create analyzer
    analyzer = XTBAnalyzer(args.input)
    
    # Run analysis
    if args.output:
        df = analyzer.analyze_all(verbose=args.verbose)
        excel_file = analyzer.save_to_excel(df, args.output)
    else:
        excel_file = analyzer.run_analysis(verbose=args.verbose)
    
    if excel_file:
        print(f"\n✓ Analysis complete: {excel_file}")
    else:
        print("\n✗ Analysis failed: No output files found")


if __name__ == '__main__':
    main()
