import sys
import argparse
from pathlib import Path
module_dir = Path(__file__).parent / 'module'
sys.path.insert(0, str(module_dir))
from pdbprocess import PDBProcessor
from topology import TopologyGenerator
from qmmm import MdinGenerator, XTBProcessor
from analysis import XTBAnalyzer

def main():
    parser = argparse.ArgumentParser(description='Process PDB files to separate protein and ligand, and generate topology')
    parser.add_argument('-i', '--input', type=str, default='samp', help='Input directory containing PDB files (default: samp)')
    parser.add_argument('-o', '--output', type=str, default='output', help='Output directory for processed files (default: output)')
    parser.add_argument('-t', '--topology', action='store_true', help='Generate topology files for ligands using acpype')
    parser.add_argument('--xtb', action='store_true', help='Run tleap, calculate QM mask, and run sander minimization')
    parser.add_argument('-d', '--distance', type=float, default=None, help='Distance cutoff for finding nearby residues (in Angstroms). If not specified, qmmask and qmcharge will not be modified')
    parser.add_argument('-p', '--parallel', action='store_true', help='Enable parallel processing for topology generation and tleap')
    parser.add_argument('-n', '--nproc', type=int, default=None, help='Number of processes for parallel processing (default: auto-detect CPU count)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('--analyze', action='store_true', help='Analyze output files and generate Excel report')
    args = parser.parse_args()
    processor = PDBProcessor()
    base_dir = Path(__file__).parent
    input_dir = base_dir / args.input if not Path(args.input).is_absolute() else Path(args.input)
    output_dir = base_dir / args.output if not Path(args.output).is_absolute() else Path(args.output)
    if not input_dir.exists():
        print(f"Error: Input directory '{input_dir}' does not exist")
        return 1
    pdb_files = list(input_dir.rglob('*.pdb'))
    if not pdb_files:
        print(f"Error: No PDB files found in '{input_dir}'")
        return 1
    results = processor.process_batch(pdb_files, output_dir, verbose=args.verbose, parallel=args.parallel, n_processes=args.nproc)
    success_count = sum((1 for r in results if r.get('status') == 'success'))
    error_count = len(results) - success_count
    if args.verbose:
        print(f'\nProcessed {len(results)} files: {success_count} successful, {error_count} errors')
        if error_count > 0:
            error_msgs = '\n'.join((f"  Error in {Path(r['input_file']).name}: {r['error_message']}" for r in results if r.get('status') == 'error'))
            print(error_msgs)
    if args.topology and success_count > 0:
        if args.verbose:
            print('\nGenerating topology files...')
        generator = TopologyGenerator()
        topo_results = generator.process_from_directory(output_dir, verbose=args.verbose, parallel=args.parallel, n_processes=args.nproc)
        topo_success = sum((1 for r in topo_results if r.get('status') == 'success'))
        topo_error = len(topo_results) - topo_success
        if args.verbose:
            print(f'\nGenerated topology for {len(topo_results)} ligands: {topo_success} successful, {topo_error} errors')
            if topo_error > 0:
                error_msgs = '\n'.join((f"  Error in {r['ligand_name']}: {r['error_message']}" for r in topo_results if r.get('status') == 'error'))
                print(error_msgs)
    if success_count > 0:
        if args.verbose:
            print('\nGenerating min.mdin files...')
        mdin_generator = MdinGenerator()
        output_subdirs = [d for d in output_dir.rglob('*') if d.is_dir() and (d / f'{d.name}_ligand.pdb').exists()]
        mdin_results = mdin_generator.generate_batch(output_subdirs)
        mdin_success = sum((1 for r in mdin_results if r.get('status') == 'success'))
        mdin_error = len(mdin_results) - mdin_success
        if args.verbose:
            print(f'\nGenerated mdin files for {len(mdin_results)} directories: {mdin_success} successful, {mdin_error} errors')
            if mdin_error > 0:
                error_msgs = '\n'.join((f"  Error in {r['output_dir']}: {r['error_message']}" for r in mdin_results if r.get('status') == 'error'))
                print(error_msgs)
    if args.xtb and success_count > 0:
        xtb_processor = XTBProcessor()
        if args.distance is not None:
            if args.verbose:
                print(f'\nRunning XTB processing (tleap and QM mask calculation with {args.distance}A cutoff)...')
            xtb_result = xtb_processor.process_directory(output_dir, distance_cutoff=args.distance, verbose=args.verbose, parallel=args.parallel, n_processes=args.nproc)
        else:
            if args.verbose:
                print(f'\nRunning tleap (QM mask calculation skipped - no distance specified)...')
            xtb_result = xtb_processor.process_directory(output_dir, distance_cutoff=None, verbose=args.verbose, parallel=args.parallel, n_processes=args.nproc)
        if xtb_result['status'] == 'success':
            if args.distance is not None:
                print(f"\nUpdated per-structure QM masks for {xtb_result['success_count']}/{xtb_result['total_count']} directories")
                if args.verbose:
                    mask_details = xtb_result.get('per_directory_masks', [])
                    for entry in mask_details:
                        qmmask_display = entry['qmmask'] or '(empty)'
                        print(f"  {entry['subdir']}: qmmask=':{qmmask_display}' | qmcharge={entry['qm_charge']:+d}")
            elif args.verbose:
                print(f"\ntleap completed for {xtb_result['total_count']} directories (QM mask not calculated)")
            if args.verbose:
                print(f'\nRunning sander minimization...')
            sander_result = xtb_processor.run_sander_batch(output_dir, verbose=args.verbose, parallel=args.parallel, n_processes=args.nproc)
            if sander_result['status'] == 'success':
                print(f"\nSander minimization completed: {sander_result['success_count']}/{sander_result['total_count']} successful")
                if sander_result['error_count'] > 0:
                    print(f"Errors: {sander_result['error_count']}")
                    if args.verbose:
                        error_msgs = '\n'.join((f"  {Path(r['output_dir']).name}: {r.get('error_message', 'Unknown error')}" for r in sander_result['results'] if r['status'] == 'error'))
                        if error_msgs:
                            print(error_msgs)
            else:
                print(f"\nSander execution failed: {sander_result.get('error_message', 'Unknown error')}")
        else:
            print(f"\nXTB processing failed: {xtb_result.get('error_message', 'Unknown error')}")
    if args.analyze:
        if args.verbose:
            print('\n' + '=' * 60)
            print('Starting analysis of output files...')
            print('=' * 60)
        try:
            analyzer = XTBAnalyzer(output_dir)
            excel_file = analyzer.run_analysis(verbose=args.verbose)
            if excel_file:
                print(f'\n✓ Analysis complete! Results saved to: {excel_file}')
            else:
                print('\n✗ Analysis failed: No output files found')
        except ImportError as e:
            print(f'\n✗ Analysis failed: Missing required package')
            print(f'  Please install: pip install pandas openpyxl')
            print(f'  Error: {e}')
        except Exception as e:
            print(f'\n✗ Analysis failed: {e}')
    return 0 if error_count == 0 else 1
if __name__ == '__main__':
    sys.exit(main())