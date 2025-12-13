# QMMM Workflow Toolkit

A streamlined toolkit for preparing, running, and analyzing QM/MM calculations with AmberTools, ACPYPE, Open Babel, XTB, and Gromacs. The workflow splits protein/ligand, builds topologies, prepares QM/MM inputs, runs Amber + XTB, and optionally converts to Gromacs for rerun analyses.

## Features
- **PDB processing**: Split protein/ligand, add hydrogens to ligands.
- **Topology generation**: ACPYPE + tleap inputs (`leap.in`, `leap_sys2.in`) with GAFF2/ff14SB.
- **MDIN generation**: `min.mdin`, `sp.mdin`, `vdwm.mdin`, plus SYS2 variants for receptor-only QM/MM.
- **QM/MM orchestration**: tleap, QM mask/charge propagation (per-directory), sander minimization/SP/VDWM, XTB QM, ligand-only XTB.
- **Gromacs rerun pipeline**: Amber→Gromacs conversion, rerun energies, Coulomb/LJ extraction.
- **Analysis**: Collect EPtot, ligand energies (Eh→kcal/mol), composite energy metrics to Excel.

## Prerequisites
- Python 3.9+ with `parmed`, `pandas`, `openpyxl`, `biopython`, `numpy`
- AmberTools (`tleap`, `sander`, `xtb` available in PATH)
- ACPYPE (`acpype` in PATH)
- Open Babel (`obabel` in PATH)
- Gromacs (`gmx` in PATH) for rerun analysis
- Optional: `openpyxl` for Excel output, `parmed` for Amber→Gromacs conversion

## Quick Start
```bash
# 1) Run the main pipeline
python main.py -i samp -o output -t --xtb -d 2.5 --parallel -n 4 -c 0 -v
#   -i: input folder with PDBs
#   -o: output folder for processed systems
#   -t: generate ligand topologies
#   --xtb: run tleap + QM mask + sander + xtb pipeline
#   -d: distance cutoff (Å) for QM region selection
#   -c: ligand net charge
#   --parallel/-n: enable and size multiprocessing
#   -v: verbose logging

# 2) Analyze results (creates analysis_results.xlsx)
python -m module.analysis.analysis_xtb -i output

# 3) Optional: Gromacs rerun and energy summaries
./run.sh output
```

## Workflow Overview
1. **Split & Hydrogenate** (`module/pdbprocess/read.py`)
   - Extract protein (`*_protein.pdb`) and ligand (`*_ligand.pdb`), add hydrogens via Open Babel.
2. **Topology** (`module/topology/topgen.py`)
   - ACPYPE generates GAFF2 parameters, tleap inputs (`leap.in`, `leap_sys2.in`) are written.
3. **MDIN Templates** (`module/qmmm/mdin.py`)
   - Emits `min.mdin`, `sp.mdin`, `vdwm.mdin`, plus `min_sys2.mdin`, `sp_sys2.mdin`.
4. **QM/MM Run** (`module/qmmm/xtb.py`)
   - tleap build, per-directory QM mask/charge, sander (min/SP/VDWM), SYS2 stages, XTB QM, ligand-only XTB with UHF fallback.
5. **Analysis** (`module/analysis/analysis_xtb.py`)
   - Extracts EPtot (sp, sp_sys2), ligand total energy (Eh→kcal/mol), and computes `E(QM/MM)_new = EPtot(sp) - EPtot(sp_sys2) - qm_ligand_kcal` plus legacy terms.
6. **Gromacs Rerun** (`run.sh`)
   - Amber→Gromacs, rerun energies, aggregate Coulomb/LJ summaries (`coulomb_summary.txt`, `lj_summary.txt`).

## Key Commands
- **Main pipeline**: `python main.py -i INPUT -o OUTPUT -t --xtb -d 2.5 -c 0 --parallel`
- **Topology only**: `python main.py -i INPUT -o OUTPUT -t`
- **QM/MM only on existing prep**: `python main.py -i INPUT -o OUTPUT --xtb -d 2.5`
- **Analysis**: `python -m module.analysis.analysis_xtb -i OUTPUT`
- **Gromacs rerun**: `./run.sh OUTPUT`

## Outputs (per system folder)
- `*_protein.pdb`, `*_ligand.pdb`
- ACPYPE dir with GAFF2 files; `leap.in`, `leap_sys2.in`
- MDIN files: `min.mdin`, `sp.mdin`, `vdwm.mdin`, `min_sys2.mdin`, `sp_sys2.mdin`
- Amber results: `wet.complex.prmtop/rst7`, `01.min.*`, `sp.*`, `vdwm.*`
- SYS2 results: `02.min.*`, `sp_sys2.*`
- XTB results: `qm.out`, `qmmm_region.{pdb,xyz}`, ligand `qmmm.xyz`, `qm_ligand.out`
- Analysis: `analysis_results.xlsx`
- Gromacs: `gromacs.top`, `gromacs.gro`, `md.tpr`, `md.edr`, `Coul_interaction.xvg`, `LJ_interaction.xvg`, summaries in output root

## Notes
- QM masks and charges are computed per directory (no global intersection). SYS2 masks exclude the ligand.
- XTB UHF fallback: QM/MM (`run_xtb_qm`) and ligand-only (`run_xtb_ligand`) retry with alternate UHF (1/2) on failure.
- Sander spin fallback: all sander stages retry with spin toggled between 1 and 2 if a run fails.
- Hydrogens are added to ligands both during initial split and again just before ACPYPE to ensure completeness.

## Troubleshooting
- **Missing executables**: ensure `acpype`, `tleap`, `sander`, `xtb`, `obabel`, `gmx` are in PATH.
- **Topology reuse**: existing ACPYPE outputs are reused; delete per-ligand `_acpype` folders to force regeneration.
- **Index creation**: `run.sh` leaves `gmx make_ndx` commented; run interactively if custom groups are needed.

## License
Internal project use. Update this section if licensing is required.
