#!/bin/bash

# QMMM Gromacs Processing Script
# This script converts Amber files to Gromacs format and runs MD analysis

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
print_message() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Check if output directory is provided
if [ $# -eq 0 ]; then
    print_error "No output directory specified!"
    echo "Usage: $0 <output_directory>"
    exit 1
fi

OUTPUT_DIR="$1"

# Check if output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    print_error "Output directory '$OUTPUT_DIR' does not exist!"
    exit 1
fi

print_message "Starting Gromacs processing in: $OUTPUT_DIR"
echo "========================================================"

# Find all subdirectories with wet.complex.prmtop and 01.min.rst7
for subdir in "$OUTPUT_DIR"/*; do
    if [ ! -d "$subdir" ]; then
        continue
    fi
    
    PRMTOP="$subdir/wet.complex.prmtop"
    RST7="$subdir/01.min.rst7"
    
    # Check if required files exist
    if [ ! -f "$PRMTOP" ] || [ ! -f "$RST7" ]; then
        print_warning "Skipping $(basename $subdir): missing wet.complex.prmtop or 01.min.rst7"
        continue
    fi
    
    print_message "Processing: $(basename $subdir)"
    cd "$subdir" || continue
    
    # Step 1: Convert Amber files to Gromacs format using ParmEd
    print_message "  Converting Amber files to Gromacs format..."
    python << EOF
import parmed as pmd
import warnings
import numpy as np

# Suppress numpy warnings
warnings.filterwarnings('ignore')

# Monkey patch for numpy 2.0 compatibility
original_array = np.array
def patched_array(*args, **kwargs):
    if 'copy' in kwargs and kwargs['copy'] is False:
        # Remove copy and other incompatible args
        kwargs.pop('copy')
        if 'subok' in kwargs:
            subok = kwargs.pop('subok')
        else:
            subok = False
        # Use asarray without incompatible arguments
        result = np.asarray(*args)
        return result
    return original_array(*args, **kwargs)
np.array = patched_array

try:
    # Load Amber topology and coordinates
    amber = pmd.load_file('wet.complex.prmtop', '01.min.rst7')
    
    # Save as Gromacs format (overwrite existing files)
    amber.save('gromacs.top', overwrite=True)
    amber.save('gromacs.gro', overwrite=True)
    
    print("Conversion complete!")
except Exception as e:
    print(f"Error during conversion: {e}")
    exit(1)
finally:
    # Restore original numpy.array
    np.array = original_array
EOF
    
    if [ $? -ne 0 ]; then
        print_error "  Failed to convert files in $(basename $subdir)"
        cd - > /dev/null
        continue
    fi
    
    # Step 2: Create md.mdp file
    print_message "  Creating md.mdp file..."
    cat > md.mdp << 'MDPEOF'
define =
integrator = md
dt         = 0.002   ; ps
nsteps     = 500000
comm-grps  = system 
comm-mode  = angular
energygrps = Protein L01
;
nstxout = 0
nstvout = 0
nstfout = 0
nstlog  = 5000
nstenergy = 1000
nstxout-compressed = 1000
compressed-x-grps  = system
;
pbc = xyz
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
vdwtype       = cut-off
rvdw          = 1.0
DispCorr      = EnerPres
;
Tcoupl  = V-rescale
tau_t   = 0.2
tc_grps = system 
ref_t   = 300
;
Pcoupl     = parrinello-rahman
pcoupltype = isotropic
tau_p = 2.0
ref_p = 1.0
compressibility = 4.5e-5
;
freezegrps  = 
freezedim   = 
constraints = hbonds
MDPEOF
    
    # Step 3: Create index file with gmx make_ndx (interactive mode)
    print_message "  Creating index file (interactive mode - follow prompts)..."
    print_warning "  You may need to define QM and MM groups manually"
    print_warning "  Type 'q' to quit and save when done"
    # gmx make_ndx -f gromacs.gro -o index.ndx
    
    if [ $? -ne 0 ]; then
        print_error "  Failed to create index file in $(basename $subdir)"
        cd - > /dev/null
        continue
    fi
    
    # Step 4: Run gmx grompp
    print_message "  Running gmx grompp..."
    gmx grompp -f md.mdp -c gromacs.gro -p gromacs.top -o md.tpr -maxwarn 10
    
    if [ $? -ne 0 ]; then
        print_error "  Failed to run grompp in $(basename $subdir)"
        cd - > /dev/null
        continue
    fi
    
    # Step 5: Run gmx mdrun in rerun mode
    print_message "  Running gmx mdrun (rerun mode)..."
    gmx mdrun -v -deffnm md -rerun gromacs.gro -ntmpi 1
    
    if [ $? -ne 0 ]; then
        print_error "  Failed to run mdrun in $(basename $subdir)"
        cd - > /dev/null
        continue
    fi
    
    # Step 6: Extract energy components
    print_message "  Extracting energy components..."
    
    echo 21 0 | gmx energy -f md.edr -o Coul_interaction.xvg
    echo 22 0 | gmx energy -f md.edr -o LJ_interaction.xvg
    
    if [ $? -eq 0 ]; then
        print_message "  ✓ Successfully completed processing for $(basename $subdir)"
    else
        print_warning "  Energy extraction completed with warnings for $(basename $subdir)"
    fi
    
    cd - > /dev/null
    echo ""
done

echo "========================================================"
print_message "Aggregating interaction energies across directories..."

coul_summary="$OUTPUT_DIR/coulomb_summary.txt"
lj_summary="$OUTPUT_DIR/lj_summary.txt"
echo "Directory\tCoulomb" > "$coul_summary"
echo "Directory\tLennard-Jones" > "$lj_summary"

while IFS= read -r -d '' xvg_file; do
    dir_name=$(basename "$(dirname "$xvg_file")")
    last_line=$(grep -hv '^[@#]' "$xvg_file" | awk 'NF{line=$0} END{print line}')
    value=$(echo "$last_line" | awk '{print $NF}')
    echo -e "$dir_name\t$value" >> "$coul_summary"
done < <(find "$OUTPUT_DIR" -type f -name 'Coul_interaction.xvg' -print0)

while IFS= read -r -d '' xvg_file; do
    dir_name=$(basename "$(dirname "$xvg_file")")
    last_line=$(grep -hv '^[@#]' "$xvg_file" | awk 'NF{line=$0} END{print line}')
    value=$(echo "$last_line" | awk '{print $NF}')
    echo -e "$dir_name\t$value" >> "$lj_summary"
done < <(find "$OUTPUT_DIR" -type f -name 'LJ_interaction.xvg' -print0)

print_message "Summaries written to:"
print_message "  $coul_summary"
print_message "  $lj_summary"
print_message "All processing complete!"
