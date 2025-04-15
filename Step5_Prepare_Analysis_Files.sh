#!/bin/bash

# Purpose: Copy selected trajectory files, extract final frames,
#          and run hbond analysis on selected trajectories.
# Depends on: CSV file listing stable replicates from Step 3.

# Exit immediately if a command fails.
set -e

# --- Configuration ---
BASE_DIR=$(pwd)
RESULT_DIR="${BASE_DIR}/MD_Result"
ANALYSIS_OUTPUT_DIR="${BASE_DIR}/Analysis_Output" # Where input CSV is
COMPARE_DIR="${BASE_DIR}/Compare_Selected"      # NEW Directory for analysis files

# Input CSV from Step 3
STABLE_REPS_CSV="${ANALYSIS_OUTPUT_DIR}/most_stable_replicates_Rg_variance.csv"

# GROMACS command
GMX="gmx"

# GROMACS selections (CHECK THESE!)
TRJCONV_DUMP_SELECT="System"        # Group for extracting PDB frame (usually System or Protein)
HBOND_GROUP_SELECT="protein protein" # Groups for hbond calculation

# --- Sanity Check ---
if [ ! -f "${STABLE_REPS_CSV}" ]; then
    echo "ERROR: Stable replicates CSV file not found: ${STABLE_REPS_CSV}"
    echo "Please run Step3_Plot_Select_Rg.R first."
    exit 1
fi

# --- Create Comparison Directory ---
echo "--- Creating directory for selected analysis files: ${COMPARE_DIR} ---"
mkdir -p "${COMPARE_DIR}"

# --- Process Selected Replicates ---
echo "--- Processing files for selected stable replicates ---"

# Read the CSV file, skipping the header line
tail -n +2 "${STABLE_REPS_CSV}" | while IFS=',' read -r temp rep_id run_id variance mean_rg; do
    # Clean up potential whitespace issues from read (though unlikely with CSV)
    temp=$(echo "$temp" | xargs)
    run_id=$(echo "$run_id" | xargs)

    echo "Processing: Temperature ${temp}, Run ID ${run_id}"

    SOURCE_DIR="${RESULT_DIR}/${temp}/${run_id}"
    TPR_FILE="md_${run_id}.tpr"
    XTС_FILE="md_${run_id}_noPBC.xtc" # Use the noPBC trajectory

    # Check if source files exist
    if [ ! -f "${SOURCE_DIR}/${TPR_FILE}" ] || [ ! -f "${SOURCE_DIR}/${XTС_FILE}" ]; then
        echo "  WARNING: Source .tpr or .xtc file missing for ${run_id}. Skipping."
        continue # Skip to the next line in the CSV
    fi

    # 1. Copy selected TPR and XTC files with simplified names
    echo "  Copying ${TPR_FILE} to ${COMPARE_DIR}/md_${temp}.tpr"
    cp "${SOURCE_DIR}/${TPR_FILE}" "${COMPARE_DIR}/md_${temp}.tpr"
    echo "  Copying ${XTС_FILE} to ${COMPARE_DIR}/traj_${temp}.xtc"
    cp "${SOURCE_DIR}/${XTС_FILE}" "${COMPARE_DIR}/traj_${temp}.xtc"

    # 2. Extract Final Frame (50 ns = 50000 ps)
    echo "  Extracting final PDB frame (50ns) to final_stable_${temp}.pdb"
    echo "${TRJCONV_DUMP_SELECT}" | ${GMX} trjconv \
        -s "${COMPARE_DIR}/md_${temp}.tpr" \
        -f "${COMPARE_DIR}/traj_${temp}.xtc" \
        -o "${COMPARE_DIR}/final_stable_${temp}.pdb" \
        -dump 50000

    # 3. Run H-bond analysis (10 ns - 50 ns)
    echo "  Running hbond analysis (10-50ns) to hbond_num_${temp}.xvg"
    echo "${HBOND_GROUP_SELECT}" | ${GMX} hbond \
        -s "${COMPARE_DIR}/md_${temp}.tpr" \
        -f "${COMPARE_DIR}/traj_${temp}.xtc" \
        -num "${COMPARE_DIR}/hbond_num_${temp}.xvg" \
        -b 10000 -e 50000 # Start at 10000 ps, end at 50000 ps

done

echo "=================================================="
echo "Preparation of analysis files completed."
echo "Selected files are in: ${COMPARE_DIR}"
echo "Next steps:"
echo " 1. Use Step5b_Plot_Hbond_Distribution.R to plot H-bond data from ${COMPARE_DIR}"
echo " 2. Use the PDB files (final_stable_*.pdb) in ${COMPARE_DIR} for PISA, DSSP, PROCHECK, and PyMOL RMSD analysis."
echo "=================================================="

exit 0