#!/bin/bash

# Exit script immediately if any command fails.
set -e

# --- Configuration ---
BASE_DIR=$(pwd)                             # Define base directory as current working directory
PDB_DIR="${BASE_DIR}/pdb_file"              # Path to PDB files
PARAM_DIR="${BASE_DIR}/MD_Parameter"        # Path to MDP parameter files
RESULT_DIR="${BASE_DIR}/MD_Result"          # Path where results will be stored

INPUT_PDB_MODELLER="${PDB_DIR}/2a5d_modeller.pdb" # Input PDB after MODELLER processing
PROTEIN_PDB="${PDB_DIR}/protein.pdb"        # Output PDB with only protein chains A+B

TEMPERATURES=("280K" "300K" "320K")         # Temperatures to simulate
REPLICATES=(1 2 3)                          # Replicate numbers for each temperature

GMX="gmx"                                   # GROMACS executable command
PYMOL="pymol"                               # PyMOL executable command

# GROMACS selections for non-interactive input (CHECK THESE FOR YOUR SETUP!)
PDB2GMX_FF_SELECT="15"                      # Force field selection (GROMOS96 54a7)
GENION_GROUP_SELECT="13"                    # Group to replace with ions (SOL)
RMSD_GROUP_SELECT="4 0"                     # Group for RMSD calculation (Backbone)
GYRATE_GROUP_SELECT="1"                     # Group for Rg calculation (Protein)
TRJCONV_CENTER_SELECT="1 0"                 # Group for centering trajectory (Protein)
# HBOND_GROUP_SELECT="protein protein"      # Groups for hbond analysis (optional, uncomment below)

# --- 1. Prepare Protein PDB (Extract Chains A+B) ---
# Use PyMOL command line to load PDB, select chains A+B, save as protein.pdb
"${PYMOL}" -c -d "load ${INPUT_PDB_MODELLER}, model; select protein, chain A+B; save ${PROTEIN_PDB}, protein; quit"

# --- 2. Create Base Result Directory ---
# Create the main directory to store all simulation results
mkdir -p "${RESULT_DIR}"

# --- 3. Loop Through Temperatures and Replicates ---
# Outer loop for each temperature
for TEMP in "${TEMPERATURES[@]}"; do
    # Inner loop for each replicate
    for REP in "${REPLICATES[@]}"; do

        # Define unique run name (e.g., 300K_1)
        RUN_NAME="${TEMP}_${REP}"
        # Define output directory path for this specific run, like 300K/300K_1
        OUTPUT_RUN_DIR="${RESULT_DIR}/${TEMP}/${RUN_NAME}"
        # Define parameter directory path for this temperature
        MDP_DIR="${PARAM_DIR}/${TEMP}"

        # Create the specific directory for this run's output
        mkdir -p "${OUTPUT_RUN_DIR}"
        # Change current directory to the run's output directory
        cd "${OUTPUT_RUN_DIR}"
        # Copy the prepared protein PDB file into the run directory
        cp "${PROTEIN_PDB}" .

        # Generate GROMACS topology (.top) and coordinate (.gro) files
        ${GMX} pdb2gmx -f protein.pdb -o protein.gro -p topol.top -water spce -ignh <<< "${PDB2GMX_FF_SELECT}"
        # Create simulation box (cubic, 1.0 nm margin)
        ${GMX} editconf -f protein.gro -o protein_box.gro -c -d 1.0 -bt cubic
        # Fill the box with water (SPCE model)
        ${GMX} solvate -cp protein_box.gro -cs spc216.gro -o system_solv.gro -p topol.top
        # Prepare run input file (.tpr) for adding ions
        ${GMX} grompp -f "${MDP_DIR}/ions.mdp" -c system_solv.gro -p topol.top -o ions.tpr -maxwarn 1
        # Add Na+/Cl- ions to neutralize the system
        ${GMX} genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral <<< "${GENION_GROUP_SELECT}"
        # Prepare run input file (.tpr) for energy minimization
        ${GMX} grompp -f "${MDP_DIR}/minim.mdp" -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
        # Execute energy minimization (using GPU for non-bonded)
        ${GMX} mdrun -v -deffnm em -nb gpu -ntomp 8 -ntmpi 2
        # Prepare run input file (.tpr) for NVT equilibration
        ${GMX} grompp -f "${MDP_DIR}/nvt.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
        # Execute NVT equilibration
        ${GMX} mdrun -deffnm nvt -nb gpu -ntomp 8 -ntmpi 2
        # Prepare run input file (.tpr) for NPT equilibration
        ${GMX} grompp -f "${MDP_DIR}/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 1
        # Execute NPT equilibration
        ${GMX} mdrun -deffnm npt -nb gpu -ntomp 8 -ntmpi 2
        # Prepare run input file (.tpr) for production simulation
        ${GMX} grompp -f "${MDP_DIR}/md.mdp" -c npt.gro -t npt.cpt -p topol.top -o md_${RUN_NAME}.tpr -maxwarn 1

        # Execute production MD simulation
        ${GMX} mdrun -deffnm md_${RUN_NAME} -nb gpu -ntmpi 2 -ntomp 8 #-cpi md_${RUN_NAME}.cpt -append


        # Remove PBC and center the trajectory on the protein
        ${GMX} trjconv -s md_${RUN_NAME}.tpr -f md_${RUN_NAME}.xtc -o md_${RUN_NAME}_noPBC.xtc -pbc mol -center <<< "${TRJCONV_CENTER_SELECT}"
        # Calculate RMSD vs starting structure (after NPT)
        ${GMX} rms -s md_${RUN_NAME}.tpr -f md_${RUN_NAME}_noPBC.xtc -o rmsd_${RUN_NAME}.xvg -tu ns <<< "${RMSD_GROUP_SELECT}"
        # Calculate RMSD vs minimized structure (closer to initial PDB)
        ${GMX} rms -s em.tpr -f md_${RUN_NAME}_noPBC.xtc -o rmsd_xtal_${RUN_NAME}.xvg -tu ns <<< "${RMSD_GROUP_SELECT}"
        # Calculate Radius of Gyration
        ${GMX} gyrate -s md_${RUN_NAME}.tpr -f md_${RUN_NAME}_noPBC.xtc -o gyrate_${RUN_NAME}.xvg <<< "${GYRATE_GROUP_SELECT}"
        # Optional: Calculate hydrogen bonds between 10ns and 50ns
        # ${GMX} hbond -s md_${RUN_NAME}.tpr -f md_${RUN_NAME}_noPBC.xtc -num hbond_num_${RUN_NAME}.xvg -b 10000 -e 50000 <<< "${HBOND_GROUP_SELECT}"
        # Optional: Extract the final frame as a PDB file
        # ${GMX} editconf -f md_${RUN_NAME}.gro -o md_${RUN_NAME}_final.pdb

        echo "Simulation ${RUN_NAME} finished processing."
        # Change directory back to the base directory for the next iteration
        cd "${BASE_DIR}"

    done # End replicate loop
done # End temperature loop

echo "All simulation commands executed."
echo "Results should be in: ${RESULT_DIR}"
exit 0