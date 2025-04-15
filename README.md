# MD Simulation Workflow: Effect of Temperature on CTA1-ARF6 Complex

This repository contains scripts for running and analyzing GROMACS Molecular Dynamics (MD) simulations to investigate the effect of temperature on the stability of the Cholera Toxin A1 subunit (CTA1) - human ARF6 protein complex (corresponding to Miniproject 7).

## Prerequisites

Running this workflow requires the following software:

*   **GROMACS:** Version 2024.3 (or other compatible version)
*   **PyMOL:** Required by `Step1_Run_All_Simulations.sh` to prepare the initial PDB file.
*   **R:** Version 4.3.1 (or other compatible version)
*   **R Packages:** `tidyverse`, `fs`, `ggpubr` (for plotting and analysis). Optionally `car`, `dunn.test` might be needed if running the statistical tests in Step 6.
*   **(Preparation Stage)** The input PDB file (`pdb_file/2a5d_modeller.pdb`) is assumed to have been pre-processed using **MODELLER** (missing residues completed).

## Repository Structure

.
├── README.md # This documentation file
├── .gitignore # Specifies files/directories for Git to ignore
│
├── pdb_file/ # Input PDB structures
│ └── 2a5d_modeller.pdb # Starting structure used by scripts (post-MODELLER)
│
├── MD_Parameter/ # GROMACS input parameter (.mdp) files
│ ├── 280K/
│ ├── 300K/
│ └── 320K/
│
├── Step1_Run_All_Simulations.sh # Step 1: Runs all 9 simulations & initial analysis
├── Step2_Monitoring_Equilibration.R # Step 2: (Optional) Plots equilibration for one run
├── Step3_Plot_Select_Rg.R # Step 3: Plots average Rg, calculates & selects most stable replicate
├── Step4_Plot_Selected_Trajectories_Manual.R # Step 4: (Optional/Manual) Plots trajectories for manually chosen replicates
├── Step5a_Prepare_Analysis_Files.sh # Step 5: Prepares files for final analysis (extracts, Hbond)
└── Step6_Plot_Hbond_Distribution.R # Step 6: Plots H-bond distribution & stats


## Workflow

Please execute the following scripts in order. Ensure Bash scripts (`.sh`) have execute permissions (`chmod +x *.sh`).

1.  **Step 1: Run All Simulations**
    *   **Command:** `bash Step1_Run_All_Simulations.sh`
    *   **Purpose:** Calls PyMOL to prepare the `protein.pdb` file, then runs all 9 GROMACS simulations (preparation, energy minimization, NVT/NPT equilibration, 50ns production simulation), and performs initial trajectory analyses (RMSD, Rg).
    *   **Output:** Creates `pdb_file/protein.pdb` and the `MD_Result/` directory. `MD_Result/` contains the complete output files for each simulation run (e.g., `280K/280K_1/` containing `.tpr`, `.xtc`, `.log`, `.edr`, `.gro`, `.xvg` files, etc.). **Note: This step is very time-consuming (potentially >10 days)!**

2.  **Step 2: (Optional) Monitor Equilibration for a Single Run**
    *   **Command:** `Rscript Step2_Monitoring_Equilibration.R`
    *   **Purpose:** Generates equilibration monitoring plots (density, potential energy, pressure) for a single simulation run hardcoded within the script (example uses `280K_2`). Primarily for debugging or checking individual runs.
    *   **Output:** Generates a `.png` image file within the corresponding `MD_Result/` subdirectory (e.g., `MD_Result/280K/280K_2/equilibration_monitor_280K_Rep_2.png`).

3.  **Step 3: Plot Average Rg & Select Stable Replicates**
    *   **Command:** `Rscript Step3_Plot_Select_Rg.R`
    *   **Purpose:** Reads the `gyrate_*.xvg` files from all 9 simulations, plots the average Rg vs. time for each temperature (with standard deviation), then calculates the Rg variance for each replicate during the production phase and outputs a CSV file identifying the most stable replicate (lowest variance) for each temperature.
    *   **Output:** Creates the `Analysis_Output/` directory. Generates `Analysis_Output/Rg_Average_vs_Time.png` plot and `Analysis_Output/most_stable_replicates_Rg_variance.csv` table. This CSV file is input for Step 5a.

4.  **Step 4: (Optional/Manual) Plot Selected Trajectories**
    *   **Command:** (First, **manually edit** the `selected_replicate_numbers` variable inside the script) Then run `Rscript Step4_Plot_Selected_Trajectories_Manual.R`
    *   **Purpose:** Plots the RMSD and Rg trajectories for the specific replicate numbers you manually choose for each temperature within the script.
    *   **Output:** Generates the `Analysis_Output/Manually_Selected_Trajectories_Plot.png` plot.

5.  **Step 5: Prepare Final Analysis Files**
    *   **Command:** `bash Step5a_Prepare_Analysis_Files.sh`
    *   **Purpose:** Reads the `most_stable_replicates_Rg_variance.csv` file generated in Step 3 to identify the most stable replicate for each temperature. Then:
        *   Creates the `Compare_Selected/` directory.
        *   Copies the selected replicate's `.tpr` and `_noPBC.xtc` files into `Compare_Selected/`, renaming them (e.g., `md_280K.tpr`, `traj_280K.xtc`).
        *   Extracts the final frame (50 ns) from the selected trajectory and saves it as a `.pdb` file (e.g., `final_stable_280K.pdb`) in `Compare_Selected/`.
        *   Runs `gmx hbond` analysis (10-50 ns) on the selected trajectory, saving the output (e.g., `hbond_num_280K.xvg`) to `Compare_Selected/`.
    *   **Output:** Creates and populates the `Compare_Selected/` directory.

6.  **Step 6: Plot H-bond Distribution**
    *   **Command:** `Rscript Step6_Plot_Hbond_Distribution.R`
    *   **Purpose:** Reads the `hbond_num_*.xvg` files generated by Step 5a from the `Compare_Selected/` directory, plots the distribution of hydrogen bond counts as a box plot, and potentially performs statistical tests and adds significance annotations.
    *   **Output:** Generates the `Analysis_Output/Hbond_Distribution_Plot.png` plot.

## Final Analysis Files

After completing all non-optional steps, the `final_stable_*.pdb` files located in the `Compare_Selected/` directory are the final structures intended for submission to PISA, DSSP, PROCHECK servers, and for RMSD comparison in PyMOL. The `hbond_num_*.xvg` files in `Compare_Selected/` contain the data used for the final H-bond box plot.
