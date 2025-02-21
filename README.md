# FoldX Simulation and Analysis for proteins

## What is FoldX?
FoldX is a protein design algorithm that uses an empirical force field to determine the energetic effect of point mutations and the interaction energy of protein complexes, including protein-DNA interactions. It can also mutate protein and DNA side chains using a probability-based rotamer library, while exploring alternative conformations of surrounding side chains. Key features of FoldX include:

   - **Mutational Analysis**: Predict the effect of mutations on protein stability.
   - **Protein Design**: Design proteins with improved stability or novel functions.
   - **Repair of Protein Structures**: Fix small structural issues in protein PDB files.

For more information, visit https://foldxsuite.crg.eu/.
<b>*"The FoldX Suite is available through academic and commercial licenses. Academic and non-profit research institutions can obtain an academic license, while companies and other users can obtain a commercial license".*</b>

## File System

```
.
├── Data
│   └── Processed_Data
│       ├── DS - All the mutations list files for the protein are stored here.
│       │   └── rbd - Name of the directory.
│       │       ├── individual_list_rbd_domain_00.txt - Mutation list file that is generated after running foldx_mutations.py.
│       │       └── split_files - The mutation list file is split into several other files, for parallel computing.
│       │           ├── individual_list_rbd_domain_00_00.txt
│       │           ├── individual_list_rbd_domain_00_01.txt
│       │           ├── individual_list_rbd_domain_00_02.txt
│       │           ├── individual_list_rbd_domain_00_03.txt
│       │           └── individual_list_rbd_domain_00_04.txt
│       ├── HM - All the protein structures are stored here in .pdb format.
│       │   └── rbd_domain.pdb
│       └── MD - The analysis and the mutated protein is stored here.
│           └── MD_rbd
│               ├── foldx_results
│               │   ├── 00 - The results are stored in seprate directories according to the number on the split file. These folders contain a lot of pdb files, only a few are mentioned here. Similar files are present in all the directories.
│               │   │   ├── Average_rbd_domain_Repair.fxout
│               │   │   ├── Dif_rbd_domain_Repair_00.fxout
│               │   │   ├── PdbList_rbd_domain_Repair.fxout
│               │   │   ├── Raw_rbd_domain_Repair.fxout
│               │   │   ├── rbd_domain_Repair_10_0.pdb
│               │   │   ├── rbd_domain_Repair_10_1.pdb
│               │   │   ├── rbd_domain_Repair_10_2.pdb
│               │   │   ├── WT_rbd_domain_Repair_10_0.pdb
│               │   │   ├── WT_rbd_domain_Repair_10_1.pdb
│               │   │   ├── WT_rbd_domain_Repair_10_2.pdb
│               │   ├── 01
│               │   │   ├── Average_rbd_domain_Repair.fxout
│               │   │   ├── Dif_rbd_domain_Repair_01.fxout
│               │   │   ├── PdbList_rbd_domain_Repair.fxout
│               │   │   ├── Raw_rbd_domain_Repair.fxout
│               │   │   ├── rbd_domain_Repair_10_0.pdb
│               │   │   ├── rbd_domain_Repair_10_1.pdb
│               │   │   ├── rbd_domain_Repair_10_2.pdb
│               │   │   ├── WT_rbd_domain_Repair_10_0.pdb
│               │   │   ├── WT_rbd_domain_Repair_10_1.pdb
│               │   │   ├── WT_rbd_domain_Repair_10_2.pdb
│               │   ├── 02
│               │   ├── 03
│               │   ├── 04
│               ├── rbd_domain_Repair.fxout - This file is created by FoldX to calculate the different types of energies assocaited with the protein. 
│               ├── rbd_domain_Repair.pdb - This file is created by FoldX that represents a clean protein structure.
│               └── temp
│                   ├── df_rbd_domain_foldx.csv - This is the csv file generated after running the foldx_analysis.py file. It cotains about the nature of stability of the protein.
│                   ├── exception_cases.txt - Since each mutations of the protein are run three times, there might be some cases that do not show the same result. These are stored in this file. Either re-run them or remove them from analysis.
│                   └── foldx_commands_file.txt - These are the command used to run FoldX simultations.
├── LICENSE
├── Program
│   ├── Bin
│   │   └── FoldX_2 - Dowload FoldX and store it in this directory.
│   └── Script
│       ├── foldx_analysis.py - Perform mutational analysis.
│       ├── foldx_mutations.py - Get the list of mutations.
│       └── foldx_simulation_file.py - Run FoldX simulation parallely.
├── README.md
├── requirements.txt

## Scripts Syntax and Execution

This document outlines the syntax and execution steps for running FoldX simulations and analysis.

**Scripts:**

*   **`foldx_mutations.py`**: Generates a list of mutations for a given PDB file.
    ```bash
    python3 /Program/Script/foldx_mutations.py -i <pdb file path> -o <output mutation list path>
    ```

*   **`foldx_simulation_file.py`**: Runs FoldX simulations for a given PDB file and mutation list.
    ```bash
    python3 /Program/Script/foldx_simulation_file.py -i <path to FoldX directory> <pdb file path> <mutation list path> -o <output results path> <output mutated structures path>
    ```

*   **`foldx_analysis.py`**: Analyzes the FoldX simulation results.
    ```bash
    python3 /Program/Script/foldx_analysis.py -i <pdb file path> <path to split mutation list files> <path to mutated structures results> -o <output temporary directory path>
    ```

**Installation and Setup:**

1.  **Uncompress FoldX:**
    ```bash
    tar -xvf foldx.tar.gz
    ```


2.  **Create Conda Environment:**
    ```bash
    cd foldx  # Navigate to the FoldX directory (where environment.yml is located)
    conda env create --name foldx_env --file=environment.yml
    ```

3.  **Activate Environment:**
    ```bash
    conda activate foldx_env
    ```

**Example Workflow:**

1.  **Generate Mutation Lists:**
    ```bash
    python3 /Program/Script/foldx_mutations.py -i /Data/Processed_Data/HM/rbd_domain.pdb -o /Data/Processed_Data/DS/rbd/
    ```
    This will create mutation list files in the `/Data/Processed_Data/DS/rbd/` directory. 

2.  **Run FoldX Simulations (using `screen` for background processes):**
    ```bash
    screen  # Create a new screen session named
    conda activate foldx_env # Activate the conda environment inside the screen session
    python3 /Program/Script/foldx_simulation_file.py -i /Program/Bin/FoldX_2 /Data/Processed_Data/HM/rbd_domain.pdb /Data/Processed_Data/DS/rbd/individual_list_rbd_domain.txt -o /Data/Processed_Data/MD/MD_rbd /Data/Processed_Data/MD/MD_rbd/foldx_results
    # Detach from the screen session: Ctrl+A then D
    ```
    *   **Managing Screen Sessions:**
        *   `screen -ls`: List active screen sessions.
        *   `screen -r <screen name here>`: Reattach to the "foldx_sim" session.
        *   `Ctrl+A then D`: Detach from a screen session (keeps the process running).

    *   **Parallelization:** To run simulations on different CPUs.

    ```bash
        python3 Program/Script/foldx_simulation_file.py -i Program/Bin/FoldX Data/Processed_Data/HM/rbd_domain.pdb Data/Processed_Data/DS/rbd/individual_list_rbd_domain.txt -o Data/Processed_Data/MD/MD_rbd  Data/Processed_Data/MD/MD_rbd/foldx_results 
    ```

3.  **Run FoldX Analysis:**
    ```bash
    python3 /Program/Script/foldx_analysis.py -i /Data/Processed_Data/HM/rbd_domain.pdb /Data/Processed_Data/DS/rbd/split_files /Data/Processed_Data/MD/MD_rbd/foldx_results -o /Data/Processed_Data/MD/MD_rbd/temp
    ```

This updated example provides a more complete workflow, including installation, environment setup, parallelization and clearer explanations of the steps involved. Remember to adapt the paths and parameters to your specific needs. Working on updates and modfying code that can make this process more robust!
