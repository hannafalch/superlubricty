#!/bin/bash
#SBATCH --job-name=ch001          # Job name
#SBATCH --output=lammps_output.log     # Output log file
#SBATCH --error=lammps_error.log       # Error log file
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks-per-node=4            # Number of tasks per node (adjust as needed)
#SBATCH --time=10:00:00                # Time limit (hh:mm:ss) adjust this based on your job's needs
#SBATCH --partition=porelab         # Partition/queue (replace with appropriate partition name)

# Load necessary modules
module load eb/use.eb
module load GCC/11.2.0
module load OpenMPI/4.1.1
module load LAMMPS/23Jun2022-kokkos

# Run the LAMMPS job
srun lmp -in run.in
