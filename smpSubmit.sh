#!/bin/bash -l

sbatch <<EOT
#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=128gb
#SBATCH --job-name="Code0"
#SBATCH --time=2:00:00
#SBATCH -p short
#SBATCH --output="priorFix.out"
#SBATCH --mail-user=rmona003@ucr.edu
#SBATCH --mail-type=ALL

# loading matlab module
module load matlab
# run matlab script
matlab -nodesktop -nosplash -r "parpool('local', 32); run_script_hpcc"
echo "----"
hostname
EOT
