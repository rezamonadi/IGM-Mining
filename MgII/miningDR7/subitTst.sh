#!/bin/bash -l

sbatch <<EOT
#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --job=test
#SBATCH --time=2:00:00
#SBATCH --output=$1.out
#SBATCH --mail-user=rmonad
#SBATCH --mail-type=ALL

# loading matlab module
# module load matlab
# run matlab script

# voigtPrep : 1 if you need making the mexa64 file 0 if you don't want to 
# voigtType: 0 --> Roman instrumental broadeing 1 --> With sigma pixel 
# maskType  0 --> masking all 1, masking most probable
# priorType 0 --> Roman, 1 --> From C13, 2 --> Singlet ~ a^n Doublet		
matlab -nodesktop -nosplash -r " tst.m;  exit"
                                
echo "----"
hostname

EOT
