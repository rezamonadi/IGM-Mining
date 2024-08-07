#!/bin/bash -l

sbatch <<EOT
#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=80gb
#SBATCH --job=$1
#SBATCH --time=200:00:00
#SBATCH -p batch
#SBATCH --output=$1.out
#SBATCH --mail-user=reza.monadi@outlook.com
#SBATCH --mail-type=ALL

# loading matlab module
module load matlab
# run matlab script

# voigtPrep : 1 if you need making the mexa64 file 0 if you don't want to 
# voigtType: 0 --> Roman instrumental broadeing 1 --> With sigma pixel 
# maskType  0 --> masking all 1, masking most probable
# priorType 0 --> Roman, 1 --> From C13, 2 --> Singlet ~ a^n Doublet	
# MaskingProb 
#    0 --> Most probable model is masked 
#    1 --> the model with P larger than the sum of rest is masked
#    5 --> model with P>0.5 is masked
matlab -nodesktop -nosplash -r "cores = 32;\
                                num_C4_samples=10000;\
                                cataloging = 0;\
                                preloading = 0;\
                                learning   = 0;\
                                sampling   = 0;\
                                plotting   = 0;\
                                processing = 1;\
                                statTesting =0;\
                                dv_mask    = 350;\
                                HPCC = 1;\
                                voigtPrep = 0;\
                                maskType = 1;\
                                priorType = 1;\
                                OccamRazor = 0;\
                                partitioning = 0;\
                                MaskingProb = 0;\
                                saving=1;\
                                run_script_hpcc_dr7;\
                                exit"
                                
echo "----"
hostname

EOT
