#!/bin/bash -l

sbatch <<EOT
#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200gb
#SBATCH --job=$1
#SBATCH --time=2:00:00
#SBATCH -p  short
#SBATCH --output=$1.out
#SBATCH --mail-user=rmona003@ucr.edu
#SBATCH --mail-type=ALL

# loading matlab module
module load matlab/R2021b
# run matlab script

# voigtPrep : 1 if you need making the mexa64 file 0 if you don't want to 
# voigtType: 0 --> Roman instrumental broadeing 1 --> With sigma pixel 
# maskType  0 --> masking all 1, masking most probable
# priorType 0 --> Roman, 1 --> From C13, 2 --> Singlet ~ a^n Doublet		
matlab -nodesktop -nosplash -r "cores = 24;\
                                num_C4_samples = 10000;\
                                cataloging = 0;\
                                preloading = 0;\
                                learning   = 0;\
                                sampling   = 0;\
                                plotting   = 0;\
                                processing = 0;\
                                merging    = 0;\
                                EWer       =0;\
                                pltP       =1;\
                                CredInt    =0;\
                                dv_mask    = 350;\
                                HPCC = 1;\
                                voigtPrep = 0;\
                                maskType = 1;\
                                priorType = 1;\
                                ind_S=1;\
                                num_quasars =5000;\
                                run_script_hpcc_dr12;\
                                exit"
                                
echo "----"
hostname

EOT
