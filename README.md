# ENTIREWORKFLOW
full start to finish workflow for basecalling, minimapping, and modkit for rna002 data


# adding a "#" sign in front makes it a comment, so anything following that would be a comment to explain. I've got no clue why Sbatch has a # in front of it, but keep and ignore it.

Starting off heres your start to any slurm submitted job. (your sbatch file)


#!/bin/bash
#SBATCH --partition=rnafold                       #where it runs (use "rnafold" or "normal")
#SBATCH --job-name=12samples_dorado                     #name your job 
#SBATCH --output=dorado1_ALLSAMPLES_output.txt              #keep _output at end but rename
#SBATCH --error=dorado1_ALLSAMPLES_debug.txt          #keep _debug at end but u can change name
#SBATCH --cpus-per-task=2   #use as written. if u need more power/gpu then do "gpus-per-node=1" 
#SBATCH --time=96:00:00           
#SBATCH --mem=64G  
