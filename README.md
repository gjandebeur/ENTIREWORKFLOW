# ENTIREWORKFLOW
full start to finish workflow for basecalling and minimapping for rna002 data


### Some important coding tips to help out, feel free to reach out w/ any questions.

### Some of these may be redundant but unsure of coding experience

**adding a "#" sign in front makes it a comment, so anything following that would be a comment to explain. I've got no clue why Sbatch has a # in front of it, but keep and ignore it.**

**If any path ends in a "/" it means its a directory and contains files, if just .txt or any format then its obv a file.**
**Another important thing is using "" around files. You can choose to use or not use, I prefer using them because it helps visually, but just stay consistent because if not things get incredibly messy.**

**Directories are very important! I recommend making a batch folder to put your scripts into, by this command
"mkdir -p "ourdisk/hpc/rnafold/ADD_UR_USERNAME/dont_archive/batch/"
You can also use "cd /" to pull yourself back into your root directory and then do "cd /ourdisk/hpc/rnafold/" for the labs data and then just (cd + whatever directory youre wanting to go into).**
*Remember though that you cant run a script without being in that specific directory, and if you dont use absolute paths then the outputs go to that specific directory too*

Starting off here's your start to any slurm submitted job. (your sbatch file). Remove comments if wanted on script but shouldn't affect code itself.


    #!/bin/bash
    #SBATCH --partition=rnafold                       #where it runs (use "rnafold" or "normal")
    #SBATCH --job-name=12samples_dorado                     #name your job 
    #SBATCH --output=dorado1_ALLSAMPLES_output.txt              #keep _output at end but rename
    #SBATCH --error=dorado1_ALLSAMPLES_debug.txt          #keep _debug at end but u can change name
    #SBATCH --cpus-per-task=2   #use as written. if u need more power/gpu then do "gpus-per-node=1" 
    #SBATCH --time=24:00:00          #this just sets to 1 day  
    #SBATCH --mem=16G               #rn its set to 16gb, i rarely go above 32gb for big jobs

**I've installed all the software you will need so just copy and paste the path I provide as given and it'll work (lmk if its prevented by permissions or anything)**

**The first step is to run dorado to basecall your data (takes fast5 and reads as A,C,T,G). The software can't write "U" so every call will say "T" instead but it means uracil.**

**For any script youll add that sbatch to top, MAKING SURE that #!/bin/bash is on line 1. I've spent hours before looking at why my code wont run to find out this is why**

**The following modules need to be loaded to run dorado, these are dependencies that are (cant think of a better analogy here) basically like modifications to a gene, necessary and won't work without them. 
Just copy paste this chunk after the sbatch part but before the dorado script.**

        module load GCC
        module load PyTorch
        module load FlexiBLAS/3.3.1-GCC-12.3.0  
        module load FFmpeg/4.4.2-GCCcore-11.3.0 
        module load HTSlib
        module load protobuf



**To run dorado youll follow this formatting 
for reference, adding "\" to the end of the line tells your system that the command isn't done being written yet. So basically if you have extra spaces or anything other than this specific thing at the end of the line, your code wont work. So for every command I only add 1 space between and change lines by pressing enter only.**


#theres a slight chance that youre gonna get an error about formatting, if so its because the model needs to be placed somewhere else, just chatgpt the chunk of code and your error and you should be able to debug it.

        "/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-0.8.2/dorado-0.8.2-    linux-x64/bin/dorado" \ 
        basecaller \ 
        "/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/dorado-0.8.2/dorado-0.8.2-linux-x64/models/rna002_70bps_hac@v3/" \
        --min-qscore 10 \
        --emit-fastq > "/ADDYOURUSERHERE/ANDOUTPUT/PATH.fastq"  

### if you ever chose to do modifications then u have to remove the "--emit-fastq" tag and change the output to be a ".bam", as fastq cant carry the mod data over. Theres a couple more fancy steps to this so lmk if you decide to do mods.
        
**I prefer to use absolute paths because its way easier to see where the outputs going. make sure to have the ">" as it says redirect this entire command into this output.**
        
**qscore of 10 only keeps calls >=90% calls, can reduce to ~6-8 if it significantly reduces your data amount.**


### The next step is minimap2, where you take the fastq output and align it to your reference. I'll attach the path for the human transcriptome/genome, but I don't have the viral reference.
    
            "/ourdisk/hpc/rnafold/gjandebeur/dont_archive/software/marginAlign/submodules/minimap2/minimap2" \
        -ax splice -uf -k14 -y --secondary=no \     #keep these settings. its importnat
        "/ourdisk/hpc/rnafold/gjandebeur/dont_archive/reference/GRCh38.fa" \ #this is your ref file (MUST BE .fasta or .fa) 
        "/this/is/your/dorado/output/that/your/mapping.fastq" > "           "/this/is/your/output/alignedsplice.sam    

### If you would rather use the transcriptome I'll attach the path to the updated version, just replace the genome part with it. its already been indexed too (if needed just do samtools index before any .fa to index it 
# path : "/ourdisk/hpc/rnafold/gjandebeur/dont_archive/reference/gencode.v47.transcripts.fa"

### The output has to stay .sam but u can change the name to whatever, I prefer naming it alignedsplice to show that its aligned and has splice variants (just improves mapping slightly but unimportant).

### Next step you have to sort your reads, just makes it easier for computer to understand and converts format to a computer-readable one only.

**output goes first in this (the -o flag), and just add the input directly on end with only a space after the output, no other characters.
output must be .bam and input is the .sam**

        samtools sort -o "/THIS/IS/OUTPUT/FIRST/SORTED.BAM" "/This/is/input/alignedsplice.sam" 

        samtools flagstat  "/THIS/IS/OUTPUT/FIRST/SORTED.BAM"  

## Flagstat shows your stats that you want! (total aligned #, mapped #, and map %). This is the info I think you'll need comparing the viral reference to using genome or transcriptome.
        
