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

**The following modules need to be loaded to run dorado, these are dependencies that are kinda like R packages you would load in (but in UNIX).
Just copy paste this chunk after the sbatch part but before the dorado script.**

        module load GCC
        module load PyTorch
        module load FlexiBLAS/3.3.1-GCC-12.3.0  
        module load FFmpeg/4.4.2-GCCcore-11.3.0 
        module load HTSlib
        module load protobuf

        module load SAMtools/1.16.1-GCC-11.3.0 #this ones for minimap but add here too



**To run dorado youll follow this formatting 
for reference, adding "\\" to the end of the line tells your system that the command isn't done being written yet. So basically if you have extra spaces or anything other than this specific thing at the end of the line, your code wont work. So for every command I only add 1 space between and change lines by pressing enter only.**


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

**If you would rather use the transcriptome I'll attach the path to the updated version, just replace the genome part with it. its already been indexed too but if needed ever just do samtools index before any .fa to index it**
### path : "/ourdisk/hpc/rnafold/gjandebeur/dont_archive/reference/gencode.v47.transcripts.fa"

### The output has to stay .sam but u can change the name to whatever, I prefer naming it alignedsplice to show that its aligned and has splice variants (just improves mapping slightly but unimportant).

### Next step you have to sort your reads, just makes it easier for computer to understand and converts format to a computer-readable one only.

**output goes first in this (the -o flag), and just add the input directly on end with only a space after the output, no other characters.
output must be .bam and input is the .sam**

        samtools sort -o "/THIS/IS/OUTPUT/FIRST/SORTED.BAM" "/This/is/input/alignedsplice.sam" 

        samtools flagstat  "/THIS/IS/OUTPUT/FIRST/SORTED.BAM"  

## Flagstat shows your stats that you want! (total aligned #, mapped #, and map %). This is the info I think you'll need comparing the viral reference to using genome or transcriptome.





# Counting Transcript Abundance and Plotting Differential Expression Volcano Plots

Now for the complicated part

after aligning to the reference, the goal is to count how many copies of each gene is expressed in each sample. This will be done using salmon, which can be installed directly into a conda environment for you (I'll add the exact code to run for this below too)

        # Install Salmon in a conda environment
        conda create -n salmon_env -c bioconda salmon
        conda activate salmon_env

So this part can be done on the login node, and add this little chunk of code to your data file, under the #SBATCH part but above your code.

        
        echo "linking the Conda environment to $ENV_PATH..."
        source /opt/oscer/software/Mamba/23.1.0-4/etc/profile.d/conda.sh
        conda activate salmon_env


Now for the salmon command itself to quantify. It's important that the transcriptome is used here as reference because it needs to be able to count isoforms. (the genome has higher accuracy in minimap step than transcript but that shouldn't change results for this software)

Salmon likes to put all results into a directory, and then the next step pulls directly from it, so run this command but change to whatever directory your using. End the file name with the sample name you're running salmon on and don't add a "/" to the very end, this way the output will separate by sample and work much better for later steps.

        #mkdir -p "/ourdisk/hpc/rnafold/UR_USER/dont_archive/salmon/BYSAMPLE"


**The following code you shouldn't need to run if you use my reference transcript file, but for the viral genome (if that ends up being used here) you'll use this.**

        salmon index -t "/ourdisk/hpc/rnafold/gjandebeur/dont_archive/reference/gencode.v47.transcripts.fa" -i "/ourdisk/hpc/rnafold/gjandebeur/dont_archive/reference/gencode.v47.transcripts_index" -k 31 


I believe the part after "-i" is saying where the script can write all the little files it creates to index the reference.

**this is the important step**

        salmon quant -t "/ourdisk/hpc/rnafold/gjandebeur/dont_archive/reference/gencode.v47.transcripts.fa" -l A -a
        "/input/ALIGNED/SORTED/data.bam" -p 8 --seqBias --gcBias -o "/this/is/output/directory/by/sample/"

Online documentation had stuff on --seqBias and --gcBias so included, but I'm not 100% sure if they're necessary

Now for running R on OSCER, just write a normal unix script, and add "Rscript (file)" in that unix script. I've heard of people shutting down nodes from running python or R directly on the command line so I would recommend against.

**No clue why but its very difficult to get the R packages working on OSCER, I use alot of ChatGPT when trying to do this and I believe there should be some modules already on oscer that have necessary packages (try module avail R).** 


## This is copy pasted from my script but theres only like 3 parts youll need to change, ive already combined the DEseq2 and volcano plot step into one so you should only have to update the parts of the following code that I have bolded.
   
    #!/usr/bin/env Rscript
    
    library(DESeq2)
    library(tximport)
    library(ggplot2)


        # --- Output directory ---
        outdir <- **"/ourdisk/hpc/rnafold/gjandebeur/dont_archive/batch"**
        dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

        cat("Loaded libraries and created output directory.\n")

        # --- Setup: Samples and files ---
        samples <- c(**"24r1a", "24r2a", "32r1a", "32r2a", "36r1a", "36r2a"**)
    salmon_dir <- "/input/directory/from/salmon/"
files <- file.path(salmon_dir, samples, "quant.sf")
names(files) <- samples

        # Sample metadata
        coldata <- data.frame(
      row.names = samples,
        **condition = factor(c("control", "CSE", "control", "CSE", "control", "CSE")** )
)

        cat("Prepared sample metadata and file paths.\n")

        # --- Import counts ---
        txi <- tximport(files, type = "salmon", txOut = TRUE)
        cat("Imported quant.sf files with tximport.\n")

    # --- Filter low count transcripts ---
    keep <- rowSums(txi$counts) > 10
    txi$counts <- txi$counts[keep, ]
    txi$abundance <- txi$abundance[keep, ]
    txi$length <- txi$length[keep, ]

    cat("Filtered low count transcripts.\n")

    #    --- Create DESeq2 dataset ---
    dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~condition)
    dds <- dds[rowSums(counts(dds)) > 0, ]  # Extra filter

    cat("Created DESeq2 dataset.\n")

    # --- Run DESeq2 ---
    dds <- DESeq(dds)
    res <- results(dds, contrast = c(**"condition", "CSE", "control")**)

        cat("Ran DESeq2 and extracted results.\n")

    # --- Handle p-values ---
    if (all(is.na(res$padj))) {
      stop("All adjusted p-values are NA. Check input or filtering.")
    }
    min_padj <- min(res$padj[res$padj > 0], na.rm = TRUE)
    res$padj[res$padj == 0] <- min_padj

    # --- Calculate -log10 adjusted p-value ---
    res$neg_log10_padj <- -log10(res$padj)

    # --- Add significance flag with relaxed cutoff ---
    pval_cutoff <- 0.075   # Change this as needed (e.g., 0.05, 0.1)
    res$sig <- ifelse(res$padj <= pval_cutoff & abs(res$log2FoldChange) >= 1, "Significant",     "Not Significant")

    # Make 'sig' a factor with explicit levels to avoid legend artifacts
    res$sig <- factor(res$sig, levels = c("Significant", "Not Significant"))

    # --- Extract gene names (6th field from pipe-separated rownames) ---
    res$gene <- sapply(strsplit(rownames(res), "\\|"), function(x) x[6])

    # --- Subset for gene labels that are not empty or NA ---
    sig_genes <- subset(res, sig == "Significant" & !is.na(gene) & gene != "")

    # --- Volcano plot ---
    p <- ggplot(as.data.frame(res), aes(x = log2FoldChange, y = neg_log10_padj, color =     sig)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_manual(
        values = c("Significant" = "firebrick", "Not Significant" = "lightgrey"),
    breaks = c("Significant", "Not Significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "blue") +
  labs(
    title = paste0("Volcano Plot of All Genes (p-adj ≤ ", pval_cutoff, ")"),
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  ) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2), limits = c(-6, 6)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +
  theme_minimal(base_size = 14) +
  geom_text(
    data = sig_genes,
    aes(label = gene),
    size = 2.5,
    vjust = -0.5,
    check_overlap = TRUE
  )

    # --- Define output directory ---
    outdir <-**"/ourdisk/hpc/rnafold/gjandebeur/dont_archive/batch"**
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    # --- Save plot ---
    ggsave(
      filename = file.path(outdir, **"nonsmokers_volcano_plot.png")**,
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

    # --- Save full results ---
    cat("Saving full results...\n")
    tryCatch({
      write.csv(as.data.frame(res), file = file.path(outdir,     **"nonsmokers_DESeq2_results.csv"**))
  cat("Full results saved.\n")
}, error = function(e) {
  cat("Error saving full results:\n")
  print(e)
})

    # --- Save significant results only ---
    cat("Saving significant results...\n")
    tryCatch({
      sig_res <- subset(res, sig == "Significant")
      write.csv(as.data.frame(sig_res), file = file.path(outdir,             **"DESeq2_results_significant.csv"**))
      cat("Significant results saved.\n")
}, error = function(e) {
  cat("Error saving significant results:\n")
  print(e)
})

    cat("All steps completed.\n")




