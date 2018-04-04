################################################################################
###  Analysis of breast cancer risk loci for allele-specific miRNA binding   ###
##            TargetScan v7.1 [Agarwal et al., 2015] Application              ##
################################################################################

# Author: Ana Jacinta-Fernandes, <a46845@ualg.pt>
# Date: Feb 2018
# The Functional Genomics of Cancer Laboratory, PI: Ana-Teresa Maia
# Centre for Biomedical Research (CBMR), University of Algarve, Portugal
# <www.maialab.org>

################################################################################
#### 00. Preparation #### 

source("./functions-gwasloci-miRNAtargets-ts.R")

check.packages()  # Require R packages 
check.files()  # Check if all necessary files exist 
file.save()  # Create folder were all new files will be generated (If you want to create a specific folder then `file.save(fullpath = "./hello-world")`)

# sink(paste0(working_path, "/R_console_output.txt"), append = TRUE) # Save all console output to .txt file. To stop sinking, use 'closeAllConnections()'



################################################################################
#### 01. Settings ####

get.miR.expression(0) # Select 0:1  #  0 = disabled | 1 = uses miRmine [Panwar et al., 2017] miRNA expression levels for normal breast (SRX513286)[Zhu et al., 2014], miRBase build 20 



################################################################################
## 02. Proxy Search

# SNAP finds proxy SNPs based on linkage disequilibrium, physical distance and/or membership in selected commercial genotyping arrays [Johnson et al., 2008] and is available as an online tool at archive.broadinstitute.org/mpg/snap/ldsearch.php
# Extract proxy SNP data using a SNP vector as input in SNAP (r^2>=0.8, dataset: 1000 Genomes Pilot 1, CEU population, distance limit= 500kb) and upload output file: SNAPResults.txt

tagproxy <- get.proxy.SNAP.data(SNAP_output = "SNAPResults.txt")



################################################################################
#### 03. Ensembl query ####

hsapiens <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

dataset <- query.snp.ensembl(snp = tagproxy)  # <---  insert here snp of interest (Ex: snp = "rs8123521")

BioRelevant_miRNAarm_3UTR <- subset(dataset,
                                    dataset$consequence_type_tv == "3_prime_UTR_variant")

refsnp_id <- unique(BioRelevant_miRNAarm_3UTR$refsnp_id)

hgnc_symbol <- unique(BioRelevant_miRNAarm_3UTR$hgnc_symbol) 
hgnc_symbol  # Genes with SNPs at their 3'UTRs



################################################################################
#### 04. Ensembl query - hg19 coordinates ####

hg19_snp = useMart(biomart = "ENSEMBL_MART_SNP", host = "grch37.ensembl.org", 
                   path = "/biomart/martservice", dataset = "hsapiens_snp")   # Access SNP database from hg19

snp_coordinates_hg19 <- snp.coordinates.hg19(refsnp_id)



################################################################################
#### 05. UTR query ####

hgnc_symbol <- hgnc_symbol[!is.na(hgnc_symbol)]  # remove NA 

writeLines("\nIf the hgnc symbol is not found in the 'UTR sequences' file try using it's aliases. Ex: KMT5A is not found but it's alias SETD8 is.\n")

UTR_sequences_list <- lapply(hgnc_symbol, grep.file, filename = "UTR_Sequences.txt")
UTR_sequences <- do.call(rbind, UTR_sequences_list)
colnames(UTR_sequences) <- c("Ensembl.ID", "Gene.ID", "Gene.Symbol", "Species.ID", "UTR.sequence")

UTR_sequences_hsa <- subset(UTR_sequences,
                            UTR_sequences$Species.ID == "9606")  # Human species = 9606

# Get UTR sequences coordinates hg19

transcript_id <- unique(UTR_sequences_hsa$Ensembl.ID)

UTR_coordinates_hg19_list <- lapply(transcript_id, grep.file, filename = "TSHuman_7_hg19_3UTRs.gff")
UTR_coordinates_hg19 <- do.call(rbind, UTR_coordinates_hg19_list)
colnames(UTR_coordinates_hg19) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

coordinates_hg19_UTR_sequences_hsa <- merge(UTR_sequences_hsa, UTR_coordinates_hg19, 
                                            by.x = "Ensembl.ID", 
                                            by.y = "attribute")



################################################################################
#### 06. Prepare dataset for TargetScan ####

snp_hgnc.symbol <- unique(BioRelevant_miRNAarm_3UTR[,c("refsnp_id","hgnc_symbol")])  # combinations of snp::gene

snp_range <- merge(snp_coordinates_hg19, snp_hgnc.symbol, 
                   by.x = "refsnp_id", 
                   by.y = "refsnp_id")

snp_range <- merge(snp_range, coordinates_hg19_UTR_sequences_hsa, 
                   by.x = "hgnc_symbol", 
                   by.y = "Gene.Symbol")

print(paste0("WARNING: ", length(unique(snp_range$refsnp_id)), " out of ", length(refsnp_id), " have an available gene sequence"), quote = FALSE)

# Get SNP position within the 3'UTR sequence alignment
snp_in_range <- get.snp.position(range_info = snp_range)

# Get allele in mRNA
split_rows_in_range <- split(snp_in_range, rownames(snp_in_range))
snp_data <- lapply(split_rows_in_range, allele.in.RNA.strand)
snp_data <- do.call("rbind", snp_data)
rownames(snp_data) <-c()

# Replace allele in UTR sequence
split_rows_snp_data <- split(snp_data, rownames(snp_data))
snp_expanded <- lapply(split_rows_snp_data, replace.allele)
snp_in_seq <- do.call("rbind", snp_expanded)
rownames(snp_in_seq) <- c()

print(paste0("WARNING: ",length(unique(snp_in_seq$refsnp_id)), " out of ", length(refsnp_id), " will be analysed!"), quote = FALSE)

write.table(snp_in_seq, 
            paste0(working_path, "/TS7_input_", length(unique(snp_data$refsnp_id)), "snp_data_", ensemblV,"_", Sys.Date(), ".txt"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)



################################################################################
#### 07. TargetScan command line ####

snps_info_df <- snp_in_seq[,c("hgnc_symbol", "Ensembl.ID", "refsnp_id", "allele_in_sequence")]

snp_in_seq$Pattern <- do.call("paste", c(snps_info_df, sep = "_"))

criterion<- unique(apply(snps_info_df[, 1:3], 1, paste0, collapse = "_"))
names(criterion) <- criterion

save(criterion, file= paste0(working_path, "/criterion.RData"))

get.command.line.TS(snps_info_df)



################################################################################
#### 08. TargetScan - 84 species ####

# Get MSA for the remaining species
UTR_sequences_species <- UTR_sequences[grep("^9606$",UTR_sequences$Species.ID, invert = TRUE),]  # Human species = 9606

# write input files for targetscan_70.pl
for (i in row(snp_in_seq)) {
  y <- snp_in_seq[i, c("Ensembl.ID", "Gene.ID", "hgnc_symbol", "Species.ID", "UTR.sequence")]
  colnames(y)[3] <- "Gene.Symbol"
  ensembl <- snp_in_seq[i, "Ensembl.ID"]
  g <- rbind(y, UTR_sequences_species)
  z <- g[grep(ensembl, g$Ensembl.ID),]
  pat <- snp_in_seq[i, "Pattern"]  
  z <- cbind(z,pat)
  write.table(z[,c(6,4,5)], paste0(working_path,"/",paste(snp_in_seq$hgnc_symbol[i],snp_in_seq$Ensembl.ID[i],snp_in_seq$refsnp_id[i], "allele",snp_in_seq$allele_in_sequence[i],"species", sep = "_"), ".txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}


for (i in command.line$command_1)    # Execute targetscan_70.pl
  system(command = i, wait = TRUE)

## Get differences without calculating context++ scores ##

# targetscan_output <- lapply(criterion, read.file, suffix = "species_output", per.allele = TRUE)  # Read output file

# Get differences between correspondent output files
# diff_targetscan_output <- lapply(targetscan_output, fun.diff, 
#                                  snp_data = snp_data, 
#                                  col = c(2,3,4,5,6,7,9), 
#                                  prefix = "diff", 
#                                  paste.col = c(8,10:14))
# diff_targetscan_output <- do.call("rbind", diff_targetscan_output)

# write.table(diff_targetscan_output, 
#             paste0(working_path,"/diff_TargetScan70_",
#                    length(unique(diff_targetscan_output$snp_id)),
#                    "snps.txt"), 
#             sep = "\t", row.names = F, quote = F)

# Read diff files according to criterion
# targetscan_output_diff <- lapply(criterion, read.file, prefix = "diff_", per.allele = FALSE) 


## ##


for (i in command.line$command_2)    # Execute targetscan_70_BL_bins.pl
  system(command = i, wait = TRUE)


for (i in command.line$command_3)    # Execute targetscan_70_BL_PCT.pl
  system(command = i, wait = TRUE)


tkey <-unique(snp_in_seq[,c("Ensembl.ID","hgnc_symbol","Pattern")])    
transcripts <- as.character(unique(tkey$Ensembl.ID))


# Get ORFs
ORF_sequences_hg19_list <- lapply(transcripts, grep.file, filename = "ORF_Sequences.txt")
ORF_sequences_hg19 <- do.call(rbind, ORF_sequences_hg19_list)
ORF_sequences_hg19 <- merge(ORF_sequences_hg19, tkey, 
                            by.x = "V1", 
                            by.y = "Ensembl.ID")

ORF_split <- split(ORF_sequences_hg19, ORF_sequences_hg19$V1)


for (i in names(ORF_split)){
  x <- ORF_split[[i]]
  y <-x[,c("Pattern","V2","V3")] 
  write.table(y, 
              paste0(working_path,"/",paste("ORF",x$hgnc_symbol[1],i, sep = "_"), ".txt"), 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE, 
              col.names = FALSE)
}

for (i in unique(command.line$command_4))    # Execute targetscan_count_8mers.pl
  system(command = i, wait = TRUE)

# From miRNA_Family_Info.txt create file containing only human miRNAs
if (dir.exists("./miR_Family_Info_hsa.txt") == FALSE){
  miR_Family_Info = read.table("./miR_Family_Info.txt", sep = "\t", header = TRUE)   # Subset miR file for human miR 
  miR_Family_Info_hsa = miR_Family_Info[grep("^9606$", miR_Family_Info$Species.ID ),]
  write.table(miR_Family_Info_hsa, "./miR_Family_Info_hsa.txt",sep = "\t", quote = FALSE, row.names = FALSE)
}

# Resize dataframe for context score application (command 5)
if (dir.exists("./miR_hsa_for_context_scores.txt") == FALSE){
  system(command = "cut -f1,3,4,5 ./miR_Family_Info_hsa.txt > ./miR_hsa_for_context_scores.txt", intern = TRUE, wait = TRUE)
}

# Rename ID in All_cell_lines.AIRs.txt to match 'Pattern', otherwise the script will not recognize them
# 'All_cell_lines.AIRs.txt' file name is hardcoded in the targetscan_70_context_scores.pl script
get.AIRs(1, tkey)  # range 0:1 || 1: Subset AIR file for used transcripts id and replace by 'Pattern' (from key) in new AIR file  | 0: Use the original AIR dataset

for (i in command.line$command_5)    # Execute targetscan_70_context_scores.pl
  system(command = i, wait = TRUE) 


targetscan_output_cs <- lapply(criterion, read.file, suffix = "output_context_score", per.allele = TRUE, read.cs = TRUE)



################################################################################
#### 09. Process TS output ####

# Get differences between correspondent output files containing context++ scores
diff_targetscan_output_cs <- lapply(targetscan_output_cs, 
                                    fun.diff, 
                                    snp_data = snp_data, 
                                    col = c(2,4:6,36), 
                                    prefix = "diff_cs", 
                                    paste.col = c(1,3,7:35,37))
diff_targetscan_output_cs <- do.call("rbind", diff_targetscan_output_cs)

# Read diff_cs files according to criterion
# targetscan_output_cs_diff <- lapply(criterion, read.file, prefix = "diff_cs_", per.allele = FALSE)

diff_targetscan_output_cs$name <- paste(diff_targetscan_output_cs$hgnc_symbol,
                                        diff_targetscan_output_cs$transcript_id,
                                        diff_targetscan_output_cs$snp_id,
                                        sep="_")

# Read diff_norm_cs files according to criterion
# targetscan_output_cs_diff <- lapply(criterion, read.file, prefix= "diff_cs_", per.allele = FALSE)

compile.tables(diff_targetscan_output_cs, 
               miRnorm = miR_normalization)  # Make a table txt file compiling all diff_cs. Also compiles miR normalization table for all. 



################################################################################
#### 10. Make plots ####

## Set scale for all
# Min
min_cs_by_diff.file <- lapply(targetscan_output_cs, 
                              get.min, 
                              by.c = "context...score")
min_cs <- unlist(min_cs_by_diff.file)
names(min_cs) <- c()
min_cs <- min(min_cs)  # Minimum context++ score for the provided dataset


# Max
max_cs_by_diff.file <- lapply(targetscan_output_cs, 
                              get.max, 
                              by.c = "context...score")
max_cs <- unlist(max_cs_by_diff.file)
names(max_cs) <- c()
max_cs <- max(max_cs)  # Maximum context++ score for the provided dataset

# Min-Max
cs_range <- c(min_cs,max_cs)  # vector of (max,min) context++ scores in dataset


# Generate plots

fun.diff.cs.plot(diff_targetscan_output_cs, miR_expression = miR_normalization, scale_limit = cs_range, c.split = "name")


# Merge plots.pdf into single pdf

Plots <- list.files(paste0(working_path, "/Plots_cs_diff"), pattern = ".pdf")
plotflow:::mergePDF(
  in.file = paste0(working_path, "/Plots_cs_diff/", Plots),
  file = paste0(working_path, "/Plots_cs_diff/Results_TS7_CS_", 
                length(unique(SNAP_data$SNP)), "tag_", 
                length(unique(snp_data$refsnp_id)), "snps.pdf"))

if (miR_normalization !=0){
  Plots_miRExpressed <- list.files(paste0(working_path, "/Plots_cs_diff/miRExpressed"), pattern = ".pdf")
  plotflow:::mergePDF(
    in.file = paste0(working_path, "/Plots_cs_diff/miRExpressed/", Plots_miRExpressed),
    file = paste0(working_path, "/Plots_cs_diff/miRExpressed/Results_TS7_CS_",
                  length(unique(SNAP_data$SNP)), "tag_", 
                  length(unique(snp_data$refsnp_id)),
                  "snps_miRExpressed", miRExp_dataset_name, ".pdf"))
}



################################################################################
#### 11. Save data  ####

file.rename(from = paste0(getwd(), "/RNAplfold_in_out"), 
            to = paste0(working_path, "/RNAplfold_in_out")) # Remove this if you want to use these results recursively
file.rename(from = paste0(getwd(), "/All_cell_lines.AIRs.txt"), 
            to = paste0(working_path, "/All_cell_lines.AIRs.sample.txt"))  # Relocate file with subset of AIR used in this analysis
file.rename(from = paste0(getwd(), "/All_cell_lines.AIRs.ORIGINAL.txt"), 
            to = paste0(getwd(), "/All_cell_lines.AIRs.txt"))  # Rename file with all AIR back to original name 


save.image(file = paste0(working_path, "/.RData"))  # Save environment
savehistory(file = paste0(working_path, "/.Rhistory"))  # Save history

sessionInfo()
