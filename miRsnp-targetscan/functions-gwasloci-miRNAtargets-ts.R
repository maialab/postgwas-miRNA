################################################################################
### f(Analysis of breast cancer risk loci for allele-specific miRNA binding) ### 
##            TargetScan v7.1 (Agarwal et al., 2015) Application              ##
################################################################################

# Author: Ana Jacinta-Fernandes, <a46845@ualg.pt>
# Date: Feb 2018
# The Functional Genomics of Cancer Laboratory, PI: Ana-Teresa Maia
# Centre for Biomedical Research (CBMR), University of Algarve, Portugal
# <www.maialab.org>

#### 00. Preparation  ####

check.packages <- function(){
  
  if(!require(biomaRt))
    stop("\nPlease install the 'biomaRt' R package in order to proceed\n")
    
  if(!require(tools))
    stop("\nPlease install the 'tools' R package in order to proceed\n")
  
  if(!require(ggplot2))
    stop("\nPlease install the 'ggplot2' R package in order to proceed\n")
  
  if(!require(devtools))
    stop("\nPlease install the 'devtools' R package in order to proceed\n")
  
  if(!require(plotflow))
    stop("\nPlease install the 'plotflow' R package in order to proceed\n") 
  
}


check.files <- function(path){
  
  if (missing(path) == TRUE)
    path = getwd()
  
  writeLines("\nTargetScan v7.1 (Agarwal et al., 2015) data and code are available at\n   www.targetscan.org/cgi-bin/targetscan/data_download.vert71.cgi\n")
  
  writeLines(paste0("Provided path:\n   ",path))
  
  # Perl scripts check
  writeLines("\n[1] Checking if Perl scripts exist...")
  
  if(!file.exists(paste0(path,"/targetscan_70.pl")))
    stop("\nPlease upload the 'targetscan_70.pl' Perl script in order to proceed\n")
  if(!file.exists(paste0(path,"/targetscan_70_BL_bins.pl")))
    stop("\nPlease upload the 'targetscan_70_BL_bins.pl' Perl script in order to proceed\n")
  if(!file.exists(paste0(path,"/targetscan_70_BL_PCT.pl")))
    stop("\nPlease upload the 'targetscan_70_BL_PCT.pl' Perl script in order to proceed\n")
  if(!file.exists(paste0(path,"/targetscan_count_8mers.pl")))
    stop("\nPlease upload the 'targetscan_count_8mers.pl' Perl script in order to proceed\n")
  if(!file.exists(paste0(path,"/targetscan_70_context_scores.pl")))
    stop("\nPlease upload the 'targetscan_70_context_scores.pl' Perl script in order to proceed\n")
  
  
  
  # Perl scripts dependencies check
  writeLines("[2] Checking if Perl script dependencies are installed...")
  
  Perl <- system("perl -v")
  if (Perl != 0)
    stop("\nPlease install the Perl programming language in order to proceed\n")
  
  Statistics<- system("perl -MStatistics::Lite -e 1")
  if (Statistics != 0)
    stop("\nPlease install the 'Statistics::Lite' Perl module in order to proceed\n")
  Tree<- system("perl -MBio::TreeIO -e 1")
  if (Tree != 0)
    stop("\nPlease install the 'Bio::TreeIO' (Part of BioPerl) Perl module in order to proceed\n")
  
  RNAplfold<- system("RNAplfold --version")
  if (RNAplfold != 0)
    warning("\nPlease install the 'RNAplfold' Application from the ViennaRNA package 2 available at www.tbi.univie.ac.at/RNA/documentation.html to calculate average pair probabilities for locally stable secondary structures\n")
  
  
  # Perl scripts file dependencies check
  writeLines("[3] Checking if Perl script data dependencies exist...")
  
  if(!file.exists(paste0(path,"/UTR_Sequences.txt")))
    stop("\nPlease upload the 'UTR_Sequences.txt' file in order to proceed\n")
  if(!file.exists(paste0(path,"/ORF_Sequences.txt")))
    stop("\nPlease upload the 'ORF_Sequences.txt' file in order to proceed\n")
  if(!file.exists(paste0(path,"/miR_Family_Info.txt")))
    stop("\nPlease upload the 'miR_Family_Info.txt' file in order to proceed\n")
  if(!file.exists(paste0(path,"/TSHuman_7_hg19_3UTRs.gff")))
    stop("\nPlease upload the 'TSHuman_7_hg19_3UTRs.gff' file in order to proceed\n")
  if(!file.exists(paste0(path,"/PCT_parameters")))
    stop("\nPlease upload the 'PCT_parameters' folder in order to proceed\n")
  if(!file.exists(paste0(path,"/Agarwal_2015_parameters.txt")))
    stop("\nPlease upload the 'Agarwal_2015_parameters.txt' file in order to proceed\n")
  if(!file.exists(paste0(path,"/All_cell_lines.AIRs.txt")))
    stop("\nPlease upload the 'All_cell_lines.AIRs.txt' file in order to proceed\n")
  if(!file.exists(paste0(path,"/TA_SPS_by_seed_region.txt")))
    stop("\nPlease upload the 'TA_SPS_by_seed_region.txt' file in order to proceed\n")
  
  writeLines("\n...Success!\n")
}


file.save <- function(fullpath){
  
  if (missing(fullpath) == TRUE)
    fullpath <- paste0(getwd(),"/Analysis_", Sys.Date())
  
  dir.create(fullpath)
  assign("working_path", fullpath, envir = .GlobalEnv)
  
  writeLines(paste0("\nAll data generated from the analysis will be stored at:\n   ", fullpath, "\n"))
}


#### 01. Settings  ####

get.miR.expression <- function(get_expression, path){
  
  # miRBase (v21) contains about 2588 mature human miRNAs
  
  assign("miR_normalization", get_expression, envir = .GlobalEnv)
  
  if (missing(path))
    path <- getwd()
  
  
  if (get_expression == 0){
    
    writeLines("\nNormalization for miRNA expression DISABLED\n")
    
    return() }
  
  else {
    if (get_expression == 1){  # Use miR expression from miRmine_v1.0 as normalization
      
      writeLines("\nUsing miRmine v1.0 (Panwar et al., 2017) miRNA expression values for breast...\nThe dataset is available at\n   guanlab.ccmb.med.umich.edu/mirmine/index.html")
      
      fullpath = paste0(path,"/miRmine - Human miRNA Expression Database.csv")
      
      if (file.exists(fullpath) == FALSE)
        stop("\nPlease upload the 'miRmine - Human miRNA Expression Database.csv' file in order to proceed\n")
      
      # sT = serum tumor || tT = tissue tumor || sN = serum Normal || tN = tissue Normal  
      
      miRmine_Breast_sT_tT_sN_tN <- read.csv(fullpath, 
                                             header = TRUE,
                                             colClasses = "character")
      miRmine_Breast_sT_tT_sN_tN[,3:6] <- lapply(miRmine_Breast_sT_tT_sN_tN[,3:6],as.numeric)
      
      # 2588 unique miRNAs (miRBase V20 || The latest version of miRBase (v21) contains about 2588 mature human miRNAs) 
      
      colnames(miRmine_Breast_sT_tT_sN_tN)[3] <- "SRX513283_serum_Tumour_Breast"
      colnames(miRmine_Breast_sT_tT_sN_tN)[4] <- "SRX513284_tissue_Tumour_Breast"
      colnames(miRmine_Breast_sT_tT_sN_tN)[5] <- "SRX513285_serum_Normal_Breast"
      colnames(miRmine_Breast_sT_tT_sN_tN)[6] <- "SRX513286_tissue_Normal_Breast"
      
      writeLines("\nAssigning 'miRmine_expression' to your .GlobalEnv containing miRNA expression levels...")
      assign("miRmine_expression",miRmine_Breast_sT_tT_sN_tN, envir = .GlobalEnv)
      
      # Select miRNAs that are expressed in normal breast tissue
      miRmine_Breast_tN_expressed <- subset(miRmine_Breast_sT_tT_sN_tN,
                                            miRmine_Breast_sT_tT_sN_tN$SRX513286_tissue_Normal_Breast > 0)
      
      normal_breast <- unique(miRmine_Breast_tN_expressed$Mature.miRNA.ID)
      
      # Get stats
      RPM_nB <- miRmine_Breast_tN_expressed$SRX513286_tissue_Normal_Breast
      
      min.RPM <- min(RPM_nB)
      max.RPM <- max(RPM_nB)
      median.RPM <- median(RPM_nB)
      
      # hist(RPM_nB,
      #      main = "Normal breast tissue miRNA expression (SRX513286)", 
      #      xlab = "RPM")
      
      writeLines(paste0("\nNormal breast tissue expression stats:\n        min = ",min.RPM ," RPM    max = ",max.RPM, " RPM    median = ",median.RPM," RPM"))
      
      writeLines("\nAssigning 'normal_breast' to your .GlobalEnv containing a character vector of expressed miRNAs in normal breast tissue...")
      assign("normal_breast",normal_breast, envir = .GlobalEnv)
      
      miR_norm_name <- "miRmine_SRX513286_normal_breast"
      
      # Extract .csv from: guanlab.ccmb.med.umich.edu/mirmine/celltissue.php?tissue=breast&cline=
      # single-end sequencing (36 bp) RNA-seq by Illumina HiSeq 2500
      # All patients had not been previously treated by chemotherapy and radiotherapy when undergoing surgery and provided informed consent to participate in the study. Fresh frozen breast cancer tumors, adjacent normal tissues, and preoperative serum from 8 patients with breast cancer and control serum sample from 8 healthy female volunteers were obtained from the Taizhou Central Hospital.
      # small RNA sequencing reads were aligned against 2578 mature miRNA sequences from miRBase build 20 using Bowtie 1.0.0 (Langmead et al., 2009) allowing at most two mismatches. The other parameters are default.
      # Values are in RPM (Reads Per Million mapped reads)              
    }
    
    assign("miRExp_dataset_name", miR_norm_name, envir = .GlobalEnv)
  }
}



#### 02. Proxy Search  ####

get.proxy.SNAP.data <- function(SNAP_output, path){
  
  if (missing(path))  # Get path 
    path = getwd()
  
  snap.txt <- paste0(path,"/",SNAP_output)  # Get file' complete path
  
  if (file.exists(snap.txt) == FALSE){
    return(
      writeLines(paste0("\nERROR: SNAP output file not found. Check directory/path\n","\n  --> Provided SNAP output file: ",snap.txt,"\n"))
    )
  }
  
  SNAP_data <- read.table(snap.txt, 
                         header = TRUE, 
                         sep = "\t", 
                         fill = TRUE, 
                         colClasses = c(rep("character",3),
                                        rep("numeric",2),
                                        rep("character",2),
                                        "character"),
                         na.strings = "")
  
  assign("SNAP_data", SNAP_data, envir = .GlobalEnv)
  writeLines("\nSNAP output can be found in your .GlobalEnv as 'SNAP_data'\n")
  
  SNAP <- SNAP_data[grep("No matching proxy snps found", SNAP_data$Distance, invert = TRUE),]
  
  warnings <- SNAP[SNAP$Proxy == "WARNING", "SNP"]
  
  writeLines(paste0("No matching proxy snps found for: ",paste(warnings,collapse=" ")))
  
  SNAP[SNAP$Proxy == "WARNING", "Proxy"] <- SNAP[SNAP$Proxy ==  "WARNING", "SNP"]
  
  tag.proxies <- unique(SNAP$Proxy)
  tag <- unique(SNAP$SNP)
  
  writeLines(paste0("\nThere are ", 
                    length(tag.proxies), 
                    " unique tag&proxy SNPs (from ",
                    length(tag),
                    " tags)"))
  
  tag.proxies
}


#### 03. Ensembl query  ####

query.snp.ensembl <- function(snp, path){
  
  if (missing(path) == TRUE){
    if (exists("working_path") == TRUE)
      path = working_path
    else path = getwd()
  }
  
  
  snp_found_in_ensembl <- getBM(attributes = "refsnp_id", 
                                filters = "snp_filter", 
                                values = snp, 
                                mart = snpmart)
  
  writeLines(paste0("\nWARNING: The Ensembl SNP mart only finds ",
                    length(unique(snp_found_in_ensembl$refsnp_id)),
                    " SNPs out of ",length(snp),
                    ", as the database excludes flagged variants"))
  
  snp_information <- getBM(attributes = c("refsnp_id","chr_name","allele","allele_1","minor_allele","chrom_strand","ensembl_gene_stable_id","ensembl_transcript_stable_id","ensembl_transcript_chrom_strand","consequence_type_tv"), 
                          filters = "snp_filter", 
                          values = snp, 
                          mart = snpmart)
  
  snp <- data.frame(snp=snp)
  
  snp_information <- merge(snp_information, 
                          snp, 
                          by.x = "refsnp_id", 
                          by.y = "snp", 
                          all = TRUE)
  
  writeLines("\nProcessing ensembl annotations...\n")
  
  #Remove "LRG" gene name.
  snp_information_processed <- snp_information[grep("LRG", snp_information$ensembl_gene_stable_id, invert = TRUE),]
  writeLines(paste0("[1] Removed ", length(grep("LRG", snp_information$ensembl_gene_stable_id)), " LRG (Locus Reference Genomic) gene names"))
  #Remove "CHR" chromossome annotations
  snp_information_processed <- snp_information_processed[grep("CHR", snp_information_processed$chr_name, invert = TRUE),]
  writeLines(paste0("[2] Removed ", length(grep("CHR", snp_information_processed$chr_name)), " CHR chromossome name annotations"))
  
  #Get hngc gene name symb
  geneID <- unique(snp_information_processed$ensembl_gene_stable_id)
  
  geneID_symbol <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), 
                         filters = "ensembl_gene_id", 
                         values = geneID, 
                         mart = hsapiens)
  
  writeLines(paste0("\nThe data includes the following ",length(unique(geneID_symbol$hgnc_symbol)) ," HGNC gene symbols for ",length(geneID)," Ensembl ID's:\n"))
  print(unique(geneID_symbol$hgnc_symbol))
  
  snp_information_processed<-  merge(snp_information_processed, 
                                     geneID_symbol, 
                                     by.x = "ensembl_gene_stable_id", 
                                     by.y = "ensembl_gene_id", all = TRUE)
  
  #reorder collums
  snp_information_processed<- snp_information_processed[,c("refsnp_id","chr_name","allele","allele_1","minor_allele","chrom_strand","ensembl_gene_stable_id","hgnc_symbol","ensembl_transcript_stable_id","ensembl_transcript_chrom_strand","consequence_type_tv")]
  
  snp_information_processed[snp_information_processed == ""] = NA    # Replace all "" with NA
  
  
  writeLines("\nThe data includes the following SNP genomic localizations:\n")
  print(unique(snp_information_processed$consequence_type_tv))
  
  #save our dataset
  
  ensemblV <- sub(" Genes ", "", listEnsembl()[1,2])  # Get Ensembl release 
  assign("ensemblV", ensemblV, envir = .GlobalEnv)
  dataset_name <- paste0(path,"/dataset_",ensemblV,"_",length(unique(snp$snp)),"snps_",Sys.Date(),".txt")
  writeLines(paste0("\nCreating .txt file with data:\n   ",dataset_name))
  
  write.table(snp_information_processed, 
              file = dataset_name, 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
  
  
  snp_information_processed
  
}



#### 04. Ensembl query - hg19 coordinates  ####

snp.coordinates.hg19 <- function(snp){
  
  snp_coordinates <- getBM(attributes = c("refsnp_id","chr_name","chrom_start","allele","allele_1","minor_allele"), 
                           filters = "snp_filter", 
                           values = snp, 
                           hg19_snp) 
  
  snp_coordinates = snp_coordinates[grep("PATCH", snp_coordinates$chr_name, invert = TRUE),]  # Remove other chromossome annotations
  snp_coordinates
}



#### 05. UTR query  ####

grep.file <- function(pattern, filename, path){  # Courtesy of Ramiro Magno
  if (missing(path) == TRUE)
    path = getwd()
  
  filename = paste0(path,"/",filename)  
  
  pat<- paste0("\"\\b",pattern,"\\b\"")
  cmd<- paste("grep -P",pat, filename, sep = " ") 
  print(cmd)
  grep.txt <- system(command = cmd, intern = TRUE)  # grep from bash --> grep directly from .txt  Returns string by row
  if (length(grep.txt) != 0){
    grep.list <- lapply(grep.txt, strsplit, split = "\t") # split string by \t and output to dataframe
    ncolumns <- length(grep.list[[1]][[1]])
    as.data.frame(matrix(unlist(grep.list),ncol = ncolumns, byrow = TRUE))
  }
  else print(paste0("Pattern ",pattern," not found (status 1)"), quote = FALSE)
}



#### 06. Prepare dataset for TargetScan  ####

check.interval <- function(intervals, c.start, c.end){
  
  c_ordered_interval <- intervals[order(intervals[,c.start]),]
  c_ordered_interval = unique(c_ordered_interval)
  
  # guarantee numeric
  c_ordered_interval[,c.start] <- as.numeric(as.character(c_ordered_interval[,c.start]))  # if is factor, then coerse to character and only then numeric
  c_ordered_interval[,c.end]<- as.numeric(as.character(c_ordered_interval[,c.end]))
  
  
  if (nrow(c_ordered_interval) >= 2){  # if 2 or more sequence coords available, see if there is a gap
    
    for (i in 1:nrow(c_ordered_interval))
      c_ordered_interval$gap[i] <- c_ordered_interval[i,c.end] - (c_ordered_interval[i+1,c.start] -1)
    
    
    for (i in 1:nrow(c_ordered_interval)) {  # the last value will be NA since there is no rows after that (see above)
      if (is.na(c_ordered_interval[i,"gap"]) == TRUE)
        c_ordered_interval[i,"gap"] = 0  # set as 0 length gap
    }
    
    intervals$gap <- sum(c_ordered_interval$gap)  
    
    if (intervals$gap[1] != 0){
      writeLines(paste0("There is a ",intervals$gap[1], " gap in the provided UTR coordinates of ",intervals$Ensembl.ID[1],"\n--> SNP position in UTR sequence will be ajusted accordingly."))
    }
    
  }
  
  else intervals$gap = 0
  
  intervals
}



get.snp.position <- function(range_info, path){
  
  # guarantee coords are numeric
  range_info$chrom_start<- as.numeric(as.character(range_info$chrom_start))
  range_info$start<- as.numeric(as.character(range_info$start))
  range_info$end <- as.numeric(as.character(range_info$end))
  
  
  c <- range_info$chrom_start  
  range_info$ranges = (c > range_info$start & c < range_info$end) # Get TRUE/FALSE if SNP is within the range of the sequence
  
  # Get min and max coords for each transcript
  t_coordinates<- range_info[,c("Ensembl.ID","start","end")]
  t_coordinates$Ensembl.ID <- as.character(t_coordinates$Ensembl.ID)
  
  # transcript_id = lapply(transcript_ranges, nrow)
  # unlist(transcript_id)
  
  transcript_ranges = split(t_coordinates, t_coordinates$Ensembl.ID)
  
  for (i in names(transcript_ranges)){
    mint <- min(transcript_ranges[i][[1]][[2]])
    maxt <- max(transcript_ranges[i][[1]][[3]])
    
    range_info[range_info$Ensembl.ID == i,"Seq_start"] <- mint
    range_info[range_info$Ensembl.ID == i,"Seq_end"] <- maxt
  }
  
  
  # Remove if it does not fall in range
  snps_in_range<- range_info[grep("TRUE", range_info$ranges),]
  
  for (i in 1:nrow(snps_in_range)){
    if (snps_in_range[i,"strand"] == "-")
      snps_in_range$position_in_sequence[i] = (snps_in_range$Seq_end[i] - snps_in_range$chrom_start[i] +1)
    else if (snps_in_range[i,"strand"] == "+")
      snps_in_range$position_in_sequence[i] = (snps_in_range$chrom_start[i] - snps_in_range$Seq_start[i] +1)
  }
  
  rownames(snps_in_range)<-c()
  
  
  # Check if provided gene ranges have any gaps, i.e., are not continuous
  
  intervals_gaps <- lapply(transcript_ranges, check.interval, c.start = "start", c.end = "end" )
  intervals_gaps <- do.call("rbind",intervals_gaps)
  intervals_gaps <- unique(intervals_gaps[,c("Ensembl.ID","gap")])
  rownames(intervals_gaps)<-c()
  
  snps_in_range <- merge(snps_in_range, intervals_gaps, by.x = "Ensembl.ID", by.y = "Ensembl.ID")
  
  # Adjust SNP position in UTR if there is a gap
  snps_in_range$position_in_sequence = snps_in_range$position_in_sequence + snps_in_range$gap
  
  snps_in_range
}


allele.in.RNA.strand <- function(SNP_info){
  
  SNP_info$allele<- as.character(SNP_info$allele)
  SNP_info$nr_alleles <- (nchar(SNP_info$allele)+1)/2  # Determine number of alleles
  
  split_alleles_gene<- strsplit(SNP_info$allele[1], "/")
  info = SNP_info[1,c("allele","strand")]
  
  info$Allele_A<- split_alleles_gene[[1]][1]
  info$Allele_B<- split_alleles_gene[[1]][2] 
  info$Allele_C<- split_alleles_gene[[1]][3] 
  info$Allele_D<- split_alleles_gene[[1]][4]
  
  # Conversion table DNA to RNA according to strand
  DNAtoRNA <- data.frame(DNA=c("A","T","G","C"), 
                         RNAplus=c("A","U","G","C"),
                         RNAminus=c("U","A","C","G"),
                         stringsAsFactors = FALSE)
  
  Allele_X <- which(!is.na(info))  # Get cols without NA
  Allele_X <- Allele_X[-c(1:2)] # Remove first 2 col #'s (allele, strand)
  
  for (i in colnames(info)[Allele_X]){
    
    alleleDNA <- as.character(info[,i])
    
    if (info$strand == "-")
      info[,i] <- DNAtoRNA[DNAtoRNA$DNA == alleleDNA, "RNAminus"]
    
    if (info$strand == "+")
      info[,i] <- DNAtoRNA[DNAtoRNA$DNA == alleleDNA, "RNAplus"]
    
  } 
  
  info<- info[,-(1:2)]
  alleles_snp_info<- cbind(SNP_info, info)
  alleles_snp_info
}



allele.in.DNA.strand <- function(SNP_info){     # used in fun.diff
  SNP_info$allele<- as.character(SNP_info$allele)
  SNP_info$strand <- as.character(SNP_info$strand)
  
  allele<- unique(SNP_info$allele)
  
  if (allele[1] %in% letters == TRUE) # in case allele is in lowercase
    SNP_info$allele <- toupper(SNP_info$allele)
  
  RNAtoDNA <- data.frame(DNA=c("A","T","G","C"), 
                         RNAplus=c("A","U","G","C"),
                         RNAminus=c("U","A","C","G"),
                         stringsAsFactors = FALSE)
  
  if (SNP_info$strand[1] == "-"){
    for (i in 1:nrow(SNP_info)){
      alleleDNA <- SNP_info[i,"allele"]
      SNP_info[i,"allele"] <- RNAtoDNA[RNAtoDNA$RNAminus == alleleDNA, "DNA"]
    }}
  
  if (SNP_info$strand[1] == "+"){
    for (i in 1:nrow(SNP_info)){
      alleleDNA <- SNP_info[i,"allele"]
      SNP_info[i,"allele"] <- RNAtoDNA[RNAtoDNA$RNAplus == alleleDNA, "DNA"]
    }}
  
  SNP_info
  
}



replace.allele<- function (info_snp){
  info_snp<- droplevels(info_snp)
  
  alleleA <- info_snp$Allele_A[1]
  alleleB <- info_snp$Allele_B[1]
  alleleC <- info_snp$Allele_C[1]
  alleleD <- info_snp$Allele_D[1]
  
  pos <- info_snp$position_in_sequence[1]
  seq0 <- as.character(info_snp$UTR.sequence)[1]
  
  seq <- strsplit(seq0, '')[[1]]    # split into letters
  
  to_replace <- seq[seq != '-'][pos]    # identify allele to replace
  to_replace<- as.character(to_replace)
  
  if (to_replace %in% letters == TRUE){   # in case allele in in lowercase (soft masking of MSA)
    alleleA <- tolower(alleleA)
    alleleB <- tolower(alleleB)
    alleleC <- tolower(alleleC)
    alleleD <- tolower(alleleD)
  }   
  
  snp_info_expanded <- info_snp[rep(row.names(info_snp), info_snp$nr_alleles),] # duplicate row by # of alleles
  rownames(snp_info_expanded) <- c()  
  snp_info_expanded$UTR.sequence<- as.character(snp_info_expanded$UTR.sequence)
  
  # assign appropriate replacement to subset
  
  if (nrow(snp_info_expanded) >= 2){
    seq1<- seq
    if (alleleA == to_replace){
      seq1[seq1 != '-'][pos] <- alleleB
      snp_info_expanded$UTR.sequence[2]<- paste(seq1, collapse = '')    # reassemble vector to string
      snp_info_expanded$allele_in_sequence[1]<- to_replace
      snp_info_expanded$allele_in_sequence[2]<- alleleB
    }
    if (alleleB == to_replace){
      seq1[seq1 != '-'][pos] <- alleleA
      snp_info_expanded$UTR.sequence[2]<- paste(seq1, collapse = '')    # reassemble vector to string
      snp_info_expanded$allele_in_sequence[1]<- to_replace
      snp_info_expanded$allele_in_sequence[2]<- alleleA
    }
  }
  
  if (nrow(snp_info_expanded) >= 3){
    seq2<- seq
    if (alleleA == to_replace){
      seq2[seq2 != '-'][pos] <-  alleleC
      snp_info_expanded$UTR.sequence[3]<- paste(seq2, collapse = '')    # reassemble vector to string
      snp_info_expanded$allele_in_sequence[3]<- alleleC
    }
    if (alleleB == to_replace){
      seq2[seq2 != '-'][pos] <-  alleleC
      snp_info_expanded$UTR.sequence[3]<- paste(seq2, collapse = '')    # reassemble vector to string
      snp_info_expanded$allele_in_sequence[3]<- alleleC
    }
    if (alleleC == to_replace){
      seq1[seq1 != '-'][pos] <-  alleleA
      seq2[seq2 != '-'][pos] <-  alleleB
      snp_info_expanded$UTR.sequence[2]<- paste(seq1, collapse = '')    # reassemble vector to string
      snp_info_expanded$UTR.sequence[3]<- paste(seq2, collapse = '')    # reassemble vector to string
      snp_info_expanded$allele_in_sequence[1]<- to_replace
      snp_info_expanded$allele_in_sequence[2]<- alleleA
      snp_info_expanded$allele_in_sequence[3]<- alleleB
    }
  }
  
  if (nrow(snp_info_expanded) == 4){
    seq3<- seq
    if (alleleA == to_replace){
      seq3[seq3 != '-'][pos] <-  alleleD
      snp_info_expanded$UTR.sequence[4]<- paste(seq3, collapse = '')    # reassemble vector to string
      snp_info_expanded$allele_in_sequence[4]<- alleleD
    }
    
    if (alleleB == to_replace){
      seq3[seq3 != '-'][pos] <-  alleleD
      snp_info_expanded$UTR.sequence[4]<- paste(seq3, collapse = '')    # reassemble vector to string
      snp_info_expanded$allele_in_sequence[4]<- alleleD
    }
    
    if (alleleC == to_replace){
      seq3[seq3 != '-'][pos] <-  alleleD
      snp_info_expanded$UTR.sequence[4]<- paste(seq3, collapse = '')    # reassemble vector to string
      snp_info_expanded$allele_in_sequence[4]<- alleleD
    }
    
    if (alleleD == to_replace){
      seq1[seq1 != '-'][pos] <-  alleleA
      seq2[seq2 != '-'][pos] <-  alleleB
      seq3[seq3 != '-'][pos] <-  alleleC
      snp_info_expanded$UTR.sequence[2]<- paste(seq1, collapse = '')    # reassemble vector to string
      snp_info_expanded$UTR.sequence[3]<- paste(seq2, collapse = '')    # reassemble vector to string
      snp_info_expanded$UTR.sequence[4]<- paste(seq3, collapse = '')    # reassemble vector to string
      snp_info_expanded$allele_in_sequence[1]<- to_replace
      snp_info_expanded$allele_in_sequence[2]<- alleleA
      snp_info_expanded$allele_in_sequence[3]<- alleleB
      snp_info_expanded$allele_in_sequence[4]<- alleleC
    }
  }
  
  print(paste0("Done: ",snp_info_expanded$refsnp_id[1]," for ref allele ", to_replace, " in ",snp_info_expanded$hgnc_symbol[1], "/",snp_info_expanded$Ensembl.ID[1]), quote = FALSE)
  
  rownames(snp_info_expanded)<- c()
  snp_info_expanded
}



#### 07. TargetScan command line  #### 

get.command.line.TS <- function(info_df, path){
  
  if (missing(path) == TRUE)
    if (exists("working_path") == TRUE)
      path = working_path
    else path = getwd()
    
    
    # Commands for targetscan_70.pl 84-species
    
    input_species<- paste0(path,"/",paste(info_df$hgnc_symbol,info_df$Ensembl.ID, info_df$refsnp_id,"allele",info_df$allele_in_sequence,"species.txt", sep = "_")) 
    command.line = as.data.frame(input_species)
    command.line$Ensembl.ID <- info_df$Ensembl.ID
    
    command.line$output_species<- paste0(file_path_sans_ext(command.line$input_species),"_output.txt")
    
    command.line$command_1<- paste("./targetscan_70.pl", "./miR_Family_Info.txt",  command.line$input_species,  command.line$output_species, sep = " ")
    
    
    # Commands for targetscan_70_BL_bins.pl 84-species
    
    command.line$output.BL_bins<- paste0(file_path_sans_ext(command.line$input_species),"_output.BL_bins.txt")
    
    command.line$command_2<- paste("./targetscan_70_BL_bins.pl", command.line$input_species,">", command.line$output.BL_bins, sep = " ")
    
    
    # Commands for targetscan_70_BL_PCT.pl 84-species
    
    command.line$output.BL_PCT<- paste0(file_path_sans_ext(command.line$input_species),"_output.BL_PCT.txt")
    
    command.line$command_3<- paste("./targetscan_70_BL_PCT.pl", "./miR_Family_Info.txt",command.line$output_species, command.line$output.BL_bins,">", command.line$output.BL_PCT, sep = " ")
    
    
    # Commands for targetscan_count_8mers.pl 84-species
    
    ORF<- unique(paste0(path,"/",paste("ORF",info_df$hgnc_symbol,info_df$Ensembl.ID, sep = "_"),".txt"))
    ORF.command.line <- as.data.frame(ORF)
    ORF.command.line$Ensembl.ID <- unique(info_df$Ensembl.ID)
    
    ORF.command.line$output.8mers<- paste0(file_path_sans_ext(ORF.command.line$ORF),"_8mer_counts.txt")
    
    ORF.command.line$command_4<- paste("./targetscan_count_8mers.pl", "./miR_Family_Info.txt",ORF.command.line$ORF,">|", ORF.command.line$output.8mers, sep = " ")
    
    ORF.command.line$output.lengths<- paste0(file_path_sans_ext(ORF.command.line$ORF),".lengths.txt")
    
    
    command.line <- merge(command.line, ORF.command.line, by.x = "Ensembl.ID", by.y = "Ensembl.ID")
    
    
    # Commands for targetscan_70_context_scores.pl 84-species
    
    command.line$output_cs<- paste0(file_path_sans_ext(command.line$input_species),"_output_context_score.txt")
    
    command.line$command_5<- paste("./targetscan_70_context_scores.pl", "./miR_hsa_for_context_scores.txt", command.line$input_species, command.line$output.BL_PCT, command.line$output.lengths, command.line$output.8mers, command.line$output_cs, sep = " ")
    
    
    assign("command.line", command.line, envir = .GlobalEnv)
    
    writeLines("Sucess! ... 'command.line' was assigned to your .GlobalEnv")
}



#### 08. TS: targetscan_70.pl ####

read.file <- function(name, prefix, suffix, per.allele, path, read.cs=FALSE) {
  if (missing(path) == TRUE)
    if (exists("working_path") == TRUE)
      path = working_path
    else path = getwd()
    
    if (missing(prefix) == TRUE) prefix <- ""
    if (missing(suffix) == TRUE) suffix <- ""
    if (per.allele == TRUE) per.allele <- ".+"
    if (per.allele == FALSE) per.allele <- ""
    
    pat <- paste0("^",prefix, name, per.allele, suffix,".txt$", collapse = "")
    #pat <- name
    filenames <- list.files(path, pattern= pat)
    fullpath_filenames <- paste0(path,"/", filenames)
    names(fullpath_filenames) <- filenames
    if (read.cs == TRUE) # in case we are reading output of targetscan_70_context_scores.pl
      lapply(X = fullpath_filenames, FUN = read.table, sep = "\t", skip = 1, col.names = c("Gene ID","Species ID","Mirbase ID","Site Type","UTR start",	"UTR end",	"Site type contribution","3' pairing contribution",	"local AU contribution",	"Min_dist contribution",	"sRNA1A contribution",	"sRNA1C contribution",	"sRNA1G contribution",	"sRNA8A contribution",	"sRNA8C contribution",	"sRNA8G contribution",	"site8A contribution",	"site8C contribution",	"site8G contribution",	"3'UTR length contribution",	"SA contribution","ORF length contribution",	"ORF 8mer contribution",	"Offset 6mer contribution",	"TA contribution",	"SPS contribution",	"PCT contribution",	"context++ score","context++ score percentile",	"AIR",	"weighted context++ score",	"weighted context++ score percentile",	"UTR region",	"UTR-miRNA pairing",	"mature miRNA sequence","miRNA family","Group #"), colClasses = c("character","numeric",rep("character",2),rep("numeric",28),rep("character",4),"numeric"))
    else lapply(X = fullpath_filenames, FUN = read.table, sep = "\t", header = TRUE)
}

get.AIRs <- function(get_AIR, key, path){
  
  if (missing(path)) path=getwd()
  if (missing(key) & get_AIR == 1) 
    stop("Please provide key with the pattern and correspondent Ensembl ID")
  if (missing(key) & get_AIR == 0) key = NULL
  
  # md5sum("All_cell_lines.AIRs.txt")  # "6d49b1951811f9cca677fe0486bf8241"
  
  AIR.txt <- paste0(path, "/All_cell_lines.AIRs.txt")
  oriAIR.txt <- paste0(path, "/All_cell_lines.AIRs.ORIGINAL.txt")
  
  if (get_AIR == 0){  # Do not want to take into account AIRs
    
    if (file.size(AIR.txt) == 8715858) 
      return("\nAIR will be set to 1 when targetscan_70_context_scores.pl is executed (Unless your sample ID/Gene ID is the ENST ID number)\n")
    
    else if (file.size(AIR.txt) != 8715858)
      warning("\nAttention: All_cell_lines.AIRs.txt is not the original file. May contain sample ID/Gene ID equal to your analysis\n")
  }
  
  
  if (get_AIR == 1){  # Take into account AIRs
    
    transcriptID <- as.character(unique(key$Ensembl.ID))
    
    # First time users:
    if (file.exists(oriAIR.txt) == FALSE  &  file.size(AIR.txt) == 8715858) {   
      
      AIR_list <-lapply(transcriptID, grep.file, filename = "All_cell_lines.AIRs.txt")
      AIR <- do.call(rbind, AIR_list)
      AIR <- merge(AIR, key, by.x = "V1", by.y = "Ensembl.ID")
      
      file.rename(AIR.txt, oriAIR.txt)  # rename original file to not overwrite
      write.table(AIR[,c(11,2:9)],  # column 11 = 'Pattern'
                  AIR.txt, 
                  sep = "\t", 
                  row.names = FALSE, 
                  col.names = FALSE, 
                  quote = FALSE)     
      return(writeLines("\nDone! ... Renamed the original file All_cell_lines.AIRs.txt to All_cell_lines.AIRs.ORIGINAL.txt\nWrote new All_cell_lines.AIRs.txt to match dataset"))
    }
    
    # Second time users:
    if (file.exists(oriAIR.txt) == FALSE  &  file.size(AIR.txt) != 8715858) {
      warning("\nAttention: All_cell_lines.AIRs.txt is not the original file and there is no backup original file. Please check this.\n") 
    }
    
    if (file.exists(oriAIR.txt) == TRUE  &  file.size(oriAIR.txt) == 8715858) {  # If this file already exists, subset to new file 
      
      AIR_list <-lapply(transcriptID,grep.file, filename = "All_cell_lines.AIRs.ORIGINAL.txt")
      AIR <- do.call(rbind, AIR_list)
      AIR <- merge(AIR, key, by.x = "V1", by.y = "Ensembl.ID")
      write.table(AIR[,c(11,2:9)], 
                  AIR.txt, 
                  sep = "\t", 
                  row.names = FALSE, 
                  col.names = FALSE, 
                  quote = FALSE)
      
      return(writeLines("\nDone! ... Queried the original file All_cell_lines.AIRs.ORIGINAL.txt and wrote new All_cell_lines.AIRs.txt to match dataset"))
    }
  }
  
}



#### 09. Process TS output ####

fun.diff <- function(list, col, snp_data, prefix, paste.col, path) {
  
  if (missing(path) == TRUE)
    if (exists("working_path") == TRUE)
      path = working_path
    else path = getwd()
    
    if (missing(paste.col))
      paste.col = NULL
    
    name <- names(list)
    x <- strsplit(name, "_")
    split_names <- as.data.frame(matrix(unlist(x), byrow = TRUE, ncol = length(x[1][[1]])))
    
    # get differences
    if (length(list) >= 2) {
      a111 <- list[[1]]  # Allele A
      a11<- a111[,c(col)]  # select columns to compare
      a1<- do.call("paste", a11)  # Make pattern
      
      a222 <- list[[2]]  # Allele B
      a22<- a222[,c(col)]  # select columns to compare
      a2<- do.call("paste", a22)  # Make pattern
      
      diff1.2 <- a11[! a1 %in% a2, ]  # what is in a1 that is not in a2?
      diff2.1 <- a22[! a2 %in% a1, ]  # what is in a2 that is not in a1?
      
      diffA <- diff1.2  # all diffs for allele A [if length(list) == 2]
      diffB <- diff2.1  # all diffs for allele B [if length(list) == 2]
    }
    
    if (length(list) >= 3){
      a333 <- list[[3]]  # Allele C
      a33<- a333[,c(col)]  # select columns to compare
      a3<- do.call("paste", a33)  # Make pattern
      
      diff1.3<- a11[! a1 %in% a3, ]  # what is in a1 that is not in a3?
      diff2.3<- a22[! a2 %in% a3, ]  # what is in a2 that is not in a3?
      diff3.1<- a33[! a3 %in% a1, ]  # what is in a3 that is not in a1?
      diff3.2<- a33[! a3 %in% a2, ]  # what is in a3 that is not in a2?
      
      diff1.23<- unique(rbind(diff1.2,diff1.3))
      diff2.13<- unique(rbind(diff2.1,diff2.3))
      diff3.12<- unique(rbind(diff3.1,diff3.2))
      
      diffA <- diff1.23  # all diffs for allele A [if length(list) == 3]
      diffB <- diff2.13  # all diffs for allele B [if length(list) == 3]
      diffC <- diff3.12  # all diffs for allele A [if length(list) == 3]
    }
    
    if (length(list) == 4){
      a444 <- list[[4]]  # Allele D
      a44<- a444[,c(col)]  # select columns to compare
      a4<- do.call("paste", a44)  # Make pattern
      
      diff1.4<- a11[! a1 %in% a4, ]  # what is in a1 that is not in a4?
      diff2.4<- a22[! a2 %in% a4, ]  # what is in a2 that is not in a4?
      diff3.4<- a33[! a3 %in% a4, ]  # what is in a3 that is not in a4?
      diff4.1<- a44[! a4 %in% a1, ]  # what is in a4 that is not in a1?
      diff4.2<- a44[! a4 %in% a2, ]  # what is in a4 that is not in a2?
      diff4.3<- a44[! a4 %in% a3, ]  # what is in a4 that is not in a2?
      
      diff1.234<- unique(rbind(diff1.23,diff1.4))
      diff2.134<- unique(rbind(diff2.13,diff2.4))
      diff3.124<- unique(rbind(diff3.12,diff3.4))
      diff4.123<- unique(rbind(diff4.1,diff4.2,diff4.3))
      
      diffA <- diff1.234  # all diffs for allele A [if length(list) == 4]
      diffB <- diff2.134  # all diffs for allele B [if length(list) == 4]
      diffC <- diff3.124  # all diffs for allele C [if length(list) == 4]
      diffD <- diff4.123  # all diffs for allele D [if length(list) == 4]
    }
    
    # Mark each df with correspondent allele and other cols of interest
    if (nrow(diffA) != 0){  
      diffA$allele<- split_names[1,5]
      diffA <- cbind(diffA,a111[rownames(diffA), paste.col])
    }
    if (nrow(diffB) != 0){
      diffB$allele<- split_names[2,5]
      diffB <- cbind(diffB,a222[rownames(diffB), paste.col])
    }  
    if (length(list) >= 3){
      if ( nrow(diffC) != 0){
        diffC$allele<- split_names[3,5]
        diffC <- cbind(diffC,a333[rownames(diffC), paste.col])
      }}  
    if (length(list) == 4){
      if (nrow(diffD) != 0){
        diffD$allele<- split_names[4,5]
        diffD <- cbind(diffD,a444[rownames(diffD), paste.col])
      }}
    
    # Get all differences in one df
    y <- rbind(diffA, diffB)
    if (length(list) == 3)  y <- rbind(y,diffC)
    if (length(list) == 4)  y <- rbind(y,diffD)
    
    
    # Go to SNP info df and extract coorespondent infos
    name_all<- unique(apply(split_names[1, 1:3], 1, paste0, collapse = "_"))  # gene_ENST00000_rs00000
    snp_data$Pattern <- paste(snp_data$hgnc_symbol,snp_data$Ensembl.ID,snp_data$refsnp_id, sep = "_")
    info_snp <- snp_data[grep(name_all, snp_data$Pattern),c("allele","allele_1","minor_allele","position_in_sequence","strand", "chr_name")]
    colnames(info_snp)[1] <- "alleles"
    
    # Add info to diff 
    y$snp_id<- split_names[1,3]  # refsnp ID
    y$hgnc_symbol<- split_names[1,1]  # gene symbol
    y$transcript_id<- split_names[1,2]  # transcript ID
    y <- cbind(y,info_snp)
    
    y <- allele.in.DNA.strand(SNP_info = y)  # convert RNA allele back to DNA
    
    table <- paste0(path,"/",paste(prefix, name_all, sep = "_"), ".txt")
    writeLines(paste0("Writing table: ",table))
    write.table(y, table, sep = "\t", row.names = FALSE, quote = FALSE)
    y
}


compile.tables <- function(df, miRnorm) {   
  
  if (miRnorm == 1){
    miR_expression = normal_breast
    miR_norm <- 'miRmine_SRX513286_normal_breast'
    miR_norm_file <- "_Normal_Breast"
  }
  
  df$Mirbase.ID <- as.character(df$Mirbase.ID)
  
  
  results <- paste0(working_path,"/targetscan70_results.txt")
  writeLines(paste0("Writing results to table: ", results))
  
  write.table(df, 
              file = results, 
              quote = FALSE,
              row.names = FALSE,
              sep = "\t")
  
  if (miRnorm != 0){  # Create another table for expressed miRNA 
    
    df_miRexp <- subset(df,
                        df$Mirbase.ID %in% miR_expression)
    df_miRexp$miR_norm <- miR_norm
    
    results_miR <- paste0(working_path,"/targetscan70_results_miRexpressed",miR_norm_file,".txt")
    writeLines(paste0("Writing results filtered by miRNA expression to table: ", results_miR))
    
    write.table(df_miRexp, 
                file = results_miR, 
                quote = FALSE,
                row.names = FALSE,
                sep = "\t")
  }
}



#### 11. Make plots ####

get.min<- function(list, by.c){
  a1<- list[[1]]
  a11<-min(a1[,by.c])
  a2<- list[[2]]
  a22<-min(a2[,by.c])
  
  if (length(list) == 2){
    z<- c(a11,a22)
    z 
  }
  
  else {
    a3<- list[[3]]
    a33<-min(a3[,by.c])
    
    if (length(list) == 3){
      z<- c(a11,a22,a33)
      z
    }
    else {
      a4<- list[[4]]
      a44<-min(a4[,by.c])
      z<- c(a11,a22,a33,a44)
      z
    }
  }
}

get.max<- function(list, by.c){
  a1<- list[[1]]
  a11<-max(a1[,by.c])
  a2<- list[[2]]
  a22<-max(a2[,by.c])
  
  if (length(list) == 2){
    z<- c(a11,a22)
    z 
  }
  
  else {
    a3<- list[[3]]
    a33<-max(a3[,by.c])
    
    if (length(list) == 3){
      z<- c(a11,a22,a33)
      z
    }
    else {
      a4<- list[[4]]
      a44<-max(a4[,by.c])
      z<- c(a11,a22,a33,a44)
      z
    }
  }
}


fun.diff.cs.plot <- function(df, miR_expression, scale_limit, c.split, path) {
  
  if (missing(path) == TRUE){
    if (exists("working_path") == TRUE)
      path = working_path
    else path = getwd() }
  
  path_cs_diff<- paste0(path,"/Plots_cs_diff/")  
  path_miRexpressed<- paste0(path,"/Plots_cs_diff/miRExpressed/")
  
  if (dir.exists(path_cs_diff) == FALSE)
    dir.create(path_cs_diff)
  
  if (dir.exists(path_miRexpressed) == FALSE)
    dir.create(path_miRexpressed)
  
  
  if (miR_expression == 1){
    miR_expression = normal_breast
    y_lab <- 'miRNA Expressed in Normal Breast'
    y_lab_file <- "_miRNA_Normal_Breast"
  }
  
  
  df_split <- split(df, df[,c.split])
  
  for (i in 1:length(df_split)){
    
    diff <- df_split[[i]]  # get df
    
    # Get constants
    crit <- diff[1,c.split]  # criterion
    pos<- as.numeric(diff[1,"position_in_sequence"])
    snp_name<- as.character(diff[1,"snp_id"])
    transcript <- as.character(diff[1,"transcript_id"])
    gene <- as.character(diff[1,"hgnc_symbol"])
    
    name_id<- paste(gene, paste0("(",transcript, ")"), sep = " ")
    
    minor_allele<- as.character(diff[1,"minor_allele"])
    allele_1<- as.character(diff[1,"allele_1"])
    
    info<- paste0("Minor Allele: ", minor_allele, "  |  ", "Ancestral Allele: ", allele_1)
    
    alleles <- diff[1,"alleles"]
    alls <- strsplit(alleles,"/")
    alls <- do.call("c",alls)  # concatenate
    
    dff<- diff[,c("Mirbase.ID","context...score","allele")]
    
    
    # missing data --> NA 
    x <- expand.grid(Mirbase.ID=unique(dff$Mirbase.ID),
                     allele = alls)  # all combinations
    x$pat <- do.call("paste",x)
    x <- as.data.frame(x$pat)
    
    dff$pat <- do.call("paste", dff[,c("Mirbase.ID","allele")])
    
    y <- merge(dff[,c("pat","context...score")], x, 
               by.x="pat",
               by.y="x$pat",
               all=TRUE)
    
    patsplit <- strsplit(y$pat, " ")  # split 'pat'
    
    for (i in 1:length(patsplit)){  # 'pat' --> 'Mirbase.ID' and 'allele'
      y$Mirbase.ID[i] <- patsplit[[i]][1]
      y$allele[i] <- patsplit[[i]][2]
    }
    
    ggplot(y, aes(x= allele, y = Mirbase.ID, fill = context...score)) + 
      geom_tile() + 
      geom_text(aes(label = round(context...score, 3)))  + 
      ylab('miRNA') + xlab('Allele') + 
      labs(title = snp_name, subtitle = name_id, caption = info) + 
      coord_equal() + 
      theme_classic() + 
      theme(plot.title = element_text(face = "bold"),
            plot.caption = element_text(hjust = 0.5, size = 7)) + 
      scale_fill_gradient(low = "red", high = "white", name = "Context++ \nScore", limits = scale_limit)
    
    
    ggsave(paste0(path_cs_diff,crit, ".png"), width = 11.7, height = 8.3)
    ggsave(paste0(path_cs_diff,crit, ".pdf"), width = 11.7, height = 8.3)
    
    if (miR_expression != 0) {
      miRexpressed <- subset(y, y$Mirbase.ID %in% miR_expression)
      if (nrow(miRexpressed) != 0){
        ggplot(miRexpressed, aes(x = allele, y = Mirbase.ID, fill = context...score)) + 
          geom_tile() + 
          geom_text(aes(label = round(context...score, 3))) + 
          scale_fill_gradient(low = "red", high = "white",name = "Context++ \nScore", limits=scale_limit) + 
          ylab(y_lab) + xlab('Allele') + 
          labs(title = snp_name, subtitle = name_id, caption = info) + 
          coord_equal() + theme_classic() + 
          theme(plot.title = element_text(face = "bold"),
                plot.caption = element_text(hjust = 0.5, size = 6)) 
        
        ggsave(paste0(path_miRexpressed,crit,y_lab_file, ".png"), width = 11.7, height = 8.3)
        ggsave(paste0(path_miRexpressed,crit,y_lab_file, ".pdf"), width = 11.7, height = 8.3)
      }
    }
  }
}
