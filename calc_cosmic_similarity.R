library(tidyverse)
library(data.table)
library(stringr)

#
# target dir
#
#lego96_dir <- "OUTPUT_lego96"
args <- commandArgs(trailingOnly = T)
lego96_dir <- args[1]

#
# fetch COSMIC v2
#
#cosmic_vs <- fread("signatures_probabilities.txt")
cosmic_v2 <- fread("https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt")
cosmic_v2 <- cosmic_v2 %>% 
  mutate(name = paste0(str_sub(`Substitution Type`,1,1),
                       str_sub(`Substitution Type`,3,3),
                       str_sub(Trinucleotide, 1, 1),
                       str_sub(Trinucleotide, 3, 3)))
cosmic_v2_mat <- cosmic_v2 %>%
  select(starts_with("Signature ")) %>% as.matrix
rownames(cosmic_v2_mat) <- cosmic_v2$name

#
# fetch COSMIC v3
#
#cosmic_v3 <- fread("sigProfiler_SBS_signatures_2019_05_22.csv")
cosmic_v3 <- fread("https://dcc.icgc.org/api/v1/download?fn=/PCAWG/mutational_signatures/Signatures/SP_Signatures/SigProfiler_reference_signatures/SigProfiler_reference_whole-genome_signatures/sigProfiler_SBS_signatures_2019_05_22.csv")
cosmic_v3 <- cosmic_v3 %>%
  mutate(name = paste0(str_sub(Type,1,1),
                       str_sub(Type,3,3),
                       str_sub(SubType, 1, 1),
                       str_sub(SubType, 3, 3)))
cosmic_v3_mat <- cosmic_v3 %>%
  select(starts_with("SBS")) %>% as.matrix
rownames(cosmic_v3_mat) <- cosmic_v3$name

#
# files that store de novo signatures
#
glob <- file.path(lego96_dir, "L1KL.lego96.*.10.MAP*.WH.RData")
denovo_files <- Sys.glob(glob)

#
# function for computing cosine similarity
#
cosine_similarity <- function(row, WH1, cosmic_mat) {
  nameA <- colnames(WH1)[row$Var1]
  nameB <- colnames(cosmic_mat)[row$Var2]
  A <- WH1[,row$Var1]
  B <- cosmic_mat[,row$Var2]
  c <- (A %*% B) / sqrt((A %*% A) * (B %*% B))
  data.frame(denovo=nameA, COSMIC=nameB, cosine_similarity=c)
}

#
# function for comparing de novo signature and COSMIC
#
compare_denovo_vs_cosmic <- function(row){
  match <- str_match(row$denovo_file, "MAP([0-9]+)\\.WH\\.RData$")
  num_sig <- match[2]

  # load de novo signature
  load(row$denovo_file)
  WH1 <- WH[[1]]
  colnames(WH1) <- paste0("W", 1:ncol(WH1))

  if(row$cosmic_version == "v2"){
    cosmic_mat <- cosmic_v2_mat
  } else {
    cosmic_mat <- cosmic_v3_mat
  }

  # align rows
  cosmic_mat <- cosmic_mat[rownames(WH1),]
  
  cs <- expand.grid(1:ncol(WH1), 1:ncol(cosmic_mat)) %>%
    rowwise %>%
    do(cosine_similarity(., WH1, cosmic_mat)) %>%
    ungroup()
  cs$Num_denovo_signatures <- num_sig
  cs$COSMIC_version = row$cosmic_version
  return(cs)
}

#
# process pairs of de novo signatures and COSMIC v2/v3
#
res <- expand.grid(denovo_file = denovo_files, cosmic_version = c("v2", "v3"), stringsAsFactors = F) %>%
  rowwise %>%
  do(compare_denovo_vs_cosmic(.)) %>%
  ungroup %>%
  select(Num_denovo_signatures, denovo, COSMIC_version, COSMIC, cosine_similarity) %>%
  write_csv("calc_cosmic_similarity.csv")
