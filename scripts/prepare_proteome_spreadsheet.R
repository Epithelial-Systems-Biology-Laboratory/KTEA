library(tidyverse)

tubule_segments <- c("S1",
                     "S2",
                     "S3",
                     "DTL1",
                     "DTL2",
                     "DTL3",
                     "ATL",
                     "mTAL",
                     "cTAL",
                     "DCT",
                     "CNT",
                     "CCD",
                     "OMCD",
                     "IMCD")

if (exists("snakemake")){
  proteome <- read_csv(snakemake@input[[1]]) %>% mutate(segment=fct_relevel(segment,tubule_segments))
  outfile <- snakemake@output[[1]]
} else {
  proteome <- read_csv("analysis/analysis_workflow/data/KTEA_proteome_processed.csv") %>% mutate(segment=fct_relevel(segment,tubule_segments))
  outfile <- "analysis/analysis_workflow/out/proteome_average_copy_number.txt"
}

proteome %>% 
  group_by(`Majority protein IDs`, `Protein names`, gene, segment) %>% 
  arrange(`Majority protein IDs`) %>% 
  summarise(geomean_copy_number = first(geomean_copy_number)) %>% 
  ungroup() %>% 
  spread(segment, geomean_copy_number) %>% 
  arrange(gene) %>% 
  write_tsv(outfile)


