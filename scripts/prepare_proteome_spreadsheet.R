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

new_sample_label <- proteome %>% 
  select(`Majority protein IDs`,`Protein names`,gene,sample,segment,copy_number) %>% 
  distinct(sample,segment) %>% 
  group_by(segment) %>% 
  mutate(n=row_number()) %>% 
  unite("new_sample_label", segment,n,sep = "_")

full_individual_data <- proteome %>% 
  arrange(segment) %>% 
  select(`Majority protein IDs`,`Protein names`,gene,sample,copy_number) %>% 
  pivot_wider(names_from = sample,values_from = copy_number) 
  
new_column_names <- full_individual_data %>% 
  colnames() %>% 
  enframe(name = NULL) %>% 
  left_join(new_sample_label, by=c("value" = "sample")) %>% 
  mutate(new_sample_label = if_else(is.na(new_sample_label),value,new_sample_label)) %>% 
  pull(new_sample_label)

colnames(full_individual_data) <- new_column_names

write_tsv(full_individual_data, snakemake@output[[2]])
