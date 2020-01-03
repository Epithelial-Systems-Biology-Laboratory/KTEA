library(tidyverse)

tubule_segments <- c("S1",
                     "S2",
                     "S3",
                     "SDL",
                     "LDLOM",
                     "LDLIM",
                     "tAL",
                     "mTAL",
                     "cTAL",
                     "DCT",
                     "CNT",
                     "CCD",
                     "OMCD",
                     "IMCD")

proteome <- read_csv(snakemake@input[[1]]) 

##### Atp8a1, Cask, Cux1, Gorasp2, Mfsd4b, Mia2, Muc4 --> these genes each have two separate protein groups. 
# These are caused by either redundancy in reference proteome or outdated gene symbol from MaxQuant.
# After manual inspection, seven rows will be excluded (keep dominant rows).
dup_prot_ids <- read_lines(snakemake@input[[2]])
proteome <- proteome %>% filter(!(`Majority protein IDs` %in% dup_prot_ids))

#When there is no gene symbol, fill with Uniprot accession (first of majority protein ids).

proteome[is.na(proteome$final_gene),]$final_gene <- proteome[is.na(proteome$final_gene),]$`Majority protein IDs` %>% 
  str_split(";", simplify = TRUE) %>% 
  magrittr::extract(,1)

proteome <- proteome %>%
  mutate(gene = final_gene) %>%
  as_tibble() %>%
  select(-contains("Sequence"), - contains("Razor")) %>%
  gather(segment, copy_number, contains("copy_number_histone")) %>%
  mutate(sample = str_remove(segment, "copy_number_histone_Intensity "),
         segment = fct_relevel(str_match(segment, " (.*)_[:digit:]")[,2] ,tubule_segments)) %>%
  # Remove protein without any quantification value
  group_by(gene) %>% 
  mutate(sum=sum(copy_number)) %>% 
  filter(sum>0) %>% 
  ungroup() %>% 
  # Calculate average
  group_by(`Majority protein IDs`, segment) %>%
  mutate(median_copy_number = median(copy_number),
         mean_copy_number= mean(copy_number),
         geomean_copy_number = exp(mean(log(copy_number+ 1))) -1) %>%
  ungroup() %>%
  group_by(sample) %>%
  # Calculate ranking
  mutate(rank = rank(geomean_copy_number, ties.method = "min"),
         percent_rank = percent_rank(na_if(geomean_copy_number, 0))) %>%
  mutate(percent_rank = replace_na(percent_rank,0)) %>% 
  ungroup() %>%
  arrange(segment)

#write_csv(proteome, snakemake@output[[1]])

new_gene_mapping <- proteome %>% 
  filter(str_detect(gene, "LOC[:digit:]|RGD[:digit:]")) %>%
  select(`Gene names`, gene) %>% 
  distinct() %>% 
  mutate(`Gene names` = str_remove_all(`Gene names`,"LOC[:digit:]*;?|RGD[:digit:]*;?" )) %>% 
  mutate(`Gene names` = na_if(`Gene names`, "")) %>% 
  unite("combined",c("Gene names", "gene"), sep = ";", remove = FALSE, na.rm = TRUE) %>% 
  mutate(new_gene = str_split(combined,";",simplify = TRUE)[,1]) %>% 
  mutate(new_gene = case_when(str_detect(new_gene, "LOC[:digit:]*;?|RGD[:digit:]*;?") ~ NA_character_, TRUE ~ new_gene)) %>% 
  drop_na(new_gene) %>% 
  select(gene, new_gene)


tmp <- proteome %>% left_join(new_gene_mapping, by = "gene") 
tmp[!is.na(tmp$new_gene),"gene"] <- tmp[!is.na(tmp$new_gene),"new_gene"]
tmp %>% 
  select(-`Gene names`, -final_gene, -id, -new_gene) %>% 
  mutate(segment = fct_recode(segment,
                              DTL1="SDL", 
                              DTL2="LDLOM", 
                              DTL3="LDLIM",
                              ATL="tAL")) %>% # Rename segment to new nomenclature https://doi.org/10.1681/ASN.2019040415
  write_csv(snakemake@output[[1]])