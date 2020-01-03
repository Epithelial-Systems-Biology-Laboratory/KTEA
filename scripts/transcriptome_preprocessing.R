library(readxl)
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


out_file <- ifelse(exists("snakemake"),snakemake@output[[1]],"data/transcriptome.txt")
tmpfile <- tempfile(fileext = ".xlsx")
download.file("https://helixweb.nih.gov/ESBL/Database/NephronRNAseq/Supplemental_Data_1.xlsx",
              tmpfile)

transcriptome <- read_excel(tmpfile,sheet = "Median for each segment", skip = 7) %>%
  dplyr::select(-c("Maximum", "Variance", "Annotation")) %>%
  mutate_if(is.numeric, ~na_if(.x,0)) %>%
  filter_if(is.numeric, any_vars(!is.na(.))) %>%
  gather(segment, rpkm, -`Gene Symbol`) %>%
  group_by(`Gene Symbol`, segment) %>%
  summarise(rpkm=sum(rpkm)) %>%
  ungroup() %>%
  mutate(rpkm = replace_na(rpkm, 0)) %>%
  mutate(log2_rpkm_plus_one = log2(rpkm + 1),
         segment=str_remove(segment," .*")) %>%
  mutate(segment = fct_relevel(segment, tubule_segments)) %>%
  rename(Gene_symbol=`Gene Symbol`) %>%
  arrange(Gene_symbol, segment) %>%
  mutate(segment = fct_recode(segment,
                              DTL1="SDL",
                              DTL2="LDLOM",
                              DTL3="LDLIM",
                              ATL="tAL")) # Rename segment to new nomenclature https://doi.org/10.1681/ASN.2019040415

# See https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/ for conversion
transcriptome <- transcriptome %>%
  group_by(segment) %>%
  mutate(tpm = 1e6*rpkm/sum(rpkm),
         log2_tpm_plus_one = log2(tpm + 1),
         rpkm = na_if(rpkm, 0),
         rank = rank(rpkm, ties.method = "min"),
         percentile = replace_na(percent_rank(rpkm),0)) %>%
  mutate(rpkm = replace_na(rpkm, 0)) %>%
  ungroup()

write_csv(transcriptome, out_file)
