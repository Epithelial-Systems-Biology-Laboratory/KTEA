library(tidyverse)

## See https://jasn.asnjournals.org/sites/default/files/digartguid.pdf
# Height < 55picas/ 9.2" / 24 cm
# Widths : 1 column 3.1"/7.8cm; 1.5 columns 4.5"/11.4cm; 2 columns 6.75"/17.1cm

small.width <- 3.1
medium.width <- 4.5
large.width <- 6.75
page.height <- 9.2


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

proteome <- read_csv(snakemake@input[['proteome']]) %>% mutate(segment=fct_relevel(segment,tubule_segments))
transcriptome <- read_csv(snakemake@input[['transcriptome']]) %>% mutate(segment=fct_relevel(segment,tubule_segments))


# protein-transcript correlation ------------------------------------------

df_prot <- proteome %>% 
  group_by(gene,segment) %>% 
  summarise(log2_protein=log2(mean(geomean_copy_number)+1)) %>% 
  mutate(protein_scaled=scale(log2_protein)) %>% 
  ungroup()

df_trans <- transcriptome %>% 
  group_by(Gene_symbol) %>% 
  mutate(tpm_scaled=scale(log2_tpm_plus_one)) %>% 
  select(Gene_symbol,segment,log2_tpm_plus_one,tpm_scaled) %>% 
  mutate(tpm_scaled = replace_na(tpm_scaled,0)) %>% 
  ungroup()

prot_trans_correlation <- left_join(df_prot,df_trans,by = c("gene"="Gene_symbol", "segment"="segment")) %>% 
  mutate(difference = protein_scaled-tpm_scaled) %>% 
  group_by(gene) %>% 
  mutate(correlation_pearson = cor(protein_scaled,tpm_scaled, method = "pearson"),
         correlation_spearman = cor(protein_scaled, tpm_scaled, method = "spearman")) %>% 
  summarise(correlation_pearson = first(correlation_pearson),
            correlation_spearman = first(correlation_spearman),
            avg_log2_prot = mean(log2_protein),
            avg_log2_tpm = mean(log2_tpm_plus_one)) %>% 
  arrange(correlation_pearson) 

dens <- density(prot_trans_correlation$correlation_spearman,na.rm = TRUE)
df <- data.frame(x=dens$x, y=dens$y)
df$quant <- factor(findInterval(df$x,0))
ggplot(df, aes(x,y)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) + 
  scale_fill_brewer(guide="none",palette = "Dark2") +
  theme_bw() +
  xlab("Gene-wise protein-mRNA Spearman correlation") +
  ylab("Density") +
  geom_vline(xintercept = median(prot_trans_correlation$correlation_spearman,na.rm = TRUE), linetype=2, color="white") + 
  geom_text(x=-0.6,y=1.1,color="black", size=10/.pt, 
            aes(label=paste0("Median: ",
                             round(median(prot_trans_correlation$correlation_spearman,na.rm = TRUE),2),
                             "\nN = ", nrow(drop_na(prot_trans_correlation))))) +
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8)) 
ggsave(snakemake@output[[1]], width = small.width, height = small.width)
dev.off()