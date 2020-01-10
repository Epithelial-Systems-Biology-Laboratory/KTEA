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
mq_proteinGroups <- read_tsv(snakemake@input[['proteinGroups']], guess_max = 8000)
cell_per_length <- read_tsv(snakemake@input[['cell_per_mm']])

#tubule_length_included <- read_tsv("experiment/tubule_length_included.txt") %>% 
tubule_length_included <- read_tsv(snakemake@input[['tubule_length_included']]) %>% 
  mutate(segment = case_when(segment == "SDL" ~ "DTL1",
                             segment == "LDLOM" ~ "DTL2",
                             segment == "LDLIM" ~ "DTL3",
                             segment == "tAL" ~ "ATL",
                             TRUE ~ segment),
         segment=fct_relevel(segment,tubule_segments))

#cell_per_length <- read_tsv("MK_analysis/cell_per_mm.txt")

n_cell <- tubule_length_included %>% 
  group_by(segment,Exp_name) %>% 
  summarise(sum_length = sum(length)) %>% 
  left_join(cell_per_length) %>% 
  drop_na() %>% 
  mutate(n_cell = sum_length*cell_per_mm)


histone <- proteome %>% 
  filter(Histone) %>% 
  pull(`Majority protein IDs`) %>% 
  unique()



histone_vs_ncell <- mq_proteinGroups %>% 
  filter(`Majority protein IDs`%in% histone) %>% 
  select(matches("Intensity ")) %>% 
  summarise_all(sum) %>% 
  gather(sample, histone_signal) %>% 
  mutate(sample = str_remove(sample,"Intensity ")) %>% 
  left_join(n_cell, by=c("sample"="Exp_name"))

lm_summary <- summary(lm(log2(histone_signal) ~ n_cell,data = histone_vs_ncell)) 

histone_vs_ncell %>% 
  drop_na() %>% 
  mutate(segment=fct_relevel(segment,tubule_segments)) %>% 
  ggplot(aes(x=n_cell,y=log10(histone_signal))) +
  geom_point(aes(col=segment)) +
  geom_smooth(method="lm") +
  geom_text(size=6/.pt,
            x=8000,y=11,
            label=paste("Adjusted~R^{2}==",round(lm_summary$adj.r.squared,2),sep="~"),
            parse = TRUE)+
  geom_text(size=6/.pt,
            x=8000,y=11.1,
            label=paste("p-value==",signif(lm_summary$coefficients[2,4],digits = 2),sep="~"),
            parse = TRUE)+
  theme_bw() +
  scale_color_viridis_d(end = 0.95,option = "C") +
  labs(x="Estimated number of cell",
       y=expression(Log[10](Histone~signal))) +
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=6),
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.position = "top", 
        legend.justification = c(0,0),
        legend.key.size = unit(0.5,"points"),
        legend.spacing = unit(1.5,"mm"),
        legend.spacing.x = unit(1.5,"mm"),
        legend.margin = margin())+
  guides(col=guide_legend(title.position = "left", nrow = 2))

ggsave(snakemake@output[[1]], width = small.width, height = small.width)
dev.off()