library(tidyverse)
library(limma)

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


proteome_unscaled <- proteome %>%
  mutate(copy_number = log2(copy_number+1)) %>%
  select(`Majority protein IDs`, gene, sample, copy_number) %>%
  group_by(`Majority protein IDs`) %>%
  spread(sample, copy_number) %>%
  as.data.frame()
rownames(proteome_unscaled) <- proteome_unscaled$gene

segment_order <- (colnames(proteome_unscaled[,-(1:2)]) %>% str_split_fixed("_",2) )[,1]
design <- model.matrix(~0+segment_order)
colnames(design) <- colnames(design) %>% str_remove_all("segment_order|\\)|\\(")

# Change segment nomenclature of descending limb
colnames(design) <- case_when(colnames(design)=="SDL"~"DTL1",
                              colnames(design)=="LDLOM"~"DTL2",
                              colnames(design)=="LDLIM"~"DTL3",
                              colnames(design)=="tAL"~"ATL",
                              TRUE ~ colnames(design))

# Compare each segment with the average of all other segments
contrast <- makeContrasts(contrasts=c("S1-(S2+S3+DTL1+DTL2+DTL3+ATL+mTAL+cTAL+DCT+CNT+CCD+OMCD+IMCD)/13",
                                      "S2-(S1+S3+DTL1+DTL2+DTL3+ATL+mTAL+cTAL+DCT+CNT+CCD+OMCD+IMCD)/13",
                                      "S3-(S1+S2+DTL1+DTL2+DTL3+ATL+mTAL+cTAL+DCT+CNT+CCD+OMCD+IMCD)/13",
                                      "DTL1-(S1+S2+S3+DTL2+DTL3+ATL+mTAL+cTAL+DCT+CNT+CCD+OMCD+IMCD)/13",
                                      "DTL2-(S1+S2+S3+DTL1+DTL3+ATL+mTAL+cTAL+DCT+CNT+CCD+OMCD+IMCD)/13",
                                      "DTL3-(S1+S2+S3+DTL1+DTL2+ATL+mTAL+cTAL+DCT+CNT+CCD+OMCD+IMCD)/13",
                                      "ATL-(S1+S2+S3+DTL1+DTL2+DTL3+mTAL+cTAL+DCT+CNT+CCD+OMCD+IMCD)/13",
                                      "mTAL-(S1+S2+S3+DTL1+DTL2+DTL3+ATL+cTAL+DCT+CNT+CCD+OMCD+IMCD)/13",
                                      "cTAL-(S1+S2+S3+DTL1+DTL2+DTL3+ATL+mTAL+DCT+CNT+CCD+OMCD+IMCD)/13",
                                      "DCT-(S1+S2+S3+DTL1+DTL2+DTL3+ATL+mTAL+cTAL+CNT+CCD+OMCD+IMCD)/13",
                                      "CNT-(S1+S2+S3+DTL1+DTL2+DTL3+ATL+mTAL+cTAL+DCT+CCD+OMCD+IMCD)/13",
                                      "CCD-(S1+S2+S3+DTL1+DTL2+DTL3+ATL+mTAL+cTAL+DCT+CNT+OMCD+IMCD)/13",
                                      "OMCD-(S1+S2+S3+DTL1+DTL2+DTL3+ATL+mTAL+cTAL+DCT+CNT+CCD+IMCD)/13",
                                      "IMCD-(S1+S2+S3+DTL1+DTL2+DTL3+ATL+mTAL+cTAL+DCT+CNT+CCD+OMCD)/13",
                                      "(S1+S2+S3)/3 - (DTL1+DTL2+DTL3+ATL+mTAL+cTAL+DCT+CNT+CCD+OMCD+IMCD)/11",
                                      "(mTAL+cTAL)/2 - (S1+S2+S3+DTL1+DTL2+DTL3+ATL+DCT+CNT+CCD+OMCD+IMCD)/12",
                                      "(CCD+OMCD+IMCD)/3 - (S1+S2+S3+DTL1+DTL2+DTL3+ATL+mTAL+cTAL+DCT+CNT)/11",
                                      "(DTL1+DTL2+DTL3+ATL+mTAL+OMCD+IMCD)/7 - (S1+S2+S3+cTAL+DCT+CNT+CCD)/7"),
                          levels=colnames(design))
colnames(contrast) <- c(tubule_segments, "Proximal_tubule", "Thick_ascending_limb", "Collecting_duct", "Medulla_VS_cortex")


fit <- lmFit(proteome_unscaled[,-(1:2)], design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

dat <- proteome %>%
  group_by(`Majority protein IDs`,gene,segment) %>%
  summarise(geomean_copy_number = log10(first(geomean_copy_number)+1)) %>%
  ungroup() %>%
  group_by(gene) %>%
  mutate(scaled_cp = scale(geomean_copy_number))    %>%
  ungroup() %>%
  select(-geomean_copy_number) %>%
  group_by(`Majority protein IDs`) %>%
  spread(segment, scaled_cp) %>%
  ungroup() %>%
  as.data.frame()
rownames(dat) <- dat$gene



## Grouped segments
lfc_threshold <- 2
n_protein_to_plot <- 20
p_threshold <- 0.01
de_genes <- map(15:18,
                ~topTable(fit2,
                          coef = .x,
                          number = 10000,
                          sort.by = "logFC",
                          resort.by = "logFC",
                          p.value = p_threshold,
                          lfc = lfc_threshold) %>%
                  rownames_to_column("gene") %>%
                  filter(logFC > 0) %>%
                  head(n_protein_to_plot) %>%
                  pull(gene))  %>%
  unlist() %>%
  unique()
selected <- dat$gene %in% de_genes

de_genes %>%
  enframe(name = NULL, value = "gene") %>%
  left_join(dat, by="gene") %>%
  gather(segment, z_score, -gene,-`Majority protein IDs`) %>%
  mutate(segment = fct_inorder(segment),
         gene = fct_inorder(gene)) %>%
  ggplot() + 
  geom_tile(aes(x=segment, y=gene, fill=z_score)) +
  scale_fill_viridis_c(option = "B") +
  theme_bw() +
  guides(fill=guide_colorbar(ticks = FALSE,
                             title.position = "top",
                             title.hjust = 0.5, 
                             title = "z-score",label.position = "top",
                             frame.colour = "black",
                             barwidth = 0.45,
                             barheight = 0.01,
                             default.unit = "npc")) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle=90, vjust = 0.5, size=7),
        axis.text.y = element_text(size = 5),
        axis.title = element_blank(),
        legend.position = "top",
        legend.box.margin = margin(t=-5,b = -10),
        legend.background = element_blank(),
        legend.text = element_text(size=5),
        legend.title = element_text(size = 7)) +
  annotate("rect", xmin=-0.5,xmax=0.5,ymin=0.5,ymax=20.5 , alpha=0.5, fill = "#fdb462") +
  annotate("text", x=-0.5,y=20.5/2,label="Proximal Tubules",angle=90, vjust=1.5, size = 8/.pt) +
  annotate("rect", xmin=-0.5,xmax=0.5,ymin=20.5,ymax=40.5 , alpha=0.5, fill = "#80b1d3") +
  annotate("text", x=-0.5,y=20.5 + 20.5/2,label="Thick ascending limb",angle=90, vjust=1.5, size = 8/.pt) +
  annotate("rect", xmin=-0.5,xmax=0.5,ymin=40.5,ymax=60.5 , alpha=0.5, fill = "#fb8072") +
  annotate("text", x=-0.5,y=40.5 + 20.5/2,label="Collecting ducts",angle=90, vjust=1.5, size = 8/.pt) + 
  annotate("rect", xmin=-0.5,xmax=0.5,ymin=60.5,ymax=80.5 , alpha=0.5, fill = "#bebada") +
  annotate("text", x=-0.5,y=60.5 + 20.5/2,label="Medulla",angle=90, vjust=1.5, size = 8/.pt) 

ggsave(snakemake@output[[1]], width = small.width)
dev.off()