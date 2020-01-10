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


transcription_factor <- read_tsv(snakemake@input[['transcription_factor_list']])
proteome <- read_csv(snakemake@input[['proteome']]) %>% mutate(segment=fct_relevel(segment,tubule_segments))

transcription_factor <- proteome %>% 
  filter(gene %in% transcription_factor$`Gene Symbol`) %>% 
  select(gene,segment,geomean_copy_number) %>% 
  group_by(gene,segment) %>% 
  summarise(copy_number=first(geomean_copy_number)) %>%
  group_by(gene) %>% 
  mutate(protein_scaled = copy_number/max(copy_number)) %>% 
  select(-copy_number) %>% 
  spread(segment, protein_scaled) %>% 
  left_join(transcription_factor, by=c("gene"="Gene Symbol")) %>% 
  ungroup()



transcription_factor %>% 
  gather(segment,value, -gene, -`TF Family`) %>% 
  mutate(segment=case_when(segment=="CTAL" ~ "cTAL",
                           segment=="MTAL" ~ "mTAL",
                           TRUE ~ segment),
         gene=fct_inorder(gene)) %>% 
  ggplot(aes(x=fct_relevel(segment,tubule_segments),fct_rev(gene),fill=value)) + 
  geom_tile(show.legend = FALSE) + 
  geom_text(aes(label=round(value,2)), size = 7/.pt) +
  geom_label(aes(label=`TF Family`, x=16.5),size=7/.pt,label.padding = unit(.2,"lines"), fill=NA) +
  theme_minimal() +
  scale_fill_gradient(low = "#FFFFFF",high = "#FFD300",aesthetics = "fill") +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = margin(r=50, t=20, b=5),
        panel.grid = element_blank(),
        axis.text = element_text(face = "bold", size=8)) +
  scale_x_discrete(position = "top") + xlab("") + ylab("") +
  annotation_custom(grid::rectGrob(gp = grid::gpar(fill="#1D6996")), 
                    xmin=0.5,xmax=3.5,ymin=24.7, ymax=25.7) +
  annotation_custom(grid::textGrob( label="Proximal Tubules",gp = grid::gpar(fontsize=8, fontface="bold",col="white")), 
                    xmin=0.5,xmax=3.5,ymin=24.7, ymax=25.7)+
  annotation_custom(grid::rectGrob(gp = grid::gpar(fill="#0F8554")), 
                    xmin=3.5,xmax=7.5,ymin=24.7, ymax=25.7) +
  annotation_custom(grid::textGrob( label="Thin limbs",gp = grid::gpar(fontsize=8, fontface="bold",col="white")), 
                    xmin=3.5,xmax=7.5,ymin=24.7, ymax=25.7)+
  annotation_custom(grid::rectGrob(gp = grid::gpar(fill="#EDAD08")), 
                    xmin=7.5,xmax=9.5,ymin=24.7, ymax=25.7) +
  annotation_custom(grid::textGrob( label="Thick limbs",gp = grid::gpar(fontsize=8, fontface="bold",col="white")), 
                    xmin=7.5,xmax=9.5,ymin=24.7, ymax=25.7)+
  annotation_custom(grid::rectGrob(gp = grid::gpar(fill="#CC503E")), 
                    xmin=9.5,xmax=10.5,ymin=24.7, ymax=25.7) +
  annotation_custom(grid::textGrob( label="DCT",gp = grid::gpar(fontsize=8, fontface="bold",col="white")), 
                    xmin=9.5,xmax=10.5,ymin=24.7, ymax=25.7)+
  annotation_custom(grid::rectGrob(gp = grid::gpar(fill="#6F4070")), 
                    xmin=10.5,xmax=14.5,ymin=24.7, ymax=25.7) +
  annotation_custom(grid::textGrob( label="Collecting duct",gp = grid::gpar(fontsize=8, fontface="bold",col="white")), 
                    xmin=10.5,xmax=14.5,ymin=24.7, ymax=25.7)

ggsave(snakemake@output[[1]], width = large.width)
dev.off()