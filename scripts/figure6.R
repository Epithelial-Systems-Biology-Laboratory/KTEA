library(tidyverse)
library(ggpubr)


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



# correlated examples -----------------------------------------------------

plot_prot_mrna <- function(gene_to_plot){
  proteome %>% 
    group_by(gene,segment) %>% 
    summarise(geomean_copy_number=first(geomean_copy_number)) %>% 
    ungroup() %>% 
    inner_join(transcriptome,by=c("gene"="Gene_symbol","segment"="segment")) %>% 
    select(gene,segment,geomean_copy_number,tpm) %>% 
    gather(key,value,-gene,-segment) %>% 
    mutate(key = case_when(key=="geomean_copy_number" ~ "Protein",
                           key=="tpm" ~ "RNA")) %>% 
    filter(gene==gene_to_plot) %>% 
    ggplot(aes(x=segment,y=value, group=key)) + geom_line() + 
    facet_wrap(~key, nrow = 2, scales = "free_y",strip.position = "right") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,size=6),
          axis.text.y = element_text(size=6),
          strip.text = element_text(size=6,face = "bold")) +
    scale_y_continuous(labels=scales::label_number_si(accuracy = 0.1)) +
    ggtitle(gene_to_plot)
}



figure6a <- map(list("Slc5a2","Slc34a1","Slc12a1","Slc12a3","Aqp2"), plot_prot_mrna)
figure6b <- map(list("Maoa","Maob","Calr","Bsg","Abcb8"), plot_prot_mrna)
fig6 <- ggarrange(ggarrange(plotlist = figure6a, ncol = 1, nrow = 5),
          ggarrange(plotlist = figure6b, ncol = 1, nrow = 5),
          labels = c("A","B"))
ggsave(filename = snakemake@output[[1]], width = large.width, height = 8)
#dev.off()