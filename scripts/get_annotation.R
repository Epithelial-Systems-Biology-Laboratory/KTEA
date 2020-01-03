
library(httr)
library(tidyverse)

# download uniprot annotation ---------------------------------------------

out_file <- ifelse(exists("snakemake"),snakemake@output[['uniprot_annotation']], "data/uniprot_annotation.txt")

proteome <- read_csv(snakemake@input[[1]])

url <- "https://www.uniprot.org/uploadlists/"
query <- paste(unique(proteome$`Majority protein IDs`), collapse = ";")

r <- POST(url,
          body=list(from="ACC",
                    to="ACC",
                    columns="id,database(InterPro),comment(SUBCELLULAR LOCATION),go(cellular component),interactor,comment(TISSUE SPECIFICITY),comment(PATHWAY),go(biological process),go(molecular function),genes(PREFERRED),protein names,keywords",
                    format="tab",
                    query=(proteome %>% 
                             distinct(`Majority protein IDs`)%>% 
                             pull(`Majority protein IDs`) %>% 
                             paste0(collapse = ";"))))

uniprot_annotation <- content(r, type = "text/tab-separated-values")
write_tsv(uniprot_annotation, out_file)

seventmr <- read_lines("https://www.uniprot.org/docs/7tmrlist.txt",
                       skip_empty_rows = TRUE) 

seventmr_uniprot <- seventmr[seventmr %>% str_detect("Human|Mouse|Rat")] %>% 
  str_split("\\(|\\)", simplify = T) %>% 
  magrittr::extract(,2) %>% 
  unique()

query <- paste(seventmr_uniprot, collapse = ";")
Sys.sleep(1)
r <- POST(url,
          body=list(from="ACC",
                    to="ACC",
                    columns="id,genes(PREFERRED),protein names,keywords",
                    format="tab",
                    query=query))

seventmr_genes <- content(r, type = "text/tab-separated-values") %>% 
  mutate(`Gene names  (primary )`=tolower(`Gene names  (primary )`)) %>% 
  distinct(`Gene names  (primary )`) %>% 
  drop_na(`Gene names  (primary )`) %>% 
  pull(`Gene names  (primary )`) %>% 
  str_split(";") %>% 
  unlist() %>% 
  str_to_sentence()

saveRDS(seventmr_genes, snakemake@output[['gpcr_genes']])

kinases <- read_lines("https://www.uniprot.org/docs/pkinfam.txt",skip_empty_rows = TRUE) 
kinases <- kinases[kinases %>% str_detect("_HUMAN|_MOUSE")] %>% 
  str_split("\\(|\\)", simplify = TRUE) %>% 
  magrittr::extract(,c(2,4)) %>% 
  str_trim() %>% 
  unique()
query <- paste(kinases, collapse = ";")
Sys.sleep(1)
r <- POST(url,
          body=list(from="ACC",
                    to="ACC",
                    columns="id,genes(PREFERRED),protein names,keywords",
                    format="tab",
                    query=query))
kinases_genes <- content(r, type = "text/tab-separated-values") %>% 
  mutate(`Gene names  (primary )`=tolower(`Gene names  (primary )`)) %>% 
  distinct(`Gene names  (primary )`) %>% 
  drop_na(`Gene names  (primary )`) %>% 
  pull(`Gene names  (primary )`) %>% 
  str_split(";") %>% 
  unlist() %>% 
  str_to_sentence()

saveRDS(kinases_genes, snakemake@output[['kinase_genes']])

query <- paste(seventmr_genes, collapse = ";")
r <- POST(url,
          body=list(from="GENENAME",
                    taxon="10116",
                    to="ACC",
                    columns="id",
                    format="tab",
                    query=query))
content(r, type = "text/tab-separated-values") %>% 
  pull(Entry) %>% 
  saveRDS(snakemake@output[['gpcr_accs']])


query <- paste(kinases_genes, collapse = ";")
Sys.sleep(1)
r <- POST(url,
          body=list(from="GENENAME",
                    taxon="10116",
                    to="ACC",
                    columns="id",
                    format="tab",
                    query=query))
content(r, type = "text/tab-separated-values") %>% 
  pull(Entry) %>% 
  saveRDS(snakemake@output[['kinase_accs']])



tmp <- uniprot_annotation %>% 
  filter(str_detect(`Gene ontology (molecular function)`, 
                    regex("transport|symporter|channel|antiporter", 
                          ignore_case = TRUE))) %>% 
  select(`Gene ontology (molecular function)` , 
         Entry,`Gene names  (primary )`) %>% 
  filter(!str_detect(`Gene ontology (molecular function)`, 
                     regex("electron", 
                           ignore_case = TRUE))) %>% 
  pull(Entry) 

transports <- msigdbr::msigdbr(species = "Rattus norvegicus", category = "C5", subcategory = "MF") %>% 
  filter(str_detect(gs_name, regex("transport|symporter|channel|antiporter", ignore_case = TRUE))) %>% 
  distinct(gene_symbol) %>% 
  pull(gene_symbol)
query <- paste(transports, collapse = ";")
Sys.sleep(1)
r <- POST(url,
          body=list(from="GENENAME",
                    taxon="10116",
                    to="ACC",
                    columns="id",
                    format="tab",
                    query=query))
transports_acc <- content(r, type = "text/tab-separated-values") %>% 
  pull(Entry) 

unique(c(tmp,transports_acc)) %>% 
  saveRDS(snakemake@output[['transport_accs']])



# msigdb_gene_to_uniprotacc -----------------------------------------------

query <- msigdbr::msigdbr(species = "Rattus norvegicus") %>% 
  distinct(gene_symbol) %>% 
  pull(gene_symbol) %>% 
  paste0(collapse = ";")
r <- POST(url,
          body=list(from="GENENAME",
                    taxon="10116",
                    to="ACC",
                    columns="id,proteome",
                    format="tab",
                    query=query))
tmp <- content(r, type = "text/tab-separated-values") %>% 
  filter(str_detect(Proteomes, "UP000002494"))
colnames(tmp) <- c("Entry", "Proteomes", "gene")
tmp <- dplyr::select(tmp, Entry, gene)
saveRDS(tmp, snakemake@output[['msigdb_gene2uniprot']])


# msigdb_GOMF -------------------------------------------------------------

msigdbr::msigdbr(species = "Rattus norvegicus", category = "C5", subcategory = "MF") %>% 
  select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(genes=paste(gene_symbol,collapse = ";")) %>% 
  mutate(gs_name=str_remove(gs_name,"GO_")) %>% 
  saveRDS(snakemake@output[['msigdb_GOMF']])


tf_ensembl<- read_tsv("http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Rattus_norvegicus_TF") %>% 
  pull(Ensembl)
query <- paste0(tf_ensembl,collapse = ";")
r <- POST(url,
          body=list(from="ENSEMBL_ID",
                    to="ACC",
                    columns="id,genes(PREFERRED),go(molecular function),proteome",
                    format="tab",
                    query=query))
tf_annotation <- content(r, type = "text/tab-separated-values") %>% 
  filter(str_detect(Proteomes, "UP000002494"))
colnames(tf_annotation) <- c("Entry","gene","GOMF", "Proteomes","ensembl", "isomap")

tf_annotation %>% 
  filter(str_detect(GOMF, "GO:0001228|GO:0003700|GO:0000981|GO:0001227")) %>% 
  pull(Entry) %>% 
  saveRDS(snakemake@output[['tf_acc']])
