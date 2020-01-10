mq_txt_folder = "data/txt"
figure_folder = "out/figures/"
figures = ["figure_1.pdf", "figure_2.pdf", "figure_3.pdf", "figure_4.pdf", "figure_5.pdf","figure_6.pdf"]
supp_figures = ["figure_s1.pdf", "figure_s2.pdf", "figure_s3.pdf"]

rule all:
     input:
        "out/transcriptome.txt",
        "out/fulldataset_individual_sample.txt",
        "out/proteome_average_copy_number.txt",
        "out/fulldataset_individual_sample.txt",
        expand("{figure_folder}{figures}", figures = figures, figure_folder=figure_folder),
        expand("{figure_folder}{supp_figures}", supp_figures = supp_figures, figure_folder=figure_folder),

rule transcriptome_preprocessing:
    output: "out/transcriptome.txt"
    conda: "envs/KTEA.yml"
    script: "scripts/transcriptome_preprocessing.R"

rule proteomic_ruler:
     input:
        mq_txt_folder,
        histone = "data/rat_histones_091918.csv",
        reference_proteome = "Reference_proteome/UP000002494_10116.fasta"
     output: "data/KTEA_proteomic_ruler.csv"
     conda: "envs/KTEA.yml"
     shell: "python scripts/proteomic_ruler.py {input[0]} {output} --histone {input.histone} --fasta {input.reference_proteome}"

rule R_preprocess:
    input: "data/KTEA_proteomic_ruler.csv", "data/duplicated_genes.txt"
    output: "out/KTEA_proteome_processed.csv"
    conda: "envs/KTEA.yml"
    script: "scripts/proteome_preprocessing.R"

rule export_proteome_spreadsheet:
    input: "out/KTEA_proteome_processed.csv"
    output: "out/proteome_average_copy_number.txt",
            "out/fulldataset_individual_sample.txt"
    conda: "envs/KTEA.yml"
    script: "scripts/prepare_proteome_spreadsheet.R"

rule get_annotation:
    input: "data/KTEA_proteomic_ruler.csv"
    output:
        uniprot_annotation="out/uniprot_annotation.txt",
        gpcr_genes="out/7tmr_gene.rds",
        kinase_genes="out/kinase_genes.rds",
        gpcr_accs="out/7tmr_acc.rds",
        kinase_accs="out/kinase_acc.rds",
        transport_accs="out/transport_acc.rds",
        msigdb_gene2uniprot="out/msigdb_gene2uniprot.rds",
        msigdb_GOMF="out/msigdb_GOMF.rds",
        tf_acc="out/tf_filtered_acc.rds"
    conda: "envs/KTEA.yml"
    script: "scripts/get_annotation.R"

rule figure1:
    input:
        proteome="out/KTEA_proteome_processed.csv",
        transporter_list="data/transporter_list.txt"
    output:
        figure_folder + figures[0]
    conda: "envs/KTEA.yml"
    script: "scripts/figure1.R"

rule figure2:
    input:
        proteome="out/KTEA_proteome_processed.csv",
        metabolic_enzyme_list="data/metabolic_enzyme_list.txt"
    output:
        figure_folder + figures[1]
    conda: "envs/KTEA.yml"
    script: "scripts/figure2.R"

rule figure3:
    input:
        proteome="out/KTEA_proteome_processed.csv",
        transcription_factor_list="data/transcription_factor_list.txt"
    output:
        figure_folder + figures[2]
    conda: "envs/KTEA.yml"
    script: "scripts/figure3.R"

rule figure4:
    input:
        proteome="out/KTEA_proteome_processed.csv"
    output:
        figure_folder + figures[3]
    conda: "envs/KTEA.yml"
    script: "scripts/figure4.R"

rule figure5:
    input:
        proteome="out/KTEA_proteome_processed.csv",
        transcriptome="out/transcriptome.txt"
    output:
        figure_folder + figures[4]
    conda: "envs/KTEA.yml"
    script: "scripts/figure5.R"

rule figure6:
    input:
        proteome="out/KTEA_proteome_processed.csv",
        transcriptome="out/transcriptome.txt"
    output:
        figure_folder + figures[5]
    conda: "envs/KTEA.yml"
    script: "scripts/figure6.R"

rule figure_s1:
    input:
        proteome="out/KTEA_proteome_processed.csv",
        proteinGroups="data/txt/proteinGroups.txt",
        cell_per_mm="data/cell_per_mm.txt",
        tubule_length_included="data/tubule_length_included.txt"
    output:
        figure_folder + supp_figures[0]
    conda: "envs/KTEA.yml"
    script: "scripts/figure_s1.R"

rule figure_s2:
    input:
        proteome="out/KTEA_proteome_processed.csv"
    output:
        figure_folder + supp_figures[1]
    conda: "envs/KTEA.yml"
    script: "scripts/figure_s2.R"

rule figure_s3:
    input:
        proteome="out/KTEA_proteome_processed.csv",
        transcriptome="out/transcriptome.txt"
    output:
        figure_folder + supp_figures[2]
    conda: "envs/KTEA.yml"
    script: "scripts/figure_s3.R"
