setwd("C:\\Users\\admin\\OneDrive - Grinnell College\\Documents\\College\\Coding\\Gene Expression Analysis Pipeline")
library(affy)
library(GEOquery)
library(tidyverse)
library(dplyr)

#get supplementary files
getGEOSuppFiles("GSE48558")

#untar files
untar("GSE48558/GSE48558_RAW.tar", exdir = 'Data/')
#list.files("Data/")

#reading .cel files
raw.data <- ReadAffy(celfile.path = "Data/")

#performing rma RNA normalization
normalized.data <- rma(raw.data)

#get expression estimates 
normalized.expr <- as.data.frame(exprs(normalized.data))

#map probe IDs to gene symbols
gse <- getGEO("GSE48558", GSEMatrix = TRUE)

#fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE48558_series_matrix.txt.gz@featureData@data
#colnames(feature.data)
#subset
#get 10th column
feature.data <- feature.data[,c(1,10)]

#delete omitted data
feature.data <- filter(feature.data, gene_assignment != "---")

#get the gene_symbol
feature.data <- feature.data%>%
  separate(
    col = gene_assignment,
    into = c("genbank", "gene_symbol", "description"),
    sep = " // ",
    remove = FALSE,
    fill = "right",
    extra = "merge"
  )
#clean up the feature.data
feature.data <- feature.data%>%
  select(ID, gene_symbol)
#change ID to character
feature.data <- feature.data%>%
  mutate(ID = as.character(ID))

#add gene_symbol to normalized.expr
normalized.expr <- normalized.expr %>%
  rownames_to_column(var = 'ID') %>%
  left_join(
    feature.data, 
    by = 'ID'
  )

#reorder gene_symbol
normalized.expr <- normalized.expr %>%
  select(ID, gene_symbol, everything())

#delete rows without gene_symbol
normalized.expr <- normalized.expr %>%
  drop_na(gene_symbol)

#remove duplicated gene_symbol
normalized.expr <- normalized.expr[!duplicated(normalized.expr$gene_symbol), ]

#get phenotype
pheno.data <- gse$GSE48558_series_matrix.txt.gz@phenoData@data
metadata <- pheno.data%>%
  select(geo_accession,title, source_name_ch1, characteristics_ch1) %>%
  mutate(
    condition = case_when(
      grepl("AML|Leukemia|Patient", characteristics_ch1, ignore.case = TRUE) ~ "leukemia",
      grepl("Normal|control", characteristics_ch1, ignore.case = TRUE) ~ "normal",
      TRUE ~ "unknown"
    )
  )

#the number of normal sample vs patients
table(metadata$condition)

#rename
metadata <- metadata %>%
  select(sample_id = geo_accession, 
         sample_name = title, 
         condition)

# Save metadata
dir.create("Processed_Data")
#metadata
write.csv(metadata,
          "Processed_Data/metadata.csv",
          row.names = FALSE)
#normalized.expr
write.csv(normalized.expr,
          "Processed_Data/normalized_expression.csv",
          row.names = TRUE)