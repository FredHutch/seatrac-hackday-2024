library(readr)
# Load count data and meta-data
# Select the data for modeling (groups and visits)
# Create DGE object and set up the model
# Run the model and volcano plot
# Does model need adjustment for covariates? e.g., sex, route, dose


# Extensions
# Use WGCNA or GSEA to identify or form gene modules from DEGs?
# How well do the genes classify protected and non-protected animals using CV-AUC?
# Are the genes also induced by BCG? At which visits?
# How do they overlap with modules in Liu et al.?

# Additional stubs
# Correlate baseline M1, M2 expression with immune responses
# Look at correlations between immune response variables and protection
# Foreman/Bromley: Identify significant CD8+ genes from protected granulomas and do the same from psuedo bulked CD8+ from Bromley control or CD4-depleted animals

data_folder <- "data/darrah_etal/"
#raw = readr::read_tsv("GSE218270_ivbcg_dosecohort_processed.txt.gz")
#rmeta = readr::read_csv("GSE218270_family.soft.csv")

raw <- readr::read_csv(pastePaths(data_folder, "liu_etal_counts.csv"))
meta <- readr::read_csv(pastePaths(data_folder, "liu_etal_metadata.csv"))
rownames(meta) = meta$sampleid

# unneccessary, since already in order
# rawordered = raw[,meta$sampleid]
# meta2 = meta[colnames(rawordered),]

dge_o = DGEList(counts=raw, genes=raw[, 1], samples=meta, group=meta[['visit']])

ncts <- calcNormFactors(raw, method = "TMM")

group <- meta$visit
design <- model.matrix(~0+group)
colnames(design) <- gsub("group","",colnames(design))
contr.matrix <- makeContrasts(wk12vPre = pre - wk12, 
                              levels = colnames(design))

voom <- voom(dge_o, design=design, plot=TRUE)

library(biomaRt)

mart <- useDataset("mmulatta_gene_ensembl", useMart("ensembl"))
genes <- voom$genes
mart <- useDataset("mmulatta_gene_ensembl", useMart("ensembl"))
gene_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)


