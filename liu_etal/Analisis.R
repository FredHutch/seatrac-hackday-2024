#input data and data preprocessing 

require(dplyr)
library(limma)
library(stringr)
library(magrittr)
library(edgeR)
library(kimma)
library(readr)
library(WGCNA)
library("Hmisc")
library(ggplot2)
library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator


GSE218157_family_soft <- read_csv("GSE218157_family.soft.csv")
GSE218157_family_soft <- as.data.frame(GSE218157_family_soft)

# Remove the first column now that it has been set as row names
GSE218157_family_soft <- GSE218157_family_soft[,-1]
counts <- read.table("GSE218157_bcg_routecohort_counts_processed.txt", 
                     sep = "\t", 
                     header = TRUE, 
                     row.names = 1)
colnames(counts) <- gsub("^X", "", colnames(counts))
meta_data <-GSE218157_family_soft
rownames(meta_data) <- meta_data$description
meta_data <- meta_data %>%
  mutate(day = recode(time_after_bcg, 
                      "pre" = "0", 
                      "d2" = "2", 
                      "wk2" = "14", 
                      "wk12" = "84"))


time_list <- list(c("0", "2"),c("0", "14"),c("0", "84"))

data_normalization_edgeR <- function(filtered_count_matrix)
{
  # creating a DGE object
  normalized_matrix <- DGEList(filtered_count_matrix)
  normalized_matrix <- calcNormFactors(normalized_matrix,method = "TMM")
  return(normalized_matrix)
  
}


# function to find DEGs using Kimma model
finding_DEGs <- function(meta_data,counts,time_list,model_parameter)
{
  sample_names = meta_data%>% 
    dplyr::filter(day %in% time_list) %>% 
    pull(description)
  
  
  mts_ = meta_data[sample_names,]
  cts_ = counts[,sample_names]
  
  
  # filtering and normalization of the dataset
  keep <- rowSums(cts_>5) >= 50
  counts_filtered= cts_[keep,]
  data1 <- DGEList(counts_filtered)
  d_ss <-  data_normalization_edgeR(data1)
  
  day <- mts_[["day"]]
  
  day <- as.character(day)
  model_matrix <-  model.matrix(~day, d_ss$samples)
  design <- model_matrix
  
  
  v <- voom(d_ss, design)

  
  
  
  v$targets = v$targets %>%
    mutate(day = mts_[["day"]]) %>%
    mutate(treatment = mts_[["vax_group"]]) %>%
    mutate(libID=colnames(cts_))%>%
    mutate(ptid=mts_[["animal_id"]])
  
  #model parameter to be used
  klm     <- kmFit(dat = v,
                   model = model_parameter,
                   run_lme = TRUE,
                   use_weights = TRUE,
                   metrics = TRUE,
                   patientID = "ptid",
                   run_contrast = TRUE,
                   contrast_var = c('day', 'treatment'),
                   processors=3)
  
  
  return(klm)
  
}

model_parameter <- "~day + treatment + (1|ptid)"
for (p in time_list) {
  nam <- paste("klm","_", paste0(p[1],"_",p[2]), sep = "")
  assign(nam, finding_DEGs(meta_data,counts,p,model_parameter))
}   


######################################WGCNA############################################################


kim_0_2_gene <- klm_0_2$lme.contrast %>% 
  filter(pval < 0.05)%>% 
  filter(abs(FDR) < 0.05) %>% 
  pull(gene)


kim_0_14_gene <- klm_0_14$lme.contrast %>% 
  filter(pval < 0.05)%>% 
  filter(abs(FDR) < 0.05) %>% 
  pull(gene)


kim_0_84_gene <- klm_0_84$lme.contrast %>% 
  filter(pval < 0.05)%>% 
  filter(abs(FDR) < 0.05) %>% 
  pull(gene)

write.csv(klm_0_84$lme.contrastE,"klm_0_84.csv")
write.csv(klm_0_14$lme.contrast,"klm_0_14.csv")
write.csv(klm_0_2$lme.contrast,"klm_0_2.csv")

gene_list = unique(c(kim_0_2_gene,kim_0_14_gene,kim_0_84_gene))   

counts_subset <- counts %>% filter(rownames(counts) %in% gene_list)



datExpr<-as.data.frame(counts_subset)
datExpr0 = as.data.frame(t(datExpr))




sampleTree = hclust(dist(datExpr0), method = "average");
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub="", xlab="", cex.lab = 2,
     cex.axis = 1.5, cex.main = 2)



# Plot a line to show the cut
abline(h = 150, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)




# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.850,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



# signed correlation for expression data #
softPower = 1
corr <- rcorr(as.matrix(datExpr))$r
adj_2<- (abs(corr)/2)^7
adjacency = adjacency(datExpr, power = 1)


# Turn adjacency into topological overlap
TOM = TOMsimilarity(adj_2);
dissTOM = 1-TOM



geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);



# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")




# Rename to moduleColors
moduleColors = dynamicColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;

a=list()

color_s <- unique( moduleColors)
for (i in color_s){
  a[[i]]<-names(datExpr)[ moduleColors == i]
}
b<-c()
for( i in 1:length(a))
{
  b[i]=max(length(a[[i]]))
}
b=max(b)


for( i in 1:length(a))
{
  a[[i]]<-c(a[[i]],rep(NA, b - length(a[[i]])))
}




picked_power = 7
temp_cor <- corr       
corr <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(datExpr,                # <= input here
                          
                          
                          # Network construction arguments: correlation options
                          
                          corType = "pearson",
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)



# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )


module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)



# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(datExpr0, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )
write.csv(MEs0,"module_eigen_genes.csv")
write.csv(mME,"module_eigen_genes2.csv")
mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
