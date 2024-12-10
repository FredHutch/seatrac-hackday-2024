library(tidyverse)
library(edgeR)
library(limma)
library(kimma)
library(ggrepel)
library(BIGpicture)
library(RNAetc)

#### Data ####
data_dir <- '~/seatrac-hackday-2024/data'

ncts <- readr::read_csv(file.path(data_dir, "liu_etal_counts.csv"))
meta <- readr::read_csv(file.path(data_dir, "liu_etal_metadata.csv"))

#### Data cleaning ####
# Focus initial analysis on high-dose IV BCG group and the pre-vaccine (week 2) time point
meta_ss = meta #%>% filter(vax_group == "IV-HD" & visit == "wk2")
keep_ids = meta_ss %>% pull(sampleid)
keep_ids = c('gene', keep_ids)

ncts_ss = ncts %>% dplyr::select(any_of(keep_ids))

# Discard genes that have low counts/prevalence
## Keep genes expressed at > 1 normalized count in at least 50% of samples
filter = rowSums(ncts_ss > 1) >= (0.5 * ncol(ncts_ss))
ncts_ss = ncts_ss[filter, ]

# Move gene ID to rownames
ncts_ss_mat <- as.matrix(ncts_ss[,-1])
rownames(ncts_ss_mat) = ncts_ss$gene

# Create the object for differential expression testing
dge_o = DGEList(counts=ncts_ss_mat,
                genes=rownames(ncts_ss_mat),
                samples=meta_ss,
                group=meta_ss[['protect_outcome']])

## Trimmed mean of means normalize
## We skip this step because counts were previously normalized by VST. If you were working with raw counts, this would need to be completed
# dge_o <- calcNormFactors(dge_o, method=c("TMM"))

# Specify the model/design matrix
design_temp = model.matrix(~protect_outcome + sex, data=dge_o$samples)

# Create the voom object and fit the model
v <- voomWithQualityWeights(dge_o, design=design_temp, plot=TRUE)
v_raw <- v

v$targets <- v$targets %>% 
  mutate(protect_outcome = factor(protect_outcome, 
                                  levels=c("not_protected","protected")),
         visit = factor(visit, levels=c("pre","d2","wk2","wk4","wk12")))

#### kimma ####
klm <- kmFit(
  dat = v,
  model = "~protect_outcome*visit + sex + (1|animalid)", #Include interaction
  run_lme = TRUE, 
  run_contrast = TRUE, contrast_var="protect_outcome:visit", #Run contrasts
  libraryID="sampleid",
  patientID="animalid",
  use_weights = TRUE,
  metrics = FALSE,
  processors=6)

distinct(klm$lme.contrast, contrast_lvl, contrast_ref) %>% View


c_key <- data.frame(
  # across time
  contrast_lvl = c("protected pre",    "protected d2",    "protected wk2",
                   "protected wk4",    "protected wk12"),
  contrast_ref = c("not_protected pre","not_protected d2","not_protected wk2",
                   "not_protected wk4","not_protected wk12")
) %>% 
  bind_rows(data.frame(
    # between outcomes
    contrast_lvl = c("protected d2", "protected wk2",
                     "protected wk4","protected wk12",
                     "not_protected d2", "not_protected wk2",
                     "not_protected wk4","not_protected wk12"),
    contrast_ref = c("protected pre","protected d2",
                     "protected wk2","protected wk4",
                     "not_protected pre","not_protected d2",
                     "not_protected wk2","not_protected wk4")
  ))

result <- klm$lme.contrast %>% 
  inner_join(c_key)

distinct(result, contrast_ref, contrast_lvl) %>% View()

dat <- collapse_voom(v, geneID = "genes", libraryID = "sampleid", include_weights = TRUE)

#### GAM ####
#Top interaction DEG
deg <- klm$lme %>% filter(variable=="protect_outcome:visit") %>% 
  slice_min(FDR, n=1) %>% pull(gene)

#GAM
library(mgcv)
dat_gam <- dat %>% filter(genes==deg) %>% 
  mutate(time = recode(as.character(visit), "pre"="0", "d2"="2", "wk2"="14", "wk4"="28", "wk12"="84"),
         time=as.numeric(time)) %>% 
  mutate(protect_ord = factor(protect_outcome, ordered = TRUE))

#gam fit
gam_deg <- mgcv::gam(formula = value ~ protect_ord + s(time, k = 5, bs = "cr") + s(time, k = 5, bs = "cr", by = protect_ord),
                     data = dat_gam)
#fit lines
dat_pred_np <- data.frame(protect_ord = c(rep("not_protected", 85)),
                          time=c(0:84))
gam_pred_np <- predict(gam_deg, dat_pred_np, type="link",se.fit=TRUE)
dat_pred_p <- data.frame(protect_ord = c(rep("protected", 85)),
                         time=c(0:84))
gam_pred_p <- predict(gam_deg, dat_pred_p, type="link",se.fit=TRUE)

gam_pred <- as.data.frame(gam_pred_np) %>% 
  rename(value=fit) %>% 
  mutate(protect_outcome = "not_protected") %>% 
  rownames_to_column() %>% 
  bind_rows(as.data.frame(gam_pred_p) %>% 
              rename(value=fit) %>% 
              mutate(protect_outcome = "protected") %>% 
              rownames_to_column()) %>% 
  mutate(time=as.numeric(rowname))

dat_gam %>% 
  filter(genes==deg) %>% 
  ggplot(aes(x=time, y=value, color=protect_outcome)) +
  geom_jitter(width=0.2, height=0) +
  geom_smooth(data=gam_pred, method="gam") +
  theme_classic() 

# Remove outlier
out <- dat_gam %>% filter(protect_outcome=="not_protected") %>% 
  filter(value<6) %>% pull(sampleid)
dat_gam2 <- dat_gam %>% 
  filter(!sampleid %in% out)

#gam fit
gam_deg <- mgcv::gam(
  formula = value ~ protect_ord + s(time, k = 4, bs = "cr") + s(time, k = 4, bs = "cr", by = protect_ord),
  data = dat_gam2)
summary(gam_deg)

gam_deg2 <- mgcv::gam(
  formula = value ~ protect_outcome + s(time, k = 4, bs = "cr") + s(time, k = 4, bs = "cr", by = protect_outcome),
  data = dat_gam2)
summary(gam_deg2)

#fit lines
dat_pred_np <- data.frame(protect_ord = c(rep("not_protected", 85)),
                          time=c(0:84))
gam_pred_np <- predict(gam_deg, dat_pred_np, type="link",se.fit=TRUE)
dat_pred_p <- data.frame(protect_ord = c(rep("protected", 85)),
                         time=c(0:84))
gam_pred_p <- predict(gam_deg, dat_pred_p, type="link",se.fit=TRUE)

gam_pred <- as.data.frame(gam_pred_np) %>% 
  rename(value=fit) %>% 
  mutate(protect_outcome = "not_protected") %>% 
  rownames_to_column() %>% 
  bind_rows(as.data.frame(gam_pred_p) %>% 
              rename(value=fit) %>% 
              mutate(protect_outcome = "protected") %>% 
              rownames_to_column()) %>% 
  mutate(time=as.numeric(rowname))

p1 <-dat_gam2 %>% 
  filter(genes==deg) %>% 
  ggplot(aes(x=time, y=value,)) +
  geom_jitter(aes(color=protect_outcome), width=0.2, height=0) +
  #fit
  geom_line(data=gam_pred %>% filter(protect_outcome!="all"), aes(color=protect_outcome)) +
  #error
  geom_ribbon(data=gam_pred %>% filter(protect_outcome!="all"), 
              aes(ymin = value-se.fit, ymax = value+se.fit, 
                  fill = protect_outcome), alpha=0.2, show.legend = FALSE) +
  theme_classic() +
  labs(x="Days", y=paste("Normalized log2 expression", deg, sep="\n"),
       color="Outcome") +
  annotate(
    "text", x=35, y=6.75, 
    label="Outcome P < 2e-16\ns(Time) P = 0.02\ns(Outcome:Time) P = 1.3e-04",
    hjust=0, size=3)

p1

ggsave("seatrac-hackday-2024/hacks/2024_DillMcFarland_hackday.png", 
       p1, width=5, height=4)
