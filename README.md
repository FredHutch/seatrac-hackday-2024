# SEATRAC TB Hackday 2024
Tuesday, December 10, 2024

9 AM to 3 PM EDT

9 AM to 3 PM PDT

Please provide feedback in this [survey](https://qfreeaccountssjc1.az1.qualtrics.com/jfe/form/SV_1zRXN9eOJjHHdgG)

## Seattle:
Fred Hutchinson Cancer Center, Steam Plant Building, O’Mack Symposium Suite, 1241 Eastlake Ave E, Seattle, WA 98102

[FH Campus Map](https://www.fredhutch.org/en/about/contact-us/campus-map.html)

## Emory University:
Emory Rollins School of Public Health, Claudia Nance Rollins Building, Room 3001, 1518 Clifton Road NE, Atlanta, GA 30322

## Virtual:

Zoom ([Zoom link](https://us02web.zoom.us/j/82494563007?from=addon)). 

Github link: [FredHutch/seatrac-hackday-2024](https://github.com/FredHutch/seatrac-hackday-2024)

## Agenda

The event will be starting in-person at Emory at 9 AM (EDT). Folks in Seattle will start in-person at 9 AM (PDT).
There will be a JOINT Seattle/Emory session at 12 PM PDT / 3 PM EDT to share hackday progress.

| Emory (EST) | Seattle (PST) | Activity |
|-------------|---------------|----------|
| 9:00 AM     |               | Alex Shalek and Johsua Bromley presenting on correlates of protection for IV BCG and Mtb immunity in NHP (live) |
| 9:45        |               | Opening remarks (Shuyi Ma, Andrew Fiore-Gartland) |
| 10:00       |               | Small group session: devise questions, hypotheses and analysis plans |
| 10:30       |               | Open hack time at Emory |
|             | 8:30 AM       | Gather with coffee and light breakfast |
|             | 9:00          | Opening remarks (Shuyi Ma, Andrew Fiore-Gartland) |
|             | 9:15          | Alex Shalek on correlates of protection for IV BCG and Mtb immunity in NHP (re-broadcast) |
|             | 10:00         | Small group session: devise questions, hypotheses and analysis plans |
|             | 10:30         | Open hack time in Seattle |
| 3:00 PM     | 12:00 PM      | **JOINT**<br>Small group readouts from Emory: Each group shares preliminary figures, analysis plans, code chunks |
|             | 3:00          | Small group readouts from Seattle |

## Objectives

Learning, teaching and collaborating! Collaborating in small groups to develop new analytical insights from published bulk and single-cell RNAseq datasets that profile NHP host gene expression during Mtb infection and BCG vaccination.

Datasets: PMIDs: 37097292, 35483355, 37267955, 37390827, 39214090

Breakfast included starting at 8:30am, lunch included at noon, happy hour snacks included starting at 3pm.

Registration: 
[Online registration](https://bit.ly/3ZdQbUq)

## Jumpstart scripts

Use these scripts/tutorials to get you going faster with the datasets:
1. [Identifying whole-blood gene expression after IV BCG in NHP that is associated with protection from Mtb challenge](https://github.com/FredHutch/seatrac-hackday-2024/blob/main/jumpstart_scripts/iv_bcg_dge.qmd)
2. [Can the adaptive responses to IV BCG be predicted from baseline or early-BCG gene expression in whole-blood?](https://github.com/FredHutch/seatrac-hackday-2024/blob/main/jumpstart_scripts/iv_bcg_gex_module_imm_resp.qmd)
3. [Exploring immune responses to i.v. BCG and protection from subsequent Mtb Challenge](https://github.com/FredHutch/seatrac-hackday-2024/blob/main/jumpstart_scripts/iv_bcg_immune_correlates.qmd)
4. [Identify granuloma associated T cell genes and validate in single-cell dataset](https://github.com/FredHutch/seatrac-hackday-2024/blob/main/jumpstart_scripts/granuloma_tcell_dge.qmd)

## Tutorials from TB Lunch & Learn sessions

1. [Data cleaning with R tidyverse](https://bigslu.github.io/workshops/2022.08.15_R.tidyverse.workshop/2_tidyverse.html)
2. [Bulk RNAseq differential expression](https://bigslu.github.io/workshops/2024_SEATRAC_series/3_linear_model_rnaseq.html)
3. [Data visualization with ggplot](https://bigslu.github.io/workshops/2022.08.15_R.tidyverse.workshop/3_ggplot.html)

# Datasets

Link to aggregated datasets on [Figshare](https://figshare.com/articles/dataset/SEATRAC_TB_Hackday_2024/27774420)

Foreman et al., 2023: bulk RNAseq from sorted cells in NHP Mtb challenge model 

Darrah et al., 2023: scRNAseq from BAL in NHP Mtb challenge (BCG route, correlates of protection) 

Liu et al., 2023: whole-blood bulk RNAseq from NHP Mtb challenge (BCG route and IV BCG dose, should match Darrah et al. study) 

Bromley et al., 2024: scRNAseq from granulomas during primary and secondary Mtb infection
 
---

## Foreman et al., 2023: bulk RNAseq from sorted cells in NHP Mtb challenge model 

Foreman TW, Nelson CE, Sallin MA, Kauffman KD, Sakai S, Otaizo-Carrasquero F, Myers TG, Barber DL. CD30 co-stimulation drives differentiation of protective T cells during Mycobacterium tuberculosis infection. J Exp Med. 2023 Aug 7;220(8):e20222090. doi: 10.1084/jem.20222090. Epub 2023 Apr 25. PMID: 37097292; PMCID: PMC10130742. 

PDF: https://rupress.org/jem/article-pdf/220/8/e20222090/1451373/jem_20222090.pdf; ([github](https://github.com/FredHutch/seatrac-hackday-2023/blob/main/foreman_etal/Foreman%20et%20al%20CD30%20drives%20differentiation%2C%202023.pdf)) 

MS: https://rupress.org/jem/article/220/8/e20222090/214054/CD30-co-stimulation-drives-differentiation-of 

DATA:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227653 (M mulatta data)

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228114 for the superset that contains both the NHP and mouse data for that study 

https://github.com/FredHutch/seatrac-hackday-2023/tree/main/foreman_etal

DATA SUMMARY: N=4 Rhesus macaques, Mtb. challenge (40–80 CFU of Mtb-Erdman-mCherry and euthanized 6–7 wk after infection); bulk RNAseq from FACS sorted T cells; [bulk RNAseq count data](https://github.com/FredHutch/seatrac-hackday-2023/blob/main/foreman_etal/GSE227653_TPM_all.csv.gz) is stored for CD4 and CD8 T cells in units of transcript per million (TPM).

---
## Darrah et al., 2023: scRNAseq from BAL in NHP Mtb challenge (BCG route, correlates of protection) 

Darrah PA, Zeppa JJ, Wang C, Irvine EB, Bucsan AN, Rodgers MA, Pokkali S, Hackney JA, Kamath M, White AG, Borish HJ, Frye LJ, Tomko J, Kracinovsky K, Lin PL, Klein E, Scanga CA, Alter G, Fortune SM, Lauffenburger DA, Flynn JL, Seder RA, Maiello P, Roederer M. Airway T cells are a correlate of i.v. Bacille Calmette-Guerin-mediated protection against tuberculosis in rhesus macaques. Cell Host Microbe. 2023 Jun 14;31(6):962-977.e8. doi: 10.1016/j.chom.2023.05.006. Epub 2023 Jun 1. PMID: 37267955; PMCID: PMC10355173. 

PDF: https://www.nature.com/articles/s41586-019-1817-8.pdf; ([github](https://github.com/FredHutch/seatrac-hackday-2023/blob/main/darrah_etal/Darrah%20et%20al.%20IV%20BCG%20correlates%20Cell%20Host%20and%20Microbe%202023.pdf))


MS: https://www.nature.com/articles/s41586-019-1817-8 

DATA: https://singlecell.broadinstitute.org/single_cell/study/SCP796/prevention-of-mycobacterium-tuberculosis-infection-and-disease-in-nonhuman-primates-following-intravenous-bcg-vaccination?scpbr=the-alexandria-project 
Single-cell data is on [Figshare](https://figshare.com/account/articles/24425053)


DATA SUMMARY: single-cell; n=15 Rhesus macaques with BCG vaccination and Mtb challenge; 3 individuals per group for AE, IDhigh, IDlow, IV, Naïve-controls; BAL collected at Weeks 13 and 25 prior to challenge at Week 26; stimulated and unstimulated conditions; 60 samples, 1000 – 5000 cells each. 

Mycobacterium tuberculosis (Mtb) is the leading cause of death from infection worldwide. Intradermal (ID) vaccination with BCG has variable efficacy against pulmonary tuberculosis, the major cause of mortality and disease transmission. Here we show that the route and dose of BCG vaccination alters circulating and lung resident T cells and subsequent protection against Mtb challenge in nonhuman primates (NHP). NHP immunized with BCG by the intravenous (IV) route induced substantially higher antigen-specific CD4 (Th1 or Th17) and CD8 responses in blood, spleen, bronchoalveolar lavage (BAL), and lung lymph nodes compared to the same BCG dose administered by ID or aerosol (AE) routes. Moreover, IV immunization was the only route that induced a high frequency of antigen-specific tissue resident T cells in lung parenchyma. Six months after BCG vaccination, NHP were challenged with virulent Mtb. Strikingly, 9 of 10 NHP that received BCG IV were highly protected, with 6 NHP showing no detectable infection as determined by PET CT imaging, mycobacterial growth, pathology, granuloma formation, or de novo immune responses to Mtb-specific antigens. The finding that BCG IV prevents or significantly limits Mtb infection in NHP has important implications for vaccine development and provides a model for determining immune correlates and mechanisms of protection against TB. 

---
## Liu et al., 2023: whole-blood bulk RNAseq from NHP Mtb challenge (BCG route and IV BCG dose) 

Liu YE, Darrah PA, Zeppa JJ, Kamath M, Laboune F, Douek DC, Maiello P, Roederer M, Flynn JL, Seder RA, Khatri P. Blood transcriptional correlates of BCG-induced protection against tuberculosis in rhesus macaques. Cell Rep Med. 2023 Jul 18;4(7):101096. doi: 10.1016/j.xcrm.2023.101096. Epub 2023 Jun 29. PMID: 37390827; PMCID: PMC10394165. 

PDF: [github](https://github.com/FredHutch/seatrac-hackday-2023/blob/main/liu_etal/Liu%20et%20al%20IV%20BCG%20NHP%20mRNA%202023.pdf)
MS: https://www.sciencedirect.com/science/article/pii/S266637912300215X?via%3Dihub#sec4.1 

DATA IV BCG DOSE STUDY (167 samples, 34 NHP): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218270 

DATA BCG ROUTE STUDY (144 samples, 36 NHP):  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218157 

DATA SUMMARY: Darah et al.  IV BCG dose and BCG route study matching; Pre BCG and Day 2, Wk 4, Wk 12 for IV BCG dose study;  Pre, Day 2, Wk2 and Wk12 for route study; bulk RNAseq from whole blood; includes protection data for each individual NHP 

---
## Bromley et al., 2024: scRNAseq from granulomas during primary and secondary Mtb infection
Bromley JD, Ganchua SKC, Nyquist SK, Maiello P, Chao M, Borish HJ, Rodgers M, Tomko J, Kracinovsky K, Mugahid D, Nguyen S, Wang QD, Rosenberg JM, Klein EC, Gideon HP, Floyd-O'Sullivan R, Berger B, Scanga CA, Lin PL, Fortune SM, Shalek AK, Flynn JL. CD4+ T cells re-wire granuloma cellularity and regulatory networks to promote immunomodulation following Mtb reinfection. Immunity. 2024 Oct 8;57(10):2380-2398.e6. doi: 10.1016/j.immuni.2024.08.002. Epub 2024 Aug 29. PMID: 39214090; PMCID: PMC11466276.

PDF: [github](https://github.com/FredHutch/seatrac-hackday-2024/blob/main/bromley_etal/Bromley%20et%20al%20CD4%20T%20cells%20re-wire%20granuloma%20cellularity%2C%20Immunity%2C%202024.pdf)
MS: [PubMed](https://pubmed.ncbi.nlm.nih.gov/39214090/) 

Bromley et al. report on an experiment with three groups of cynomolgus macaques: (1) anti-CD4 treated, Mtb-exposed (n=7), (2) IgG control, Mtb-exposed (n=6), (3) No treatment, Mtb-naive (n=6). Groups (2) and (3) were infected with Mtb, then given an anti-CD4 antibody to deplete CD4+ T cells or an isotype control (IgG) antibody, and finally challenged with a secondary Mtb infection. Group 3 only received a primary Mtb infection. Granulomas were then analyzed using single-cell RNAseq, with 3 NHP from (1) and 2 NHP from (2) and (3) each. We have created pseudo-bulk datasets using two different clustering of the cells, offering two levels of cell type granularity for analysis: `bromley_X_pseudobulk_counts.csv` where `X` is `coarse` or `subclustering`, with clusters defined by the authors of the manuscript. By analyzing data from primary infection, reinfection, and reinfection-CD4+ T cell-depleted granulomas, they found that the presence of CD4+ T cells during reinfection resulted in a less inflammatory lung milieu characterized by reprogrammed CD8+ T cells, reduced neutrophilia, and blunted type 1 immune signaling among myeloid cells.
