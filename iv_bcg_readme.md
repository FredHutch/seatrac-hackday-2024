# README for datasets relating to IV BCG in NHP

## Manuscripts

Hacks related to IV BCG in NHP will build on these studies and the published datasets:

## Cleaned datasets and data dictionaries

### Bulk transcriptomics

Bulk RNA sequencing from whole blood was published in Liu et al. Data has been aggregated at the gene level, with files specified below. Values are log2-counts.

[meta-data](link)
 - filepath: `[DATA PATH]/liu_etal_metadata.csv`
 - columns: `sampleid` [OTHERS]
 - rows: one row per animal
 - combines both the Darrah et al. "route" study and "IV dose" studies (2020 and 2023)

[counts table (long-form)](link)
 - filepath: `[DATA PATH]/liu_etal_counts_long.csv`
 - columns: `sampleid`, `gene` symbol (macca mulatta), and log2 `count`
 - rows: one row per sample, per gene
 - can be joined on sampleid with the meta-data CSV
 - combines "route" and "dose" Darrah et al. studies

[counts table (wide-form)](link)
 - filepath: `[DATA PATH]/liu_etal_counts.csv`
 - columns: `gene` and `sampleid`s (one column per sampleid)
 - rows: one row per gene
 - this format is suitable for reading straight into differential expression analysis, with `sampleid` order matching that of the `sampleid` order in the meta-data CSV
  - combines "route" and "dose" Darrah et al. studies

### BCG-induced immune responses

Both Darrah et al. studies meausured cellular and humoral immune responses before and after BCG vaccination from lung (BAL) and blood.

[immune_responses](link)
 - filepath: `[DATA PATH]/darrah_etal_immune_resp.csv`
 - columns: `sampleid`, `tissue` (blood or BAL), `bioid` (unique identifier of the immune response), [OTHER COLUMS DESCRIBING IMMUNE RESPONSE], `Value`, `Unit`
 - rows: one row per animal
 - combines both the Darrah et al. "route" study and "IV dose" studies (2020 and 2023)
 - importantly, the `sampleid` match those from the Liu et al. meta-data file above to allow for joining these datasets (e.g., `0NR_65` from "route" study and `OBC_wk12Pre` from "dose" study)

[KMB: here I'm thinking that the immune response metadata describing which assay etc will be pre-joined withe the assay values themselves. Let me know if you think it should be kept separate]


## TODO:
 - Clean and prep each file
 - Write startup scripts
 - Organize repo
 - Test on Rstudio servers