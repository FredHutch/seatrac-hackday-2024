import pandas as pd
from os.path import join as opj
from fg_shared import _fg_data

"""EXAMPLE CHUNK of SOFT file:
FROM GSE218157

^SAMPLE = GSM6735570
!Sample_title = 0NR IV pre
!Sample_geo_accession = GSM6735570
!Sample_status = Public on Jun 19 2023
!Sample_submission_date = Nov 16 2022
!Sample_last_update_date = Jun 19 2023
!Sample_type = SRA
!Sample_channel_count = 1
!Sample_source_name_ch1 = Whole Blood
!Sample_organism_ch1 = Macaca mulatta
!Sample_taxid_ch1 = 9544
!Sample_characteristics_ch1 = tissue: Whole Blood
!Sample_characteristics_ch1 = animal_id: 0NR
!Sample_characteristics_ch1 = pitt_id: 13417
!Sample_characteristics_ch1 = vax_group: IV
!Sample_characteristics_ch1 = vax_dose: 3.70E+07
!Sample_characteristics_ch1 = dose_group: high
!Sample_characteristics_ch1 = time after bcg: pre
!Sample_characteristics_ch1 = total_mtb_cfu: 30
!Sample_characteristics_ch1 = log_mtb_cfu: 1.491
!Sample_characteristics_ch1 = grans_nx: 1
!Sample_characteristics_ch1 = protect_outcome: protected
!Sample_molecule_ch1 = total RNA
!Sample_extract_protocol_ch1 = Whole blood was collected in PAXgene blood RNA tubes (Qiagen) and stored at -80C for batch processing at end of study. RNA was extracted using the PAXgene Blood RNA Tube kit (PreAnalytiX) as instructed. Globin mRNA was removed using GLOBINclear Kit (Life Technologies) and remaining mRNA concentration and quality was measured on an Agilent Bioanalyzer using an Agilent nano 6000 kit.
!Sample_extract_protocol_ch1 = Illumina-ready libraries were generated using NEBNext Ultra II RNA Preparation reagents (New England BioLabs).
!Sample_description = 0NR_65

FROM GSE218270

^SAMPLE = GSM6738257
!Sample_title = 13N022 high dose IV BCG pre
!Sample_geo_accession = GSM6738257
!Sample_status = Public on Jun 19 2023
!Sample_submission_date = Nov 17 2022
!Sample_last_update_date = Jun 19 2023
!Sample_type = SRA
!Sample_channel_count = 1
!Sample_source_name_ch1 = Whole Blood
!Sample_organism_ch1 = Macaca mulatta
!Sample_taxid_ch1 = 9544
!Sample_characteristics_ch1 = tissue: Whole Blood
!Sample_characteristics_ch1 = animal id: 13N022
!Sample_characteristics_ch1 = pitt id: 3218
!Sample_characteristics_ch1 = bcg dose_log10: 6.39
!Sample_characteristics_ch1 = dose group: high
!Sample_characteristics_ch1 = timeafterbcg: pre
!Sample_characteristics_ch1 = total mtb_cfu: 0
!Sample_characteristics_ch1 = log mtb_cfu: 0
!Sample_characteristics_ch1 = grans nx: 0
!Sample_characteristics_ch1 = log grans_nx: 0
!Sample_characteristics_ch1 = protect outcome: protected
!Sample_molecule_ch1 = total RNA
!Sample_extract_protocol_ch1 = Whole blood was collected in PAXgene blood RNA tubes (Qiagen) and stored at -80C for batch processing at end of study. RNA was extracted using the PAXgene Blood RNA Tube kit (PreAnalytiX) as instructed. Globin mRNA was removed using GLOBINclear Kit (Life Technologies) and remaining mRNA concentration and quality was measured on an Agilent Bioanalyzer using an Agilent nano 6000 kit.
!Sample_extract_protocol_ch1 = Illumina-ready libraries were generated using NEBNext Ultra II RNA Preparation reagents (New England BioLabs).
!Sample_description = 13N022_wkm4Pre

"""

"""FROM NCBI methods: We trimmed reads with Cutadapt and performed quality control checks with FastQC before aligning to the Macaca mulatta genome (mmul_10) using kallisto. We removed genes in the globin family, genes that were unannotated, and genes with fewer than 5 reads in over two-thirds of all samples. We applied DESeq2’s variance stabilizing transformation to generate normalized expression values."""


def parse_file(fn):
    # fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/liu_etal/GSE218157_family.soft')
    # fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/liu_etal/GSE218270_family.soft')
    s = open(fn, 'r')

    keepers = ['Sample_title',
               'Sample_geo_accession',
               'Sample_source_name_ch1',
               'Sample_organism_ch1',
               'Sample_taxid_ch1',
               'Sample_characteristics_ch1',
               'Sample_description']

    convert2float = ['bcg dose_log10', 'total mtb_cfu', 'log mtb_cfu',
                     'total_mtb_cfu', 'log_mtb_cfu']

    map_cols = {'animal id':'subjid'}

    tab = []
    started = False
    for line in s.readlines():
        if line[:7] == '^SAMPLE':
            """Find each row that starts with "^SAMPLE" and then parse the next 17 rows if they match the
            keeper texts above"""
            if started:
                tab.append(tmp)
            started = True
            tmp = {}
            continue
        if started:
            trimmed = line[1:].split('=')[0].strip()
            if trimmed in keepers:
                value = line.split('=')[1].strip()
                if trimmed == 'Sample_characteristics_ch1':
                    key = value.split(':')[0].strip()
                    value = value.split(':')[1].strip()
                    if key in convert2float:
                        value = float(value)
                else:
                    key = trimmed.replace('_ch1', '').replace('Sample_', '')

                key = map_cols.get(key, key).replace(' ', '_')
                tmp.update({key : value})
    tab.append(tmp)
    out = pd.DataFrame(tab)
    return out

def clean_data():
    out_folder = opj(_fg_data, 'SEATRAC/TB_hackday_2024/data/liu_etal')
    meta_fn = opj(out_folder, 'liu_etal_metadata.csv')
    cts_long_fn = opj(out_folder, 'liu_etal_counts_long.csv')
    cts_fn = opj(out_folder, 'liu_etal_counts.csv')

    """Load and merge meta-data into a pd DataFrame"""
    route_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2024/data/liu_etal/GSE218157_bcg_routecohort_processed.txt.gz')
    dose_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2024/data/liu_etal/GSE218270_ivbcg_dosecohort_processed.txt.gz')

    route_meta_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2024/data/liu_etal/GSE218157_family.soft.csv')
    dose_meta_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2024/data/liu_etal/GSE218270_family.soft.csv')
    more_meta_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2024/data/darrah_etal_2023/Darrah et al supplement/Project_Code/Frozen_Dataset', '220329b Correlates build.xlsx')

    route = pd.read_csv(route_fn, sep='\t')
    dose = pd.read_csv(dose_fn, sep='\t')

    mmulatta_genes_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/seatrac-hackday-2023/mmulatta_to_human.tsv')
    mmulatta = pd.read_csv(mmulatta_genes_fn, sep='\t')
    """Some macaque genes have 2 human homologs which is messing things up"""
    mmulatta = mmulatta.drop_duplicates(subset=['Gene stable ID', 'Gene name'])
    """Some gene names are NA, so use ID"""
    na_ind = mmulatta['Gene name'].isnull()
    mmulatta.loc[na_ind, 'Gene name'] = mmulatta.loc[na_ind, 'Gene stable ID']
    
    route_ids = route.columns[1:]
    route_cts = pd.merge(route.rename({'Unnamed: 0': 'geneid'}, axis=1).set_index('geneid').stack().reset_index(),
                         mmulatta.rename({'Gene stable ID':'geneid', 'Gene name':'gene'}, axis=1),
                         on='geneid', how='left')
    route_cts = route_cts.rename({'level_1':'sampleid',
                                  0:'count'}, axis=1)

    dose_ids = dose.columns[1:]
    dose_cts = pd.merge(dose.rename({'Unnamed: 0': 'geneid'}, axis=1).set_index('geneid').stack().reset_index(),
                         mmulatta.rename({'Gene stable ID':'geneid', 'Gene name':'gene'}, axis=1),
                         on='geneid', how='left')
    dose_cts = dose_cts.rename({'level_1':'sampleid',
                                  0:'count'}, axis=1)

    ct_cols = ['gene', 'sampleid', 'count']
    cts = pd.concat((route_cts[ct_cols], dose_cts[ct_cols]), axis=0)
    """required since there are multiple ENSMMUG IDs that have the same Gene name so aggregating across these"""
    cts = cts.groupby(['sampleid', 'gene'])['count'].sum().reset_index()
    cts.to_csv(cts_long_fn, index=False)
    cts = cts.set_index(['gene', 'sampleid'])['count'].unstack('sampleid')
    sample_order = cts.columns
    cts = cts.fillna(1).reset_index()
    cts.to_csv(cts_fn, index=False)

    """Process the meta-data"""
    """more is keyed on animalid and VRC ID"""
    more = pd.read_excel(more_meta_fn)
    more_cols = ['VRC ID', 'Pitt ID', 'Cohort', 'Study', 'Use for Corr?', 'Sterile?',
                   'Acquired Immunity?', 'BCG Route', 'Log IV BCG Dose', 'IV BCG Dose',
                   'IV Dose Group', 'Days','Gender', 'Age at Vaccination (yrs)']
    more = more.rename(columns={'VRC ID':'animalid',
                                'Pitt ID':'pitt_id',
                                'Gender':'sex',
                                'Age at Vaccination (yrs)':'age_yrs_at_bcg'})
    more = more.assign(studyid=more['Study'].str.lower(),
                       use_for_corr=more['Use for Corr?'].str.split(' ', expand=True).iloc[:, -1],
                       sterile=more['Sterile?'].str.split(' ', expand=True).iloc[:, -1],
                       acquired_immunity=more['Acquired Immunity?'].str.split(' ', expand=True).iloc[:, -1])
    keep_more_cols = ['animalid', 'use_for_corr', 'sterile', 'acquired_immunity', 'sex', 'age_yrs_at_bcg']
    """ONLY animalid = 0NR and 0TL were NOT in the immune response meta-data "more" and there were no additional animals in "more" 
    To make sure we have complete sex data I have sex-typed 0NR and 0TL by hand looking at expression of y-linked genes
    'DDX3Y', 'UTY', 'ZFY'. They are both female with levels comparable to other females in the cohort"""
    more = more[keep_more_col]
    more = pd.concat((more, pd.DataFrame({'animalid':['0NR','0TL'], 'sex':['F', 'F']})))

    route_meta = pd.read_csv(route_meta_fn)
    dose_meta = pd.read_csv(dose_meta_fn)
    def _parse_dose(s):
        return np.log10(float(s.split('/')[0]))
    def _try_float_as_str(v):
        try:
            return f'{v:1.0f}'
        except ValueError:
            return str(v)
    route_meta = route_meta.assign(studyid='route',
                                   route=route_meta['vax_group'].str.replace('HD-',''),
                                   bcg_dose_log10=route_meta['vax_dose'].map(_parse_dose),
                                   grans_nx=route_meta['grans_nx'].map(_try_float_as_str))
    route_meta = route_meta.assign(vax_group=route_meta.apply(lambda r: r['route'] + '-' + {'low':'LD', 'high':'HD'}[r['dose_group']], axis=1))
    dose_meta = dose_meta.assign(studyid='dose',
                                 vax_group=dose_meta['dose_group'].map(lambda s: 'IV-' + {'low':'LD', 'high':'HD'}[s]),
                                 route='IV',
                                 vax_dose='NA')

    dose_cols = ['sampleid','subjid', 'pitt_id',
                  'vax_group','bcg_dose_log10', 'dose_group',
                   'timeafterbcg', 'total_mtb_cfu', 'log_mtb_cfu',
                   'grans_nx', 'protect_outcome']

    dose_meta = dose_meta.rename(columns={'timeafterbcg':'visit',
                                            'subjid':'animalid',
                                            'description':'sampleid'})

    route_cols = ['sampleid', 'animal_id', 'pitt_id',
                    'vax_group', 'bcg_dose_log10', 'vax_dose','dose_group',
                    'time_after_bcg', 'total_mtb_cfu', 'log_mtb_cfu',
                   'grans_nx', 'protect_outcome']
    route_meta = route_meta.rename(columns={'time_after_bcg':'visit',
                                            'animal_id':'animalid',
                                            'description':'sampleid'})

    keep_cols = ['studyid','sampleid', 'animalid', 'visit', 'pitt_id',
                   'route', 'vax_group', 'bcg_dose_log10', 'vax_dose','dose_group',
                   'total_mtb_cfu', 'log_mtb_cfu',
                   'grans_nx', 'protect_outcome']
    meta = pd.concat((route_meta[keep_cols], dose_meta[keep_cols]), axis=0)
    meta = pd.merge(meta, more, on='animalid', how='left')

    """re-order to match counts"""
    meta = meta.set_index('sampleid').loc[sample_order].reset_index()
    meta.to_csv(meta_fn, index=False)
    
if __name__ == '__main__':
    from fg_shared import _fg_data
    from os.path import join as opj

    files = ['GSE218157_family.soft', 'GSE218270_family.soft']
    for fn in files:
        full_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/liu_etal', fn)
        out = parse_file(full_fn)
        out.to_csv(full_fn + '.csv')