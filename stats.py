import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from patsy import dmatrices

f = "/Users/mitchjo/Analyses/DmPA_Lucas/LucasPellets_051823_09_07_2023/filtered_feature_tables/log_transformed_normalized_lucas_pellets_for_SL_Feature_table.tsv"
data = pd.read_csv(f,sep='\t')
data.head(3)

samples6 = [x.replace('cell_', 'DmPA_Lucas_cellpellets_') for x in ['cell_IV7', 'cell_IV8', 'cell_IV9',
               'cell_IV13', 'cell_IV14', 'cell_IV15',
                'cell_IV28', 'cell_IV29', 'cell_IV30', 
                'cell_IV40', 'cell_IV41', 'cell_IV42',
               ]]
data6hr = data[samples6]

df = pd.DataFrame(
    {'genotype': pd.Series(['tim3wt']*6 + ['tim3ko']*6, index=samples6),
    'treatment': pd.Series(['Normoxia']*3 + ['Hypoxia']*3 + ['Normoxia']*3 + ['Hypoxia']*3, index=samples6),}
)

def get_2anova_pvals(d, df):
    import scipy.stats as stats
    df['data'] = d
    y, X = dmatrices('data ~ genotype + treatment', data=df, return_type='dataframe')
    mod = sm.OLS(y, X)
    res = mod.fit()

    groups = {}
    for value, genotype, treatment in zip(df['data'], df['genotype'], df['treatment']):
        group = genotype + "_" + treatment
        if group not in groups:
            groups[group] = []
        groups[group].append(value)

    p_tim3_hypoxia = stats.ttest_ind(groups['tim3ko_Hypoxia'], groups['tim3wt_Hypoxia'], equal_var=False).pvalue
    p_tim3_normoxia = stats.ttest_ind(groups['tim3ko_Normoxia'], groups['tim3wt_Normoxia'], equal_var=False).pvalue
    p_o2_tim3wt = stats.ttest_ind(groups['tim3wt_Hypoxia'], groups['tim3wt_Normoxia'], equal_var=False).pvalue
    p_o2_tim3ko = stats.ttest_ind(groups['tim3ko_Hypoxia'], groups['tim3ko_Normoxia'], equal_var=False).pvalue
    p_tim3_hypoxia = 1 if not p_tim3_hypoxia else p_tim3_hypoxia
    p_tim3_normoxia = 1 if not p_tim3_normoxia else p_tim3_normoxia
    p_o2_tim3ko = 1 if not p_o2_tim3ko else p_o2_tim3ko
    p_o2_tim3wt = 1 if not p_o2_tim3wt else p_o2_tim3wt

    return res.pvalues['genotype[T.tim3wt]'], res.pvalues['treatment[T.Normoxia]'], p_tim3_hypoxia, p_tim3_normoxia, p_o2_tim3wt, p_o2_tim3ko
  

p_genotype, p_treatment = [], []
p_tim3_hypoxia = []
p_tim3_normoxia = []
p_o2_tim3wt = []
p_o2_tim3ko = []
for ii in range(data6hr.shape[0]):
    p1, p2, p3, p4, p5, p6 = get_2anova_pvals(data6hr.iloc[ii, :], df)
    p_genotype.append(p1)
    p_treatment.append(p2)
    p_tim3_hypoxia.append(p3)
    p_tim3_normoxia.append(p4)
    p_o2_tim3wt.append(p5)
    p_o2_tim3ko.append(p6)

data6hr[['id_number', 'mz', 'rtime', 'snr',]] = data[
    ['id_number', 'mz', 'rtime', 'snr',]]

data6hr['p_genotype'] = p_genotype
data6hr['p_treatment'] = p_treatment
data6hr['p_tim3_hypoxia'] = p_tim3_hypoxia
data6hr['p_tim3_normoxia'] = p_tim3_normoxia
data6hr['p_o2_tim3wt'] = p_o2_tim3wt
data6hr['p_o2_tim3ko'] = p_o2_tim3ko

data6hr.to_csv('./cells_data6hr_stats_09_08_23.tsv', index=True, sep="\t")

samples3hr = [x.replace('cell_', 'DmPA_Lucas_cellpellets_') for x in [
    'cell_IV4', 'cell_IV5', 'cell_IV6', 'cell_IV10', 'cell_IV11',
       'cell_IV12',
    'cell_IV22', 'cell_IV23', 'cell_IV24', 'cell_IV34', 'cell_IV35',
       'cell_IV36',
]]
data3hr = data[samples3hr]
df3 = pd.DataFrame(
    {'genotype': pd.Series(['tim3wt']*6 + ['tim3ko']*6, index=samples3hr),
    'treatment': pd.Series(['Normoxia']*3 + ['Hypoxia']*3 + ['Normoxia']*3 + ['Hypoxia']*3, index=samples3hr),}
)

p_genotype, p_treatment = [], []
p_tim3_hypoxia = []
p_tim3_normoxia = []
p_o2_tim3wt = []
p_o2_tim3ko = []
for ii in range(data3hr.shape[0]):
    p1, p2, p3, p4, p5, p6 = get_2anova_pvals(data3hr.iloc[ii, :], df3)
    p_genotype.append(p1)
    p_treatment.append(p2)
    p_tim3_hypoxia.append(p3)
    p_tim3_normoxia.append(p4)
    p_o2_tim3wt.append(p5)
    p_o2_tim3ko.append(p6)

data3hr[['id_number', 'mz', 'rtime', 'snr',]] = data[
    ['id_number', 'mz', 'rtime', 'snr',]]

data3hr['p_genotype'] = p_genotype
data3hr['p_treatment'] = p_treatment
data3hr['p_tim3_hypoxia'] = p_tim3_hypoxia
data3hr['p_tim3_normoxia'] = p_tim3_normoxia
data3hr['p_o2_tim3wt'] = p_o2_tim3wt
data3hr['p_o2_tim3ko'] = p_o2_tim3ko

data3hr.to_csv('./cells_data3hr_stats_09_08_23.tsv', index=True, sep="\t")

f = '/Users/mitchjo/Analyses/DmPA_Lucas/LucasMedia_051523_09_07_2023/filtered_feature_tables/log_transformed_normalized_lucas_pellets_for_SL_Feature_table.tsv'
# Assuming typo above
data = pd.read_csv(f, header=0, sep='\t')
data.head(3)

samples6 = [x.replace('cell_', 'DmPA_Lucas_cellmedia_') for x in ['cell_IV7', 'cell_IV8', 'cell_IV9',
               'cell_IV13', 'cell_IV14', 'cell_IV15',
                'cell_IV28', 'cell_IV29', 'cell_IV30', 
                'cell_IV40', 'cell_IV41', 'cell_IV42',
               ]]
data6hr = data[samples6]

df = pd.DataFrame(
    {'genotype': pd.Series(['tim3wt']*6 + ['tim3ko']*6, index=samples6),
    'treatment': pd.Series(['Normoxia']*3 + ['Hypoxia']*3 + ['Normoxia']*3 + ['Hypoxia']*3, index=samples6),}
)

p_genotype, p_treatment = [], []
p_tim3_hypoxia = []
p_tim3_normoxia = []
p_o2_tim3wt = []
p_o2_tim3ko = []
for ii in range(data6hr.shape[0]):
    p1, p2, p3, p4, p5, p6 = get_2anova_pvals(data6hr.iloc[ii, :], df)
    p_genotype.append(p1)
    p_treatment.append(p2)
    p_tim3_hypoxia.append(p3)
    p_tim3_normoxia.append(p4)
    p_o2_tim3wt.append(p5)
    p_o2_tim3ko.append(p6)

data6hr[['id_number', 'mz', 'rtime', 'snr',]] = data[
    ['id_number', 'mz', 'rtime', 'snr',]]

data6hr['p_genotype'] = p_genotype
data6hr['p_treatment'] = p_treatment
data6hr['p_tim3_hypoxia'] = p_tim3_hypoxia
data6hr['p_tim3_normoxia'] = p_tim3_normoxia
data6hr['p_o2_tim3wt'] = p_o2_tim3wt
data6hr['p_o2_tim3ko'] = p_o2_tim3ko
data6hr.to_csv('./supernatant_data6hr_stats_09_08_23.tsv', index=True, sep="\t")

samples3hr = [x.replace('cell_', 'DmPA_Lucas_cellmedia_') for x in [
    'cell_IV4', 'cell_IV5', 'cell_IV6', 'cell_IV10', 'cell_IV11',
       'cell_IV12',
    'cell_IV22', 'cell_IV23', 'cell_IV24', 'cell_IV34', 'cell_IV35',
       'cell_IV36',
]]
data3hr = data[samples3hr]
df3 = pd.DataFrame(
    {'genotype': pd.Series(['tim3wt']*6 + ['tim3ko']*6, index=samples3hr),
    'treatment': pd.Series(['Normoxia']*3 + ['Hypoxia']*3 + ['Normoxia']*3 + ['Hypoxia']*3, index=samples3hr),}
)

p_genotype, p_treatment = [], []
p_tim3_hypoxia = []
p_tim3_normoxia = []
p_o2_tim3wt = []
p_o2_tim3ko = []
for ii in range(data3hr.shape[0]):
    p1, p2, p3, p4, p5, p6 = get_2anova_pvals(data3hr.iloc[ii, :], df3)
    p_genotype.append(p1)
    p_treatment.append(p2)
    p_tim3_hypoxia.append(p3)
    p_tim3_normoxia.append(p4)
    p_o2_tim3wt.append(p5)
    p_o2_tim3ko.append(p6)

data3hr[['id_number', 'mz', 'rtime', 'snr',]] = data[
    ['id_number', 'mz', 'rtime', 'snr',]]

data3hr['p_genotype'] = p_genotype
data3hr['p_treatment'] = p_treatment
data3hr['p_tim3_hypoxia'] = p_tim3_hypoxia
data3hr['p_tim3_normoxia'] = p_tim3_normoxia
data3hr['p_o2_tim3wt'] = p_o2_tim3wt
data3hr['p_o2_tim3ko'] = p_o2_tim3ko

data3hr.to_csv('./supernatant_data3hr_stats_09_08_23.tsv', index=True, sep="\t")