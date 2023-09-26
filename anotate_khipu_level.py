import json
import numpy as np
from mass2chem.search import build_centurion_tree, find_all_matches_centurion_indexed_list
from jms.io import read_table_to_peaks
from khipu.epdsConstructor import epdsConstructor
from mass2chem.formula import PROTON, calculate_formula_mass, atom_mass_dict
from khipu.utils import isotope_search_patterns, adduct_search_patterns, extended_adducts

#adduct_search_patterns = [(PROTON, "H+"), (atom_mass_dict["K"], "K+"), (atom_mass_dict["Na"], "Na+"), (calculate_formula_mass("NH4"), "NH4+"),(calculate_formula_mass("CH3CN"), "ACN+")]

dmpa = calculate_formula_mass("C10H11O1N1")
proton = PROTON
RT_TOL = 2
carboxy_numbers = json.load(open("./stds_carboxy.json"))

for k, v in carboxy_numbers.items():
    if v is None:
        carboxy_numbers[k] = 0

import mummichog.JSON_metabolicModels as mcg
hmfn = mcg.metabolicModels['human_model_mfn']
# Copy cpds from above
Glycolysis = []
for c in ['C00068', 'C00027', 'C00149', 'C00125', 'C00126', 'C00024', 'C00288', 'C00035', 'C00036', 'C00900', 'C05378', 'C00103', 'C05125', 'C00022', 'C00579', 'C00668', 'C00092', 'C00020', 'C01342', 'C00236', 'C00010', 'C15973', 'C15972', 'C00111', 'C00033', 'C00031', 'C00354', 'C00118', 'C00074', 'C00197', 'C00085', 'C01136', 'C00665', 'C05345', 'C00227', 'C01231', 'C00084', 'C00631', 'C00221', 'C01159', 'C00248', 'C01172', 'C00267', 'C00469', 'C00661', 'C00044', 'C00186', 'C05993', 'C00169'
]:
    if c in hmfn['dict_cpds_mass']:
        # named tuples in future for better readability
        Glycolysis.append( (c, hmfn['dict_cpds_mass'].get(c, 0), hmfn['dict_cpds_def'].get(c, '')) )
    
print(len(Glycolysis), Glycolysis[:3])

TCA = []
for c in ['C00068', 'C00122', 'C00149', 'C00024', 'C00028', 'C00026', 'C00104', 'C05379', 'C01169', 'C00399', 'C00016', 'C00390', 'C00579', 'C00311', 'C00010', 'C00091', 'C00036', 'C00158', 'C15973', 'C15972', 'C00033', 'C00030', 'C00074', 'C00081', 'C00248', 'C00035', 'C00042', 'C01352', 'C00417', 'C00044', 'C05381'
]:
    if c in hmfn['dict_cpds_mass']:
        TCA.append( (c, hmfn['dict_cpds_mass'].get(c, 0), hmfn['dict_cpds_def'].get(c, '')) )
    
print(len(TCA), TCA[:3])

isotope_search_patterns = isotope_search_patterns[:6]

infile = '/Users/mitchjo/Analyses/DmPA_Lucas/LucasPellets_051823_09_07_2023/asari_results/_asari_project_97194854/preferred_Feature_table.tsv'

features = read_table_to_peaks(infile, feature_id=0, mz_col=1, rtime_col=2, intensity=(6, 7), full_extract=True)

for f in features:
    f['id'] = f['id_number']

ECON = epdsConstructor(features, mode='pos')

khipu_dict = ECON.peaks_to_epdDict(
    isotope_search_patterns = isotope_search_patterns,
    adduct_search_patterns = adduct_search_patterns,
    extended_adducts = extended_adducts,
    mz_tolerance_ppm=5, 
    rt_tolerance= RT_TOL,             # cleaner if RT window smaller
    charges=[1, 2, 3],
    has_parent_masstrack=True,
 )
epds = khipu_dict

neutrals = []
for k,v in epds.items():
    if v['neutral_formula_mass']:
        p = {}
        p['id'] = k    #  k may differ from v['interim_id']
        p['mz'] = v['neutral_formula_mass']
        p['rtime'] = np.mean([x['rtime'] for x in v['MS1_pseudo_Spectra']])
        neutrals.append(p)
    else:
        for x in v['MS1_pseudo_Spectra']:
            p = {}
            p['id'] = x['id_number']
            p['mz'] = x['mz'] - 1.0073
            p['rtime'] = x['rtime']
            neutrals.append(p)
   
print(len(neutrals))

mztree = build_centurion_tree(neutrals)

# Glycolysis

# non-derivatized
gmatched = []
for x in Glycolysis: # neutral
    match = find_all_matches_centurion_indexed_list(x[1], mztree, 5)
    if match:
        print(x, match)
        gmatched.append((x, match))
        
# single-derivatized
gmatched1 = []
for x in Glycolysis:
    cid = x[0]
    num_carboxy = carboxy_numbers[cid]
    if num_carboxy >= 1:
        mz = x[1] + dmpa 
        name = x[2] + "_dmpa*1"
        z = (cid, mz, name)
        match = find_all_matches_centurion_indexed_list(mz, mztree, 5)
        if match:
            gmatched1.append((z, match))

# double-derivatized
gmatched2 = []
for x in Glycolysis:
    cid = x[0]
    num_carboxy = carboxy_numbers[cid]
    if num_carboxy >= 2:
        mz = x[1] + 2*dmpa 
        name = x[2] + "_dmpa*1"
        z = (cid, mz, name)
        match = find_all_matches_centurion_indexed_list(mz, mztree, 5)
        if match:
            gmatched2.append((z, match))
        

print("Number of matches: ", [len(x) for x in [gmatched, gmatched1, gmatched2]])
    
# compile list of features

gwanted = []
for x in gmatched + gmatched1 + gmatched2:
    new = []
    for E in x[1]:
        if E['id'] in epds:
            for F in epds[E['id']]['MS1_pseudo_Spectra']:
                gwanted.append([x[0][0], x[0][1], x[0][2], E['id'], F['id'], F.get('ion_relation', '')])
        else:
            # gwanted.append( E )
            print("Error: ", x, E)
    
print(len(gwanted), gwanted[:3])

# TCA
# [('C00068', 425.045, 'Thiamine diphosphate'), ...]

# non-derivatized

tmatched = []
for x in TCA: # neutral
    match = find_all_matches_centurion_indexed_list(x[1], mztree, 5)
    if match:
        tmatched.append((x, match))

# single-derivatized

tmatched1 = []
for x in TCA:
    cid = x[0]
    num_carboxy = carboxy_numbers[cid]
    if num_carboxy >= 1:
        mz = x[1] + dmpa 
        name = x[2] + "_dmpa*1"
        z = (cid, mz, name)
        match = find_all_matches_centurion_indexed_list(mz, mztree, 5)
        if match:
            tmatched1.append((z, match))

# double-derivatized

tmatched2 = []
for x in TCA:
    cid = x[0]
    num_carboxy = carboxy_numbers[cid]
    if num_carboxy >= 2:
        mz = x[1] + 2* dmpa 
        name = x[2] + "_dmpa*2"
        z = (cid, mz, name)
        match = find_all_matches_centurion_indexed_list(mz, mztree, 5)
        if match:
            tmatched2.append((z, match))
            

print("Number of matches: ", [len(x) for x in [tmatched, tmatched1, tmatched2]])
    
# compile list of features

twanted = []
for x in tmatched + tmatched1 + tmatched2:
    new = []
    for E in x[1]:
        if E['id'] in epds:
            for F in epds[E['id']]['MS1_pseudo_Spectra']:
                twanted.append([x[0][0], x[0][1], x[0][2], E['id'], F['id'], F.get('ion_relation', '')])
        else:
            print("Error: ", x, E)
    
print(len(twanted), twanted[:3])

def get_dict_features(infile, sep='\t'):
    '''Returns feature dict and header.
    '''
    features = open(infile).readlines()
    header = features[0]
    idx = header.split(sep).index('id_number')
    dict_features = {}
    for f in features[1:]:
        a = f.split(sep)
        dict_features[a[idx]] = f

    print(len(dict_features))
    return dict_features, header

dict_cell3hr, header3hr = get_dict_features('/Users/mitchjo/Analyses/DmPA_Lucas/cells_data3hr_stats_09_08_23.tsv')
s = 'metabolite\tmz_expected\tname\tempCpd\tFeature_ID\tion_relation\t' + header3hr
for x in gwanted:
    x[1] = str(x[1])
    s += '\t'.join(x) + '\t' + dict_cell3hr.get(x[4], "\n")
    
with open('./Glycolysis_cells_data3hr_stats_0908_khipu_annot.tsv', 'w') as O:
    O.write(s)

dict_cell6hr, header6hr = get_dict_features('/Users/mitchjo/Analyses/DmPA_Lucas/cells_data6hr_stats_09_08_23.tsv')
s = 'metabolite\tmz_expected\tname\tempCpd\tFeature_ID\tion_relation\t' + header6hr
for x in gwanted:
    x[1] = str(x[1])
    s += '\t'.join(x) + '\t' + dict_cell6hr.get(x[4], "\n")
    
with open('./Glycolysis_cells_data6hr_stats_0908_khipu_annot.tsv', 'w') as O:
    O.write(s)

# TCA

s = 'metabolite\tmz_expected\tname\tempCpd\tFeature_ID\tion_relation\t' + header3hr
for x in twanted:
    x[1] = str(x[1])
    s += '\t'.join(x) + '\t' + dict_cell3hr.get(x[4], "\n")
    
with open('./TCA_cells_data3hr_stats_0908_khipu_annot.tsv', 'w') as O:
    O.write(s)

# TCA

s = 'metabolite\tmz_expected\tname\tempCpd\tFeature_ID\tion_relation\t' + header6hr
for x in twanted:
    x[1] = str(x[1])
    s += '\t'.join(x) + '\t' + dict_cell6hr.get(x[4], "\n")
    
with open('./TCA_cells_data6hr_stats_0908_khipu_annot.tsv', 'w') as O:
    O.write(s)

infile = '/Users/mitchjo/Analyses/DmPA_Lucas/LucasMedia_051523_09_07_2023/asari_results/_asari_project_9720159/preferred_Feature_table.tsv'

features = read_table_to_peaks(infile, feature_id=0,
                                  mz_col=1, rtime_col=2, intensity=(6, 7), full_extract=True)
for f in features:
    f['id'] = f['id_number']

ECON = epdsConstructor(features, mode='pos')

khipu_dict = ECON.peaks_to_epdDict(
    isotope_search_patterns = isotope_search_patterns,
    adduct_search_patterns = adduct_search_patterns,
    extended_adducts = extended_adducts,
    mz_tolerance_ppm=5, 
    rt_tolerance= RT_TOL,             # cleaner if RT window smaller
    charges=[1, 2, 3],
    has_parent_masstrack=True,
 )

epds = khipu_dict
neutrals = []
for k,v in epds.items():
    if v['neutral_formula_mass']:
        p = {}
        p['id'] = k    #  k may differ from v['interim_id']
        p['mz'] = v['neutral_formula_mass']
        p['rtime'] = np.mean([x['rtime'] for x in v['MS1_pseudo_Spectra']])
        neutrals.append(p)
    else:
        for x in v['MS1_pseudo_Spectra']:
            p = {}
            p['id'] = x['id_number']
            p['mz'] = x['mz'] - 1.0073
            p['rtime'] = x['rtime']
            neutrals.append(p)
   
print(len(neutrals))
mztree = build_centurion_tree(neutrals)

# Glycolysis

# non-derivatized
gmatched = []
for x in Glycolysis: # neutral
    match = find_all_matches_centurion_indexed_list(x[1], mztree, 5)
    if match:
        print(x, match)
        gmatched.append((x, match))
        
# single-derivatized
gmatched1 = []
for x in Glycolysis:
    cid = x[0]
    num_carboxy = carboxy_numbers[cid]
    if num_carboxy >= 1:
        mz = x[1] + dmpa 
        name = x[2] + "_dmpa*1"
        z = (cid, mz, name)
        match = find_all_matches_centurion_indexed_list(mz, mztree, 5)
        if match:
            gmatched1.append((z, match))

# double-derivatized
gmatched2 = []
for x in Glycolysis:
    cid = x[0]
    num_carboxy = carboxy_numbers[cid]
    if num_carboxy >= 2:
        mz = x[1] + 2*dmpa 
        name = x[2] + "_dmpa*1"
        z = (cid, mz, name)
        match = find_all_matches_centurion_indexed_list(mz, mztree, 5)
        if match:
            gmatched2.append((z, match))
        

print("Number of matches: ", [len(x) for x in [gmatched, gmatched1, gmatched2]])
    
# compile list of features

gwanted = []
for x in gmatched + gmatched1 + gmatched2:
    new = []
    for E in x[1]:
        if E['id'] in epds:
            for F in epds[E['id']]['MS1_pseudo_Spectra']:
                gwanted.append([x[0][0], x[0][1], x[0][2], E['id'], F['id'], F.get('ion_relation', '')])
        else:
            # gwanted.append( E )
            print("Error: ", x, E)
    
print(len(gwanted), gwanted[:3])
# TCA
# [('C00068', 425.045, 'Thiamine diphosphate'), ...]

# non-derivatized

tmatched = []
for x in TCA: # neutral
    match = find_all_matches_centurion_indexed_list(x[1], mztree, 5)
    if match:
        tmatched.append((x, match))

# single-derivatized

tmatched1 = []
for x in TCA:
    cid = x[0]
    num_carboxy = carboxy_numbers[cid]
    if num_carboxy >= 1:
        mz = x[1] + dmpa 
        name = x[2] + "_dmpa*1"
        z = (cid, mz, name)
        match = find_all_matches_centurion_indexed_list(mz, mztree, 5)
        if match:
            tmatched1.append((z, match))

# double-derivatized

tmatched2 = []
for x in TCA:
    cid = x[0]
    num_carboxy = carboxy_numbers[cid]
    if num_carboxy >= 2:
        mz = x[1] + 2* dmpa 
        name = x[2] + "_dmpa*2"
        z = (cid, mz, name)
        match = find_all_matches_centurion_indexed_list(mz, mztree, 5)
        if match:
            tmatched2.append((z, match))
        

print("Number of matches: ", [len(x) for x in [tmatched, tmatched1, tmatched2]])
    
# compile list of features

twanted = []
for x in tmatched + tmatched1 + tmatched2:
    new = []
    for E in x[1]:
        if E['id'] in epds:
            for F in epds[E['id']]['MS1_pseudo_Spectra']:
                twanted.append([x[0][0], x[0][1], x[0][2], E['id'], F['id'], F.get('ion_relation', '')])
        else:
            print("Error: ", x, E)
    
print(len(twanted), twanted[:3])

dict_sup3hr, header3hr = get_dict_features('/Users/mitchjo/Analyses/DmPA_Lucas/supernatant_data3hr_stats_09_08_23.tsv')
dict_sup6hr, header6hr = get_dict_features('/Users/mitchjo/Analyses/DmPA_Lucas/supernatant_data6hr_stats_09_08_23.tsv')

s = 'metabolite\tmz_expected\tname\tempCpd\tFeature_ID\tion_relation\t' + header3hr
for x in gwanted:
    x[1] = str(x[1])
    s += '\t'.join(x) + '\t' + dict_sup3hr.get(x[4], "\n")
    
with open('./Glycolysis_supernatant_data3hr_stats_0908_khipu_annot.tsv', 'w') as O:
    O.write(s)

s = 'metabolite\tmz_expected\tname\tempCpd\tFeature_ID\tion_relation\t' + header6hr
for x in gwanted:
    x[1] = str(x[1])
    s += '\t'.join(x) + '\t' + dict_sup6hr.get(x[4], "\n")
    
with open('./Glycolysis_supernatant_data6hr_stats_0908_khipu_annot.tsv', 'w') as O:
    O.write(s)

# TCA

s = 'metabolite\tmz_neutral\tname\tempCpd\tFeature_ID\tion_relation\t' + header3hr
for x in twanted:
    x[1] = str(x[1])
    s += '\t'.join(x) + '\t' + dict_sup3hr.get(x[4], "\n")
    
with open('./TCA_supernatant_data3hr_stats_0908_khipu_annot.tsv', 'w') as O:
    O.write(s)

# TCA

s = 'metabolite\tmz_neutral\tname\tempCpd\tFeature_ID\tion_relation\t' + header6hr
for x in twanted:
    x[1] = str(x[1])
    s += '\t'.join(x) + '\t' + dict_sup6hr.get(x[4], "\n")
    
with open('./TCA_supernatant_data6hr_stats_0908_khipu_annot.tsv', 'w') as O:
    O.write(s)