import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.stats import binom

from mass2chem.formula import parse_chemformula_dict

RT_TOL = 0
STD_MODE = "confirmed"

#path = "/Users/mitchjo/Analyses/DmPA_Lucas/LucasPellets_051823_09_07_2023/asari_results/_asari_project_97194854/preferred_Feature_tablefeature_annot.tsv"
#path = "/Users/mitchjo/Analyses/DmPA_Lucas/LucasPellets_051823_09_07_2023/filtered_feature_tables/log_transformed_normalized_lucas_pellets_for_SL_Feature_tablefeature_annot.tsv"
#path = "/Users/mitchjo/Analyses/DmPA_Lucas/LucasPellets_051823_09_07_2023/filtered_feature_tables/normalized_90_samples_drop_02_Feature_tablefeature_annot.tsv"
path = "/Users/mitchjo/Analyses/DmPA_Lucas/LucasMedia_051523_09_07_2023/filtered_feature_tables/normalized_90_samples_drop_02_Feature_tablefeature_annot.tsv"
cell_data = pd.read_csv(path, sep="\t")

def build_sample_map():
    sample_mapping = {}
    for s in ["IV1", "IV2", "IV3"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3+", "Normoxia", 0)
    for s in ["IV4", "IV5", "IV6"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3+", "Normoxia", 3)
    for s in ["IV7", "IV8", "IV9"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3+", "Normoxia", 6)
    for s in ["IV10", "IV11", "IV12"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3+", "Hypoxia", 3)
    for s in ["IV13", "IV14", "IV15"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3+", "Hypoxia", 6)
    for s in ["IV16", "IV17", "IV18"]: #, "IV19", "IV20", "IV21"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3-", "Normoxia", 0)
    for s in ["IV22", "IV23", "IV24"]: #, "IV25", "IV26", "IV27"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3-", "Normoxia", 3)
    for s in ["IV28", "IV29", "IV30"]: #, "IV31", "IV32", "IV33"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3-", "Normoxia", 6)
    for s in ["IV34", "IV35", "IV36"]: #, "IV37", "IV38", "IV39"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3-", "Hypoxia", 3)
    for s in ["IV40", "IV41", "IV42"]: #, "IV43", "IV44", "IV45"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3-", "Hypoxia", 6)
    
    to_delete = []
    for k in sample_mapping.keys():
        if k.endswith("M"):
            to_delete.append(k)

    for d in to_delete:
        del sample_mapping[d]

    reverse_mapping = {}
    for k,v in sample_mapping.items():
        if v not in reverse_mapping:
            reverse_mapping[v] = []
        reverse_mapping[v].append(k)
    
    return sample_mapping, reverse_mapping

def build_color_map():
    color_map = {
        ("Tim3+", None, 0): ('/', "blue"),
        ("Tim3-", None, 0): ('*', "blue"),
        ("Tim3+", "Normoxia", 3): ('/', "palegreen",),
        ("Tim3+", "Normoxia", 6): ('/', "green",),
        ("Tim3+", "Hypoxia", 3) : ('/', "lightcoral"),
        ("Tim3+", "Hypoxia", 6) : ('/', "red"),
        ("Tim3-", "Normoxia", 3): ('*', "palegreen"),
        ("Tim3-", "Normoxia", 6): ('*', "green"),
        ("Tim3-", "Hypoxia", 3) : ('*', "lightcoral"),
        ("Tim3-", "Hypoxia", 6) : ('*', "red")
    }
    return color_map

def columns_to_ids(sample_map, reverse_map, table):
    _d = {}
    _r = {}
    for column in table.columns:
        possible_id = column.split("_")[-1]
        if possible_id != column:
            _d[column] = possible_id
            _r[possible_id] = column
    return _d, _r

def build_rtime_map():
    rtime_map = {}
    MODE = "feature"
    if MODE == 'feature':
        for k, v in json.load(open("rtime_mapping_g1_feature.json")).items():
            rtime_map[k] = v
        for k, v in json.load(open("rtime_mapping_g2_feature.json")).items():
            rtime_map[k] = v
    elif MODE == 'khipu':
        for k, v in json.load(open("rtime_mapping_g1_khipu.json")).items():
            rtime_map[k] = v
        for k, v in json.load(open("rtime_mapping_g2_khipu.json")).items():
            rtime_map[k] = v
    return rtime_map

def extract(row, columns):
    return {column: row[column] for column in columns}        

def get_num_carbons(ft):
    return {n: parse_chemformula_dict(f)["C"] for (n,f) in zip(ft["name"], ft["formula"])}

def C_number_to_iso(iterable):
    return ["m+13C*" + str(num_C) if num_C != 0 else "M0" for num_C in iterable]

def plot(datas, adduct, cpd, means, id_to_sample, selection=[(None, None, None)], mode="Mean", log_transform=False, abundance=True, correct=True, drop_neg=True):
    print("plotting")
    _sm, _rm = build_sample_map()
    matching = []
    for k in _rm.keys():
        for selected in selection:
            if np.all([x == y or x is None for x,y in zip(selected, k)]):
                matching.append(k)

    if drop_neg:
        min_val = -np.inf
    else:
        min_val = 0

    new_df = []
    ISOs = []
    rv = binom(len(means), .01078)
    if correct:
        corrections = [rv.pmf(i) for i in range(len(means))]
    else:
        corrections = [0 for i in range(len(means))]
    if mode == "Mean":
        sums = {}
        for mean in means:
            for key in matching:
                if key not in sums:
                    sums[key] = 0
                if key in mean:
                    sums[key] += mean[key]
        plot = False
        values = {k: [] for k in matching}
        for i, mean in enumerate(means):
            iso = C_number_to_iso([i])[0]
            ISOs.append(iso)
            if mean:
                plot = True
                for key in matching:
                    #print(mean[key])
                    if log_transform:
                        values[key].append(np.log2(mean[key]))
                    elif abundance:
                        values[key].append(max(min_val, (mean[key] / sums[key]) - corrections[i]))
                    else:
                        if i > 0 and correct:
                            #print(mean[key])
                            values[key].append(mean[key] - corrections[i] * mean[key])
                        else:
                            values[key].append(mean[key])
            else:
                for key in matching:
                    values[key].append(0)

    else:
        plot = False
        values = {}
        for k in matching:
            for id in _rm[k]:
                sample = id_to_sample[id]
                values[sample] = []

        for i , data in enumerate(datas):
            iso = C_number_to_iso([i])[0]
            ISOs.append(iso)
            if data:
                plot = True
                for k in matching:
                    for id in _rm[k]:
                        sample = id_to_sample[id]
                        if log_transform:
                            values[sample].append(np.log2(data[sample][0]))
                        else:
                            values[sample].append(data[sample][0])
            
            else:
                for k in matching:
                    for id in _rm[k]:
                        sample = id_to_sample[id]
                        values[sample].append(0)
    if plot:
        width = 0.05
        multiplier = 0
        x = np.arange(len(ISOs))
        fig, ax = plt.subplots(layout='constrained')
        for attribute, measurement in values.items():
            print(measurement)
            offset = width * multiplier
            rects = ax.bar(x + offset, measurement, width, label='_'.join([str(x) for x in attribute]))
            ax.bar_label(rects, padding=3)
            multiplier += 1
        ax.set_xticks(x+width, ISOs)
        #ax.legend(loc='lower right', ncols=3)
        ax.set_title(cpd + " " + adduct)
        manager = plt.get_current_fig_manager()
        #manager.full_screen_toggle()
        #plt.show()
        plt.savefig("./figures/supernatant_" + cpd + "_adduct_no_legend.png")

def compare_rts(df, rtime_map, mode="confirmed_retention_times"):
    matches_rt = []
    for name, rtime in zip(df["name"], df['rtime']):
        name = name.lower()
        if name in rtime_map:
            matched_rt = False
            for interval in rtime_map[name][mode]:
                if interval[0] - RT_TOL < rtime < interval[1] + RT_TOL:
                    matched_rt = True
                    break
            matches_rt.append(matched_rt)
        else:
            if STD_MODE == "confirmed":
                matches_rt.append(False)
            elif STD_MODE == "permissive":
                matches_rt.append(True)
    df["rt_match"] = matches_rt
    df = df[df["rt_match"].isin([None, True])]
    return df

rtime_map = build_rtime_map()

cell_data = compare_rts(cell_data, rtime_map)
#cell_data = cell_data[cell_data["rt_match"] == True] 


num_carbons = get_num_carbons(cell_data)
sample_mapping, reverse_mapping = build_sample_map()
sample_to_id, id_to_sample = columns_to_ids(sample_mapping, reverse_mapping, cell_data)

all_cpds = set([x for x in cell_data["name"]])
for cpd in all_cpds:
    selected = cell_data[cell_data["name"] == cpd]
    #print(cpd, selected.shape)
    for adduct in set(list(cell_data["ion_relation"])):
        a_selected = selected[selected["ion_relation"] == adduct]
        #print("\t", adduct, a_selected.shape)
        datas = []
        means = []
        for iso in C_number_to_iso(range(num_carbons[cpd])):
            data = {}
            mean = {}
            i_a_selected = a_selected[a_selected["isotopologue"] == iso]
            #print("\t\t", iso, i_a_selected.shape)
            for entry in i_a_selected.apply(extract, axis=1, args=(i_a_selected.columns,)):
                for k, v in reverse_mapping.items():
                    if k not in mean:
                        mean[k] = []
                    for s_i in v:
                        sample = id_to_sample[s_i]
                        if sample not in data:
                            data[sample] = []
                        data[sample].append(entry[sample])
                        mean[k].append(entry[sample])
                    mean[k] = np.mean(mean[k], axis=0)
                break
            means.append(mean)
            datas.append(data)
        plot(datas, adduct, cpd, means, id_to_sample)




