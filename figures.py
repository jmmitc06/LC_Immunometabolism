import pandas as pd
import matplotlib.pyplot as plt
import os 
import matplotlib.patches as mpatches
import json
import sys
import seaborn as sns

ALT_COLORS = True
RT_TOL = 0
MODE = sys.argv[1]
STD_MODE = sys.argv[2]
SIG_ONLY = sys.argv[3]
print(sys.argv)
if SIG_ONLY == "1":
    SIG_ONLY = True
else:
    SIG_ONLY = False

RESULT_DIR = "./results_09_08_23"
if MODE == "feature":
    RESULT_DIR += "_feature_level"
elif MODE == "khipu":
    RESULT_DIR += "_khipu_level"

if STD_MODE == "confirmed":
    RESULT_DIR += "_confirmed"
elif STD_MODE == "permissive":
    RESULT_DIR += "_pemissive"

if SIG_ONLY:
    RESULT_DIR += "_sig_only"

RESULT_DIR += "_" + str(RT_TOL)
RESULT_DIR += "/"

USED_CPDS=set()
TRANSLATION = [
    ['Lactic acid', '(S)-Lactate'],
    ['Lactic acid', 'l-lactate'],
    ['Malate', '(S)-Malate'],
    ['Malate', 'l-malate'],
    ['D-Glucose', 'alpha-D-Glucose'],
    ['D-Glucose', 'beta-D-Glucose'],
    ['Citric acid', 'Citrate'],
    ['DL-alpha-Lipoamide', 'Lipoamide'],
    ['Fructose 6-phosphate (F6P)', 'd-fructose 6-phosphate'],
    ['Fructose 6-phosphate (F6P)', 'beta-d-fructose 6-phosphate'],
    ['Fructose 6-phosphate (F6P)', 'alpha-d-fructose 6-phosphate'],
    ['Glucose 6-phosphate (G6P)', 'd-glucose 6-phosphate'],
    ['Glucose 6-phosphate (G6P)', 'beta-d-glucose 6-phosphate'],
    ['Glucose 6-phosphate (G6P)', 'alpha-d-glucose 6-phosphate'],
    ['Coenzyme Q10', 'ubiquinone'],
    ['Thiamine pyrophosphate hydrochloride', 'thiamin diphosphate'],
    ['3-Phosphoglycerate (3PG)', 'dl-glyceraldehyde 3-phosphate'],
    ['Isocitric acid', 'isocitrate'],
    ['beta-nicotinamide adenine dinucleotide hydrate', 'amp'],
    ['beta-nicotinamide adenine dinucleotide hydrate', 'amp'],
    ['Glucose 1-phosphate (G1P)', 'd-glucose 1-phosphate'],
    ['Glucose 1-phosphate (G1P)', 'beta-d-glucose 1-phosphate'],
    ['3-Phosphoglycerate (3PG)', '3-phospho-d-glycerate'],
    ['3-Phosphoglycerate (3PG)', 'd-glyceraldehyde 3-phosphate'],
    ['2-Phosphoglycerate (2PG)', '2-phospho-d-glycerate'],
    ["Coenzyme Q10", "ubiquinol"],
    ["Oxaloacetic acid", "oxaloacetate"],
    ["Inosine 5-triphosphate trisodium salt", 'itp'],
    ["Fructose 1,6-bisphosphate (FBP)", "beta-d-fructose 1,6-bisphosphate"],
    ["Fructose 1,6-bisphosphate (FBP)", "d-fructose 1,6-bisphosphate"],
    ['phosphoenolpyruvate (PEP)', 'phosphoenolpyruvate'],  
]

def build_rtime_map():
    rtime_map = {}
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
    for s in ["IV16", "IV17", "IV18", "IV19", "IV20", "IV21"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3-", "Normoxia", 0)
    for s in ["IV22", "IV23", "IV24", "IV25", "IV26", "IV27"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3-", "Normoxia", 3)
    for s in ["IV28", "IV29", "IV30", "IV31", "IV32", "IV33"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3-", "Normoxia", 6)
    for s in ["IV34", "IV35", "IV36", "IV37", "IV38", "IV39"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3-", "Hypoxia", 3)
    for s in ["IV40", "IV41", "IV42", "IV43", "IV44", "IV45"]:
        sample_mapping[s] = sample_mapping[s+"M"] = ("Tim3-", "Hypoxia", 6)
    return sample_mapping

def build_color_map():
    color_map = {
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

def build_patches(color_map):
    if not ALT_COLORS:
        patches3hr = [mpatches.Patch(color=v[1], hatch=v[0], alpha=0.25, label='_'.join([str(x) for x in k]) + 'hrs') for k,v in color_map.items() if k[-1] == 3]
        patches6hr = [mpatches.Patch(color=v[1], hatch=v[0], alpha=0.25, label='_'.join([str(x) for x in k]) + 'hrs') for k,v in color_map.items() if k[-1] == 6]
    else:
        patches3hr = [mpatches.Patch(color=v[1], fill=True if v[0] == '/' else False ,label='_'.join([str(x) for x in k]) + 'hrs') for k,v in color_map.items() if k[-1] == 3]
        patches6hr = [mpatches.Patch(color=v[1], fill=True if v[0] == '/' else False , label='_'.join([str(x) for x in k]) + 'hrs') for k,v in color_map.items() if k[-1] == 6]
    return patches3hr, patches6hr

def check_all_interpolated(row, samples):
    values = set()
    for sample in samples:
        values.add(row[sample])
    if len(values) == 1:
        return True
    else:
        return False

def make_results_dir(result_dir):
    try:
        os.makedirs(result_dir)
    except:
        import shutil
        shutil.rmtree(result_dir)
        os.makedirs(result_dir)

def find_annotation_files(search_dir, extra=''):
    annotation_dfs = []
    base_names = []
    sup_titles = []
    hours = []
    for path in os.listdir(os.path.abspath(search_dir)):
        sample_type = None
        hour = None
        annot_type = None
        if '3hr' in path:
            hour = '3hr'
        elif '6hr' in path:
            hour = '6hr'

        if 'cell' in path:
            sample_type = 'cell'
        elif 'supernatant' in path:
            sample_type = 'supernatant'
        
        if MODE == 'khipu':
            if 'khipu' in path:
                annot_type = 'khipu'
        elif MODE == 'feature':
            if 'feature' in path:
                annot_type = 'feature'
        
        if hour and sample_type and annot_type and path.endswith(".tsv") and extra in path:
            basename = '_'.join([sample_type, hour]) + '_'
            annotation_df = pd.read_csv(os.path.join(os.path.abspath("."), path), sep="\t", index_col=0)
            sup_titles.append(sample_type + "_" + hour)
            base_names.append(basename)
            annotation_dfs.append(annotation_df)
            hours.append(hour)
    return sup_titles, base_names, annotation_dfs, hours

def filter_df(df):
    if "isotopologue" in df.columns:
        df["ion_relation"] = df["isotopologue"] + "," + df["ion_relation"]
    df = df[df['rtime'].notna()].copy()
    df["drop"] = df.apply(check_all_interpolated, axis=1, args=(get_sample_names(df),))
    df = df[df['drop'] == False].copy()
    df["name"] = [x.lower() for x in df["name"]]
    if SIG_ONLY:
        df = df[(df["p_genotype"] < .05) | (df["p_treatment"] < .05) | (df["p_tim3_hypoxia"] < .05) | (df['p_tim3_normoxia'] < .05) | (df['p_o2_tim3wt'] < .05) | (df['p_o2_tim3ko'] < .05)]
    return df

def get_sample_names(df):
    return [column for column in df.columns if column.startswith("DmPA")]

def compare_rts(df, rtime_map, mode="confirmed_retention_times"):
    matches_rt = []
    for name, rtime in zip(df["name"], df['rtime']):
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

def extract_cpds(df):
    return set(list(df["name"]))

def plot_cpd(cpd_name, df, title, hour, sort="genotype"):
    sup_title = title
    for translation in TRANSLATION:
        if translation[1].lower() in cpd_name.lower():
            cpd_name.replace(translation[1].lower(), translation[0].lower())
    if cpd_name in USED_CPDS:
        return
    USED_CPDS.add(cpd_name)
    sort_dict = {"genotype": ["Tim3+", "Tim3-"], "treatment": ["Normoxia", "Hypoxia"], "time": [3, 6]}
    sort_index = {"genotype": 0, "treatment": 1, "time": 2}
    df = df[df["name"] == cpd_name]
    working_df = df.copy()
    sort_values = []
    for ion_rel in zip(working_df["ion_relation"]):
        ion_rel = ion_rel[0]
        iso = ion_rel.split(",")[0]
        adduct = ion_rel.split(",")[1]
        val = iso[-1]
        if val == "C":
            val = 1
        sort_values.append(adduct + str(val))
    working_df["sort"] = sort_values
    if not df.empty:
        working_df.sort_values("sort", inplace=True, ascending=True)
        fig, ax = plt.subplots()
        #num = df.shape[0]
        sorted = []
        for s in get_sample_names(working_df):
            for val in sort_dict[sort]:
                if sample_mapping[s.split("_")[-1]][sort_index[sort]] == val:
                    sorted.append(s)
        g_hatch = [color_map[sample_mapping[x.split("_")[-1]]][0] for x in sorted]
        g_color = [color_map[sample_mapping[x.split("_")[-1]]][1] for x in sorted]
        sub_title = cpd_name + ", " + MODE + ", " + STD_MODE
        if SIG_ONLY:
            sub_title += ", Significant Only"

        #g = working_df.boxplot(column=sorted)

        g = working_df.plot(x="sort", y=sorted, kind="bar", title=sub_title, ax=ax, legend=False)        
        bars = ax.patches
        ax.set_xticklabels([z + "\nr=" + str(r) + "\np_C=" + str(round(y,2)) + "\np_t=" + str(round(y2,2)) + "\n\np_tim_h=" + str(round(a, 2)) +"\np_tim_n=" + str(round(b, 2)) + "\np_o2_t+=" + str(round(c, 2)) + "\np_o2_t-=" + str(round(d, 2)) for z, y, y2, r, id, m, a, b, c, d in 
                            zip(working_df['ion_relation'], 
                                working_df['p_genotype'], 
                                working_df['p_treatment'], 
                                working_df['rtime'], 
                                working_df['id_number'],
                                working_df['mz'],
                                working_df['p_tim3_hypoxia'],
                                working_df['p_tim3_normoxia'],
                                working_df['p_o2_tim3wt'],
                                working_df['p_o2_tim3ko'])], rotation=0)
        box = ax.get_position()
        plt.ylabel("Log2 peak intensity")
        ax.set_position([box.x0, box.y0, box.width * 1, box.height * 1.1])
        if hour == '3hr':
            ax.legend(handles=patches3hr, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
        elif hour == '6hr':
            ax.legend(handles=patches6hr, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
        hatches = []
        colors = []
        for c, h in zip(g_color, g_hatch):
            for i in range(int(len(bars) / len(g_hatch))):
                hatches.append(h)
                colors.append(c)
        for bar, h, c in zip(bars, hatches, colors):
            bar.set_fill(False)
            bar.set_color(c)
            if not ALT_COLORS:
                bar.set_hatch(h)
            else:
                if h == '/':
                    bar.set_fill(True)
                else:
                    bar.set_fill(False)
        plt.tight_layout()
        plt.xticks(fontsize=6)
        plt.suptitle(sup_title)
        plt.savefig(os.path.join(os.path.abspath(RESULT_DIR), basename + cpd_name + ".png"))
        plt.close()

rtime_map = build_rtime_map()
sample_mapping = build_sample_map()
color_map = build_color_map()
patches3hr, patches6hr = build_patches(color_map)
make_results_dir(RESULT_DIR)

import tqdm
for title, basename, df, hour in list(zip(*find_annotation_files("."))):    
    print(title)
    USED_CPDS = set()
    df = filter_df(df)
    df = compare_rts(df, rtime_map)
    cpd_names = extract_cpds(df)
    for cpd_name in tqdm.tqdm(cpd_names):
        plot_cpd(cpd_name, df, title, hour)








