import pandas as pd
import sys
from mass2chem.formula import PROTON, ELECTRON, calculate_formula_mass, atom_mass_dict, parse_chemformula_dict, calculate_mass
from khipu.epdsConstructor import epdsConstructor
from khipu.extended import *
from jms.io import read_table_to_peaks
from intervaltree import IntervalTree
from copy import deepcopy
import numpy as np
import re
import matplotlib.pyplot as plt
import json

RT_ERR_WITHIN = .5
RT_ERR_BETWEEN = 2.5
MZ_ERR = 5

NUM_DMPA = 3
DMPA_TOL_RATIO = 1/2
DMPA_TOL = [DMPA_TOL_RATIO, 1/DMPA_TOL_RATIO] #[.67, 1.5]
MIN_GOOD_MATCHES = 3
DMPA_dict = parse_chemformula_dict("C10H11O1N1")
DMPA_mass = calculate_formula_mass("C10H11O1N1")
delta_13C = atom_mass_dict["[13C]"] - atom_mass_dict["[12C]"]
print(delta_13C)

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

standards = pd.read_excel(sys.argv[1])
recalc_masses = []
for formula in standards["Formula"]:
    recalc_masses.append(calculate_formula_mass(formula))
standards["Monoisotopic Molecular Weight"] = recalc_masses

uft = pd.read_csv(sys.argv[2], sep="\t")
dft = pd.read_csv(sys.argv[3], sep="\t")

def overlaps(a, b, ERR):
    # per The Unfun Cat stackoverflow
    # https://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
    return min(a[1] + ERR, b[1] + ERR) - max(a[0] - ERR, b[0] - ERR)

def build_groups(ft):
    groups = {}
    for column in ft.columns:
        m = re.search("Group(\d+)", column)
        if m:
            group_num = int(m.group(1))
            if group_num not in groups:
                groups[group_num] = []
            groups[group_num].append(column)
        else:
            if "Process_Blank" in column:
                if "Blanks" not in groups:
                    groups["Blanks"] = []
                groups["Blanks"].append(column)
    return groups

def build_derivatives(standards, dft_groups, uft_groups):
    underivatized = {}
    derivatized = {}
    all = {}
    for cpd_name, cpd_group, cpd_formula, num_COOH in zip(standards['Name'], standards['Group'], standards['Formula'], standards['NumCarboxyl']):
        cpd_name = cpd_name.strip().lower()
        if cpd_group in dft_groups and cpd_group in uft_groups:
            cpd = {}
            cpd["formula_dict"] = parse_chemformula_dict(cpd_formula)
            cpd["group"] = cpd_group
            for n in range(num_COOH + 1):
                if n != 0:
                    pos_dmpa_13C = [n*2]
                else:
                    pos_dmpa_13C = []
                for num_dmpa_13C in [0] + pos_dmpa_13C:
                    for adduct in [(PROTON, "H+"),
                                (0, "+"), 
                                (atom_mass_dict["K"], "K+"), 
                                (atom_mass_dict["Na"], "Na+"), 
                                (calculate_formula_mass("NH4"), "NH4+"),
                                (calculate_formula_mass("CH3CN"), "ACN+")]:
                        derivative = dict()
                        derivative["group"] = cpd_group
                        derivative["parent"] = cpd_name
                        derivative["formula_dict"] = deepcopy(cpd["formula_dict"])
                        derivative["name"] = cpd_name
                        derivative["DMPA"] = n
                        if n != 0:
                            derivative["name"] = derivative["name"] + "_dmpa*" + str(n)
                            for _ in range(n):
                                for ele, count in DMPA_dict.items():
                                    if ele not in derivative["formula_dict"]:
                                        derivative["formula_dict"][ele] = 0
                                    derivative["formula_dict"][ele] += count
                        if num_dmpa_13C != 0:
                            derivative["isotopologue"] = "m+13C*" + str(num_dmpa_13C)
                            derivative["formula_dict"]["C"] -= num_dmpa_13C
                            derivative["formula_dict"]["[13C]"] = num_dmpa_13C
                        else:
                            derivative["isotopologue"] = "M0"
                        derivative["adduct"] = adduct[1]
                        derivative["mass"] = calculate_mass(derivative["formula_dict"])
                        derivative["adduct_mass"] = derivative["mass"] + adduct[0] - ELECTRON
                        derivative["confirmed_retention_times"] = {}
                        derivative["isotope_confirmed_retention_times"] = []
                        derivative["low_confidence_retention_times"] = []
                        if derivative["DMPA"]:
                            id = 'd' + str(len(derivatized))
                            derivatized[id] = derivative
                        else:
                            id = 'u' + str(len(underivatized))
                            underivatized[id] = derivative
                        all[id] = derivative
                        #print(json.dumps(derivative, indent=4))
    return underivatized, derivatized, all

def build_std_tree(stds):
    tree = IntervalTree()
    for c_id, cpd in stds.items():
        mass = cpd["adduct_mass"]
        mass_err = mass / 1e6 * MZ_ERR
        tree.addi(mass - mass_err, mass + mass_err, c_id)
    return tree

def build_ft_trees(ft):
    mz_tree = IntervalTree()
    rt_tree = IntervalTree()
    for f_id, mz, rtime in zip(ft['id_number'], ft['mz'], ft['rtime']):
        mz_err = mz / 1e6 * MZ_ERR
        mz_tree.addi(mz - mz_err, mz + mz_err, f_id)
        rt_tree.addi(rtime - RT_ERR_BETWEEN, rtime + RT_ERR_BETWEEN, f_id)
    return mz_tree, rt_tree

def organize_derivatives(derivatives):
    partners = {}
    for id1, derivative_1 in derivatives.items():
        for id2, derivative_2 in derivatives.items():
            if derivative_1['name'] == derivative_2['name'] and derivative_1['adduct'] == derivative_2['adduct']:
                if derivative_1['isotopologue'] != derivative_2['isotopologue']:
                    if derivative_1['isotopologue'] == 'M0':
                        partners[id1] = id2
    return partners

dft_groups = build_groups(dft)
uft_groups = build_groups(uft)

underivatized_stds, derivatized_stds, all_stds = build_derivatives(standards, dft_groups, uft_groups)
u_tree = build_std_tree(underivatized_stds)
d_tree = build_std_tree(derivatized_stds)
a_tree = build_std_tree(all_stds)

dft_mz_tree, dft_rt_tree = build_ft_trees(dft) 
uft_mz_tree, uft_rt_tree = build_ft_trees(uft)
partners = organize_derivatives(derivatized_stds)

# FIND UNDERIVATIZED STANDARDS
for mz, rtime, left_rtime, right_rtime in zip(dft['mz'], dft['rtime'], dft['rtime_left_base'], dft['rtime_right_base']):
    #print(mz)
    for std_match in u_tree.at(mz):
        #print("\t", std_match.data)
        u_std = underivatized_stds[std_match.data]
        for uft_pos_m0 in uft_mz_tree.at(mz):
            #print("\t\t", uft_pos_m0)
            uft_pos_m0_rtime_interval = uft[uft['id_number'] == uft_pos_m0.data][['rtime_left_base', 'rtime_right_base']].values
            if overlaps([left_rtime, right_rtime], uft_pos_m0_rtime_interval[0], RT_ERR_BETWEEN):
                #print("\t\t\t overlaps ", uft_pos_m0.data)
                for dft_pos_m13c in dft_mz_tree.at(mz + delta_13C):
                    dft_pos_m13c_rtime_interval = dft[dft['id_number'] == dft_pos_m13c.data][['rtime_left_base', 'rtime_right_base']].values
                    if overlaps([left_rtime, right_rtime], dft_pos_m13c_rtime_interval[0], RT_ERR_WITHIN):
                        #print("\t\t\t\t overlaps iso dfs ", dft_pos_m13c.data)
                        for uft_pos_m13c in uft_mz_tree.at(mz + delta_13C):
                            uft_pos_m13c_rtime_interval = uft[uft['id_number'] == uft_pos_m13c.data][['rtime_left_base', 'rtime_right_base']].values
                            if overlaps([left_rtime, right_rtime],uft_pos_m13c_rtime_interval[0], RT_ERR_BETWEEN):
                                #print("\t\t\t\t\t overlaps iso uft ", uft_pos_m13c.data)
                                min_left_rtime = min(left_rtime, uft_pos_m0_rtime_interval[0][0], dft_pos_m13c_rtime_interval[0][0], uft_pos_m13c_rtime_interval[0][0])
                                max_right_rtime = max(right_rtime, uft_pos_m0_rtime_interval[0][1], dft_pos_m13c_rtime_interval[0][1], uft_pos_m13c_rtime_interval[0][1])
                                u_std['confirmed_retention_times'].append([min_left_rtime, max_right_rtime])

for mz, rtime, left_rtime, right_rtime in zip(uft['mz'], uft['rtime'], uft['rtime_left_base'], uft['rtime_right_base']):
    for pos_m0 in u_tree.at(mz):
        print( pos_m0.data, underivatized_stds[pos_m0.data]["name"])
        for uft_pos_m13c in uft_mz_tree.at(mz+delta_13C):
            print("\t", uft_pos_m13c.data)
            pos_m13c_rtime_interval = uft[uft['id_number'] == uft_pos_m13c.data][['rtime_left_base', 'rtime_right_base']].values
            if left_rtime < rtime < right_rtime:
                if type(underivatized_stds[pos_m0.data]["confirmed_retention_times"]) is not list:
                    underivatized_stds[pos_m0.data]["confirmed_retention_times"] = []
                underivatized_stds[pos_m0.data]["confirmed_retention_times"].append([left_rtime, right_rtime])
    
# FIND DERIVATIZED STANDARDS
for d_std_id, std_m0 in derivatized_stds.items():
    if d_std_id in partners:
        std_m13c = derivatized_stds[partners[d_std_id]]
        dft_m0_mz_matches = dft_mz_tree.at(std_m0['adduct_mass'])
        dft_m13c_mz_matches = dft_mz_tree.at(std_m13c['adduct_mass'])
        if dft_m0_mz_matches and dft_m13c_mz_matches:
            for pos_m0 in dft_m0_mz_matches:
                pos_m0_fid = pos_m0.data
                pos_m0_rtime_interval = dft[dft['id_number'] == pos_m0_fid][['rtime_left_base', 'rtime_right_base']].values
                for pos_m13c in dft_m13c_mz_matches:
                    pos_m13c_fid = pos_m13c.data
                    pos_m13c_rtime_interval = dft[dft['id_number'] == pos_m13c_fid][['rtime_left_base', 'rtime_right_base']].values
                    if overlaps(pos_m0_rtime_interval[0], pos_m13c_rtime_interval[0], RT_ERR_WITHIN):
                        pos_m0_values = dft[dft['id_number'] == pos_m0_fid][dft_groups[std_m0["group"]]].values[0]
                        pos_m13c_values = dft[dft['id_number'] == pos_m13c_fid][dft_groups[std_m13c["group"]]].values[0]
                        good_matches = 0
                        ratios = []
                        for pos_m0_val, pos_m13c_val, sample_name in zip(pos_m0_values, pos_m13c_values, dft_groups[std_m0["group"]]):
                            ratio = pos_m0_val / pos_m13c_val if pos_m13c_val else np.inf
                            intensities = []
                            if DMPA_TOL[0] < ratio < DMPA_TOL[1]:
                                good_matches += 1
                                ratios.append(ratio)
                                intensities.append(pos_m0_val)
                        by_adduct = False
                        if good_matches > MIN_GOOD_MATCHES and by_adduct:
                            if std_m0["adduct"] not in std_m0["confirmed_retention_times"]:
                                std_m0['confirmed_retention_times'][std_m0['adduct']] = [None, 0]
                            if np.mean(intensities) > std_m0['confirmed_retention_times'][std_m0['adduct']][1]:
                                min_left_rtime = min(pos_m0_rtime_interval[0][0], pos_m13c_rtime_interval[0][0])
                                max_right_rtime = max(pos_m0_rtime_interval[0][1], pos_m13c_rtime_interval[0][1])
                                std_m0['confirmed_retention_times'][std_m0['adduct']] = [[min_left_rtime, max_right_rtime, np.mean(ratios), np.mean(intensities), std_m0['adduct']], np.mean(intensities)]
                        elif good_matches > MIN_GOOD_MATCHES and not by_adduct:
                            if "zzz" not in std_m0["confirmed_retention_times"]:
                                std_m0['confirmed_retention_times']["zzz"] = [None, 0]
                            if np.mean(intensities) > std_m0['confirmed_retention_times']["zzz"][1]:
                                min_left_rtime = min(pos_m0_rtime_interval[0][0], pos_m13c_rtime_interval[0][0])
                                max_right_rtime = max(pos_m0_rtime_interval[0][1], pos_m13c_rtime_interval[0][1])
                                std_m0['confirmed_retention_times']["zzz"] = [[min_left_rtime, max_right_rtime, np.mean(ratios), np.mean(intensities), std_m0['adduct']], np.mean(intensities)]

for std in derivatized_stds.values():
    z = []
    for adduct, values in std['confirmed_retention_times'].items():
        if values[0] is not None:
            z.append(values[0])
    std['confirmed_retention_times'] = z 


combined = {}
for std in [derivatized_stds, underivatized_stds]:
    for k, v in std.items():
        if v['low_confidence_retention_times'] or v['confirmed_retention_times'] or v['isotope_confirmed_retention_times']:
            for translate in TRANSLATION:
                if translate[0].lower() in v["name"].lower():
                    working_name = v["name"].lower().replace(translate[0].lower(), translate[1].lower())
                    if working_name  not in combined:
                        combined[working_name] = {
                            "low_confidence_retention_times": [],
                            "confirmed_retention_times": [],
                            "isotope_confirmed_retention_times": [],
                            "all_retention_times": [],
                            "search_mass": []
                        }
                    combined[working_name]['low_confidence_retention_times'].extend(v['low_confidence_retention_times'])
                    combined[working_name]['confirmed_retention_times'].extend(v['confirmed_retention_times'])
                    combined[working_name]['isotope_confirmed_retention_times'].extend(v['isotope_confirmed_retention_times'])
                    combined[working_name]['all_retention_times'].extend(v['low_confidence_retention_times'])
                    combined[working_name]['all_retention_times'].extend(v['confirmed_retention_times'])
                    combined[working_name]['all_retention_times'].extend(v['isotope_confirmed_retention_times'])
                    combined[working_name]['search_mass'].append((v["adduct_mass"], v["adduct"]))

            if v["name"] not in combined:
                combined[v["name"]] = {
                    "low_confidence_retention_times": [],
                    "confirmed_retention_times": [],
                    "isotope_confirmed_retention_times": [],
                    "all_retention_times": [],
                    "search_mass": []
                }
            combined[v["name"]]['low_confidence_retention_times'].extend(v['low_confidence_retention_times'])
            combined[v["name"]]['confirmed_retention_times'].extend(v['confirmed_retention_times'])
            combined[v["name"]]['isotope_confirmed_retention_times'].extend(v['isotope_confirmed_retention_times'])
            combined[v["name"]]['all_retention_times'].extend(v['low_confidence_retention_times'])
            combined[v["name"]]['all_retention_times'].extend(v['confirmed_retention_times'])
            combined[v["name"]]['all_retention_times'].extend(v['isotope_confirmed_retention_times'])
            combined[v["name"]]['search_mass'].append((v["adduct_mass"], v["adduct"]))


if not by_adduct:
    for k, v in combined.items():
        highest_intensity = 0
        filtered = None
        for x in v["confirmed_retention_times"]:
            print(x)
            if x[-2] > highest_intensity:
                highest_intensity = x[-2]
                filtered = x
        if filtered is not None:
            v['confirmed_retention_times'] = [filtered]
        else:
            v['confirmed_retention_times'] = []


if 1 in uft_groups and uft_groups[1]:
    json.dump(combined, open("rtime_mapping_g1_feature.json", "w+"), indent=4)
else:
    json.dump(combined, open("rtime_mapping_g2_feature.json", "w+"), indent=4)
#print(json.dumps(combined, indent=4))