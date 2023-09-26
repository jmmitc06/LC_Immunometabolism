import pandas as pd
import sys
from mass2chem.formula import PROTON, ELECTRON, calculate_formula_mass, atom_mass_dict, parse_chemformula_dict, calculate_mass
from mass2chem.search import build_centurion_tree, find_all_matches_centurion_indexed_list

from khipu.epdsConstructor import epdsConstructor
from khipu.extended import isotope_search_patterns, adduct_search_patterns, extended_adducts
from jms.io import read_table_to_peaks
from intervaltree import IntervalTree
from copy import deepcopy
import numpy as np
import re
import matplotlib.pyplot as plt
import json

#adduct_search_patterns = [(PROTON, "H+"), (0, "+"),  (atom_mass_dict["K"], "K+"),  (atom_mass_dict["Na"], "Na+"),  (calculate_formula_mass("NH4"), "NH4+"), (calculate_formula_mass("CH3CN"), "ACN+")]

RT_TOL_WITHIN = 2
RT_ERR_WITHIN = .5
RT_ERR_BETWEEN = 5
MZ_ERR = 5

NUM_DMPA = 3
DMPA_TOL_RATIO = 1/2
DMPA_TOL = [DMPA_TOL_RATIO, 1/DMPA_TOL_RATIO] #[.67, 1.5]
MIN_GOOD_MATCHES = 3
DMPA_dict = parse_chemformula_dict("C10H11O1N1")
DMPA_mass = calculate_formula_mass("C10H11O1N1")
delta_13C = atom_mass_dict["[13C]"] - atom_mass_dict["[12C]"]

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

def build_unadducted_derivatives(standards, dft_groups, uft_groups):
    underivatized = {}
    derivatized = {}
    all = {}
    for cpd_name, cpd_group, cpd_formula in zip(standards['Name'], standards['Group'], standards['Formula']):
        cpd_name = cpd_name.strip().lower()
        if cpd_group in dft_groups and cpd_group in uft_groups:
            cpd = {}
            cpd["formula_dict"] = parse_chemformula_dict(cpd_formula)
            cpd["group"] = cpd_group
            for n in range(NUM_DMPA + 1):
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
                derivative["isotopologue"] = "M0"
                derivative["mass"] = calculate_mass(derivative["formula_dict"])
                derivative["confirmed_retention_times"] = []
                derivative["isotope_confirmed_retention_times"] = []
                derivative["low_confidence_retention_times"] = []
                if derivative["DMPA"]:
                    id = 'd' + str(len(derivatized))
                    derivatized[id] = derivative
                else:
                    id = 'u' + str(len(underivatized))
                    underivatized[id] = derivative
                all[id] = derivative
    return underivatized, derivatized, all

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
                        derivative["confirmed_retention_times"] = []
                        derivative["isotope_confirmed_retention_times"] = []
                        derivative["low_confidence_retention_times"] = []
                        if derivative["DMPA"]:
                            id = 'd' + str(len(derivatized))
                            derivatized[id] = derivative
                        else:
                            id = 'u' + str(len(underivatized))
                            underivatized[id] = derivative
                        all[id] = derivative
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

def build_khipu(ft_path):
    features = read_table_to_peaks(ft_path, feature_id=0, mz_col=1, rtime_col=2, intensity=(6,7), full_extract=True)
    for f in features:
        f['id'] = f['id_number']
    ECON = epdsConstructor(features, mode='pos')
    khipu_dict = ECON.peaks_to_epdDict(
        isotope_search_patterns = isotope_search_patterns[:6],
        adduct_search_patterns = adduct_search_patterns, #[(PROTON, "H+"), (atom_mass_dict["K"], "K+"), (atom_mass_dict["Na"], "Na+"), (calculate_formula_mass("NH4"), "NH4+"),(calculate_formula_mass("CH3CN"), "ACN+")],
        extended_adducts = extended_adducts,
        mz_tolerance_ppm=5, 
        rt_tolerance= RT_TOL_WITHIN,             # cleaner if RT window smaller
        charges=[1, 2, 3],
        has_parent_masstrack=True,
    )
    return khipu_dict

def build_neutral_tree(epds):
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
    return build_centurion_tree(neutrals)

def organize_empcpd(empcpd):
    organized_peaks = {}
    for peak in empcpd["MS1_pseudo_Spectra"]:
        modification = peak["modification"]
        if modification not in organized_peaks:
            organized_peaks[modification] = {}
        isotope = peak["isotope"]
        organized_peaks[modification][isotope] = peak
    empcpd["organized_peaks"] = organized_peaks
    return empcpd

def calc_rt_interval(empcpd):
    max_rtime_right = -np.inf
    min_rtime_left = np.inf
    for peak in empcpd["MS1_pseudo_Spectra"]:
        rtime_left_base = float(peak["rtime_left_base"])
        rtime_right_base = float(peak["rtime_right_base"])
        max_rtime_right = max(max_rtime_right, rtime_right_base)
        min_rtime_left = min(min_rtime_left, rtime_left_base)
    return [min_rtime_left, max_rtime_right]

dft_groups = build_groups(dft)
uft_groups = build_groups(uft)

underivatized_stds, derivatized_stds, all_stds = build_unadducted_derivatives(standards, dft_groups, uft_groups)
#u_tree = build_std_tree(underivatized_stds)
#d_tree = build_std_tree(derivatized_stds)
#a_tree = build_std_tree(all_stds)

#dft_mz_tree, dft_rt_tree = build_ft_trees(dft) 
#uft_mz_tree, uft_rt_tree = build_ft_trees(uft)
#partners = organize_derivatives(derivatized_stds)

uft_empcpd = build_khipu(sys.argv[2])
dft_empcpd = build_khipu(sys.argv[3])

uft_centurion_tree = build_neutral_tree(uft_empcpd)
dft_centurion_tree = build_neutral_tree(dft_empcpd)

for std in underivatized_stds.values():
    if std['isotopologue'] == 'M0':
        mz = std["mass"]
        uft_matches = find_all_matches_centurion_indexed_list(mz, uft_centurion_tree, MZ_ERR)
        for match in uft_matches:
            matching_empcpd = uft_empcpd[match['id']]
            matching_empcpd = organize_empcpd(matching_empcpd)
            for modification in matching_empcpd["organized_peaks"].keys():
                m0 = "M0"
                m13c = "13C/12C"
                if m0 in matching_empcpd["organized_peaks"][modification] and m13c in matching_empcpd["organized_peaks"][modification]:
                    std["isotope_confirmed_retention_times"].append(calc_rt_interval(matching_empcpd))
            std["low_confidence_retention_times"].append(calc_rt_interval(matching_empcpd))
        dft_matches = find_all_matches_centurion_indexed_list(mz, dft_centurion_tree, MZ_ERR)
        for uft_match in uft_matches:
            uft_rtime_interval = calc_rt_interval(uft_empcpd[uft_match['id']])
            for dft_match in dft_matches:
                dft_rtime_interval = calc_rt_interval(dft_empcpd[dft_match['id']])
                if overlaps(uft_rtime_interval, dft_rtime_interval, RT_ERR_BETWEEN):
                    combine_interval = [min(uft_rtime_interval[0], dft_rtime_interval[0]), max(uft_rtime_interval[1], dft_rtime_interval[1])]
                    std["confirmed_retention_times"].append(combine_interval)

for std in derivatized_stds.values():
    if std['isotopologue'] == 'M0':
        mz = std["mass"]
        matches = find_all_matches_centurion_indexed_list(mz, dft_centurion_tree, MZ_ERR)
        for match in matches:
            matching_empcpd = dft_empcpd[match['id']]
            matching_empcpd = organize_empcpd(matching_empcpd)
            num_DMPA = std["DMPA"]
            if num_DMPA != 0:
                for modification in matching_empcpd["organized_peaks"].keys():
                    m0 = "M0"
                    iso = "13C/12C*" + str(num_DMPA)
                    if m0 in matching_empcpd["organized_peaks"][modification] and iso in matching_empcpd["organized_peaks"][modification]:
                        num_good_matches = 0
                        for sample in dft_groups[std["group"]]:
                            print(sample)
                            m0_intensity = matching_empcpd["organized_peaks"][modification][m0][sample]
                            iso_intensity = matching_empcpd["organized_peaks"][modification][m0][sample]
                            if iso_intensity:
                                ratio = float(m0_intensity) / float(iso_intensity)
                                if DMPA_TOL[0] < ratio < DMPA_TOL[1]:
                                    num_good_matches += 1
                        if num_good_matches > MIN_GOOD_MATCHES:
                            if "confirmed_retention_times" not in std:
                                std["confirmed_retention_times"] = []
                            std["confirmed_retention_times"].append(calc_rt_interval(matching_empcpd))
                            
combined = {}
for v in list(derivatized_stds.values()) + list(underivatized_stds.values()):
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
                                "search_mass": None
                            }
                        combined[working_name]['low_confidence_retention_times'].extend(v['low_confidence_retention_times'])
                        combined[working_name]['confirmed_retention_times'].extend(v['confirmed_retention_times'])
                        combined[working_name]['isotope_confirmed_retention_times'].extend(v['isotope_confirmed_retention_times'])
                        combined[working_name]['all_retention_times'].extend(v['low_confidence_retention_times'])
                        combined[working_name]['all_retention_times'].extend(v['confirmed_retention_times'])
                        combined[working_name]['all_retention_times'].extend(v['isotope_confirmed_retention_times'])
                        combined[working_name]['search_mass'] = v["mass"]

                if v["name"] not in combined:
                    combined[v["name"]] = {
                        "low_confidence_retention_times": [],
                        "confirmed_retention_times": [],
                        "isotope_confirmed_retention_times": [],
                        "all_retention_times": []
                    }
                combined[v["name"]]['low_confidence_retention_times'].extend(v['low_confidence_retention_times'])
                combined[v["name"]]['confirmed_retention_times'].extend(v['confirmed_retention_times'])
                combined[v["name"]]['isotope_confirmed_retention_times'].extend(v['isotope_confirmed_retention_times'])
                combined[v["name"]]['all_retention_times'].extend(v['low_confidence_retention_times'])
                combined[v["name"]]['all_retention_times'].extend(v['confirmed_retention_times'])
                combined[v["name"]]['all_retention_times'].extend(v['isotope_confirmed_retention_times'])
                combined[v["name"]]['search_mass'] = v["mass"]

if 1 in uft_groups and uft_groups[1]:
    json.dump(combined, open("rtime_mapping_g1_khipu.json", "w+"), indent=4)
else:
    json.dump(combined, open("rtime_mapping_g2_khipu.json", "w+"), indent=4)
#print(json.dumps(combined, indent=4))