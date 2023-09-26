import json
import numpy as np
from mass2chem.search import build_centurion_tree, find_all_matches_centurion_indexed_list
from jms.io import read_table_to_peaks
from khipu.epdsConstructor import epdsConstructor
from mass2chem.formula import PROTON, ELECTRON, calculate_formula_mass, parse_chemformula_dict, atom_mass_dict, calculate_mass
import mummichog.JSON_metabolicModels as mcg
from intervaltree import IntervalTree
from copy import deepcopy
import sys
import pandas as pd

dmPA_mass = calculate_formula_mass("C10H11O1N1")
dmPA_dict = parse_chemformula_dict("C10H11O1N1")
hmfn = mcg.metabolicModels['human_model_mfn']
num_DMPA = 3

Glycolysis = []
for c in ['C00068', 'C00027', 'C00149', 'C00125', 'C00126', 'C00024', 'C00288', 'C00035', 'C00036', 'C00900', 
          'C05378', 'C00103', 'C05125', 'C00022', 'C00579', 'C00668', 'C00092', 'C00020', 'C01342', 'C00236', 
          'C00010', 'C15973', 'C15972', 'C00111', 'C00033', 'C00031', 'C00354', 'C00118', 'C00074', 'C00197', 
          'C00085', 'C01136', 'C00665', 'C05345', 'C00227', 'C01231', 'C00084', 'C00631', 'C00221', 'C01159', 
          'C00248', 'C01172', 'C00267', 'C00469', 'C00661', 'C00044', 'C00186', 'C05993', 'C00169']:
    if c in hmfn["Compounds"]:
        Glycolysis.append((c, hmfn["Compounds"][c]))

TCA = []
for c in ['C00068', 'C00122', 'C00149', 'C00024', 'C00028', 'C00026', 'C00104', 'C05379', 'C01169', 'C00399', 
          'C00016', 'C00390', 'C00579', 'C00311', 'C00010', 'C00091', 'C00036', 'C00158', 'C15973', 'C15972', 
          'C00033', 'C00030', 'C00074', 'C00081', 'C00248', 'C00035', 'C00042', 'C01352', 'C00417', 'C00044', 
          'C05381']:
    if c in hmfn["Compounds"]:
        TCA.append((c, hmfn["Compounds"][c]))   

carboxy_numbers = json.load(open("./stds_carboxy.json"))

derivatives = {}
observed = set()
for cpd_id, cpd in Glycolysis + TCA:
    if cpd_id not in observed:
        if cpd["name"] == "alpha-D-Glucose 6-phosphate":
            cpd["formula"] = "C6H13O9P"
        if cpd["name"] == "D-Glyceraldehyde 3-phosphate":
            cpd["formula"] = "C3H7O6P"
        if cpd["name"] == "Ubiquinol":
            cpd["formula"] = "C59H92O4"
        observed.add(cpd_id)
        cpd["formula_dict"] = parse_chemformula_dict(cpd["formula"])
        if cpd["formula"]:
            num_carboxy = carboxy_numbers[cpd_id]
            if num_carboxy is None:
                num_carboxy = num_DMPA
            for n in range(num_carboxy + 1):
                num_C = cpd["formula_dict"]["C"] if "C" in cpd["formula_dict"] else 0
                for c in range(num_C + 1):
                    for adduct in [(PROTON, "H+"),
                                (0, "+"), 
                                (atom_mass_dict["K"], "K+"), 
                                (atom_mass_dict["Na"], "Na+"), 
                                (calculate_formula_mass("NH4"), "NH4+"),
                                (calculate_formula_mass("CH3CN"), "ACN+")]:
                        derivative = dict()
                        derivative["parent_formula"] = deepcopy(cpd["formula"])
                        derivative["parent"] = cpd_id
                        derivative["formula_dict"] = deepcopy(cpd["formula_dict"])
                        derivative["name"] = cpd["name"].split(";")[0]
                        if n != 0:
                            derivative["name"] = derivative["name"] + "_dmpa*" + str(n)
                            for _ in range(n):
                                for ele, count in dmPA_dict.items():
                                    if ele not in derivative["formula_dict"]:
                                        derivative["formula_dict"][ele] = 0
                                    derivative["formula_dict"][ele] += count
                        if c != 0:
                            derivative["isotopologue"] = "m+13C*" + str(c)
                            derivative["formula_dict"]["C"] -= c
                            derivative["formula_dict"]["[13C]"] = c
                        else:
                            derivative["isotopologue"] = "M0"
                        derivative["adduct"] = adduct[1]
                        derivative["mass"] = calculate_mass(derivative["formula_dict"])
                        derivative["adduct_mass"] = derivative["mass"] + adduct[0] - ELECTRON
                        derivatives['d' + str(len(derivatives))] = derivative

annotation_tree = IntervalTree()
for k, v in derivatives.items():
    try:
        mass = v["adduct_mass"]
        mass_error = mass / 1e6 * 5
        annotation_tree.addi(mass - mass_error, mass + mass_error, k)
    except:
        pass

in_tables = [
    #"/Users/mitchjo/Analyses/DmPA_Lucas/LucasPellets_051823_09_07_2023/asari_results/_asari_project_97194854/preferred_Feature_table.tsv",
    "/Users/mitchjo/Analyses/DmPA_Lucas/LucasPellets_051823_09_07_2023/filtered_feature_tables/log_transformed_normalized_lucas_pellets_for_SL_Feature_table.tsv",
    "/Users/mitchjo/Analyses/DmPA_Lucas/LucasPellets_051823_09_07_2023/filtered_feature_tables/normalized_90_samples_drop_02_Feature_table.tsv",
    "/Users/mitchjo/Analyses/DmPA_Lucas/LucasMedia_051523_09_07_2023/filtered_feature_tables/normalized_90_samples_drop_02_Feature_table.tsv",
    #"/Users/mitchjo/Analyses/DmPA_Lucas/cells_data3hr_stats_09_08_23.tsv",
    #"/Users/mitchjo/Analyses/DmPA_Lucas/cells_data6hr_stats_09_08_23.tsv",
    #"/Users/mitchjo/Analyses/DmPA_Lucas/supernatant_data3hr_stats_09_08_23.tsv",
    #"/Users/mitchjo/Analyses/DmPA_Lucas/supernatant_data6hr_stats_09_08_23.tsv"
]

for table_path in in_tables:
    ftable = pd.read_csv(table_path, sep="\t")
    annotations = []
    ion_relations = []
    isotopologues = []
    annotated_ids = []
    out_ftable = []
    for id_no, mz, rtime in zip(ftable['id_number'], ftable['mz'], ftable['rtime']):
        match = annotation_tree.at(float(mz))
        if match:
            annotated_ids.append(id_no)
            for m in match:
                deriv_id = m.data
                annotations.append(derivatives[deriv_id]['name'])
                ion_relations.append(derivatives[deriv_id]['adduct'])
                isotopologues.append(derivatives[deriv_id]['isotopologue'])
                out_dict_row = dict(zip(ftable.columns, ftable[ftable['id_number'] == id_no][ftable.columns].values[0]))
                out_dict_row["name"] = derivatives[deriv_id]['name']
                out_dict_row["ion_relation"] = derivatives[deriv_id]['adduct']
                out_dict_row["isotopologue"] = derivatives[deriv_id]["isotopologue"]
                out_dict_row["formula"] = derivatives[deriv_id]["parent_formula"]
                out_ftable.append(out_dict_row)
    report_table = pd.DataFrame(out_ftable)
    report_table.to_csv(table_path.replace(".tsv", "feature_annot.tsv"), sep="\t")
