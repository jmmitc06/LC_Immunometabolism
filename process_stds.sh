#pcpfm assemble_experiment_from_CSV ./Lucas_Deriv_Stds_G1_09_08_2023 /Users/mitchjo/Analyses/DmPA_Lucas/Standards_Deriv/08292023_Lucas_glycolysis_TCA_authenticstandards_derivatization/08292023_Lucas_glycolysis_TCA_authenticstandards_derivatization.csv auto --filter='{"File Name": {"lacks": ["Group2"]}}'
#pcpfm convert_to_mzML ./Lucas_Deriv_Stds_G1_09_08_2023 $(which mono) ~/Projects/PythonCentricPipelineForMetabolomics-1/ThermoRawFileParser/ThermoRawFileParser.exe
#pcpfm asari_full_processing ./Lucas_Deriv_Stds_G1_09_08_2023 
#pcpfm feature_QCQA ./Lucas_Deriv_Stds_G1_09_08_2023 preferred --all --save_figs
pcpfm blank_masking ./Lucas_Deriv_Stds_G1_09_08_2023 preferred --new_table_moniker=preferred_blank_masked --blank_intensity_ratio=3 --blank_type="solblank"
pcpfm blank_masking ./Lucas_Deriv_Stds_G1_09_08_2023 preferred_blank_masked --new_table_moniker=preferred_blank_masked --blank_intensity_ratio=3 --blank_type="blank"
pcpfm drop_samples ./Lucas_Deriv_Stds_G1_09_08_2023 preferred_blank_masked --new_table_moniker=preferred_blank_masked --field_value="solblank"
pcpfm drop_samples ./Lucas_Deriv_Stds_G1_09_08_2023 preferred_blank_masked --new_table_moniker=preferred_blank_masked --field_value="blank"


#pcpfm assemble_experiment_from_CSV ./Lucas_Deriv_Stds_G2_09_08_2023 /Users/mitchjo/Analyses/DmPA_Lucas/Standards_Deriv/08292023_Lucas_glycolysis_TCA_authenticstandards_derivatization/08292023_Lucas_glycolysis_TCA_authenticstandards_derivatization.csv auto --filter='{"File Name": {"lacks": ["Group1"]}}'
#pcpfm convert_to_mzML ./Lucas_Deriv_Stds_G2_09_08_2023 $(which mono) ~/Projects/PythonCentricPipelineForMetabolomics-1/ThermoRawFileParser/ThermoRawFileParser.exe
#pcpfm asari_full_processing ./Lucas_Deriv_Stds_G2_09_08_2023 
#pcpfm feature_QCQA ./Lucas_Deriv_Stds_G2_09_08_2023 preferred --all --save_figs
pcpfm blank_masking ./Lucas_Deriv_Stds_G2_09_08_2023 preferred --new_table_moniker=preferred_blank_masked --blank_intensity_ratio=3 --blank_type="solblank"
pcpfm blank_masking ./Lucas_Deriv_Stds_G2_09_08_2023 preferred_blank_masked --new_table_moniker=preferred_blank_masked --blank_intensity_ratio=3 --blank_type="blank"
pcpfm drop_samples ./Lucas_Deriv_Stds_G2_09_08_2023 preferred_blank_masked --new_table_moniker=preferred_blank_masked --field_value="solblank"
pcpfm drop_samples ./Lucas_Deriv_Stds_G2_09_08_2023 preferred_blank_masked --new_table_moniker=preferred_blank_masked --field_value="blank"



#pcpfm assemble_experiment_from_CSV ./Lucas_UnDeriv_Stds_G1_09_08_2023 /Users/mitchjo/Analyses/DmPA_Lucas/Standards_NoDeriv/08282023_lucas_glycolysis_TCA_authenticstandards_underivatization/08282023_Lucas_glycolysis_TCA_authenticstandards_underivatization.csv auto --filter='{"File Name": {"lacks": ["Group2"]}}'
#pcpfm convert_to_mzML ./Lucas_UnDeriv_Stds_G1_09_08_2023 $(which mono) ~/Projects/PythonCentricPipelineForMetabolomics-1/ThermoRawFileParser/ThermoRawFileParser.exe
#pcpfm asari_full_processing ./Lucas_UnDeriv_Stds_G1_09_08_2023 
#pcpfm feature_QCQA ./Lucas_UnDeriv_Stds_G1_09_08_2023 preferred --all --save_figs
pcpfm blank_masking ./Lucas_UnDeriv_Stds_G1_09_08_2023 preferred --new_table_moniker=preferred_blank_masked --blank_intensity_ratio=3 --blank_type="solblank"
pcpfm blank_masking ./Lucas_UnDeriv_Stds_G1_09_08_2023 preferred_blank_masked --new_table_moniker=preferred_blank_masked --blank_intensity_ratio=3 --blank_type="blank"
pcpfm drop_samples ./Lucas_UnDeriv_Stds_G1_09_08_2023 preferred_blank_masked --new_table_moniker=preferred_blank_masked --field_value="solblank"
pcpfm drop_samples ./Lucas_UnDeriv_Stds_G1_09_08_2023 preferred_blank_masked --new_table_moniker=preferred_blank_masked --field_value="blank"

#pcpfm assemble_experiment_from_CSV ./Lucas_UnDeriv_Stds_G2_09_08_2023 /Users/mitchjo/Analyses/DmPA_Lucas/Standards_NoDeriv/08282023_lucas_glycolysis_TCA_authenticstandards_underivatization/08282023_Lucas_glycolysis_TCA_authenticstandards_underivatization.csv auto --filter='{"File Name": {"lacks": ["Group1"]}}'
#pcpfm convert_to_mzML ./Lucas_UnDeriv_Stds_G2_09_08_2023 $(which mono) ~/Projects/PythonCentricPipelineForMetabolomics-1/ThermoRawFileParser/ThermoRawFileParser.exe
#pcpfm asari_full_processing ./Lucas_UnDeriv_Stds_G2_09_08_2023 
#pcpfm feature_QCQA ./Lucas_UnDeriv_Stds_G2_09_08_2023 preferred --all --save_figs
pcpfm blank_masking ./Lucas_UnDeriv_Stds_G2_09_08_2023 preferred --new_table_moniker=preferred_blank_masked --blank_intensity_ratio=3 --blank_type="solblank"
pcpfm blank_masking ./Lucas_UnDeriv_Stds_G2_09_08_2023 preferred_blank_masked --new_table_moniker=preferred_blank_masked --blank_intensity_ratio=3 --blank_type="blank"
pcpfm drop_samples ./Lucas_UnDeriv_Stds_G2_09_08_2023 preferred_blank_masked --new_table_moniker=preferred_blank_masked --field_value="solblank"
pcpfm drop_samples ./Lucas_UnDeriv_Stds_G2_09_08_2023 preferred_blank_masked --new_table_moniker=preferred_blank_masked --field_value="blank"
