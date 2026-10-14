sbatch auc_plots_2.sh adeno-controls_crc/ Controls-Adeno CRC 207 104
sbatch auc_plots_2.sh adeno_crc/ Adenoma CRC 105 104
sbatch auc_plots_2.sh controls_adeno/ Controls Adenoma 102 105
sbatch auc_plots_2.sh crc-adeno_controls/ Controls CRC-Adenoma 102 209


###
sbatch shap_plots_2.sh Metadata/config_file_crc_controls.yaml crc_controls/ nested_random-forest_crc_controls_*joblib  CORSA_nestedCV figures/
sbatch shap_plots_2.sh Metadata/config_file_adeno-controls_crc.yaml adeno-controls_crc/ nested_random-forest*joblib  CORSA_nestedCV figures/
sbatch shap_plots_2.sh Metadata/config_file_adeno_crc.yaml adeno_crc/ nested_random-forest*joblib  CORSA_nestedCV figures/
sbatch shap_plots_2.sh Metadata/config_file_controls_adeno.yaml controls_adeno/ nested_random-forest*joblib  CORSA_nestedCV figures/
sbatch shap_plots_2.sh Metadata/config_file_crc-adeno_controls.yaml crc-adeno_controls/ nested_random-forest*joblib  CORSA_nestedCV figures/
