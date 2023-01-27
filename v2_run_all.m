clc; clear; close all;

import matlab.internal.liveeditor.LiveEditorUtilities

run('set_training_settings.m')
run('set_generalization_settings.m')
run_save_export_MLX_v1('define_curves_for_generalization.mlx')

%% Training
run_save_export_MLX_v1('v2_train_versionLyudmila.mlx')

%% Generalizing to Fig8

run_save_export_MLX_v1('v2_generalize.mlx')
run('v3_compare_ppl.m')

%% Generalizing to new figures

run('v2_generalize_new_TRAJ.mlx')