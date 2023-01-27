# ELM_learns_tracing_tasks
Tracing curves in the plane: geometric-invariant learning from human demonstrations 

# Abstract
The empirical laws governing human-curvilinear movements have been studied using various relationships, including minimum jerk, the 2/3 power law, and the piecewise power law. These laws quantify the speed-curvature relationships of human movements during curve tracing using critical speed and curvature as regressors. In this work, we provide a reservoir computing-based framework that can learn and reproduce human-like movements. Specifically, the geometric invariance of the observations, i.e., lateral distance from the closest point on the curve, instantaneous velocity, and curvature, when viewed from the moving frame of reference, are exploited to train the reservoir system. The artificially produced movements are evaluated using the power law to assess whether they are indistinguishable from their human counterparts. The generalisation capabilities of the trained reservoir to curves that have not been used during training are also shown.

# Order of running the files in this repository

- set_training_settings.m
- set_generalization_settings.m
- define_curves_for_generalization.mlx
- v2_train_versionLyudmila.mlx
- v2_generalize.mlx
- v3_compare_ppl.m
- v2_generalize_new_TRAJ.mlx
