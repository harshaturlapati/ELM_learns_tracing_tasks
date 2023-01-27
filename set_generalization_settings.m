function set_generalization_settings
clc; clear; close all;
curvature_threshold = 1;                %% select minorant for curvature
generalization_duration = 50; %s - number of seconds to run the ELM for
beta = 0.05;
cutoff_frequency = nan;   %% Hz
save('generalization_settings.mat');
end