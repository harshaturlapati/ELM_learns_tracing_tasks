
clc; clear; close all;

ppl_pathname = [pwd '/datafiles/'];
ppl_file_list = dir([ppl_pathname '*.ELM_self_looping.mat'])

 real_K = []; %- 1xn
 real_V = []; 
 real_K_CI = []; 
 real_V_CI = []; 
 ELM_K = []; 
 ELM_V = []; 
 ELM_K_CI = []; 
 ELM_V_CI = []; 
  lambda_ridge_collection = []; 
  length(ppl_file_list)

%% Figs 1, 2 and 3
for ppl_j_files = 1 : length(ppl_file_list)
 
%% Beginning of main for loop

   
    clearvars -except ppl_file_list ppl_j_files ppl_pathname curvature_threshold real_K real_V real_K_CI real_V_CI ELM_K ELM_V ELM_K_CI ELM_V_CI lambda_ridge_collection elm_pathname
        
    ppl_filename = ppl_file_list(ppl_j_files).name
    
    subject = load( strcat(ppl_pathname,ppl_filename) );
    
    subject_name = strsplit(ppl_filename,'.ELM_self_looping.mat'); subject_name = subject_name(1);
    training_filename = strcat(subject_name,'.trained_elm.mat');
    training_data = load(char(strcat(ppl_pathname,training_filename)));

        real_K = [real_K subject.VK_real(2)]; %- 1xn
        real_V = [real_V subject.VK_real(1)]; 
        real_K_CI = [real_K_CI subject.VK_CI_real(2)]; 
        real_V_CI = [real_V_CI subject.VK_CI_real(1)]; 
        ELM_K = [ELM_K subject.VK_ELM(2)]; 
        ELM_V = [ELM_V subject.VK_ELM(1)]; 
        ELM_K_CI = [ELM_K_CI subject.VK_CI_ELM(2)]; 
        ELM_V_CI = [ELM_V_CI subject.VK_CI_ELM(1)];
        lambda_ridge_collection = [lambda_ridge_collection training_data.ridgeLambda];

end

total_subjects = length(real_K);
h_ppl = errorbar_plot_v2(total_subjects,ppl_pathname,real_K,real_V,real_K_CI,real_V_CI,ELM_K,ELM_V,ELM_K_CI,ELM_V_CI);


%% Exporting the table for real ppl and elm ppl

clear input;
fprintf('\n\nExample 3: using an array as data input\n\n');

% numeric values you want to tabulate:
input.data = [real_V' real_K' ELM_V' ELM_K'];

% Optional fields:

% Set column labels (use empty string for no label):
input.tableColLabels = {'$v^{REAL}$ (m/s)','$\kappa^{REAL} (1/m)$','$v^{ELM} (m/s)$','$\kappa^{ELM} (1/m)$'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'S1','S10','S2','S3','S4','S5','S6','S7','S8','S9'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used

%
no_cols = size(input.data,2);
%input.dataFormat = {'%.3f',2,'%.1f',1}; % three digits precision for first two columns, one digit for the last
input.dataFormat = {'%.3f',no_cols}; % three digits precision for all columns

% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';

% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';

% Switch table borders on/off (borders are enabled by default):
input.tableBorders = 1;

% LaTex table caption:
input.tableCaption = 'Piece-wise power-law quantities evaluated from tracing data of seven subjects.';

% LaTex table label:
input.tableLabel = 'tab:power_law_real_ELM';

% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 1;

% call latexTable:
latex = latexTable(input);
writecell(latex,'./datafiles/ELM_v_REAL_PPL_comparison_LATEX_TAB.txt')

exportgraphics(h_ppl,strcat('./datafiles/', 'all_subjects.power_law_comparison.fig8','.pdf'),'BackgroundColor','none','ContentType','vector');
