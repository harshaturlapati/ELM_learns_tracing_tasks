clc;clear;close all;


pathname = [pwd '/datafiles/'];
file_list = dir([pathname '*.fig']);
length(file_list)

%% Figs 1, 2 and 3
for j_files = 1 : length(file_list)
    
    clearvars -except file_list j_files pathname
        
    filename = file_list(j_files).name;
    
    curr_fig = openfig( strcat(pathname,filename) );
    
    subject_name = strsplit(filename,'.fig'); subject_name = subject_name(1);
    
    jpeg_filename = strcat(pathname,char(subject_name),'.jpg');

    saveas(curr_fig, jpeg_filename);
end

clc;clear;close all;