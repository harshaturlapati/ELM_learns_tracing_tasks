function run_save_export_MLX_v1(filename)
filePath = fullfile(pwd);
% set the options according to the use case
options = struct('format','pdf','outputDir',filePath);
fileNames = dir(filePath);


% specify path to the live scripts
for j = 1:size(fileNames,1)
    if ~endsWith(fileNames(j).name,filename)
       continue;
    end
    matlab.internal.liveeditor.executeAndSave(fullfile(fileNames(j).folder,fileNames(j).name));
end

pdf_filename = strsplit(filename,'.mlx');
pdf_filename = char(strcat(pdf_filename(1),'.pdf'));

matlab.internal.liveeditor.openAndConvert (filename, pdf_filename);

end