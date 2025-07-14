function [outputArg1] = deletebashwriteforcontinue(inputArg1,N_strip,savetime_file)
%inputArg1 = 'fluids'
% N_strip = 2;
fileFolder=fullfile([inputArg1,'/']); 
dirOutput=dir(fullfile(fileFolder));

filename = {dirOutput.name};
filename = filename(3:3+N_strip-1);

run_copy = [];
run_copy{1} = '#!/bin/bash';
run_copy{2} = 'cd fluids';

for i = 1:N_strip

name_temp = [inputArg1,'/',filename{i}];
fileFolder_temp = fullfile(name_temp);
dirOutput_temp = dir(fullfile(fileFolder_temp));
filename_temp = {dirOutput_temp.name};

run_copy{3+3*(i-1)} = ['cd',32,filename{i}];

delete_name = [];
k=1;
for j = 1:length(filename_temp)
    if str2num(filename_temp{j}) > savetime_file
    delete_name  =[delete_name,32,filename_temp{j}];
    k=k+1;
    end
end

temp_delete_bash = ['rm -rf',delete_name];

run_copy{4+3*(i-1)} = temp_delete_bash;
run_copy{5+3*(i-1)} = ['cd ..'];
end

fileID_new=fopen('delete_bash_forcontinue','w+');
[M,N]=size(run_copy);
for q=1:N
	fprintf(fileID_new,'%s\n',run_copy{q});
end
fclose(fileID_new);

outputArg1=1;


end