function [outputArg1] = deletebashwrite(inputArg1,N_strip,savetime,all_step,rightnowtime,steptime,outputtime)
%inputArg1 = 'fluids'
% N_strip = 2;
fileFolder=fullfile([inputArg1,'/']); 
dirOutput=dir(fullfile(fileFolder));

filename = {dirOutput.name};
filename = filename(3:3+N_strip-1);

run_copy = [];
run_copy{1} = '#!/bin/bash';
run_copy{2} = 'cd fluids';
pause(10);
for i = 1:N_strip

name_temp = [inputArg1,'/',filename{i}];
fileFolder_temp = fullfile(name_temp);
dirOutput_temp = dir(fullfile(fileFolder_temp));
filename_temp = {dirOutput_temp.name};

run_copy{3+3*(i-1)} = ['cd',32,filename{i}];


num_zero = find(strcmp(filename_temp,num2str(rightnowtime)));

deletenum = rightnowtime - steptime*2;
num_delete = find(strcmp(filename_temp,num2str(deletenum)));

temp_delete_bash = [];
% for j = 1 : all_step
%     if mod(j,savetime)  == 0
%     else
%     temp_delete_bash =[temp_delete_bash,filename_temp{num_zero+j},32]
%     end
% end
temp_delete_bash = filename_temp{num_delete};
temp_delete_bash = ['rm -rf',32,temp_delete_bash];

run_copy{4+3*(i-1)} = temp_delete_bash;
run_copy{5+3*(i-1)} = ['cd ..'];
end

fileID_new=fopen('delete_bash','w+');
[M,N]=size(run_copy);
for q=1:N
	fprintf(fileID_new,'%s\n',run_copy{q});
end
fclose(fileID_new);

if mod(deletenum,outputtime) == 0 
outputArg1 = 0;
else
outputArg1 = 1;
end



end