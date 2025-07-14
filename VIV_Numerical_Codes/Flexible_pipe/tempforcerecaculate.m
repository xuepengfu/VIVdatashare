function [outputArg1] = tempforcerecaculate(inputArg1,runningstep,N_strip,rightnowtime,steptime,stripnum)
%inputArg1 = 'fluids'
%
fileFolder=fullfile([inputArg1,'/']); 
dirOutput=dir(fullfile(fileFolder));

filename = {dirOutput.name};
filename = filename(3:3+N_strip-1);

run_delete = [];
run_delete{1} = '#!/bin/bash';
run_delete{2} = 'cd fluids';


name_temp = [inputArg1,'/',filename{stripnum}];
fileFolder_temp = fullfile(name_temp);
dirOutput_temp = dir(fullfile(fileFolder_temp));
filename_temp = {dirOutput_temp.name};

run_delete{3} = ['cd',32,filename{stripnum}];

num_zero_file = find(strcmp(filename_temp,num2str(rightnowtime)));

temp_delete_bash = [];
temp_delete_bash = filename_temp{num_zero_file};
temp_delete_bash = ['rm -rf',32,temp_delete_bash];

run_delete{4} = temp_delete_bash;
run_delete{5} = ['wait'];
run_delete{6} = ['cd postProcessing/forceCoeffs'];


name_temp_force = [inputArg1,'/',filename{stripnum},'/postProcessing/forceCoeffs'];
fileFolder_temp = fullfile(name_temp_force);
dirOutput_temp = dir(fullfile(fileFolder_temp));
filename_temp = {dirOutput_temp.name};

num_zero_force = find(strcmp(filename_temp,num2str(rightnowtime-steptime)));

temp_delete_bash = [];
temp_delete_bash = filename_temp{num_zero_force};
temp_delete_bash = ['rm -rf',32,temp_delete_bash];

run_delete{7} = temp_delete_bash;
run_delete{8} = ['wait'];
run_delete{9} = ['cd ..'];
run_delete{10} = ['cd ..'];


fileID_new=fopen('forcerecaculate_bash','w+');
[M,N]=size(run_delete);
for q=1:N
	fprintf(fileID_new,'%s\n',run_delete{q});
end
fclose(fileID_new);

run_delete{11} = ['pimpleFoam'];

fileID_new=fopen('forcerecaculate_bashv1','w+');
[M,N]=size(run_delete);
for q=1:N
	fprintf(fileID_new,'%s\n',run_delete{q});
end
fclose(fileID_new);

outputArg1 = 1;

end