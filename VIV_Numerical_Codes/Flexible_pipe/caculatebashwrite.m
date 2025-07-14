function [outputArg1] = caculatebashwrite(inputArg1,N_strip)
%inputArg1 = 'fluids'
fileFolder=fullfile([inputArg1,'/']); 
dirOutput=dir(fullfile(fileFolder));

filename = {dirOutput.name};
filename = filename(3:3+N_strip-1);

run_copy = [];
run_copy{1} = '#!/bin/bash';
run_copy{2} = 'cd fluids';

for i = 1:N_strip-1
run_copy{3+3*(i-1)} = ['(cd',32,filename{i},';\'];
run_copy{4+3*(i-1)} = ['pimpleFoam >log',';\'];
run_copy{5+3*(i-1)} = ['cd ..',')&\'];
end

for i = N_strip
run_copy{3+3*(i-1)} = ['(cd',32,filename{i},';\'];
run_copy{4+3*(i-1)} = ['pimpleFoam >log',';\'];
run_copy{5+3*(i-1)} = ['cd ..',')'];
run_copy{6+3*(i-1)} = ['cd ..'];
run_copy{7+3*(i-1)} = ['wait'];
run_copy{8+3*(i-1)} = ['echo ''finish openfoam bash caculation'''];
end

fileID_new=fopen('runpimpleFoam_bash','w+');
[M,N]=size(run_copy);
for q=1:N
	fprintf(fileID_new,'%s\n',run_copy{q});
end
fclose(fileID_new);

outputArg1 = 1;

end