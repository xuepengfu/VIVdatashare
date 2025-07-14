function [outputArg1] = initialdis(inputArg1,N_strip)
%inputArg1 = 'fluids'

fileFolder=fullfile([inputArg1,'/']); 
dirOutput=dir(fullfile(fileFolder));

filename = {dirOutput.name};
filename = filename(3:3+N_strip-1);

for j = 1:length(filename)
    name_temp = [inputArg1,'/',filename{j}];
    fileFolder_temp = fullfile(name_temp);
    dirOutput_temp = dir(fullfile(fileFolder_temp));
    filename_temp = {dirOutput_temp.name};

    num_zero = find(strcmp(filename_temp,'0'));

    file_old = [];
    p = 1;
    fid=fopen([name_temp,'/',filename_temp{num_zero},'/pointDisplacement']);
    while(feof(fid) ~= 1)
        line = fgetl(fid);
        file_old{p} = line;
        p = p+1;
    end

    start = strfind(file_old,'CYLINDER');

    start_em = cellfun('isempty',start);
    num_cylinder = find(start_em==0);

    file_new = file_old;
    temp = ['        value           uniform (',num2str(0),32, ...
        num2str(0),32,'0);'];
    file_new{num_cylinder+3} = temp;

    fileID_new=fopen([name_temp,'/',filename_temp{num_zero},'/pointDisplacement'],'w+');
    [M,N]=size(file_new);
    for q=1:N
	fprintf(fileID_new,'%s\n',file_new{q});
    end
    fclose(fileID_new); 


end


outputArg1 = 1;

end