function [CL_all,CD_all] = forceread(inputArg1,runningstep,N_strip,rightnowtime,steptime)
%inputArg1 = 'fluids'
%runningstep = 2

fileFolder=fullfile([inputArg1,'/']); 
dirOutput=dir(fullfile(fileFolder));

filename = {dirOutput.name};
filename = filename(3:3+N_strip-1);

CL_all = [];
CD_all = [];
for j = 1:length(filename)
    stripnum = j;
    
    name_temp = [inputArg1,'/',filename{j},'/postProcessing/forceCoeffs'];
    fileFolder_temp = fullfile(name_temp);
    dirOutput_temp = dir(fullfile(fileFolder_temp));
    filename_temp = {dirOutput_temp.name};

    num_zero = find(strcmp(filename_temp,num2str(rightnowtime-steptime)));

    fid=fopen([inputArg1,'/',filename{j},'/postProcessing/forceCoeffs/'...
        ,filename_temp{num_zero},'/forceCoeffs.dat']);

    ForceallArray = textscan(fid,'%f %f %f %f %f %f', ...
        'Delimiter','\t', 'HeaderLines',9, 'ReturnOnError', false);

    temle = length(ForceallArray{4});

    if temle == 1
    disp(['exsiting error here, do recaculation']);

    [outputArg1] = tempforcerecaculate(inputArg1,runningstep,N_strip,rightnowtime,steptime,stripnum)
    system('source forcerecaculate_bash','-echo') %run OF in every strip for one strip

    pause(5);
    system('source forcerecaculate_bashv1','-echo') %run OF in every strip for one strip

    disp(['reopen ',inputArg1,'/',filename{j},'/postProcessing/forceCoeffs/'...
        ,filename_temp{num_zero},'/forceCoeffs.dat now']);

    fid=fopen([inputArg1,'/',filename{j},'/postProcessing/forceCoeffs/'...
        ,filename_temp{num_zero},'/forceCoeffs.dat']);

    ForceallArray = textscan(fid,'%f %f %f %f %f %f', ...
        'Delimiter','\t', 'HeaderLines',9, 'ReturnOnError', false);

	CL_temp=ForceallArray{4}(2);
	CD_temp=ForceallArray{3}(2);

    else
	CL_temp=ForceallArray{4}(2);
	CD_temp=ForceallArray{3}(2);
    end


    CL_all(j) = CL_temp;
    CD_all(j) = CD_temp;

end