%% find the close time of breaktime
dis_breakall = load("displacement.mat");
dis_break = dis_breakall.dis;
[mm1,nn1] = find(dis_break(5,:)~=0,1,'last');%remove the additional zero

breaktime = steptime*(nn1-1);
savetime_file = steptime*(nn1-2);

%%savetime is the savetime_file-1 due to some problem
savetime_file = steptime*(nn1-3);

%% delete file beyond the time
[outputArg1] = deletebashwriteforcontinue('fluids',N_strip,savetime_file);
system('source delete_bash_forcontinue'); %run OF in every strip for one strip

[outputArg1] = deleteforcebashwriteforcontinue('fluids',N_strip,savetime_file);
system('source delete_force_bash_forcontinue'); %run OF in every strip for one strip

%% load the data from breaking point
dis_breakall = load("displacement.mat");

dis = dis_breakall.dis;
acc = dis_breakall.acc;
vel = dis_breakall.vel;

load("forcedata.mat");

breaki = nn1-2;

disp(['finish break point process'])




