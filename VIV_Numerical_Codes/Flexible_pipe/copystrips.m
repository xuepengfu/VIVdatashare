function [outputArg1] = copystrips(N_strip)

run_copy = [];
run_copy{1} = '#!/bin/bash';
run_copy{2} = 'cd fluids';

for i = 2:N_strip
    if i < 10
    run_copy{i+1} = ['cp -r strip01 strip0',num2str(i),''];
    else
    run_copy{i+1} = ['cp -r strip01 strip',num2str(i),''];
    end
end

fileID_new=fopen('copy_bash','w+');
[M,N]=size(run_copy);
for q=1:N
	fprintf(fileID_new,'%s\n',run_copy{q});
end
fclose(fileID_new);


system('source copy_bash','-echo');

end