clear 
clc

%% FEM model parameter
nel=100;                % FEM model node number

d=0.028;                 %outter diameter
d0=0;                    %inner diameter
area=pi*(d^2-d0^2)/4;    %area
I=pi*(d^4-d0^4)/64;      %moment of inertia
Ip=pi*(d^4-d0^4)/32;     %polar moment of inertia 
J=Ip;                    
EI = 58.6;               %bending stiffness
em=EI/I;                 %stiffness module
miu=0.3;                 %Poisson's Ratio
G=em/2/(1+miu);          %Shear Modulus
L=3.88;                  %total length
l=L/nel;                 %unit length
mbar=1.24;               %mass per meter

dalpha=0.07;             %damping model coefficients
dbata=0.0001;            %damping model coefficients
fre_number=10;           % mode number
Tension=550;             %tension %change tension

steptime=0.0005;          %time step
totaltime=41;             %running time
savestep = 500;          % save per step
savetime = 500;

bata=0.25;               %newmarkbeta
gama=0.5;                %newmarkbeta
nnel=2;                 % node per unit
ndof=6;                 % node DOF
nnode=(nnel-1)*nel+1;    % total nodes number
sdof=nnode*ndof;         % total DOF number
rho=mbar/(pi*(d/2)^2);   %density

all_step = floor(totaltime/steptime+0.01)-1;
cord_FEM = linspace(0,L,nnode);
%% CFD model parameter
N_strip = 11; %number of strips
D = d; %outer diameter
rho_flu = 998;% density of fluid
N_sim = floor(totaltime/steptime+0.01)-1; %simulation step number
L_mesh = 0.032311; % the z length of mesh in the OpenFOAM (differ with case!!!)
vel_begin = 0;
vel_end = 0.8;
Vel_matrix = linspace(vel_begin,vel_end,N_strip); %define velocity matrix (uniform,linearly sheard flow and BSF)
Vel_matrix_FEM = linspace(vel_begin,vel_end,nnode); 

Cord_CFD = linspace(0,L,N_strip);

%% FEM initial
x=zeros(nnode,1);
y=zeros(nnode,1);
z=zeros(nnode,1);
for i=1:nnode
x(i,1)=(i-1)*L/nel;
y(i,1)=0;
z(i,1)=0;
end
massmatrixass;
K=zeros(nnode*ndof);
stiffnessmatrix;
U=zeros(nnode*ndof,1);
%pin-pin end
BC1=zeros(nnode+3+2,3);   
BC1(1:3,:)=[  1    1    0;
              1    2    0;
              1    3    0;]; 
for yn=1:nnode             
    BC1(yn+3,1)=yn;
    BC1(yn+3,2)=4;
    BC1(yn+3,3)=0;
end
BC1(nnode+3+1,:)=[nnode   2   0];    
BC1(nnode+3+2,:)=[nnode   3   0];
    
[bc1_number,dummy] = size( BC1 ) ;
 w2max = max( diag(K)./diag(M) ) ;
    for ibc=1:1:bc1_number
        n = BC1(ibc, 1 ) ;
        dd = BC1(ibc, 2 ) ;
        mm = (n-1)*ndof + dd ;
        K(:,mm) = zeros( nnode*ndof, 1 ) ;
        Kb(:,mm) = zeros( nnode*ndof, 1 );
        K(mm,:) = zeros( 1, nnode*ndof ) ;
        Kb(mm,:) = zeros( 1, nnode*ndof ) ;
        K(mm,mm) = 1;
        Kb(mm,mm) = 1;
        M(:,mm) = zeros( nnode*ndof, 1 ) ;
        M(mm,:) = zeros( 1, nnode*ndof ) ;
        M(mm,mm) = K(mm,mm)/w2max/1e10 ;
    end

for i=1:nnode*ndof
    for j=i:nnode*ndof
        K(j,i) = K(i,j) ;
        Kb(j,i) = Kb(i,j) ;
        M(j,i) = M(i,j) ;
    end
end
    

%[EigVector, EigValue] = eig(inv(M)*K);
[EigVector, EigValue] = eigs(K, M, fre_number, 'SM' );

for ibc=1:1:bc1_number
    n = BC1(ibc, 1 ) ;
    dd = BC1(ibc, 2 ) ;
    mm = (n-1)*ndof + dd ;
    EigVector(mm,:) = BC1(ibc,3) ;
end
    
frequency = zeros(fre_number,3);
for j=fre_number:-1:1
    %frequency( fre_number-j+1,:) = [EigValue(j,j), sqrt(EigValue(j,j)), sqrt(EigValue(j,j))/2/pi];
    frequency(j,:) = [EigValue(j,j), sqrt(EigValue(j,j)), sqrt(EigValue(j,j))/2/pi];        
end
frequency_y = zeros(fre_number/2,3);
for jjj=1:fre_number/2
    %frequency_y( fre_number/2-jjj+1,:) = [EigValue(2*jjj,2*jjj), sqrt(EigValue(2*jjj,2*jjj)), sqrt(EigValue(2*jjj,2*jjj))/2/pi];
    frequency_y(jjj,:) = [EigValue(2*jjj,2*jjj), sqrt(EigValue(2*jjj,2*jjj)), sqrt(EigValue(2*jjj,2*jjj))/2/pi];          
end

fid=fopen('frequency.txt','w');
fprintf( fid,'Eigenvalue List \n' ) ;
fprintf( fid,'   mode       Eigenvaluey        f_rad(Hz)     f_hz(Hz)  \n' ) ;
for i=1:fre_number
    fprintf( fid,'%6d   %15.7e   %15.7e   %15.7e\n', i, ...
        frequency(i,1), frequency(i,2), frequency(i,3) ) ;
end
fclose(fid);
fid=fopen('frequency_y.txt','w');
fprintf( fid,'Eigenvalue List in y direction \n' ) ;
fprintf( fid,'   mode       Eigenvaluey        fy_rad(Hz)     fy_hz(Hz)  \n' ) ;
for i=1:fre_number/2
    fprintf( fid,'%6d   %15.7e   %15.7e   %15.7e\n', i, ...
        frequency_y(i,1), frequency_y(i,2), frequency_y(i,3) ) ;
end
fclose(fid);


modalshape=zeros(nnode,fre_number);
for j=1:fre_number
modalshape(:,fre_number-j+1) = EigVector(2:ndof:nnode*ndof, j );
factor = max( abs(modalshape(:,fre_number-j+1)) );
modalshape(:,fre_number-j+1) = modalshape(:,fre_number-j+1) / factor;
end

dampingmatrix;

F=zeros(ndof*nnode,totaltime/steptime+1);

dis=zeros(ndof*nnode,totaltime/steptime+1);
vel=zeros(ndof*nnode,totaltime/steptime+1);
acc=zeros(ndof*nnode,totaltime/steptime+1);

a0=1/(bata*steptime^2);
a1=gama/(bata*steptime);
a2=1/(bata*steptime);
a3=1/(2*bata)-1;
a4=gama/(bata)-1;
a5=steptime/2*(gama/bata-2);
a6=steptime*(1-gama);
a7=steptime*gama;

%% copy the strips 
copystrips(N_strip); %copy the strips

%% initial the displace ment
initialdis('fluids',N_strip);
%% initial the velocity
initialdvelxocity('fluids',N_strip,Vel_matrix);
%% begin caculate; transfer fluids force; caculate dis; transfer to openfoam
%% write cacalulation bash
[outputArg1] = caculatebashwrite('fluids',N_strip); %output the bash for running of
%% caculate
CL_all_rare = []; %rate CL data from of
CD_all_rare = [];
CL_all_true = []; %true CL data from of
CD_all_true = [];
CL_all_FEM = []; % CL of FEM model
CD_all_FEM = [];
F_CL_FEM = [];% CL of FEM model 
F_CD_FEM = [];

outputtime=0.5; % save the data per outputtime

parpool('local',20,'IdleTimeout',Inf)

%Nsim floor(totaltime/steptime+0.01)-1
for i = 1 : floor(totaltime/steptime+0.01)-1
%display i
t00 = tic;
disp(['current i is ',num2str(i)])
% begin caculate
system('source runpimpleFoam_bash'); %run OF in every strip for one strip
pause(0.5);
disp(['finish openfoam caculation'])
% read  rare cl/cd
runningstep = i;
rightnowtime = steptime*runningstep;

disp(['current timestep is ',num2str(rightnowtime)])


countt0 = toc(t00);
t11 = tic;

[CL_all_temp,CD_all_temp] = forceread('fluids',runningstep,N_strip,rightnowtime,steptime);%rare data from of
CL_all_rare(i,:) = CL_all_temp;
CD_all_rare(i,:) = CD_all_temp;

% process for the CL CD true and input F
CL_all_true(i,:) = CL_all_temp./(D*L_mesh*Vel_matrix.^2);
CD_all_true(i,:) = CD_all_temp./(D*L_mesh*Vel_matrix.^2);

CL_all_true(isnan(CL_all_true))=0;
CL_all_true(isinf(CL_all_true))=0;

CD_all_true(isnan(CD_all_true))=0;
CD_all_true(isinf(CD_all_true))=0;

CL_all_FEM(i,:) = spline(Cord_CFD,CL_all_true(i,:),cord_FEM);
CD_all_FEM(i,:) = spline(Cord_CFD,CD_all_true(i,:),cord_FEM);

F_CL_FEM(i,:) = 0.5 * rho_flu * D * l * Vel_matrix_FEM.^2.*CL_all_FEM(i,:);
F_CL_FEM(i,1) = F_CL_FEM(i,1) / 2;
F_CL_FEM(i,nnode) = F_CL_FEM(i,nnode) / 2;

F_CD_FEM(i,:) = 0.5 * rho_flu * D * l * Vel_matrix_FEM.^2.*CD_all_FEM(i,:);
F_CD_FEM(i,1) = F_CD_FEM(i,1) / 2;
F_CD_FEM(i,nnode) = F_CD_FEM(i,nnode) / 2;
 
for  j=2:1:nnode-1   
 F(ndof*(j-1)+2,runningstep+1) = F_CL_FEM(i,j);   %Y
 F(ndof*(j-1)+3,runningstep+1) = F_CD_FEM(i,j);   %x
end  
% caculate the displacement of strips based FEM
Q(:,runningstep+1)=F(:,runningstep+1)+M*(a0*dis(:,runningstep)...
    +a2*vel(:,runningstep)+a3*acc(:,runningstep))+C*(a1*dis(:,runningstep)...
    +a4*vel(:,runningstep)+a5*acc(:,runningstep));

K1=a0*M+a1*C+K;
dis(:,runningstep+1)=K1\Q(:,runningstep+1); % 2:CF. 3:IL
acc(:,runningstep+1)=a0*(dis(:,runningstep+1)-dis(:,runningstep))-a2*vel(:,runningstep)-a3*acc(:,runningstep);
vel(:,runningstep+1)=vel(:,runningstep)+(1-gama)*acc(:,runningstep)*steptime...
    +gama*acc(:,runningstep+1)*steptime;

tempdis = [];
for  j=1:1:nnode    
 tempdis(j,1) = dis(ndof*(j-1)+3,runningstep+1);   %x il
 tempdis(j,2) = dis(ndof*(j-1)+2,runningstep+1);   %y cf
end  

tempdis_CFD = [];
tempdis_CFD (:,1) = spline(cord_FEM,tempdis(:,1),Cord_CFD);
tempdis_CFD (:,2) = spline(cord_FEM,tempdis(:,2),Cord_CFD);

countt11 = toc(t11);
t22 = tic;

% write the displacement for of
[outputArg1] = transferdis('fluids',runningstep,tempdis_CFD,N_strip,rightnowtime,steptime);
disp(['finish displacement transfer'])

countt22 = toc(t22);
t33 = tic;

if i >= 3 && mod(steptime*i,outputtime)~=0 
 disp(['begin deleting timestep file'])
 [outputArg1] = deletebashwrite('fluids',N_strip,savestep,all_step,rightnowtime,steptime,outputtime);
if outputArg1 == 0
disp(['No delete in this time step'])
else if outputArg1 == 1
system('source delete_bash','-echo'); %run OF in every strip for one strip
end
end

 [outputArg1] = deleteforcebashwrite('fluids',N_strip,savestep,all_step,rightnowtime,steptime,outputtime);
if outputArg1 == 1
system('source delete_force_bash','-echo'); %run OF in every strip for one strip
disp(['finish deleting timestep file'])
end
end

if i >= 3 && mod(steptime*i,5)==0 
save('forcedata.mat','CL_all_rare','CD_all_rare','CL_all_true',...
    'CD_all_true','CL_all_FEM','CD_all_FEM','F_CL_FEM',...
    'F_CD_FEM','-v7.3')

save('displacement.mat','dis','vel','acc','-v7.3')
end

save('time.mat','rightnowtime','steptime','-v7.3');

ST = fclose('all');
countt33 = toc(t33);

fprintf(['Openfoam cacu takes', num2str(countt0), ' seconds.\n'])
fprintf(['DisplacementFEM cacu takes', num2str(countt11), ' seconds.\n'])
fprintf(['Transfer dis', num2str(countt22), ' seconds.\n'])
fprintf(['delete takes', num2str(countt33), ' seconds.\n\n'])
end

%    xlswrite('displacement.xlsx',dis);
%    xlswrite('velocity.xlsx',vel);
%    xlswrite('acceleration.xlsx',acc);
%    save('displacement.mat','dis','vel','acc','-v7.3')
