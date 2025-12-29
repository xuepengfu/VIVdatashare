clear;
clc;

fileFolder = fullfile('I:\YOUR_DATA_FOLDER\water_rotation');
dirOutput  = dir(fullfile(fileFolder,'*.csv'));
fileNames  = {dirOutput.name}';
filenumber = length(fileNames);

L  = 7.64;
D  = 0.02841;
D1 = 0.025;
R  = D/2;
EA = 9.4E5;
EA1 = 800*10E6*pi*(0.025^2-(0.025-0.00168*2)^2);
EI = 58.6;

Fs = 250;
dt = 1/250;

flow  = 0.01;
fhigh = 20;

CFFBG = [121 186 251.25 316.35 381.45 446.45 511.45 576.7 641.7];
CFFBG = CFFBG/100;

ILFBG = [88.5 133.6 178.7 224.05 269.2 314.35 359.5 404.65 449.55 ...
         495.05 540.3 585 630.5 675.4];
ILFBG = ILFBG/100;

for i = 62

    filename = fileNames{i};
    data = [];
    name = [];
    [data,name] = xlsread(fullfile(fileFolder,filename));

    name(1,:) = [];
    name(10) = {'encoder'};

    [~,n] = size(name);
    for j = 11:n
        name(j) = cellstr(name{j}(end-5:end));
    end

    [m1,~] = size(data);
    time = (0:dt:dt*(m1-1))';
    disp([filename,' start']);

    for j = 1:4
        chname1 = 'CF1_4';
        value1  = data(:,find(strcmp([chname1 num2str(j)],name)));
        value1  = value1 - mean(value1);

        chname2 = 'CF2_4';
        value2  = data(:,find(strcmp([chname2 num2str(j)],name)));
        value2  = value2 - mean(value2);

        value3 = (value1-value2)/2;
        eval(['VIVCF',num2str(j),'=','value3',';']);

        value4 = (value1+value2)/2;
        eval(['TENCF',num2str(j),'=','value4',';']);
    end

    for j = 1:5
        k = 4+j;

        chname1 = 'CF1_5';
        value1  = data(:,find(strcmp([chname1 num2str(j)],name)));
        value1  = value1 - mean(value1);

        chname2 = 'CF2_5';
        value2  = data(:,find(strcmp([chname2 num2str(j)],name)));
        value2  = value2 - mean(value2);

        value3 = (value1-value2)/2;
        eval(['VIVCF',num2str(k),'=','value3',';']);

        value4 = (value1+value2)/2;
        eval(['TENCF',num2str(k),'=','value4',';']);
    end

    for j = 1:6
        chname1 = 'IL1_6';
        value1  = data(:,strcmp([chname1 num2str(j)],name));
        value1  = value1 - mean(value1(1:1000));

        chname2 = 'IL2_6';
        value2  = data(:,find(strcmp([chname2 num2str(j)],name)));
        value2  = value2 - mean(value2(1:1000));

        value3 = (value1-value2)/2;
        eval(['VIVIL',num2str(j),'=','value3',';']);

        value4 = (value1+value2)/2;
        eval(['TENIL',num2str(j),'=','value4',';']);
    end

    for j = 1:8
        k = 6+j;

        chname1 = 'IL1_8';
        value1  = data(:,find(strcmp([chname1 num2str(j)],name)));
        value1  = value1 - mean(value1(1:1000));

        chname2 = 'IL2_8';
        value2  = data(:,find(strcmp([chname2 num2str(j)],name)));
        value2  = value2 - mean(value2(1:1000));

        value3 = (value1-value2)/2;
        eval(['VIVIL',num2str(k),'=','value3',';']);

        value4 = (value1+value2)/2;
        eval(['TENIL',num2str(k),'=','value4',';']);
    end

    inIL = [];
    figure
    plot(VIVIL5,'b')
    [range1,~] = ginput;
    [range2,~] = ginput;
    range = floor(range1):floor(range2);

    for j = 1:14
        eval(['inIL(j)=mean(VIVIL',num2str(j),'(range));']);
    end
    inIL = inIL*(D/D1);

    CFStrain = [];
    for j = 1:9
        eval(['CFStrain=[CFStrain; VIVCF',num2str(j),'''];']);
    end

    ILStrain = [];
    for j = 1:14
        eval(['ILStrain=[ILStrain; VIVIL',num2str(j),'''];']);
    end

    CFStrainmax = max(abs(CFStrain'));
    ILStrainmax = max(abs(ILStrain'));

    ppp = 1;
    huaidianCF = [];
    huaidianIL = [];

    for j = 1:length(CFStrainmax)
        if (CFStrainmax(j)>2000)
            disp([filename,' CF spike'])
            huaidianCF(ppp)=j;
            ppp=ppp+1;
        elseif (j==length(CFStrainmax))
            disp([filename,' CF no spike'])
        end
    end

    ppp = 1;
    for j = 1:length(ILStrainmax)
        if (ILStrainmax(j)>2000)
            disp([filename,' IL spike'])
            huaidianIL(ppp)=j;
            ppp=ppp+1;
        elseif (j==length(ILStrainmax))
            disp([filename,' IL no spike'])
        end
    end

    if ~isempty(huaidianCF)
        disp('fix CF spikes')
        for p = 1:length(huaidianCF)
            eval(['huaidiandata=VIVCF',num2str(huaidianCF(p)),';'])
            huaidiandata = huaidiandata - mean(huaidiandata(1:100));
            for q = 1:length(huaidiandata)
                if abs(huaidiandata(q))>2000
                    huaidiandata(q) = huaidiandata(q-1) + huaidiandata(q-1) - huaidiandata(q-2);
                end
            end
            eval(['VIVCF',num2str(huaidianCF(p)),'=huaidiandata;'])
        end
        disp('CF fixed')
    end

    if ~isempty(huaidianIL)
        disp('fix IL spikes')
        for p = 1:length(huaidianIL)
            eval(['huaidiandata=VIVIL',num2str(huaidianIL(p)),';'])
            huaidiandata = huaidiandata - mean(huaidiandata(1:100));
            for q = 1:length(huaidiandata)
                if abs(huaidiandata(q))>2000
                    huaidiandata(q) = huaidiandata(q-1) + huaidiandata(q-1) - huaidiandata(q-2);
                end
            end
            eval(['VIVIL',num2str(huaidianIL(p)),'=huaidiandata;'])
        end
        disp('IL fixed')
    end

    for j = 1:9
        eval(['VIVCF',num2str(j),'=bpass(VIVCF',num2str(j),',dt,flow,fhigh);']);
    end

    for j = 1:14
        eval(['VIVIL',num2str(j),'=bpass(VIVIL',num2str(j),',dt,flow,fhigh);']);
    end

    CFStrain = [];
    for j = 1:9
        eval(['CFStrain=[CFStrain; VIVCF',num2str(j),'''];']);
    end
    CFStrain = CFStrain/(10^6);
    CFStrain = CFStrain*(D/D1);

    ILStrain = [];
    for j = 1:14
        eval(['ILStrain=[ILStrain; VIVIL',num2str(j),'''];']);
    end
    ILStrain = ILStrain/(10^6);
    ILStrain = ILStrain*(D/D1);

    CFTransModeShape = [];
    for j = 1:9
        for k = 1:8
            CFTransModeShape(j,k) = -R*(k*pi/L)^2*sin(k*pi*CFFBG(j)/L);
        end
    end

    ILTransModeShape = [];
    for j = 1:14
        for k = 1:13
            ILTransModeShape(j,k) = -R*(k*pi/L)^2*sin(k*pi*ILFBG(j)/L);
        end
    end

    CF_DisModeWeight = inv(CFTransModeShape'*CFTransModeShape)*CFTransModeShape'*CFStrain;
    IL_DisModeWeight = inv(ILTransModeShape'*ILTransModeShape)*ILTransModeShape'*ILStrain;

    ell = 0:L/100:L;

    ModeShapeForCFDis = zeros(length(ell),8);
    for j = 1:length(ell)
        for k = 1:8
            ModeShapeForCFDis(j,k) = sin(k*pi/L*ell(j));
        end
    end

    ModeShapeForILDis = zeros(length(ell),13);
    for j = 1:length(ell)
        for k = 1:13
            ModeShapeForILDis(j,k) = sin(k*pi/L*ell(j));
        end
    end

    CFDIS = ModeShapeForCFDis*CF_DisModeWeight;
    ILDIS = ModeShapeForILDis*IL_DisModeWeight;

    strain_by_diff = diff(CFDIS,2,1)/((L/100)^2)*R*10^6;

    inIL = inIL/(1E6);
    IL_IniDisModeWeight = inv(ILTransModeShape'*ILTransModeShape)*ILTransModeShape'*inIL';
    ILIniDIS  = ModeShapeForILDis*IL_IniDisModeWeight;
    ILIniDIS  = ILIniDIS/D;
    ILIniDISm = ILIniDIS*D;

    Fz = data(:,find(strcmp('TF2_Fz',name)));
    Fz = Fz/54.94505495*56.17977528;
    Fz = Fz*9.8;
    Fz_filter = bpass(Fz,dt,flow,fhigh);

    Fx = data(:,find(strcmp('TF2_Fx',name)));
    Fx = Fx*9.8;

    Fy = data(:,find(strcmp('TF2_Fy',name)));
    Fy = Fy*9.8;

    Fzini = mean(Fz(1:1000));

    Fz_strain = (TENCF1+TENCF2+TENCF3+TENCF4+TENCF5+TENCF6+TENCF7+TENCF8+TENCF9+ ...
                 TENIL1+TENIL2+TENIL3+TENIL4+TENIL5+TENIL6+TENIL7+TENIL8+TENIL9+TENIL10+ ...
                 TENIL11+TENIL12+TENIL13+TENIL14)/(9+14)*(10^-6)*EA1;

    Fz_strainfilter = bpass(Fz_strain,dt,flow,fhigh);

    disp([filename,' done']);
    close all
end

