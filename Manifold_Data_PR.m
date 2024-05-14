clear all
close all

%Change Session Name
% databasemaker2_2M('394.txt')
% EgoGEN('394_db.txt','00394_db.xls')
% databasemaker2_2M('503.txt')
% EgoGEN('503_db.txt','00503_db.xls')
% databasemaker2_2M('504.txt')
% EgoGEN('504_db.txt','00504_db.xls')
addpath D:\GPS\cal\pippin-master\MEC_Final\xcorr


% HighDimData_Name='HighDimData_Nosepoke_Rcount_1';

% HighDimData_Name={'Nosepoke_Rcount_1','Nosepoke_Rcount_2','Nosepoke_Rcount_3','Trayon_pump','Trayon_nonpump'};
% HighDimData_Name={'Nose_Poke_1'};
HighDimData_Name={'Nosepoke_PRstate_4'};
for i =1:length(HighDimData_Name)

    clear all；

ManifoldData={};
ManifoldDataNumber=1;
BinWidth = 50; %bin 50ms
FR_Fixed_Count=3;
binwidth=0.1;
WindowWidth=200;
Li=-2; % the event interval between -4s to 4s
Ui=2;
state_start=2; % PR state
state_end=8;
max_state=10;
% data_save_name='HighDimData_PR_np2';  % 2 9 200
data_save_name='HighDimData_PR_np_2_8_200';  % 2 9 200
data_save_name='HighDimData_PR_np_2_8_400';  % 2 9 200
data_save_name='HighDimData_PR_np_2_8_2_2_400';  % 2 9 200
data_save_name='HighDimData_PR_np_2_8_2_2_800';  % 2 9 200
data_save_name='HighDimData_PR_np_2_8_2_2_200';  % 2 9 200
fid = fopen('warning.txt', 'a+');



%% Data Load

ExcelPath='E:\GPS\cal\Motivation\Motivation_gps\Code\ABET\';
DataPath='E:\GPS\cal\Motivation\Motivation_gps\Data\';

addpath(ExcelPath)
addpath(DataPath)

ExceLgg1 = ls([ExcelPath '*csv']);
ExceLhh1 = cellstr(ExceLgg1);
ExceLnumber1 = length(ExceLhh1);
ExceLff1 = ExceLhh1;

for ii=1:ExceLnumber1  %%%session number
    fnumIRT =[ExceLff1{ii}];

    if ~isempty(findstr(fnumIRT,'IRT')), continue, end


    Egoggg = ls([DataPath fnumIRT(1:end-8) '*mat']);
    Egohhh = cellstr(Egoggg);
    Egonumber = length(Egohhh);

    if isempty(Egoggg), continue, end

    try
        [timeSeriesData,Event]=preData_event([ExcelPath fnumIRT],FR_Fixed_Count);

    catch exception
               
        fprintf(fid, 'Warning: %s\n', ['cannot read ' ExcelPath fnumIRT ' ' exception.message]);
        warning(['cannot read ' ExcelPath fnumIRT ' ' exception.message])
        continue 
    end

        if any(contains(fieldnames(timeSeriesData), 'FR_Count'))
        continue
    end

%     nosepokeon = timeSeriesData.Nose_Poke_1;
    trayonon = timeSeriesData.Tray_1;
    nosepokeon = timeSeriesData.Nose_Poke_1;
%     trayonon = timeSeriesData.Pump_1;  



    %  histogram(ton,'BinWidth',20);   %3600 second 3600/20=180
    [N1,nedges] = histcounts(nosepokeon,'BinWidth',BinWidth);
    [N2,nedges] = histcounts(trayonon,'BinWidth',BinWidth);


        for kk=1:Egonumber
            close all
            fnumspk =[Egohhh{kk}];
            load(fnumspk);
            Egoname=fnumspk(1:19);

            %                 figure
            cellTS = ego.cellTS;


%% Plot



            [FiringRate,sedges] = histcounts(cellTS,'BinWidth',BinWidth); % 200s

            fr=FiringRate/100;

            numS = length(FiringRate);
            % FiringRate=FiringRate/max(FiringRate);
            FiringRate=(FiringRate-min(FiringRate))/(max(FiringRate)-min(FiringRate));
            ts = 0:BinWidth:BinWidth*(numS-1);

            num = length(N1);
            if num < numS
                NoseRate = zeros(numS);
                NoseRate(1:num) = N1;
                t = ts;
            elseif num > numS
                FiringRate1=FiringRate;
                FiringRate = zeros(num,1);
                FiringRate(1:numS) = FiringRate1;
                ts = t;

            else
                NoseRate = N1;
                t = 0:BinWidth:BinWidth*(num-1);
            end

            % If the time diff between behavior and electrophysiology
            % recording is so large ,then reject it
            if abs((numS-num)/numS)>0.5
                continue
            end

            %                 [NTrayon,sedges] = histcounts(N2,'BinWidth',BinWidth); % 200s
            %                 numS = length(N2);
            num = length(N2);
            if num < numS
                Trayonrate = zeros(numS);
                Trayonrate(1:num) = N2;
                t = ts;
            else
                Trayonrate = N2;
                t = 0:BinWidth:BinWidth*(num-1);
            end
            Trayonrate=(Trayonrate-min(Trayonrate))/(max(Trayonrate)-min(Trayonrate));

            if abs((numS-num)/numS)>0.5
                continue
            end

            %                 NoseRate=NoseRate/max(NoseRate);
            NoseRate=(NoseRate-min(NoseRate))/(max(NoseRate)-min(NoseRate));
            Trayonrate=(Trayonrate-min(Trayonrate))/(max(Trayonrate)-min(Trayonrate));

         
            EndTimePoint=1000;


%             corrValue = corrcoef(NoseRate,FiringRate);
            corrValue = corrcoef(Trayonrate,FiringRate);
            CORR = corrValue(1,2);
            b=['Correlation = ',num2str(CORR)];
            %                 text(2400,0.75, b,'FontSize',12,'Color','red')
            title(b)
            % filename = strcat(fnum(1:end-4),'_all');

%% 
% 


%             if abs(CORR)<0.2
        ManifoldData{ManifoldDataNumber,1}=fnumspk;
for i=2:max_state
        
    try
        % Note!!!: start from state 2
                [SpkWindow_Nosepoke]  = WindowPick(Event.(['Nosepoke_PRstate_' num2str(i)]),cellTS,WindowWidth,Ui,Li);
    catch exception
        fprintf(fid, 'Warning: %s\n', [fnumspk ' ' exception.message]);
warning([fnumspk ' ' exception.message])
        continue
    end
              

                ManifoldData{ManifoldDataNumber,i}=SpkWindow_Nosepoke;
%                 ManifoldData{ManifoldDataNumber,3}=SpkWindow_Reward;
%                 ManifoldData{ManifoldDataNumber,4}=CORR;
   

                
end
%             end


     ManifoldDataNumber=ManifoldDataNumber+1;

            close all
        end
    end



end






Size_NosepokeWiondow=size(ManifoldData{1,2});
% Size_RewardWiondow=size(ManifoldData{1,3});

HighDimData=zeros(Size_NosepokeWiondow(2),length(ManifoldData));





state=1;
for iii=state_start:state_end
for iiii=1:length(ManifoldData)
% HighDimData(i)
HighDimData((state-1)*WindowWidth+1:state*WindowWidth,iiii)=mean(ManifoldData{iiii,iii});

% HighDimData(Size_NosepokeWiondow(2)+1 : Size_NosepokeWiondow(2)+Size_RewardWiondow(2),i)...
%     =mean(ManifoldData{i,3});

% HighDimData=vertcat(HighDimData,mean(ManifoldData{i,2}));
% HighDimData=vertcat(HighDimData,mean(ManifoldData{i,3}));


end
state=state+1;
end
Timelag_per_bin=linspace(Li,Ui,WindowWidth)';
% Timelag_per_bin=[Timelag_per_bin];

% 初始化一个逻辑向量，用于标记每一列是否全为0
is_zero_column = all(HighDimData == 0, 1);
HighDimData = HighDimData(:, ~is_zero_column);


% save(['HighDimData_PR_np',HighDimData_Name{1,i}],'HighDimData','Timelag_per_bin')
save(data_save_name,'HighDimData','Timelag_per_bin')




%%
function sMap = hdFlatWindowSmoothing(map, numSmoothingBins)

% Number of bins in the map
N = length(map);

% Allocate memory for the smoothed map
sMap = zeros(1, N);

% Make sure the number of smoothing bins is a odd number
if mod(numSmoothingBins, 2) == 0
    numSmoothingBins = numSmoothingBins + 1;
end

% Number of bins to each side of the current bin when smoothing
d = (numSmoothingBins-1) / 2;

for ii = 1:N
    if ii-d <= 0 || ii+d > N
        if ii-d <= 0
            sumRate = sum(map(1:ii+d)) + sum(map(N-(d-ii):N));
            sMap(ii) = sumRate / numSmoothingBins;
        end
        if ii+d > N
            sumRate = sum(map(ii-d:N)) + sum(map(1:(ii+d-N)));
            sMap(ii) = sumRate / numSmoothingBins;
        end
    else
        sMap(ii) = nanmean(map(ii-d:ii+d));
    end
end
end