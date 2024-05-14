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

HighDimData_Name={'Nosepoke_slice1','Nosepoke_slice2','Trayon_slice1','Trayon_slice2'};

for i =3:length(HighDimData_Name)

    clear all；

ManifoldData={};
ManifoldDataNumber=1;
BinWidth = 50; %bin 50ms
FR_Fixed_Count=3;
binwidth=0.1;
%% Data Load

ExcelPath='E:\GPS\cal\Motivation\Motivation_gps\Code\ABET\'
ExceLgg1 = ls([ExcelPath '*csv']);
ExceLhh1 = cellstr(ExceLgg1);
ExceLnumber1 = length(ExceLhh1);
ExceLff1 = ExceLhh1;

for ii=1:ExceLnumber1  %%%session number
    fnumIRT =[ExceLff1{ii}];

    if ~isempty(findstr(fnumIRT,'IRT')), continue, end

    try
        [timeSeriesData,Event]=preData_event([ExcelPath fnumIRT],FR_Fixed_Count);

    catch
        continue
        warning(['cannot read ' ExcelPath fnumIRT ])
    end

%     nosepokeon = timeSeriesData.Nose_Poke_1;
%     rewardon = timeSeriesData.Tray_1;
    nosepokeon = Event.(HighDimData_Name{1,i});


    
    trayonon = Event.Trayon_pump;  



    %  histogram(ton,'BinWidth',20);   %3600 second 3600/20=180
    [N1,nedges] = histcounts(nosepokeon,'BinWidth',BinWidth);
    [N2,nedges] = histcounts(trayonon,'BinWidth',BinWidth);

    %         ggg = ls('00504_1*mat');
    Egoggg = ls([fnumIRT(1:end-8) '*mat']);
    Egohhh = cellstr(Egoggg);
    Egonumber = length(Egohhh);
    Egofff = Egohhh;


    if ~isempty(Egoggg)
        for kk=1:Egonumber
            close all
            fnumspk =[Egofff{kk}];
            % fnumspk = '00394_12082301_T1C3_ego.mat';
            % %             fnumspk = '00397_03082301_T4C1_ego.mat';
            %             indexspk = strcmp(fnumIRT(1:end-8),fnumspk(1:end-13));
            % indexspk=exist (fnumspk,'file');
            %             if indexspk~=0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            load(fnumspk);
            Egoname=fnumspk(1:19);

            %                 figure
            cellTS = ego.cellTS;


%% Plot



            [FiringRate,sedges] = histcounts(cellTS,'BinWidth',BinWidth); % 200s

            fr=FiringRate/200;

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


            if abs(CORR)<0.2
                [SpkWindow_Nosepoke]  = WindowPick(nosepokeon,cellTS,200,2,-4);
                [SpkWindow_Reward]  = WindowPick(trayonon,cellTS,200,2,-4);
                [SpkWindow_Nosepoke]  = WindowPick(nosepokeon,cellTS,200,2,-4);
                [SpkWindow_Reward]  = WindowPick(trayonon,cellTS,200,2,-4);

                ManifoldData{ManifoldDataNumber,1}=fnumspk;
                ManifoldData{ManifoldDataNumber,2}=SpkWindow_Nosepoke;
                ManifoldData{ManifoldDataNumber,3}=SpkWindow_Reward;
                ManifoldData{ManifoldDataNumber,4}=CORR;
                % if CORR>0
                %     ManifoldData{ManifoldDataNumber,4}=1;
                % else
                %     ManifoldData{ManifoldDataNumber,4}=-1;
                % end

                ManifoldDataNumber=ManifoldDataNumber+1;
                
            end


     

            close all
        end
    end



end



            Size_NosepokeWiondow=size(ManifoldData{1,2});
Size_RewardWiondow=size(ManifoldData{1,3});

HighDimData=zeros(Size_NosepokeWiondow(2),length(ManifoldData));

for iiii=1:length(ManifoldData)
% HighDimData(i)
HighDimData(1:Size_NosepokeWiondow(2),iiii)=mean(ManifoldData{iiii,2});
% HighDimData(Size_NosepokeWiondow(2)+1 : Size_NosepokeWiondow(2)+Size_RewardWiondow(2),i)...
%     =mean(ManifoldData{i,3});

% HighDimData=vertcat(HighDimData,mean(ManifoldData{i,2}));
% HighDimData=vertcat(HighDimData,mean(ManifoldData{i,3}));


end

Timelag_per_bin=linspace(-4,2,200)';
Timelag_per_bin=[Timelag_per_bin];




save(['HighDimData_',HighDimData_Name{1,i}],'HighDimData','Timelag_per_bin')


end


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