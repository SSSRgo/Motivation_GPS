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


%% Parameter Setup
ManifoldData={};
ManifoldDataNumber=1;
BinWidth = 50; %bin 50ms
FR_Fixed_Count=3;
binwidth=0.1

%% Data Load

ExcelPath='E:\GPS\cal\Motivation\Motivation_SR\ABET\'
ExceLgg1 = ls([ExcelPath '*csv']);
ExceLhh1 = cellstr(ExceLgg1);
ExceLnumber1 = length(ExceLhh1);
ExceLff1 = ExceLhh1;

for ii=1:ExceLnumber1  %%%session number
    fnumIRT =[ExceLff1{ii}];

    if ~isempty(findstr(fnumIRT,'IRT')), continue, end

    try
        fnumIRT='00397_08092301.csv';
        [timeSeriesData,Event]=preData_event([ExcelPath fnumIRT],FR_Fixed_Count);

    catch
        continue
        warning(['cannot read ' ExcelPath fnumIRT ])
    end

%     nosepokeon = timeSeriesData.Nose_Poke_1;
%     rewardon = timeSeriesData.Tray_1;
%     nosepokeon = Event.Nosepoke_Rcount_1;  
%     rewardon = Event.Trayon_pump;  
% 
% 
% 
%     %  histogram(ton,'BinWidth',20);   %3600 second 3600/20=180
%     [N1,nedges] = histcounts(nosepokeon,'BinWidth',BinWidth);
%     [N2,nedges] = histcounts(rewardon,'BinWidth',BinWidth);

    %         ggg = ls('00504_1*mat');
    Egoggg = ls([fnumIRT(1:end-8) '*mat']);
    Egohhh = cellstr(Egoggg);
    Egonumber = length(Egohhh);
    Egofff = Egohhh;


    if ~isempty(Egoggg)
        for kk=1:Egonumber
            close all
            fnumspk =[Egofff{kk}];
            fnumspk='00397_08092301_T3C2_ego.mat';
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


% 
%             [FiringRate,sedges] = histcounts(cellTS,'BinWidth',BinWidth); % 200s
% 
%             fr=FiringRate/200;
% 
%             numS = length(FiringRate);
%             % FiringRate=FiringRate/max(FiringRate);
%             FiringRate=(FiringRate-min(FiringRate))/(max(FiringRate)-min(FiringRate));
%             ts = 0:BinWidth:BinWidth*(numS-1);
% 
%             num = length(N1);
%             if num < numS
%                 NoseRate = zeros(numS);
%                 NoseRate(1:num) = N1;
%                 t = ts;
%             elseif num > numS
%                 FiringRate1=FiringRate;
%                 FiringRate = zeros(num,1);
%                 FiringRate(1:numS) = FiringRate1;
%                 ts = t;
% 
%             else
%                 NoseRate = N1;
%                 t = 0:BinWidth:BinWidth*(num-1);
%             end
% 
%             % If the time diff between behavior and electrophysiology
%             % recording is so large ,then reject it
%             if abs((numS-num)/numS)>0.5
%                 continue
%             end
% 
%             %                 [NTrayon,sedges] = histcounts(N2,'BinWidth',BinWidth); % 200s
%             %                 numS = length(N2);
%             num = length(N2);
%             if num < numS
%                 Trayonrate = zeros(numS);
%                 Trayonrate(1:num) = N2;
%                 t = ts;
%             else
%                 Trayonrate = N2;
%                 t = 0:BinWidth:BinWidth*(num-1);
%             end
%             Trayonrate=(Trayonrate-min(Trayonrate))/(max(Trayonrate)-min(Trayonrate));
% 
%             if abs((numS-num)/numS)>0.5
%                 continue
%             end
% 
%             %                 NoseRate=NoseRate/max(NoseRate);
%             NoseRate=(NoseRate-min(NoseRate))/(max(NoseRate)-min(NoseRate));
% 
%             f_corr=figure;
%             plot(ts,FiringRate,'LineWidth',1.5)
%             hold on
%             plot(ts,Trayonrate,'LineWidth',1.5)
%             plot(ts,NoseRate,'LineWidth',1.5)
%             %                 scatter(N2,PeakInternalTrayon)
%             legend('Firing Rate','Trayon Rate','Nosepoke Rate','Peak Gain')
%             set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.35])
%             EndTimePoint=1000;
% 
% 
%             corrValue = corrcoef(NoseRate,FiringRate);
%             CORR = corrValue(1,2);
%             b=['Correlation = ',num2str(CORR)];
%             %                 text(2400,0.75, b,'FontSize',12,'Color','red')
%             title(b)
%             % filename = strcat(fnum(1:end-4),'_all');

%% 
% 

% 
%             if abs(CORR)>0.3
%                 continue
%                 [SpkWindow_Nosepoke]  = WindowPick(nosepokeon,cellTS,200,2,-4);
%                 [SpkWindow_Reward]  = WindowPick(rewardon,cellTS,200,2,-4);
%                 [SpkWindow_Nosepoke]  = WindowPick(nosepokeon,cellTS,200,2,-4);
%                 [SpkWindow_Reward]  = WindowPick(rewardon,cellTS,200,2,-4);
% 
%                 ManifoldData{ManifoldDataNumber,1}=fnumspk;
%                 ManifoldData{ManifoldDataNumber,2}=SpkWindow_Nosepoke;
%                 ManifoldData{ManifoldDataNumber,3}=SpkWindow_Reward;
%                 ManifoldData{ManifoldDataNumber,4}=CORR;
%                 % if CORR>0
%                 %     ManifoldData{ManifoldDataNumber,4}=1;
%                 % else
%                 %     ManifoldData{ManifoldDataNumber,4}=-1;
%                 % end
% 
%                 ManifoldDataNumber=ManifoldDataNumber+1;
%             end
%% 
% 
%% Analysis of Peak Gain



            %%% Time of Reward
            %                 gg = ls('E:\GPS\cal\Motivation\Motivation_SR\ABET\*csv');
            %                 hh = cellstr(gg);
            %                 number = length(hh);
            %                 ff = hh;
%             jj = 1;
%             jjj = 1;
% 
%             clear rewardtime trayontime
%             k=0;
% 
%             index = findstr(fnumIRT,'IRT');
%             fnum=fnumIRT(1:index-2);
%             fnum=[fnum '.csv'];
%             rawdata = RAW_Read(fnum);
%             pump = rawdata(:,4);
%             event = rawdata(:,3);
%             time = rawdata(:,1);
%             for ii = 1:length(pump)
% 
%                 if strcmp(pump{ii},'Pump #1')
%                     rewardtime(jj,1) = time{ii};
%                     jj = jj+1;
%                     k=1;
%                 end
% 
%                 if strcmp(pump{ii},'Tray #1')&&strcmp(event{ii},'Input Transition On Event')&&k==1
%                     trayontime(jjj,1) = time{ii};
%                     jjj = jjj+1;
%                     k=0;
%                 end
%             end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% trayontime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             trayontime=rewardon;
%             [NTrayon,sedges] = histcounts(trayontime,'BinWidth',BinWidth); % 200s
% 
%             num = length(NTrayon);
%             if num < numS
%                 Trayonrate = zeros(numS);
%                 Trayonrate(1:num) = NTrayon;
%                 t = ts;
%             else
%                 Trayonrate = NTrayon;
%                 t = 0:BinWidth:BinWidth*(num-1);
%             end
%             Trayonrate=(Trayonrate-min(Trayonrate))/(max(Trayonrate)-min(Trayonrate));
% 
% 
%             %            CorrelationComponent=corrcoef([NoseRate,Trayonrate,NoseRate+Trayonrate,FiringRate']);
%             %%%%%%%%%% tray on or Nose %%%%%%%%%%
%             CorrInWindow=zeros(length(20:20),length(5:40));


%% Event Response Tray on


            %%%%%%%%%%%%%%%%%%% 1: All Tray on
            Li=-3;
            Inti=8;
            LeftBound=Li; RightBound=Li+Inti;
            TWindow_trayon=0:binwidth:(RightBound-LeftBound);
            [NWindow_trayon,trayon_starpoint] = EventResonse(timeSeriesData.Tray_1,cellTS,Li,Inti,binwidth);

            f_peak_trayon=figure;
            set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.35])
            plot(linspace(LeftBound,RightBound,length(TWindow_trayon)-1),NWindow_trayon);
            xlabel('Time (s)')
            ylabel('Firing Rate (Hz)')
            title([Egoname ' Trayontime'], 'Interpreter', 'none')

            % peak definition
            PeakInternalTrayon=mean(NWindow_trayon(14:18,:));
            % aviod
            if sum(abs(PeakInternalTrayon))==0, continue
            end
            PeakInternalTrayon=(PeakInternalTrayon-min(PeakInternalTrayon))/(max(PeakInternalTrayon)-min(PeakInternalTrayon));

            %%%%%%%%%%%%%%%%%%% 2: Tray on


%% Event Response Nose Poke


            Li=-5;
            Inti=15;
            
            LeftBound=Li; RightBound=Li+Inti;
            TWindow_nosepoke=0:binwidth:(RightBound-LeftBound);
            [NWindow_nosepoke,nosepoke_starpoint] = EventResonse(timeSeriesData.Nose_Poke_1,cellTS,Li,Inti,binwidth);

            f_peak_nosepoke=figure;
            set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.35])
            plot(linspace(LeftBound,RightBound,length(TWindow_nosepoke)-1),NWindow_nosepoke);
            xlabel('Time (s)')
            ylabel('Firing Rate (Hz)')
            title([Egoname ' Nosepokeon'], 'Interpreter', 'none')

%% Mean Event Response Nose Poke of FR count 1-3

            Li=-2;
            Inti=6;
            LeftBound=Li; RightBound=Li+Inti;
            TWindow_nosepoke=0:binwidth:(RightBound-LeftBound);
            [NWindow_nosepoke_Rcount_1,nosepoke_starpoint_Rcount_1] = EventResonse(Event.Nosepoke_Rcount_1,cellTS,Li,Inti,binwidth);
            [NWindow_nosepoke_Rcount_2,nosepoke_starpoint_Rcount_2] = EventResonse(Event.Nosepoke_Rcount_2,cellTS,Li,Inti,binwidth);
            [NWindow_nosepoke_Rcount_3,nosepoke_starpoint_Rcount_3] = EventResonse(Event.Nosepoke_Rcount_3,cellTS,Li,Inti,binwidth);

            f_MeanEventResponse_Nosepoke=figure;
%             set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.5])
            shadedErrorBar(linspace(LeftBound,RightBound,length(TWindow_nosepoke)-1),NWindow_nosepoke_Rcount_1',{@mean,@std},'lineprops','-r','transparent',true,'patchSaturation',0.4); 
            hold on 
            shadedErrorBar(linspace(LeftBound,RightBound,length(TWindow_nosepoke)-1),NWindow_nosepoke_Rcount_2',{@mean,@std},'lineprops','-b','transparent',true,'patchSaturation',0.4); 
            shadedErrorBar(linspace(LeftBound,RightBound,length(TWindow_nosepoke)-1),NWindow_nosepoke_Rcount_3',{@mean,@std},'lineprops','-g','transparent',true,'patchSaturation',0.4); 

%             plot(linspace(LeftBound,RightBound,length(TWindow_nosepoke)-1),NWindow_nosepoke);
            xlabel('Time (s)')
            ylabel('Firing Rate (Hz)')
            title([Egoname ' Event Response '], 'Interpreter', 'none')
            legend(['#1';'#2';'#3'],'location','northwest')
%% Mean Event Response Tray on Comparing Pump and Non-pump

            Li=-2;
            Inti=6;
            LeftBound=Li; RightBound=Li+Inti;
            TWindow_trayon=0:binwidth:(RightBound-LeftBound);



   
            TWindow_nosepoke=0:binwidth:(RightBound-LeftBound);
            [NWindow_Trayon_pump,] = EventResonse(Event.Trayon_pump,cellTS,Li,Inti,binwidth);
            [NWindow_Trayon_nonpump,] = EventResonse(Event.Trayon_nonpump,cellTS,Li,Inti,binwidth);
        

            f_MeanEventResponse_trayon_pump=figure;
%             set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.5])
            shadedErrorBar(linspace(LeftBound,RightBound,length(TWindow_trayon)-1),NWindow_Trayon_pump',{@mean,@std},'lineprops','-r','transparent',true,'patchSaturation',0.4); 
            hold on 
            shadedErrorBar(linspace(LeftBound,RightBound,length(TWindow_trayon)-1),NWindow_Trayon_nonpump',{@mean,@std},'lineprops','-b','transparent',true,'patchSaturation',0.4); 
           

%             plot(linspace(LeftBound,RightBound,length(TWindow_nosepoke)-1),NWindow_nosepoke);
            xlabel('Time (s)')
            ylabel('Firing Rate (Hz)')
            title([Egoname ' Event Response '], 'Interpreter', 'none')
            legend('pump','nonpump','location','northwest')



%% Plot Peak Gain

% 
%             figure
%             plot(ts,FiringRate,'LineWidth',1.5)
%             hold on
%             plot(ts,Trayonrate,'LineWidth',1.5)
%             plot(ts,NoseRate,'LineWidth',1.5)
%             scatter(trayontime(trayon_starpoint:end),PeakInternalTrayon)
%             legend('Firing Rate','Trayon Rate','Nosepoke Rate','Peak Gain')
%             set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.35])
%             EndTimePoint=1000;
% 
% 
%             figure
%             %                 f_peak_gain=
%             plot(ts,FiringRate,'LineWidth',1.5)
%             hold on
%             plot(ts(1:EndTimePoint/BinWidth),Trayonrate(1:EndTimePoint/BinWidth),'LineWidth',1.5)
%             plot(ts(1:EndTimePoint/BinWidth),NoseRate(1:EndTimePoint/BinWidth),'LineWidth',1.5)
%             legend('Firing Rate','Trayon Rate','Nosepoke Rate')
% 
% 
%             [fitresult, gof,f_peak_gain] = createFit2(trayontime(trayon_starpoint:end), PeakInternalTrayon);
%             hold on
%             plot(ts,Trayonrate,'LineWidth',1.5)
%             plot(ts,NoseRate,'LineWidth',1.5)
%             %                 scatter(trayontime,PeakInternalTrayon)
%             legend('Peak Gain','Fitting Curve','Trayon Rate','Nosepoke Rate')
%             set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.35])
%% Save figure

%             exportgraphics(f_corr,[Egoname '_1' 'Trayon_pump&Nosepoke_Rcount_1' '.jpg'],'Resolution',300)
%             exportgraphics(f_peak_trayon,[Egoname '_2' '_Trayon_pump' '.jpg'],'Resolution',300)
%             exportgraphics(f_peak_nosepoke,[Egoname '_3' '_Nosepoke_Rcount_1' '.jpg'],'Resolution',300)
             exportgraphics(f_peak_nosepoke_count,[Egoname '_5' '_peak_nosepoke_count' '.jpg'],'Resolution',300)
             exportgraphics(f_peak_trayon,[Egoname '_5' '_peak_nosepoke_count' '.jpg'],'Resolution',300)
             exportgraphics(f_MeanEventResponse_trayon_pump,[Egoname '_5' '_MeanEventResponse_trayon_pump' '.jpg'],'Resolution',300)
             exportgraphics(f_MeanEventResponse_Nosepoke,[Egoname '_5' '_MeanEventResponse_nosepoke_count' '.jpg'],'Resolution',300)

%             exportgraphics(f_peak_gain,[Egoname '_4' '.jpg'],'Resolution',300)
            % exportgraphics(f5,[Egoname '_5' '.jpg'],'Resolution',300)
            % exportgraphics(f6,[Egoname '_6' '.jpg'],'Resolution',300)
            % exportgraphics(f7,[Egoname '_7' '.jpg'],'Resolution',300)
%% 
% 

            close all
        end
    end



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