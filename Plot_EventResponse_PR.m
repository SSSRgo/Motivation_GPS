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
binwidth=0.1;

Li_Trayon=-4;
Inti_Trayon=8;

Li_nosepoke=-4;
Inti_nosepoke=8;

Li_state_np=-4;
Inti_state_np=8;

Li_pumpOrnopump=-4;
Inti_pumpOrnopump=8;

ExcelPath='E:\GPS\cal\Motivation\Motivation_SR\ABET\';
DataPath='E:\GPS\cal\Motivation\Motivation_SR\Data\';
figure_save_path='E:\GPS\cal\Motivation\Motivation_SR\Figure_PR\';

addpath(ExcelPath)
addpath(DataPath)

%% Data Load



ExceLgg1 = ls([ExcelPath '*csv']);
ExceLhh1 = cellstr(ExceLgg1);
ExceLnumber1 = length(ExceLhh1);
ExceLff1 = ExceLhh1;

for session_number=663:ExceLnumber1  %%%session number
    fnumIRT =[ExceLff1{session_number}];

    if ~isempty(findstr(fnumIRT,'IRT')), continue, end

    try
        %         fnumIRT='00397_08092301.csv';
        [timeSeriesData,Event]=preData_event([ExcelPath fnumIRT],FR_Fixed_Count);

    catch
        continue
        warning(['cannot read ' ExcelPath fnumIRT ])
    end

    %%%%%%%%%%% if it is
    
    if any(contains(fieldnames(timeSeriesData), 'FR_Count'))
        continue
    end



    %         ggg = ls('00504_1*mat');
    Egoggg = ls([DataPath fnumIRT(1:end-8) '*mat']);
    Egohhh = cellstr(Egoggg);
    Egonumber = length(Egohhh);
    Egofff = Egohhh;


    if ~isempty(Egoggg)
        for kk=1:Egonumber
            close all
            fnumspk =[Egofff{kk}];

            load(fnumspk);

            %             save(['a',fnumspk],ego)
            Egoname=fnumspk(1:19);

            %                 figure
            cellTS = ego.cellTS;


            %% Firing rate (correlation) Plot


            %  histogram(ton,'BinWidth',20);   %3600 second 3600/20=180
            [N1,nedges] = histcounts(timeSeriesData.Nose_Poke_1,'BinWidth',BinWidth);
            [N2,nedges] = histcounts(timeSeriesData.Tray_1,'BinWidth',BinWidth);

            [FiringRate,sedges] = histcounts(cellTS,'BinWidth',BinWidth); % 200s

            fr=FiringRate/200;
            FiringRate=(FiringRate-min(FiringRate))/(max(FiringRate)-min(FiringRate));

            numS = length(FiringRate);
            num = length(N1);
            % FiringRate=FiringRate/max(FiringRate);
            
            ts = 0:BinWidth:BinWidth*(numS-1);
            t = 0:BinWidth:BinWidth*(num-1);

            
            if num < numS
                NoseRate = zeros(numS);
                NoseRate(1:num) = N1;
                t = ts;
            else
                FiringRate1=FiringRate;
                FiringRate = zeros(num,1);
                FiringRate(1:numS) = FiringRate1;
                ts = t;
                NoseRate = N1;

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

            f_corr=figure('Visible','off');
            
            plot(ts,FiringRate,'LineWidth',1.5)
            hold on
            plot(ts,Trayonrate,'LineWidth',1.5)
            plot(ts,NoseRate,'LineWidth',1.5)
            %                 scatter(N2,PeakInternalTrayon)
            legend('Firing Rate','Trayon Rate','Nosepoke Rate','Peak Gain')
            set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.35])

            corrValue = corrcoef(NoseRate,FiringRate);
            CORR = corrValue(1,2);
            b=['Correlation = ',num2str(CORR)];
            %                 text(2400,0.75, b,'FontSize',12,'Color','red')
            title(b)
            % filename = strcat(fnum(1:end-4),'_all');


            %% Event Response Tray on


            %%%%%%%%%%%%%%%%%%% 1: All Tray on

            LeftBound=Li_Trayon; RightBound=Li_Trayon+Inti_Trayon;
            TWindow_trayon=0:binwidth:(RightBound-LeftBound);
            [NWindow_trayon,trayon_starpoint] = EventResonse(timeSeriesData.Tray_1,cellTS,Li_Trayon,Inti_Trayon,binwidth);

            f_peak_trayon=figure('Visible','off');
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

            %% Event Response Nose Poke

            LeftBound=Li_nosepoke; RightBound=Li_nosepoke+Inti_nosepoke;
            TWindow_nosepoke=0:binwidth:(RightBound-LeftBound);
            [NWindow_nosepoke,nosepoke_starpoint] = EventResonse(timeSeriesData.Nose_Poke_1,cellTS,Li_nosepoke,Inti_nosepoke,binwidth);

            f_peak_nosepoke=figure('Visible','off');
            set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.35])
            plot(linspace(LeftBound,RightBound,length(TWindow_nosepoke)-1),NWindow_nosepoke);
            xlabel('Time (s)')
            ylabel('Firing Rate (Hz)')
            title([Egoname ' Nosepokeon'], 'Interpreter', 'none')



            %% Mean Event Response Tray on Comparing Pump and Non-pump

            LeftBound=Li_pumpOrnopump; RightBound=Li_pumpOrnopump+Inti_pumpOrnopump;
            TWindow_trayon=0:binwidth:(RightBound-LeftBound);


            TWindow_nosepoke=0:binwidth:(RightBound-LeftBound);
            [NWindow_Trayon_pump,] = EventResonse(Event.Trayon_pump,cellTS,Li_pumpOrnopump,Inti_pumpOrnopump,binwidth);
            [NWindow_Trayon_nonpump,] = EventResonse(Event.Trayon_nonpump,cellTS,Li_pumpOrnopump,Inti_pumpOrnopump,binwidth);


            f_MeanEventResponse_trayon_pump=figure('Visible','off');
            %             set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.5])
            shadedErrorBar(linspace(LeftBound,RightBound,length(TWindow_trayon)-1),NWindow_Trayon_pump',{@mean,@std},'lineprops','-r','transparent',true,'patchSaturation',0.4);
            hold on
            shadedErrorBar(linspace(LeftBound,RightBound,length(TWindow_trayon)-1),NWindow_Trayon_nonpump',{@mean,@std},'lineprops','-b','transparent',true,'patchSaturation',0.4);


            %             plot(linspace(LeftBound,RightBound,length(TWindow_nosepoke)-1),NWindow_nosepoke);
            xlabel('Time (s)')
            ylabel('Firing Rate (Hz)')
            title([Egoname ' Event Response '], 'Interpreter', 'none')
            legend('Pump','Non-pump','location','northwest')

            %% PR: Mean Event Response of Nose Poke in different states (states 1-end)


            LeftBound=Li_state_np; RightBound=Li_state_np+Inti_state_np;
            TWindow_nosepoke=0:binwidth:(RightBound-LeftBound);
            [NWindow_nosepoke_Rcount_1,nosepoke_starpoint_Rcount_1] = EventResonse(Event.Nosepoke_PRstate_4,cellTS,Li_state_np,Inti_state_np,binwidth);
            f_MeanEventResponse_Nosepoke=figure('Visible','off');
            set(gcf, 'Position', [100, 100, 900, 600]);
            hold on
            ii=1;
            win_act_state=[];
            %             cmap22=colormap(othercolor('YlOrBr7',length(timeSeriesData.Pump_1)));
            cmap22=parula(length(timeSeriesData.Pump_1));
            for i = 1:length(timeSeriesData.Pump_1)



                if i==length(timeSeriesData.Pump_1)

                    if length(Event.('Nosepoke_PRstate_last'))>1
                        [NWindow_nosepoke_Rcount,nosepoke_starpoint_Rcount_1] = EventResonse(Event.('Nosepoke_PRstate_last'),cellTS,Li_state_np,Inti_state_np,binwidth);
                       
                    else
                        continue
                    end

                else

                    if length(Event.(['Nosepoke_PRstate_',num2str(i)]))<=1
                        continue
                    end

                    %                 [NWindow_nosepoke_Rcount,nosepoke_starpoint_Rcount_1] = EventResonse(Event.(['Nosepoke_PRstate_',num2str(i)]),cellTS,Li_state_np,Inti_state_np,binwidth);
                    [NWindow_nosepoke_Rcount,nosepoke_starpoint_Rcount_1] = EventResonse(Event.(['Nosepoke_PRstate_',num2str(i)]),cellTS,Li_state_np,Inti_state_np,binwidth);
                end

                 win_act_state(ii,:)=mean(NWindow_nosepoke_Rcount,2);

                if nosepoke_starpoint_Rcount_1~=1
                    
                    continue
                end


                %                 h =shadedErrorBar(linspace(LeftBound,RightBound,length(TWindow_nosepoke)-1),NWindow_nosepoke_Rcount',{@mean,@std},'lineprops','-r','transparent',true,'patchSaturation',0.4);
                p=shadedErrorBar(linspace(LeftBound,RightBound,length(TWindow_nosepoke)-1),NWindow_nosepoke_Rcount',{@mean,@std},'lineprops','-','transparent',true,'patchSaturation',0.1);
                p.patch.FaceColor=cmap22(i,:);
                p.mainLine.Color=cmap22(i,:);
                p.mainLine.LineWidth =2;
                p.edge(1).Color='none';
                p.edge(2).Color='none';
                p.edge=[];
                if i==length(timeSeriesData.Pump_1)
                    legend_labels{ii} = ['State ', num2str(i),', Num: ',num2str(length(Event.(['Nosepoke_PRstate_last'])))];
                else
                    legend_labels{ii} = ['State ', num2str(i),', Num: ',num2str(length(Event.(['Nosepoke_PRstate_',num2str(i)])))]; % 设置每个图例的标签
                end
                ii=ii+1;
            end

            legend(legend_labels, 'Location', 'best'); % 添加图例

            %             set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.5])

            %             plot(linspace(LeftBound,RightBound,length(TWindow_nosepoke)-1),NWindow_nosepoke);
            ylim([0 inf])
            xlabel('Time (s)')
            ylabel('Firing Rate (Hz)')
            title([Egoname ' Event Response '], 'Interpreter', 'none')
            

            %% Oder heatmap
            win_act_state=(win_act_state-min(min(win_act_state)))./(max(max(win_act_state))-min(min(win_act_state)));
            win_act_state = smoothdata(win_act_state,2);
            f_orderheatmap=figure('Visible','off');
            imagesc(win_act_state);
            
            xlabel('Time (s)')
            ylabel('State')
            colormap("hot")
%             colormap(othercolor('Reds8'))
h=colorbar;
         set(get(h,'Title'),'string','fr');
% colormap("hot")
set(gca,'YDIR','reverse');

daspect([2,1,1])
ax = gca; % current axes
% ax.FontSize = 12;
% ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
ax.XTick=[];
ax.YTick=[1 size(win_act_state,1)];
title([Egoname ' Event Response '], 'Interpreter', 'none');

%            colormap(othercolor('Oranges7'))
            



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
            file_path=[figure_save_path,Egoname(1:5),'\'];
%             exportgraphics(f_corr,[file_path,Egoname '_1' 'Trayon_pump&Nosepoke_Rcount_1' '.jpg'],'Resolution',300)
%             exportgraphics(f_peak_trayon,[file_path,Egoname '_2' '_Trayon_pump' '.jpg'],'Resolution',300)
%             exportgraphics(f_peak_nosepoke,[file_path,Egoname '_3' '_Nosepoke_Rcount_1' '.jpg'],'Resolution',300)
            %              exportgraphics(f_peak_nosepoke_count,[file_path,Egoname '_5' '_peak_nosepoke_count' '.jpg'],'Resolution',300)
            %              exportgraphics(f_peak_trayon,[file_path,Egoname '_5' '_peak_nosepoke_count' '.jpg'],'Resolution',300)
            %              exportgraphics(f_MeanEventResponse_trayon_pump,[file_path,Egoname '_5' '_MeanEventResponse_trayon_pump' '.jpg'],'Resolution',300)
%             exportgraphics(f_MeanEventResponse_Nosepoke,[file_path,Egoname '_5_PR' '_MeanEventResponse_nosepoke_state' '.jpg'],'Resolution',300)
            exportgraphics(f_orderheatmap,[file_path,Egoname '_6_PR' '_MeanEventResponse_heat_nosepoke_state' '.jpg'],'Resolution',300)

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