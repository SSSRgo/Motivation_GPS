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

ExceLgg1 = ls('E:\GPS\cal\Motivation\Motivation_gps\Code\ABET\*csv');
ExceLhh1 = cellstr(ExceLgg1);
ExceLnumber1 = length(ExceLhh1);
ExceLff1 = ExceLhh1;
        ManifoldData={};
ManifoldDataNumber=1;


for ii=1:ExceLnumber1  %%%session number
    fnumIRT =[ExceLff1{ii}];
    index = findstr(fnumIRT,'IRT');

    if index
        IRT = IRT_Read(fnumIRT);
        count = IRT(10:end,1);
        nosepoke = IRT(10:end,2);
        reward = IRT(10:end,3);

        M = length(count);
        nosepokeon = zeros(M,1);
        rewardon = zeros(M,1);

        % ton(1,1) = arrayIRI(1,2);
        % toff(1,1) = arrayIRI(1,3);
        for i = 2:M
            nosepokeon(i) = nosepokeon(i-1)+nosepoke(i);
            rewardon(i) = rewardon(i-1)+reward(i);
        end

        BinWidth = 50;
        %  histogram(ton,'BinWidth',20);   %3600 second 3600/20=180
        [N1,nedges] = histcounts(nosepokeon,'BinWidth',BinWidth);
        [N2,nedges] = histcounts(rewardon,'BinWidth',BinWidth);

        %         ggg = ls('00504_1*mat');
        Egoggg = ls([fnumIRT(1:end-8) '*mat']);
        Egohhh = cellstr(Egoggg);
        Egonumber = length(Egohhh);
        Egofff = Egohhh;


        if ~isempty(Egoggg)
        for kk=1:Egonumber
            
            fnumspk =[Egofff{kk}];
            % fnumspk = '00394_12082301_T1C3_ego.mat';
            % %             fnumspk = '00397_03082301_T4C1_ego.mat';
%             indexspk = strcmp(fnumIRT(1:end-8),fnumspk(1:end-13));
% indexspk=exist (fnumspk,'file');
%             if indexspk~=0
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                load(fnumspk);
                figure
                cellTS = ego.cellTS;
                [FiringRate,sedges] = histcounts(cellTS,'BinWidth',BinWidth); % 200s
                
                fr=FiringRate/200;

                numS = length(FiringRate);
                % FiringRate=FiringRate/max(FiringRate);
                FiringRate=(FiringRate-min(FiringRate))/(max(FiringRate)-min(FiringRate));
                ts = 0:BinWidth:BinWidth*(numS-1);

                % legend('Normalized Firing Rate')

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


                %                 NoseRate=NoseRate/max(NoseRate);
                NoseRate=(NoseRate-min(NoseRate))/(max(NoseRate)-min(NoseRate));
                plot(ts,FiringRate,'LineWidth',1.5)
                hold on
                plot(t,NoseRate,'LineWidth',1.5)


                corrValue = corrcoef(NoseRate,FiringRate);
                CORR = corrValue(1,2)
                b=['Correlation = ',num2str(CORR)];
                %                 text(2400,0.75, b,'FontSize',12,'Color','red')
                title(b)
                % filename = strcat(fnum(1:end-4),'_all');
                if abs(CORR)>0.3
[SpkWindow_Nosepoke,SpkWindow_Reward]  = WindowPick(nosepokeon,rewardon,cellTS,100,1,-2);

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

                %%% Time of Reward
                %                 gg = ls('E:\GPS\cal\Motivation\Motivation_gps\Code\ABET\*csv');
                %                 hh = cellstr(gg);
                %                 number = length(hh);
                %                 ff = hh;
                jj = 1;
                jjj = 1;

                clear rewardtime trayontime
                k=0;

                index = findstr(fnumIRT,'IRT');
                fnum=fnumIRT(1:index-2);
                fnum=[fnum '.csv'];
                rawdata = RAW_Read(fnum);
                pump = rawdata(:,4);
                event = rawdata(:,3);
                time = rawdata(:,1);
                for ii = 1:length(pump)

                    if strcmp(pump{ii},'Pump #1')
                        rewardtime(jj,1) = time{ii};
                        jj = jj+1;
                        k=1;
                    end

                    if strcmp(pump{ii},'Tray #1')&&strcmp(event{ii},'Input Transition On Event')&&k==1
                        trayontime(jjj,1) = time{ii};
                        jjj = jjj+1;
                        k=0;
                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%% trayontime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                [NTrayon,sedges] = histcounts(trayontime,'BinWidth',BinWidth); % 200s

                num = length(NTrayon);
                if num < numS
                    Trayonrate = zeros(numS);
                    Trayonrate(1:num) = NTrayon;
                    t = ts;
                else
                    Trayonrate = NTrayon;
                    t = 0:BinWidth:BinWidth*(num-1);
                end
                Trayonrate=(Trayonrate-min(Trayonrate))/(max(Trayonrate)-min(Trayonrate));


                %            CorrelationComponent=corrcoef([NoseRate,Trayonrate,NoseRate+Trayonrate,FiringRate']);
                %%%%%%%%%% tray on or Nose %%%%%%%%%%
                CorrInWindow=zeros(length(20:20),length(5:40));
                figure
                set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.35])
                for Li=-3:-3
                    for Inti=8:8
                        LeftBound=Li;
                        RightBound=Li+Inti;
                        CurveN=1;
                        clear NWindow
                        for i=1:length(trayontime) % nosepokeon trayontime
                            timestemp=trayontime(i);
                            if timestemp<-LeftBound
                                continue
                            end
                            WindowCellTS=cellTS(cellTS>(timestemp+LeftBound)&cellTS<(timestemp+RightBound));
                            WindowCellTS=WindowCellTS-timestemp-LeftBound;
                            TWindow=0:0.2:(RightBound-LeftBound);
                            [WindowCurve,sedges]=histcounts(WindowCellTS,TWindow);
                            NWindow(:,CurveN)=hdFlatWindowSmoothing(WindowCurve,4);
                            CurveN=CurveN+1;

                            %
                            %                     aaa=corrcoef(NWindow);
                            %                      a=triu(aaa);
                            %                      a(isnan(a))=0;
                            %                      a(a==1)=0;
                            %                      CorrInWindow(Li+21,Inti-4)=sum(a(:))/nnz(a);
                        end

                        plot(linspace(LeftBound,RightBound,length(TWindow)-1),NWindow);
                        xlabel('Time (s)')
                        ylabel('Firing Rate (Hz)')
                        title('trayontime')


                    end
                end

                PeakInternalTrayon=mean(NWindow(14:18,:));
                PeakInternalTrayon=(PeakInternalTrayon-min(PeakInternalTrayon))/(max(PeakInternalTrayon)-min(PeakInternalTrayon));
                %                 scatter(trayontime,PeakInternalTrayon)
                %                 hold on
                %                 plot(ts(1:ceil(trayontime(end)/BinWidth)),Trayonrate(1:ceil(trayontime(end)/BinWidth)))



                figure
                set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.35])
                for Li=-5:-5
                    for Inti=15:15
                        LeftBound=Li;
                        RightBound=Li+Inti;
                        CurveN=1;
                        clear NWindow
                        for i=1:length(nosepokeon) % nosepokeon trayontime
                            timestemp=nosepokeon(i);
                            if timestemp<-LeftBound
                                continue
                            end
                            WindowCellTS=cellTS(cellTS>(timestemp+LeftBound)&cellTS<(timestemp+RightBound));
                            WindowCellTS=WindowCellTS-timestemp-LeftBound;
                            TWindow=0:0.2:(RightBound-LeftBound);
                            [WindowCurve,sedges]=histcounts(WindowCellTS,TWindow);
                            NWindow(:,CurveN)=hdFlatWindowSmoothing(WindowCurve,4);
                            CurveN=CurveN+1;

                            %
                            %                     aaa=corrcoef(NWindow);
                            %                      a=triu(aaa);
                            %                      a(isnan(a))=0;
                            %                      a(a==1)=0;
                            %                      CorrInWindow(Li+21,Inti-4)=sum(a(:))/nnz(a);
                        end

                        plot(linspace(LeftBound,RightBound,length(TWindow)-1),NWindow);
                        xlabel('Time (s)')
                        ylabel('Firing Rate (Hz)')
                        title('nosepokeon')


                    end
                end

                figure
                plot(ts,FiringRate,'LineWidth',1.5)
                hold on
                plot(ts,Trayonrate,'LineWidth',1.5)
                plot(ts,NoseRate,'LineWidth',1.5)
                scatter(trayontime,PeakInternalTrayon)
                legend('Firing Rate','Trayon Rate','Nosepoke Rate','Peak Gain')
                set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.35])
                EndTimePoint=1000;

                figure
                plot(ts,FiringRate,'LineWidth',1.5)
                hold on
                plot(ts(1:EndTimePoint/BinWidth),Trayonrate(1:EndTimePoint/BinWidth),'LineWidth',1.5)
                plot(ts(1:EndTimePoint/BinWidth),NoseRate(1:EndTimePoint/BinWidth),'LineWidth',1.5)
                legend('Firing Rate','Trayon Rate','Nosepoke Rate')


                [fitresult, gof] = createFit2(trayontime, PeakInternalTrayon)
                hold on
                plot(ts,Trayonrate,'LineWidth',1.5)
                plot(ts,NoseRate,'LineWidth',1.5)
                %                 scatter(trayontime,PeakInternalTrayon)
                legend('Peak Gain','Fitting Curve','Trayon Rate','Nosepoke Rate')
                set(gcf,'Units','normalized','position',[0.2,0.2,0.6,0.35])





% [SpkWindow_Nosepoke,SpkWindow_Reward]  = WindowPick(nosepokeon,rewardon,cellTS,100,1,-2);




%             end


            close all
        end
        end
    end


end














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
