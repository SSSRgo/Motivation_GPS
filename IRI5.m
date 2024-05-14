clear all
close all
% filename = '00397_18122301';
% inputf = strcat(filename,'.txt');
% inputfile = strcat(filename,'_db','.txt');
% outputfile = strcat(filename,'_db','.xls');
% databasemaker2_2M(inputf)
% EgoGEN(inputfile,outputfile)

gg1 = ls('E:\GPS\cal\Motivation\Motivation_SR\ABET\*csv');
hh1 = cellstr(gg1);
number1 = length(hh1);
ff1 = hh1;
for ii=1:number1  %%%session number
    fnumcsv =[ff1{ii}];
    index = findstr(fnumcsv,'IRT');
    if index
        IRT = IRT_Read(fnumcsv);
        count = IRT(10:end,1);
        nosepoke = IRT(10:end,2);
        reward = IRT(10:end,3);
        
        M = length(count);
        nosepokeon = zeros(M,1);
        rewardon = zeros(M,1);
        ton(1) = nosepoke(1);
        toff(1) = reward(1);
        % ton(1,1) = arrayIRI(1,2);
        % toff(1,1) = arrayIRI(1,3);
        for i = 2:M
            nosepokeon(i) = nosepokeon(i-1)+nosepoke(i);
            rewardon(i) = rewardon(i-1)+reward(i);
        end
        
        BinWidth = 200;
        %  histogram(ton,'BinWidth',20);   %3600 second 3600/20=180
        [N1,nedges] = histcounts(nosepokeon,'BinWidth',BinWidth);
        [N2,nedges] = histcounts(rewardon,'BinWidth',BinWidth);
        
        
        ggg = ls('*mat');
        hhh = cellstr(ggg);
        number = length(hhh);
        fff = hhh;
        
        
        
        
        for kk=1:number
            fnumspk =[fff{kk}];
            indexspk = strcmp(fnumcsv(1:end-8),fnumspk(1:end-13));
            if indexspk
                
                load(fnumspk);
                figure(1)
                cellTS = ego.cellTS;
                [S,sedges] = histcounts(cellTS,'BinWidth',BinWidth);
                
                numS = length(S);
                S=S/max(S);
                ts = 0:BinWidth:BinWidth*(numS-1);
                
                
                plot(ts,S,'LineWidth',1.5)
                % legend('Normalized Firing Rate')
                
                hold on
                
                num = length(N1); % nose poke length
                if num <= numS
                    NN = zeros(numS); %add zero to nose poke
                    NN(1:num) = N1;
                    t = ts;
                else
                    %                     NN = N1;
                    %                     t = 0:BinWidth:BinWidth*(num-1);
                    
                    NN = zeros(numS);
                    NN(1:numS) = N1(1:numS);
                    %                     t = 0:BinWidth:BinWidth*(numS-1);
                    t = ts;
                end
                
                NN=NN/max(NN);
                plot(t,NN,'LineWidth',1.5)
                hold on
                corrValue = corrcoef(NN,S);
                CORR(1) = corrValue(1,2);
                %                 b=['Correlation = ',num2str(CORR)];
                %                 text(2400,0.75, b,'FontSize',12,'Color','red')
                %                 title(b)
                % filename = strcat(fnum(1:end-4),'_all');
                
                
                %%% Time of Reward
                gg = ls('*csv');
                hh = cellstr(gg);
                number = length(hh);
                ff = hh;
                jj = 1;
                for ii=1:number
                    fnum =[ff{ii}];
                    index = findstr(fnum,'IRT');
                    if isempty(index)
                        rawdata = RAW_Read(fnum);
                        pump = rawdata(:,4);
                        time = rawdata(:,1);
                        for ii = 1:length(pump)
                            if strcmp(pump{ii},'Pump #1')
                                rewardtime(jj) = time{ii};
                                jj = jj+1;
                            end
                        end
                    end
                end
                
                amp1 = 0.2;
                amp2 = 0.3;
                yy = [amp1 amp2];
                %                 hold on
                [N3,nedges] = histcounts(rewardtime,'BinWidth',BinWidth);
                
                num = length(N3);
                if num <= numS
                    NN2 = zeros(numS);
                    NN2(1:num) = N3;
                    t = ts;
                else
                    %                     NN2 = N3;
                    %                     t = 0:BinWidth:BinWidth*(num-1);
                    
                    
                    NN2 = zeros(numS);
                    NN2(1:numS) = N3(1:numS);
                    %                     t = 0:BinWidth:BinWidth*(numS-1);
                    t = ts;
                end
                
                
                
                NN2=NN2/max(NN2);
                plot(t,NN2,'LineWidth',1.5)
                
                corrValue = corrcoef(NN2,S);
                CORR(2) = corrValue(1,2);
                b=['Nose-Reward Correlation = ',num2str(CORR)];
                %                 text(2400,0.75, b,'FontSize',12,'Color','red')
                title(b)
                % %                 yy = 0.7*ones(length(rewardtime),1);
                %                 for hh = 1:length(rewardtime)
                %                     hold on
                %                     line([rewardtime(hh)' rewardtime(hh)'],yy,'Color','black');
                %                 end
                % %                     amp = amp+0.05;
                % %                 end
                legend('Firing Rate','Nosepoke Rate','Reward Rate')
                
                filename = strcat(fnumspk(1:end-4),'_Nosepoke');
                imageStore4b(figure(1),2,filename,300);
                close
            end
        end
    end
end