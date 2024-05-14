%% Check if the data have both .mat and .xsl files.
%% Rui 2024.04.14


clear all
close all

addpath D:\GPS\cal\pippin-master\MEC_Final\xcorr


%% Parameter Setup
ManifoldData={};
ManifoldDataNumber=1;
BinWidth = 50; %bin 50ms
FR_Fixed_Count=3;
binwidth=0.1;

ExcelPath='E:\GPS\cal\Motivation\Motivation_gps\Code\ABET\';
DataPath='E:\GPS\cal\Motivation\Motivation_gps\Data\';


%% Data Load

ExceLgg1 = ls([ExcelPath '*csv']);
ExceLhh1 = cellstr(ExceLgg1);
ExceLnumber1 = length(ExceLhh1);
ExceLff1 = ExceLhh1;

kkk=1;

for ii=1:ExceLnumber1  %%%session number
    fnumIRT =[ExceLff1{ii}];

    if ~isempty(findstr(fnumIRT,'IRT')), continue, end

%     try
%         %         fnumIRT='00397_08092301.csv';
%         [timeSeriesData,Event]=preData_event([ExcelPath fnumIRT],FR_Fixed_Count);
% 
%     catch
%         continue
%         warning(['cannot read ' ExcelPath fnumIRT ])
%     end

%     if any(contains(fieldnames(timeSeriesData), 'FR_Count'))
%         continue
%     end

    
    Egoggg = ls([DataPath fnumIRT(1:end-8) '*mat']);
    Egohhh = cellstr(Egoggg);
    Egonumber = length(Egohhh);
    Egofff = Egohhh;

    if ~isempty(Egoggg)
        for kk=1:Egonumber
            close all
            fnumspk =[Egofff{kk}];

            Data_ego{kkk}=fnumspk;
            kkk=kkk+1;
            Egoname=fnumspk(1:19);


% % 取出a中有而b中没有的元素，即a的补集
% complement_a = setdiff(Egohhh, Data_ego);


        end
    end



end


Allegodata = ls([DataPath '*mat']);
complement_a = setdiff(Allegodata, Data_ego)











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