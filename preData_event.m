
function [timeSeriesData,Event]=preData_event(file,FR_Fixed_Count,varargin)

% params = inputParser;
% params.parse(varargin{:});

if any(cellfun(@(x) strcmp(x,'state_compress'),varargin))
state_compress=varargin{cellfun(@(x) strcmp(x,'state_compress'),varargin),2};
end


opts = delimitedTextImportOptions("NumVariables", 17);

% 指定范围和分隔符
opts.DataLines = [18, Inf];
opts.Delimiter = ",";

% 指定列名称和类型
opts.VariableNames = ["Evnt_Time", "Evnt_ID", "Event_Name", "Item_Name", "Alias_Name", "Group_ID", "Num_Args", "Arg1_Name", "Arg1_Value", "Arg2_Name", "Arg2_Value", "Arg3_Name", "Arg3_Value", "Arg4_Name", "Arg4_Value", "Arg5_Name", "Arg5_Value"];
opts.VariableTypes = ["double", "double", "char", "char", "char", "double", "double", "char", "double", "char", "double", "char", "double", "char", "double", "char", "double"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 指定变量属性
opts = setvaropts(opts, ["Event_Name", "Item_Name", "Alias_Name", "Arg1_Name", "Arg2_Name", "Arg3_Name", "Arg4_Name", "Arg5_Name"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Event_Name", "Item_Name", "Alias_Name", "Arg1_Name", "Arg2_Name", "Arg3_Name", "Arg4_Name", "Arg5_Name"], "EmptyFieldRule", "auto");

% 导入数据
Untitled = readtable(file, opts);


%% Parameter
if ~exist("FR_Fixed_Count")
    FR_Fixed_Count=3; end
% PR_threshold=40;
% FR_Fixed_Count=3;
FR_Fixed_Ratio_Count=1:FR_Fixed_Count;

%% Data

% % Define the file path
% filePath = 'E:\GPS\cal\Motivation\Motivation_gps\Code\ABET\00394_01082301.csv';
%
% % Read the CSV file
% dataTable = readtable(filePath);

Item_Name=Untitled.Item_Name;
Event_Name=Untitled.Event_Name;
% Extract unique event names
eventNames = unique(Untitled.Item_Name);

% Initialize a structure to hold time series data
timeSeriesData = struct();
Arg2_ValueSeriesData= struct();
parfor i = 1:length(Item_Name)
    eventName = Untitled.Item_Name{i};
    eventName = strrep(eventName, ' ', '_'); % 将空格替换为下划线
    eventName = regexprep(eventName, '[^\w\s]', ''); % 移除特殊字符
    Item_Name{i}=eventName;
end
Untitled.Item_Name=Item_Name;

% Loop through each event name to extract its time series
for i = 1:length(eventNames)
    eventName = eventNames{i};
    eventName = strrep(eventName, ' ', '_'); % 将空格替换为下划线
    eventName = regexprep(eventName, '[^\w\s]', ''); % 移除特殊字符
    % Find rows corresponding to the current event
    if strcmp(eventName,'Tray_1') || strcmp(eventName,'Nose_Poke_1')
        eventRows = Untitled(strcmp(Event_Name, 'Input Transition On Event')...
            &strcmp(Item_Name, eventName), :);
    else
        eventRows = Untitled(strcmp(Item_Name, eventName), :);
    end



    % Extract event times and store in the structure
    try
        timeSeriesData.(eventName) = eventRows.Evnt_Time;
        Arg2_ValueSeriesData.(eventName)=eventRows.Arg2_Value;
    catch
        continue
    end
end

% Display the time series for each event
% disp(timeSeriesData);

% timeSeriesData.Tray_1
% timeSeriesData.Pump_1
%% Divide different types fo Tray on 


a_sequence=timeSeriesData.Tray_1;
b_sequence=timeSeriesData.Pump_1;
c_sequence=timeSeriesData.Nose_Poke_1;

Event.Trayon_pump = [];
for i = 1:length(a_sequence)
    a_time = a_sequence(i);

    % 找到上一个 b 发生的时间点
    last_b_idx = find(b_sequence < a_time, 1, 'last');

    if ~isempty(last_b_idx)
        last_b_time = b_sequence(last_b_idx);

        % 检查在 a 和上一个 b 之间是否有 c 发生
        if isempty(find(c_sequence > last_b_time & c_sequence < a_time, 1))
            Event.Trayon_pump = [Event.Trayon_pump, a_time];
        end
    end
end
Event.Trayon_pump=Event.Trayon_pump';
Event.Trayon_nonpump = setdiff(a_sequence, Event.Trayon_pump);

%%%%%%%%%%%%%%%%%%% PR %%%%%%%%%%%%%%%%%%%

%% Divide different types fo Nose poke FR or PR
%%%%%%%%%%%%%%%%%%% FR %%%%%%%%%%%%%%%%%%%
if any(contains(eventNames, 'FR Count'))
    a_sequence=timeSeriesData.Pump_1;
    b_sequence=timeSeriesData.FR_Count;
    c_sequence=timeSeriesData.Nose_Poke_1;

FR_Count_sequence=zeros(size(timeSeriesData.FR_Count));
FR_State_sequence=zeros(size(timeSeriesData.FR_Count));

state_number=1;

       state_compress_num=1;
for i = 1:length(a_sequence)
    a_time = a_sequence(i);

    % 找到上一个 b 发生的时间点
    last_b_idx = find(b_sequence <= a_time, FR_Fixed_Count, 'last');
    if ~isempty(last_b_idx)
        last_b_time = b_sequence(last_b_idx);
        FR_Count_sequence(last_b_idx)=FR_Fixed_Ratio_Count;
        FR_State_sequence(last_b_idx)=state_number;

        
        if mod(state_compress_num,state_compress)==0
        state_number=state_number+1;
        end
        state_compress_num=state_compress_num+1;
    end

end
% FR_Count_np_idx=FR_Count_sequence(discretize(timeSeriesData.Nose_Poke_1,timeSeriesData.FR_Count));

nosepoke_FR_count=struct();
for i=FR_Fixed_Ratio_Count
  Event.(['Nosepoke_Rcount_' num2str(i)])=timeSeriesData.Nose_Poke_1(find(FR_Count_sequence==i));
end

for i=1:state_number
  Event.(['Nosepoke_FRstate_' num2str(i)])=timeSeriesData.Nose_Poke_1(find(FR_State_sequence==i));
end

end

Event.(['Nosepoke_FRstate_last'])=[];




%%%%%%%%%%%%%%%%%%% PR %%%%%%%%%%%%%%%%%%%
if any(contains(eventNames, 'PR_count'))

a=1;

    a_sequence=timeSeriesData.Pump_1(2:end);
    b_sequence=timeSeriesData.PR_Count;


PR_Count_sequence=zeros(size(timeSeriesData.PR_Count));

for i = 1:length(a_sequence)
    a_time = a_sequence(i);

    % 找到上一个 b 发生的时间点
    last_b_idx = find(b_sequence <= a_time, Arg2_ValueSeriesData.PR_Value(end-length(a_sequence)+i-1), 'last');
    if ~isempty(last_b_idx)
        last_b_time = b_sequence(last_b_idx);
        PR_Count_sequence(last_b_idx)=i;
    end

end

for i=1:length(a_sequence)
  Event.(['Nosepoke_PRstate_' num2str(i)])=b_sequence(find(PR_Count_sequence==i));
end

last_b_idx = find(b_sequence > a_sequence(end));
PR_Count_sequence(last_b_idx)=0;
Event.(['Nosepoke_PRstate_last'])=b_sequence(find(PR_Count_sequence==0));

end




%% Divide half of np trayon
%%%%%%%%%%%%%%%%%%% FR %%%%%%%%%%%%%%%%%%%
if any(contains(eventNames, 'FR Count'))
Event.Nosepoke_slice1=timeSeriesData.FR_Count(1:round(length(timeSeriesData.FR_Count)/2));
Event.Nosepoke_slice2=timeSeriesData.FR_Count(round(length(timeSeriesData.FR_Count)/2)+1:end);
end

%%%%%%%%%%%%%%%%%%% PR %%%%%%%%%%%%%%%%%%%
if any(contains(eventNames, 'PR Count'))
Event.Nosepoke_slice1=timeSeriesData.PR_Count(1:round(length(timeSeriesData.PR_Count)/2));
Event.Nosepoke_slice2=timeSeriesData.PR_Count(round(length(timeSeriesData.PR_Count)/2)+1:end);

end

Event.Trayon_slice1=timeSeriesData.Tray_1(1:round(length(timeSeriesData.Tray_1)/2));
Event.Trayon_slice2=timeSeriesData.Tray_1(round(length(timeSeriesData.Tray_1)/2)+1:end);
%% Divide different types fo Nose poke



% disp(Trayon_pump);
% disp(Event.);
% Trayon_nonpump = setdiff(a_sequence, Trayon_pump);


%% Summarize Evnet

% Event.nosepoke_FR_count=nosepoke_FR_count;
% Event.Trayon_pump=Trayon_pump;
% Event.Trayon_nonpump=Trayon_pump;

%% 转换为输出类型

% Untitled = table2cell(Untitled);
% numIdx = cellfun(@(x) ~isnan(str2double(x)), Untitled);
% Untitled(numIdx) = cellfun(@(x) {str2double(x)}, Untitled(numIdx))

%% 清除临时变量
 
% clear numIdx opts