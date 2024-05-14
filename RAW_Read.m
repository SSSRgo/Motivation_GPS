%% 导入文本文件中的数据。
% 用于从以下文本文件导入数据的脚本:
%
%    C:\Users\Sheng-Jia\Desktop\IRI\00397_28072301.csv
%
% 要将代码扩展到其他选定数据或其他文本文件，请生成函数来代替脚本。

% 由 MATLAB 自动生成于 2023/08/04 12:14:29
function rawdata = RAW_Read(filename)
%% 初始化变量。
filename = filename;

% filename = '00397_28072301.csv';  
delimiter = ',';

%% 将数据列作为文本读取:
% 有关详细信息，请参阅 TEXTSCAN 文档。
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% 打开文本文件。
fileID = fopen(filename,'r');

%% 根据格式读取数据列。
% 该调用基于生成此代码所用的文件的结构。如果其他文件出现错误，请尝试通过导入工具重新生成代码。
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

%% 关闭文本文件。
fclose(fileID);

%% 将包含数值文本的列内容转换为数值。
% 将非数值文本替换为 NaN。
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,6,7,9,11,13,15,17]
    % 将输入元胞数组中的文本转换为数值。已将非数值文本替换为 NaN。
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % 创建正则表达式以检测并删除非数值前缀和后缀。
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % 在非千位位置中检测到逗号。
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % 将数值文本转换为数值。
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% 将数据拆分为数值和字符串列。
rawNumericColumns = raw(:, [1,2,6,7,9,11,13,15,17]);
rawStringColumns = string(raw(:, [3,4,5,8,10,12,14,16]));


%% 将非数值元胞替换为 NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % 查找非数值元胞
rawNumericColumns(R) = {NaN}; % 替换非数值元胞

%% 确保包含 <undefined> 的任何文本都已正确转换为 <undefined> 分类值
for catIdx = [1,2,3,4]
    idx = (rawStringColumns(:, catIdx) == "<undefined>");
    rawStringColumns(idx, catIdx) = "";
end

%% 创建输出变量
Untitled = table;
Untitled.Evnt_Time = cell2mat(rawNumericColumns(:, 1));
Untitled.Evnt_ID = cell2mat(rawNumericColumns(:, 2));
Untitled.Event_Name = categorical(rawStringColumns(:, 1));
Untitled.Item_Name = categorical(rawStringColumns(:, 2));
Untitled.Alias_Name = categorical(rawStringColumns(:, 3));
Untitled.Group_ID = cell2mat(rawNumericColumns(:, 3));
Untitled.Num_Args = cell2mat(rawNumericColumns(:, 4));
Untitled.Arg1_Name = categorical(rawStringColumns(:, 4));
Untitled.Arg1_Value = cell2mat(rawNumericColumns(:, 5));
Untitled.Arg2_Name = rawStringColumns(:, 5);
Untitled.Arg2_Value = cell2mat(rawNumericColumns(:, 6));
Untitled.Arg3_Name = rawStringColumns(:, 6);
Untitled.Arg3_Value = cell2mat(rawNumericColumns(:, 7));
Untitled.Arg4_Name = rawStringColumns(:, 7);
Untitled.Arg4_Value = cell2mat(rawNumericColumns(:, 8));
Untitled.Arg5_Name = rawStringColumns(:, 8);
Untitled.Arg5_Value = cell2mat(rawNumericColumns(:, 9));

rawdata = raw;
%% 清除临时变量
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R catIdx idx;