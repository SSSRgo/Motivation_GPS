%% �����ı��ļ��е����ݡ�
% ���ڴ������ı��ļ��������ݵĽű�:
%
%    C:\Users\Sheng-Jia\Desktop\IRI\00397_28072301.csv
%
% Ҫ��������չ������ѡ�����ݻ������ı��ļ��������ɺ���������ű���

% �� MATLAB �Զ������� 2023/08/04 12:14:29
function rawdata = RAW_Read(filename)
%% ��ʼ��������
filename = filename;

% filename = '00397_28072301.csv';  
delimiter = ',';

%% ����������Ϊ�ı���ȡ:
% �й���ϸ��Ϣ������� TEXTSCAN �ĵ���
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% ���ı��ļ���
fileID = fopen(filename,'r');

%% ���ݸ�ʽ��ȡ�����С�
% �õ��û������ɴ˴������õ��ļ��Ľṹ����������ļ����ִ����볢��ͨ�����빤���������ɴ��롣
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

%% �ر��ı��ļ���
fclose(fileID);

%% ��������ֵ�ı���������ת��Ϊ��ֵ��
% ������ֵ�ı��滻Ϊ NaN��
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,6,7,9,11,13,15,17]
    % ������Ԫ�������е��ı�ת��Ϊ��ֵ���ѽ�����ֵ�ı��滻Ϊ NaN��
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % ����������ʽ�Լ�Ⲣɾ������ֵǰ׺�ͺ�׺��
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % �ڷ�ǧλλ���м�⵽���š�
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % ����ֵ�ı�ת��Ϊ��ֵ��
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


%% �����ݲ��Ϊ��ֵ���ַ����С�
rawNumericColumns = raw(:, [1,2,6,7,9,11,13,15,17]);
rawStringColumns = string(raw(:, [3,4,5,8,10,12,14,16]));


%% ������ֵԪ���滻Ϊ NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % ���ҷ���ֵԪ��
rawNumericColumns(R) = {NaN}; % �滻����ֵԪ��

%% ȷ������ <undefined> ���κ��ı�������ȷת��Ϊ <undefined> ����ֵ
for catIdx = [1,2,3,4]
    idx = (rawStringColumns(:, catIdx) == "<undefined>");
    rawStringColumns(idx, catIdx) = "";
end

%% �����������
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
%% �����ʱ����
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R catIdx idx;