function [featureCount, data, allFeatureFlag]=ClassificationVoicedUnvoiced(allData, allFeatureFlag)
%Function Summary: Identify voiced and unvoiced frames and mark them as
%voiced an unvoiced
%   Detailed explanation goes here

% open .scp file`
file = fopen('/home/jitu/btp/data/diarization_data/gt.scp/IS1001d.scp', 'r');
tline = fgetl(file);
i = 1;
while ischar(tline)
    %disp(tline);    
    line = regexp(tline, '[_=]', 'split');
    %disp(line(2));
    %disp(line(3));
    for j = str2double(line(2)):str2double(line(3))
        allFeatureFlag(j, 1) = 1;
        allFeatureFlag(j, 2) = j;
    end
    tline = fgetl(file);    
    i = i+1;
end

siz = size(allData);
data = zeros(siz(1, 1), 19);
i = 1;
for j=1:siz(1,1)
    if allFeatureFlag(j, 1) == 1
        data(i, :) = allData(j, :);
        allFeatureFlag(j, 2) = i;
        i = i+1;
        disp(i);
    end            
end
featureCount = i-1;
end

