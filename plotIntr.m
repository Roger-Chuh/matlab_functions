function plotIntr()
inputDir = 'G:\matlab\data\direct\gt\D2_011\4\intr';

paraInfo = dir(fullfile(inputDir,'1*'));

% paraInfo = paraInfo(1:3);
% paraInfo = paraInfo(4:6);
% paraInfo = paraInfo(7:10);
% paraInfo = paraInfo([7 10]);
% paraInfo = paraInfo([8 9]);
% paraInfo = paraInfo([8]);
% paraInfo = paraInfo([3 6 8]);
% paraInfo = paraInfo([11 12 13]);
% paraInfo = paraInfo([1 4]);
% paraInfo = paraInfo([6 8]);

% paraInfo = paraInfo([16 17]);
% paraInfo = paraInfo([9 18]);
paraInfo = paraInfo([1 19]);

AngCell = cell(length(paraInfo), 10);
LineCell = cell(length(paraInfo), 10);
for para_id = 1 : length(paraInfo)
    angInfo = dir(fullfile(inputDir, paraInfo(para_id).name,'camera_*_angle.txt'));
    lineInfo = dir(fullfile(inputDir, paraInfo(para_id).name, 'camera_*_line.txt'));
    
    cam_num = length(angInfo);
    for i = 1 : cam_num
        a = load(fullfile(inputDir, paraInfo(para_id).name, angInfo(i).name));
        AngCell{para_id, i} = a;
    end
    
    for j = 1 : cam_num
        b = load(fullfile(inputDir, paraInfo(para_id).name, lineInfo(j).name));
        LineCell{para_id, j} = b;
    end
    
end

AngCell = AngCell(:,1:cam_num);
LineCell = LineCell(:,1:cam_num);


figure(20),clf;
for i = 1 : cam_num
    subplot(2,2,i);hold on;
    name = '';
    Names = cell(size(AngCell(:,i),1),1);
    for j = 1:size(AngCell(:,i),1)
%         histogram(AngCell{j,i},200, 'Normalization','probability','BinWidth',0.01);grid on;
%         histogram(AngCell{j,i},'Normalization','probability','BinWidth',0.008);grid on;
        histogram(AngCell{j,i},'BinWidth',0.005);grid on;
        if j == 1
            name = paraInfo(j).name;
        else
            name = strcat(name,strcat(',',paraInfo(j).name));
        end
        %         legend(paraInfo(j).name);
        Names{j,1} = paraInfo(j).name;
    end
    %     legend(name);
    legend(Names,'Location','NorthWest');
    title(sprintf('camera %d',i-1));
end


figure(21),clf;
for i = 1 : cam_num
    subplot(2,2,i);hold on;
    name = '';
    Names = cell(size(LineCell(:,i),1),1);
    for j = 1:size(LineCell(:,i),1)
%         histogram(AngCell{j,i},200, 'Normalization','probability','BinWidth',0.01);grid on;
%         histogram(LineCell{j,i},'Normalization','pdf','BinWidth',0.008);grid on;
        histogram(LineCell{j,i},'BinWidth',0.005);grid on;
        if j == 1
            name = paraInfo(j).name;
        else
            name = strcat(name,strcat(',',paraInfo(j).name));
        end
        %         legend(paraInfo(j).name);
        Names{j,1} = paraInfo(j).name;
    end
    %     legend(name);
    legend(Names,'Location','NorthEast');
    title(sprintf('camera %d',i-1));
end

end