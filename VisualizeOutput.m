function VisualizeOutput()
% close all;

data = load('G:\matlab\data\direct\gt\D2_011\4\tbc\output.txt');

data = data(:,2:8);

poseMat = [];
Twc_stack = {};
for i = 1 : size(data,1)
    data1 = data(i,:);
    xyzw = data1(4:7);
    trans = data1(1:3);
    R = quat2rotm(xyzw([4 1 2 3]));
    poseMat = [poseMat; [reshape(R,1,9), trans]];
    Twc_stack{i,1} = [R trans';0 0 0 1];
end


figure,plotPath(poseMat(1:1:end,:));
figure,plot(poseMat(:,10:12))
norm(poseMat(486,10:12) - poseMat(908,10:12)) * 1000 - 6000
norm(poseMat(676,10:12) - poseMat(1111,10:12)) * 1000 - 6000
end