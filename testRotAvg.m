function testRotAvg()
a = load('G:\matlab\data\direct\gt\D2_011\4\tbc\rotavg_poses_single_0.txt');

Data = a(:,6:9);
poseMat = [];
poseMat = [];
for i = 1 : size(Data,1)
    data = Data(i,:);
    xyzw = data(1:4);
    trans = [0 0 10 * i];
    
    %     R2 = quatern2rotMat(xyzw([4 1 2 3]));
    R = quat2rotm(xyzw([4 1 2 3]));
    poseMat = [poseMat; [reshape(R,1,9), trans]];
    Rwc{i,1} = R;
    %     rotm2quat(R)
    %     rotMat2quatern(R2)
end
figure,plotPath(poseMat(1:1:end,:), 0.2);

id = [3 : 10 : 100]';

for i = 1 : length(id)-3
id1 = id(i);
id2 = id(i + 1);
id3 = id(i + 2);

Rwc1 = Rwc{id1};
Rwc2 = Rwc{id2};
Rwc3 = Rwc{id3};

Rc1c2 = inv(Rwc1) * Rwc2;
Rc2c3 = inv(Rwc2) * Rwc3;
Rwc3_check = Rwc1 * Rc1c2 * Rc2c3;

err = rad2deg(norm(rodrigues(Rwc3_check' * Rwc3)));

    
    
end
end