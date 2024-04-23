% function [Trans, trans_cat, Rot_diff, rot_cat] = testOnlineCalibration(inputDir, draw)
function [Trans_diff, trans_cat, Rot_diff, rot_cat] = testOnlineCalibration(inputDir, draw)
% function [Trans, trans_cat, Rot, rot_cat] = testOnlineCalibration(inputDir, draw)


% inputDir = 'G:\matlab\data\vio\with_epiplane\oc_time_010';
% inputDir = 'G:\matlab\data\vio\without_epiplane\oc_time_010';

dirInfo = dir(fullfile(inputDir,'online_calibration_time_consume*.txt'));
a = load(fullfile(inputDir,dirInfo(1).name));
% figure,subplot(1,2,1);plot(a(:,2));grid on;subplot(1,2,2),hist(a(:,2),50);title(sprintf('mean time: %f ms', mean(a(a(:,2)>1,2))));grid on;

 fprintf(sprintf('mean consume time: %f ms\n',mean(a(a(:,2)>1,2))));

Trans = [];
Rot = [];
for i = 1 : size(a,1)
   ex = a(i,3:end);
   Ex_ = reshape(ex, 4,3,4);
   Ex = permute(Ex_, [2 1 3]);
   trans = permute(Ex(:,4,:), [2 1 3]);
   rot_ =[];
   for j = 1 : 4
       rot_(:,j) = rad2deg(rodrigues(Ex(1:3,1:3,j)));
   end
   rot = reshape(rot_,1,3,4);
   Trans = [Trans; trans];
   Rot = [Rot; rot];
end

a = permute(Trans,[1 3 2]);
trans_cat = reshape(a,[],3);

b = permute(Rot,[1 3 2]);
rot_cat = reshape(b,[],3);

Trans_diff = 1000.*diff(Trans);
Rot_diff = rad2deg(diff(Rot));

if draw
    if 1
        figure,subplot(2,4,1);plot(Trans(:,:,1));title('cam0 trans');grid on;
        subplot(2,4,2);plot(Trans(:,:,2));title('cam1 trans');grid on;
        subplot(2,4,3);plot(Trans(:,:,3));title('cam2 trans');grid on;
        subplot(2,4,4);plot(Trans(:,:,4));title('cam3 trans');grid on;
        subplot(2,4,5);plot(Rot(:,:,1));title('cam0 rot');grid on;
        subplot(2,4,6);plot(Rot(:,:,2));title('cam1 rot');grid on;
        subplot(2,4,7);plot(Rot(:,:,3));title('cam2 rot');grid on;
        subplot(2,4,8);plot(Rot(:,:,4));title('cam3 rot');grid on;
    else
        figure;subplot(1,2,1);plot(trans_cat);grid on;  title('trans');
               subplot(1,2,2);plot(rot_cat);grid on;  title('rot');
    end
end

end