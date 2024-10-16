function delta_ang_trans = AnalyseVI(inputDir)



if 0
    
    delta_ang_trans1_1 = AnalyseVI('G:\matlab\data\traj\12_ext\yvr\4');delta_ang_trans1 = delta_ang_trans1_1;
    delta_ang_trans2_1 = AnalyseVI('G:\matlab\data\traj\12_ext\qvr\4');delta_ang_trans2 = delta_ang_trans2_1;
    
    close all
    delta_ang_trans1_1 = AnalyseVI('\\192.168.9.225\original_hand_dataset\tmp\1\data');delta_ang_trans1 = delta_ang_trans1_1;
    delta_ang_trans2_1 = AnalyseVI('\\192.168.9.225\original_hand_dataset\tmp\2\data');delta_ang_trans2 = delta_ang_trans2_1;
    mean(delta_ang_trans1)-mean(delta_ang_trans2)
    median(delta_ang_trans1)-median(delta_ang_trans2)
    num = min([size(delta_ang_trans1,1) size(delta_ang_trans2,1)]);
    figure,hist([delta_ang_trans1(1:num,1) delta_ang_trans2(1:num,1)],100); title('rotation error (deg)'); legend('yvr','qualcomm');grid on;
    figure,hist([delta_ang_trans1(1:num,2) delta_ang_trans2(1:num,2)],100); title('translation error (m)'); legend('yvr','qualcomm');grid on;
    figure,hist([delta_ang_trans1(1:num,3) delta_ang_trans2(1:num,3)],100); title('reprojection error (pixel)'); legend('yvr','qualcomm');grid on;
    figure,hist([delta_ang_trans1(1:num,1) - delta_ang_trans2(1:num,1)], 100); title('rotation error (deg)'); legend('yvr - qualcomm');grid on;
    figure,hist([delta_ang_trans1(1:num,2) - delta_ang_trans2(1:num,2)], 100); title('translation error (m)'); legend('yvr - qualcomm');grid on;
    figure,hist([delta_ang_trans1(1:num,3) - delta_ang_trans2(1:num,3)], 100); title('reprojection error (pixel)'); legend('yvr - qualcomm');grid on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    delta_ang_trans1_1 = AnalyseVI('\\192.168.9.225\original_hand_dataset\tmp\device_calibration_ucmrtp_pre\4\data');delta_ang_trans1 = delta_ang_trans1_1;
    delta_ang_trans2_1 = AnalyseVI('\\192.168.9.225\original_hand_dataset\tmp\device_calibration_ucmrtp_new\4\data');delta_ang_trans2 = delta_ang_trans2_1;
    mean(delta_ang_trans1)-mean(delta_ang_trans2)
    median(delta_ang_trans1)-median(delta_ang_trans2)
    num = min([size(delta_ang_trans1,1) size(delta_ang_trans2,1)]);
    figure,hist([delta_ang_trans1(1:num,1) delta_ang_trans2(1:num,1)],100); title('rotation error (deg)'); legend('old','new');grid on;
    figure,hist([delta_ang_trans1(1:num,2) delta_ang_trans2(1:num,2)],100); title('translation error (m)'); legend('old','new');grid on;
    figure,hist([delta_ang_trans1(1:num,3) delta_ang_trans2(1:num,3)],100); title('reprojection error (pixel)'); legend('old','new');grid on;
    figure,hist([delta_ang_trans1(1:num,1) - delta_ang_trans2(1:num,1)], 100); title('rotation error (deg)'); legend('old - new');grid on;
    figure,hist([delta_ang_trans1(1:num,2) - delta_ang_trans2(1:num,2)], 100); title('translation error (m)'); legend('old - new');grid on;
    figure,hist([delta_ang_trans1(1:num,3) - delta_ang_trans2(1:num,3)], 100); title('reprojection error (pixel)'); legend('old - new');grid on;
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    delta_ang_trans1 = AnalyseVI;
    delta_ang_trans2 = AnalyseVI;
    delta_ang_trans3 = AnalyseVI;
    
    figure,hist([delta_ang_trans1(:,1) delta_ang_trans2(:,1)],100); title('rotation error (deg)'); legend('yvr','qualcomm');
    figure,hist([delta_ang_trans1(:,2) delta_ang_trans2(:,2)],100); title('translation error (m)'); legend('yvr','qualcomm');
    figure,hist([delta_ang_trans1(:,3) delta_ang_trans2(:,3)],100); title('reprojection error (pixel)'); legend('yvr','qualcomm');
    
    mean(delta_ang_trans1)-mean(delta_ang_trans2)
    median(delta_ang_trans1)-median(delta_ang_trans2)
    
    
    acc = load('G:\matlab\data\traj\12_ext\acc_file.txt');
    gyro = load('G:\matlab\data\traj\12_ext\gyro_file.txt');
    a = load('G:\matlab\data\traj\12_ext\reProject_error.txt');
    figure,plot(a(:,1), a(:,2),'.r');axis equal;grid on;
    figure,plot([acc(:,2)  acc(:,5)])
    figure,plot([acc(:,2+1)  acc(:,5+1)])
    figure,plot([acc(:,2+2)  acc(:,5+2)])
    
    
    figure,plot(acc(:,2:4));
    figure,plot(gyro(:,2:4));
    
    mean(a)
    [~,normerr] = NormalizeVector(a);
    figure,hist(normerr,1000)
    
    
    
    c = load('G:\matlab\data\traj\12_ext\accelerometer.xml');
    figure,plot(c(:,1:3));
    
end



if 0
    inputDir = 'G:\matlab\data\traj\1';
    inputDir = 'G:\matlab\func\data';
    inputDir = 'G:\matlab\data\traj\2';
    inputDir = 'G:\matlab\data\traj\3';
    inputDir = 'G:\matlab\data\traj\4';
    inputDir = 'G:\matlab\data\traj\5';
    inputDir = 'G:\matlab\data\traj\6';
    inputDir = 'G:\matlab\data\traj\7_real';
    inputDir = 'G:\matlab\data\traj\8_ext';
    inputDir = 'G:\matlab\data\traj\9_ext';
    inputDir = 'G:\matlab\data\traj\10_ext';
    inputDir = 'G:\matlab\data\traj\11_ext';
    inputDir = 'G:\matlab\data\traj\12_ext';
    inputDir = 'G:\matlab\data\traj\yvr_jitter_1';
    inputDir = 'G:\matlab\data\traj\qualcomm_jitter_1';
    
    inputDir = 'G:\matlab\data\traj\yvr_smooth_onecam_1';
    inputDir = 'G:\matlab\data\traj\qualcomm_smooth_onecam_1';
    
    inputDir = 'G:\matlab\data\traj\yvr_smooth_twocam_1';
    inputDir = 'G:\matlab\data\traj\qualcomm_smooth_twocam_1';
    
    inputDir = 'G:\matlab\data\traj\yvr_jitter_four_cam_new_1';
    inputDir = 'G:\matlab\data\traj\qualcomm_jitter_four_cam_new_1';
    inputDir = 'G:\matlab\data\traj\yvr_jitter_four_cam_new_-td_1';
    inputDir = 'G:\matlab\data\traj\qualcomm_jitter_four_cam_new_-td_1';
    
    inputDir = 'G:\matlab\data\traj\y_jitter_4_-td_200_1';
    inputDir = 'G:\matlab\data\traj\q_jitter_4_-td_200_1';
    
    inputDir = 'G:\matlab\data\traj\12_ext';
end


inputCamPathFile = fullfile(inputDir,'camerapose_trajectory.csv');
inputImuPathFile =  fullfile(inputDir,'camimucalib_trajectory.csv');
inputExtrFile = fullfile(inputDir,'camimu_calib_extrinsic.csv');
inputRepFile = fullfile(inputDir, 'reprojection_error.csv');


data_camera_traj = csvread(inputCamPathFile);
data_camera_traj_imu_smoothed = csvread(inputImuPathFile);
vi = csvread(inputExtrFile);
rep = csvread(inputRepFile);

rep = rep(2:end);

I_T_C = ([ 0.173801    0.973502    0.148619   0.0947449;
    -0.767464    0.228463   -0.599003  -0.0184703;
    -0.617084 -0.00995257    0.786834   0.0764879;
    0           0           0           1]);


% I_T_C = [ -0.0228711  -0.0223319    0.999489   0.0819141
%    0.999735  0.00196504   0.0229206     0.13145
% -0.00247589    0.999749   0.0222811  -0.0715627
%           0           0           0           1];


I_T_C = [ 0.206576  0.960373  0.187111 0.0961678;
    -0.750003  0.278238 -0.600066 -0.004979;
    -0.628348 -0.016375   0.77776 0.0417894;
    0         0         0         1];

G_T_I0 = [0.963576    0.264303     -0.0408           0;
    2.77556e-17   -0.152561   -0.988294           0;
    -0.267434    0.952297   -0.147004           0;
    0           0           0           1];


firstId = 2;

camPath = [];
camPathComp = [];
camPath_ = [];
camPathComp_ = [];
delta_ang_trans = [];
for i = firstId : size(data_camera_traj,1)
    
    rot_cam =  quat_2_Rot(data_camera_traj(i,1:4)')';
    rot_imu =  quat_2_Rot(data_camera_traj_imu_smoothed(i,1:4)')';
    trans_cam = data_camera_traj(i,5:7);
    trans_imu = data_camera_traj_imu_smoothed(i,5:7);
    
    rot_vi =  quat_2_Rot(vi(i,1:4)')';
    trans_vi = vi(i,5:7);
    
    I_T_C = [rot_vi trans_vi'; 0 0 0 1];
    
    if i == firstId
        T_cam_base = [rot_cam trans_cam'; 0 0 0 1];
        T_imu_base = [rot_imu trans_imu'; 0 0 0 1];
    end
    
    
    T_cam_cur = [rot_cam trans_cam'; 0 0 0 1];
    T_imu_cur = [rot_imu trans_imu'; 0 0 0 1];
    
    
    %     T_imu_cur =inv(G_T_I0) * T_imu_cur;
    
    T_cam_w2c = inv(inv(T_cam_base)*T_cam_cur);
    T_imu_w2c = inv(inv(T_imu_base)*T_imu_cur);
    T_cam_w2c_comp = inv(I_T_C) * T_imu_w2c * I_T_C;
    
    comp = [T_cam_w2c T_cam_w2c_comp];
    
    deltaPose = inv(T_cam_w2c) * T_cam_w2c_comp;
    delta_ang_trans = [delta_ang_trans; [rad2deg(norm(rodrigues(deltaPose(1:3,1:3)))) 1*norm(deltaPose(1:3,4))]];
    
    T_cam_c2w = inv(T_cam_w2c);
    T_cam_c2w_comp = inv(T_cam_w2c_comp);
    camPath = [camPath; reshape(T_cam_c2w(1:3,1:3),1,9) T_cam_c2w(1:3,4)'];
    camPathComp = [camPathComp; reshape(T_cam_c2w_comp(1:3,1:3),1,9) T_cam_c2w_comp(1:3,4)'];
    
    camPath_ = [camPath_; reshape(T_cam_cur(1:3,1:3),1,9) T_cam_cur(1:3,4)'];
    camPathComp_ = [camPathComp_; reshape(T_imu_cur(1:3,1:3),1,9) T_imu_cur(1:3,4)'];
end
delta_ang_trans = [delta_ang_trans rep];
figure,subplot(2,2,2);plot(delta_ang_trans);legend('deg','m','pixel');subplot(2,2,3),plot(vi(:,5:7));legend('x','y','z');
subplot(2,2,4);hist(delta_ang_trans, 100);legend('deg','m','pixel'); subplot(2,2,1);plotPath(camPath);hold on;plotPath(camPathComp);

% figure,plotPath(camPath_);hold on;plotPath(camPathComp_);








end
function Rot =  quat_2_Rot(q)
q_x = SkewSymMat(q(1:3));
Rot = (2 * q(4)^2 - 1) * eye(3) - 2 * q(4) * q_x +  2 * q(1:3) * q(1:3)';
end