function testScaleOpt()
global use_dist_error use_5dof pose_size opt_3dof add_noise right_perturbation include_first_pose offset fixed_pose_num use_reproj_factor use_epiplane_factor...
    fix_trans_norm use_xyz fix_rot

close all;

add_noise = true;
% Tbcx = [rodrigues(rand(3, 1)) rand(3, 1); 0 0 0 1];
% Tbcy = [rodrigues(rand(3, 1)) rand(3, 1); 0 0 0 1];
% Twb1 = [rodrigues(rand(3, 1)) rand(3, 1); 0 0 0 1];
% Twb2 = [rodrigues(rand(3, 1)) rand(3, 1); 0 0 0 1];
%
% n = rand(3,1);
% n = n./norm(n);
% n = sign(n(3)).*n;
% lambda = 0.6;
% scale = 1.3;
% Tb2b1 = inv(Twb2) * Twb1;
%
% Tb2b1_scaled = Tb2b1;
% Tb2b1_scaled(1:3,4) = scale * Tb2b1(1:3,4);
%
%
% Tth = inv(Tbcy) * Tb2b1_scaled * Tbcx;
%
%
% X1 = Tth(1:3,1:3) * n + lambda * Tth(1:3,4);
%
% R_th = Tbcy(1:3,1:3)' * Twb2(1:3,1:3)' * Twb1(1:3,1:3) * Tbcx(1:3,1:3);
% t_th = Tbcy(1:3,1:3)' * (Twb2(1:3,1:3)' * (Twb1(1:3,1:3) * Tbcx(1:3,4) + scale * (Twb1(1:3,4) - Twb2(1:3,4))) - Tbcy(1:3,4));
% X2 = R_th * n + lambda * t_th;
%
%
% X1 - X2
K = [300 0 320; 0 300 240; 0 0 1];

xyz = [10 200 1000;200 -20 1002;-10 500 1003;-100 -50 990;300 -50 990;-100 -250 990;-100 -150 990];

xyz = [xyz; xyz + [-10 20 3]; xyz + [10 50 3]; xyz + [-70 -20 3]; xyz + [-10 20 30]; xyz + [-50 15 30]; xyz + [-10 27 -25]; xyz + [75 -20 -30];xyz + [14 -9 31]];

xyz = xyz(1:61,:);

fixed_pose_num = 3;
pose_num = 10;

each_host = ceil(size(xyz,1) / pose_num);

host_ids = [];
target_ids = [];
for(i = 1: size(xyz,1))
    
    host_ids = [host_ids; mod(i, pose_num-1)+1];
    target_ids = [target_ids; mod(i, pose_num-1)+2];
    
end


cam_num = 3;

Tbc = {};
for i = 1 : cam_num
    if i == 1
        Tbc{i, 1} = eye(4);
    else
        Tbc{i, 1} = [rodrigues(0.1 * (rand(3,1) - 0.5)) 0.1 * (rand(3,1) - 0.5); 0 0 0 1];
    end
end


host_pids_each_fid = round(size(xyz,1) / pose_num);

pid_range = 1 : host_pids_each_fid : size(xyz,1);
pid_range = [pid_range size(xyz,1)];

pid_range = unique(pid_range);

fid_to_pid_range = {};
pid_to_host_fid = [];
for i = 1 : length(pid_range)-1
    if i ~= length(pid_range)-1
        fid_to_pid_range{i,1} = [pid_range(i) : pid_range(i+1)-1];
    else
        fid_to_pid_range{i,1} = [pid_range(i) : pid_range(i+1)];
    end
    pid_to_host_fid = [pid_to_host_fid;i * ones(length(fid_to_pid_range{i,1}),1)];
end

[~,depths] = NormalizeVector(xyz);

% idps = 1./depths;

idps_gt = {};
idps = {};
pose_wb{1,1} = eye(4);
trans_ang_err{1,1} = [0 0];
for i = 1 : pose_num %size(xyz,1)
    if i == 1
        pose_wb{i,1} = eye(4);
    else
        pose_wb{i,1} = [rodrigues(0.2*(rand(3,1)-0.5)) 0.5*(rand(3,1)-0.5); 0 0 0 1];
    end
    
    for j = 1 : cam_num
        pose_cw = inv(pose_wb{i,1} * Tbc{j});
        bearing = pose_cw(1:3,1:3)*xyz' + repmat(pose_cw(1:3,4), 1, size(xyz,1));
        [bearing_, dist] = NormalizeVector(bearing');
        bearings{j,i} = bearing_;
        if j == 1
            idps{1, i} = 1./dist + 0.0001;
            idps_gt{1, i} = 1./dist;
            if ~add_noise
                idps{1, i} = 1./dist;
                bearings{j,i} .* repmat(dist,1,3);
            end
        end
    end
end


opt_fids = fixed_pose_num+1:pose_num;%size(xyz,1);

pose_Twb_gt = pose_wb;
pose_Twb = pose_wb;
for k = 1 : fixed_pose_num
    pose_err{k,1} = eye(4);
    pose_Twb{k,1} = (pose_wb{k,1});
end
for i = opt_fids
    pose_err{i,1} = [rodrigues(0.1*(rand(3,1)-0.5)) 0.05*(rand(3,1)-0.5); 0 0 0 1];
    if add_noise
        pose_Twb{i,1} = (pose_wb{i,1}) * pose_err{i,1};
    else
        pose_Twb{i,1} = (pose_wb{i,1});
    end
end


[xMat, yMat] = meshgrid(1 : pose_num, 1 :pose_num);
pairs = [xMat(:) yMat(:)];

scale = 2;

pose_size = 6;
start_idp_idx = pose_size*(pose_num-fixed_pose_num);
start_scale_idx = pose_size*(pose_num-fixed_pose_num) + size(xyz,1);
loss_vec = [];
for iter = 1 : 30
    
    H = zeros(pose_size*(pose_num-fixed_pose_num) + size(xyz,1) + 1, pose_size*(pose_num-fixed_pose_num) + size(xyz,1) + 1);
    b = zeros(pose_size*(pose_num-fixed_pose_num) + size(xyz,1) + 1, 1);
    err_sum = 0;
    err_count = 0;
    for pid = 1 : size(xyz,1)
        pt = xyz(pid,:);
        host_fid = pid_to_host_fid(pid);
        host_idp = idps{1, host_fid}(pid);
        host_bearing = bearings{1,host_fid}(pid,:);
        pose_wc_host = pose_Twb{host_fid,1} * Tbc{1};
        for target_fid = 1 : pose_num
            if target_fid == host_fid
                continue;
            end
            for cid = 1 : cam_num
                if cid == 1
                    %                     continue;
                end
                target_bearing = bearings{cid, target_fid}(pid,:);
                pose_cw_target = inv(pose_Twb{target_fid,1} * Tbc{cid});
                Tth = pose_cw_target * pose_wc_host;
                reproj = Tth(1:3,1:3)*(host_bearing'/host_idp) +Tth(1:3,4);
                if ~add_noise
                    assert(norm(target_bearing' - reproj./norm(reproj)) < 0.000001);
                end
                [err, d_err_d_Twb1, d_err_d_Twb2, d_err_d_scale, d_err_d_rho] = computeScaledReprojFactor(pose_Twb{host_fid,1}, pose_Twb{target_fid,1}, Tbc{1}, Tbc{cid}, host_idp, scale, host_bearing', target_bearing');
                err_sum = err_sum + 235 * norm(err);
                err_count = err_count + 1;
                
                host_fid_valid = host_fid > fixed_pose_num;
                target_fid_valid = target_fid > fixed_pose_num;
                
                if host_fid_valid
                    H = FillMatrix(H, (host_fid - fixed_pose_num - 1) * pose_size + 1, (host_fid - fixed_pose_num - 1) * pose_size + 1, pose_size, pose_size, d_err_d_Twb1' * d_err_d_Twb1);
                    b = FillMatrix(b, (host_fid - fixed_pose_num - 1) * pose_size + 1, 1, pose_size, 1, -d_err_d_Twb1' * err);
                    
                    H = FillMatrix(H, (host_fid - fixed_pose_num - 1) * pose_size + 1, start_idp_idx + pid, pose_size, 1, d_err_d_Twb1' * d_err_d_rho);
                    H = FillMatrix(H, start_idp_idx + pid, (host_fid - fixed_pose_num - 1) * pose_size + 1, 1, pose_size, d_err_d_rho' * d_err_d_Twb1);
                    
                    H = FillMatrix(H, (host_fid - fixed_pose_num - 1) * pose_size + 1, start_scale_idx + 1, pose_size, 1, d_err_d_Twb1' * d_err_d_scale);
                    H = FillMatrix(H, start_scale_idx + 1, (host_fid - fixed_pose_num - 1) * pose_size + 1, 1, pose_size, d_err_d_scale' * d_err_d_Twb1);
                end
                if target_fid_valid
                    H = FillMatrix(H, (target_fid - fixed_pose_num - 1) * pose_size + 1, (target_fid - fixed_pose_num - 1) * pose_size + 1, pose_size, pose_size, d_err_d_Twb2' * d_err_d_Twb2);
                    b = FillMatrix(b, (target_fid - fixed_pose_num - 1) * pose_size + 1, 1, pose_size, 1, -d_err_d_Twb2' * err);
                    
                    H = FillMatrix(H, (target_fid - fixed_pose_num - 1) * pose_size + 1, start_idp_idx + pid, pose_size, 1, d_err_d_Twb2' * d_err_d_rho);
                    H = FillMatrix(H, start_idp_idx + pid, (target_fid - fixed_pose_num - 1) * pose_size + 1, 1, pose_size, d_err_d_rho' * d_err_d_Twb2);
                    
                    H = FillMatrix(H, (target_fid - fixed_pose_num - 1) * pose_size + 1, start_scale_idx + 1, pose_size, 1, d_err_d_Twb2' * d_err_d_scale);
                    H = FillMatrix(H, start_scale_idx + 1, (target_fid - fixed_pose_num - 1) * pose_size + 1, 1, pose_size, d_err_d_scale' * d_err_d_Twb2);
                end
                H = FillMatrix(H, start_idp_idx + pid, start_idp_idx + pid, 1, 1, d_err_d_rho' * d_err_d_rho);
                b = FillMatrix(b, start_idp_idx + pid, 1, 1, 1, -d_err_d_rho' * err);
                H = FillMatrix(H, start_scale_idx + 1, start_scale_idx + 1, 1, 1, d_err_d_scale' * d_err_d_scale);
                b = FillMatrix(b, start_scale_idx + 1, 1, 1, 1, -d_err_d_scale' * err);
                
                H = FillMatrix(H, start_idp_idx + pid, start_scale_idx + 1, 1, 1, d_err_d_rho' * d_err_d_scale);
                H = FillMatrix(H, start_scale_idx + 1, start_idp_idx + pid, 1, 1, d_err_d_scale' * d_err_d_rho);
                if host_fid_valid && target_fid_valid
                    H = FillMatrix(H, (host_fid - fixed_pose_num - 1) * pose_size + 1, (target_fid - fixed_pose_num - 1) * pose_size + 1, pose_size, pose_size, d_err_d_Twb1' * d_err_d_Twb2);
                    H = FillMatrix(H, (target_fid - fixed_pose_num - 1) * pose_size + 1, (host_fid - fixed_pose_num - 1) * pose_size + 1, pose_size, pose_size, d_err_d_Twb2' * d_err_d_Twb1);
                end
            end
        end
    end
    loss_vec = [loss_vec; err_sum];
    dx = inv(H) * b;
    fprintf(sprintf('iter: %d, err_sum: %.15f, err_count: %d, mean_err: %.15f, scale: %.15f\n', iter, err_sum, err_count, err_sum / err_count, scale));
    dp = reshape(dx(1:pose_size * (pose_num - fixed_pose_num)), 6, []);
    for id = 1 : size(dp,2)
        pose_id = id + fixed_pose_num;
        pose_Twb{pose_id,1} = Exp(dp(1:3,id), dp(4:6,id)) * pose_Twb{pose_id,1};
    end
    
    for point_id = 1 : size(xyz,1)
        host_frame_id = pid_to_host_fid(point_id);
        idps{1, host_frame_id}(point_id) = idps{1, host_frame_id}(point_id) + dx(start_idp_idx + point_id);
    end
    
    scale = scale + dx(start_scale_idx + 1);
end

figure,plot(diff(loss_vec));

end
function result =  Exp(w, v)


omega = w;
theta_sq = sum(omega.^2);

if (theta_sq < 1e-10)
    theta = 0;
else
    theta = sqrt(theta_sq);
end
so3 = rodrigues(omega);
Omega = SkewSymMat(omega);
Omega_sq = Omega * Omega;
V = zeros(3, 3);

if (theta < 1e-5)
    V = so3;
    
else
    theta_sq = theta * theta;
    V = (eye(3) + ((1) - cos(theta)) / (theta_sq)*Omega + (theta - sin(theta)) / (theta_sq * theta) * Omega_sq);
end

tran = V * v;
result = eye(4);

result(1:3,1:3) = so3;
result(1:3,4) = tran;

end
function H = FillMatrix(H, start_row, start_col, row_size, col_size, value)

if (start_row < 1 || start_col < 1)
    return
end

H(start_row:start_row+row_size-1, start_col:start_col+col_size-1) = H(start_row:start_row+row_size-1, start_col:start_col+col_size-1) + value;
end
function [err, d_err_d_Twb1, d_err_d_Twb2, d_err_d_scale, d_err_d_rho] = computeScaledReprojFactor(Twb1, Twb2, Tbcx, Tbcy, rho, scale, host_bearing, target_bearing)
err = []; d_err_d_Twb1 = []; d_err_d_Twb2 = []; d_err_d_scale = []; d_err_d_rho = [];

Rwb1 = Twb1(1:3,1:3);
twb1 = Twb1(1:3,4);
Rwb2 = Twb2(1:3,1:3);
twb2 = Twb2(1:3,4);
Rbcx = Tbcx(1:3,1:3);
tbcx = Tbcx(1:3,4);
Rbcy = Tbcy(1:3,1:3);
tbcy = Tbcy(1:3,4);

Tb2b1 = inv(Twb2) * Twb1;

Tb2b1_scaled = Tb2b1;
Tb2b1_scaled(1:3,4) = scale * Tb2b1(1:3,4);

Tth = inv(Tbcy) * Tb2b1_scaled * Tbcx;
X = Tth(1:3,1:3) * host_bearing + rho * Tth(1:3,4);
err = X./norm(X) - target_bearing;

d_bearing_d_X = compute_d_bearing_d_pt_jac(X);

d_X_d_Rwb1 = -Rbcy' * Rwb2' * SkewSymMat(Rwb1 * Rbcx * host_bearing) - rho * Rbcy' * (Rwb2' * SkewSymMat(Rwb1 * tbcx + scale * twb1));
d_X_d_Rwb2 =  Rbcy' * Rwb2' * SkewSymMat(Rwb1 * Rbcx * host_bearing) + rho * Rbcy' * (Rwb2' * SkewSymMat(Rwb1 * tbcx + scale * twb1));

d_X_d_twb1 = scale * rho * Rbcy' * Rwb2';
d_X_d_twb2 = -scale * rho * Rbcy' * Rwb2';

d_X_d_rho = Rbcy' * (Rwb2' * (Rwb1 * tbcx + scale * (twb1 - twb2)) - tbcy);

d_X_d_scale = rho * Rbcy' * Rwb2' * (twb1 - twb2);

d_X_d_Twb1 = [d_X_d_Rwb1 d_X_d_twb1];
d_X_d_Twb2 = [d_X_d_Rwb2 d_X_d_twb2];


d_err_d_Twb1 = d_bearing_d_X * d_X_d_Twb1;
d_err_d_Twb2 = d_bearing_d_X * d_X_d_Twb2;
d_err_d_scale = d_bearing_d_X * d_X_d_scale;
d_err_d_rho = d_bearing_d_X * d_X_d_rho;


end
function d_bearing_d_pt = compute_d_bearing_d_pt_jac(pt)
d_bearing_d_pt = (norm(pt).*eye(3) - pt * pt'./norm(pt))./(norm(pt)^2);
end