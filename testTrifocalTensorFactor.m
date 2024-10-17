function testTrifocalTensorFactor()

close all;

data = load('G:\matlab\data\direct\gt\D2_011\4\output.txt');

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


% figure,plotPath(poseMat(1:10:100,:));


intrMat = [500 0 320; 0 500 240; 0 0 1];

host_uv = [100 10; 110 15; 95 8; 300 30; 350 25; 200 8];
host_rho = [0.5; 0.51; 0.52; 0.49; 0.48; 0.47];

pose_num = 20;

T_wc_host = [reshape(poseMat(1, 1:9), 3, 3) poseMat(1, 10:12)';0 0 0 1];

use_bearing = true;false;
is_5dof = true; true;
is_left_5dof = true;
is_left_5dof_right_update = true;
point_trace = pextend(host_uv')';
% point_trace_gt = [host_uv 1];

rho_diff = zeros(pose_num, 1);
rho_noise = 0.0 * (rand(pose_num, 1) - 0.5);
rho_diff(3:end) = rho_noise(3:end);
T_cw_stack{1, 1} = inv(T_wc_host);
% xyz_err = 1./(host_rho + rho_diff) - 1/host_rho;

use_line = true;

if ~use_line
    if 1
        is_5dof = false;
        is_left_5dof = false;
    end
    add_noise = true;
    uv_stack = {};
    for id = 1 : pose_num
        T_wc_cur = Twc_stack{id, 1};
        T_cw_stack{id, 1} = inv(T_wc_cur);
        T_th = inv(T_wc_cur) * T_wc_host;
        xyz_in_host = inv(intrMat) * pextend(host_uv')./repmat(host_rho', 3, 1);
        xyz_in_cur = T_th(1:3,1:3) * xyz_in_host + repmat(T_th(1:3,4), 1, size(host_uv,1));
        [target_braring, ~] = NormalizeVector(xyz_in_cur');
        uv_stack{id, 1} = target_braring ;
    end
    vm_num = size(target_braring, 1);
else
    if 1
        is_5dof = false;
        use_dist_2_plane = true;
    end
    add_noise = true;
    combinations = nchoosek(1:length(host_rho), 2);
    vm_num = size(combinations, 1);
    uv_stack = {};
    for id = 1 : pose_num
        T_wc_cur = Twc_stack{id, 1};
        T_cw_stack{id, 1} = inv(T_wc_cur);
        T_th = inv(T_wc_cur) * T_wc_host;
        uv_stack{id, 1} = [];
        for lid = 1 : size(combinations, 1)
            id1 = combinations(lid,1);
            id2 = combinations(lid,2);
            xyz_in_host1 = inv(intrMat) * pextend(host_uv(id1,:)')./repmat(host_rho(id1)', 3, 1);
            xyz_in_host2 = inv(intrMat) * pextend(host_uv(id2,:)')./repmat(host_rho(id2)', 3, 1);
            xyz_in_cur1 = T_th(1:3,1:3) * xyz_in_host1 + T_th(1:3,4);
            xyz_in_cur2 = T_th(1:3,1:3) * xyz_in_host2 + T_th(1:3,4);
            plane = cross(xyz_in_cur1, xyz_in_cur2);
            uv_stack{id, 1} = [ uv_stack{id, 1}; [plane'./norm(plane) xyz_in_cur1'./norm(xyz_in_cur1) xyz_in_cur2'./norm(xyz_in_cur2)]] ;
        end
    end
end
Twc_stack_gt = Twc_stack;

Twc_stack_noise = Twc_stack_gt;
fix_pose_num = 2;
if add_noise
    for i = 1 : pose_num
        if i > fix_pose_num
            if is_5dof
                Twc_stack_noise{i,1} = Twc_stack{i,1} * [rodrigues(0.01 * (rand(3,1)-0.5)) zeros(3, 1);0 0 0 1];
            else
                if ~use_line
                    Twc_stack_noise{i,1} = Twc_stack{i,1} * [rodrigues(0.01 * (rand(3,1)-0.5)) 0.01 * (rand(3,1)-0.5);0 0 0 1];
                else
                    if ~use_dist_2_plane
                        Twc_stack_noise{i,1} = Twc_stack{i,1} * [rodrigues(0.005 * (rand(3,1)-0.5)) 0.005 * (rand(3,1)-0.5);0 0 0 1];
                    else
                        Twc_stack_noise{i,1} = Twc_stack{i,1} * [rodrigues(0.005 * (rand(3,1)-0.5)) 0.005 * (rand(3,1)-0.5);0 0 0 1];
                    end
                end
            end
        end
    end
end

iter_max = 10;
for iter = 1 : iter_max
    if is_5dof
        pose_size = 5;
    else
        pose_size = 6;
    end
    H = zeros(pose_size * (pose_num), pose_size * (pose_num));
    b = zeros(pose_size * (pose_num), 1);
    error = 0;
    for id = 1 : pose_num
        T_wc_cur = Twc_stack_noise{id, 1};
        T_cw_stack{id, 1} = inv(T_wc_cur);
        if id >= fix_pose_num + 1 % 3
            target_braring_predict = zeros(vm_num, 3);
            for pid = 1 : vm_num
                if ~use_line
                    [err, target_braring_predict(pid,:), d_err_d_Twc1, d_err_d_Twc2, d_err_d_Twc3] = trifocalTransfer(T_cw_stack{id-2, 1}, T_cw_stack{id-1, 1}, T_cw_stack{id, 1}, eye(3), uv_stack{id-2}(pid,:), uv_stack{id-1}(pid,:), uv_stack{id}(pid,:), use_bearing, is_5dof, is_left_5dof, is_left_5dof_right_update);
                else
                    [err, target_braring_predict(pid,:), d_err_d_Twc1, d_err_d_Twc2, d_err_d_Twc3] = trifocalTransferLine(inv(T_cw_stack{id-2, 1}), inv(T_cw_stack{id-1, 1}), inv(T_cw_stack{id, 1}), uv_stack{id-2}(pid,:)', uv_stack{id-1}(pid,:)', uv_stack{id}(pid,:)', is_5dof, use_dist_2_plane);
%                     continue;
                end
                pose_1_start_idx = pose_size * ((id-2)-1) + 1;
                pose_2_start_idx = pose_size * ((id-1)-1) + 1;
                pose_3_start_idx = pose_size * ((id)-1) + 1;
                if (id-2 <= fix_pose_num)
                    d_err_d_Twc1 = zeros(size(d_err_d_Twc1));
                end
                if (id-1 <= fix_pose_num)
                    d_err_d_Twc2 = zeros(size(d_err_d_Twc2));
                end
                if (id <= fix_pose_num)
                    d_err_d_Twc3 = zeros(size(d_err_d_Twc3));
                end
                H = FillMatrix(H, pose_size, pose_size, pose_1_start_idx, pose_1_start_idx, d_err_d_Twc1' * d_err_d_Twc1);
                H = FillMatrix(H, pose_size, pose_size, pose_2_start_idx, pose_2_start_idx, d_err_d_Twc2' * d_err_d_Twc2);
                H = FillMatrix(H, pose_size, pose_size, pose_3_start_idx, pose_3_start_idx, d_err_d_Twc3' * d_err_d_Twc3);
                
                H = FillMatrix(H, pose_size, pose_size, pose_2_start_idx, pose_1_start_idx, d_err_d_Twc2' * d_err_d_Twc1);
                H = FillMatrix(H, pose_size, pose_size, pose_3_start_idx, pose_1_start_idx, d_err_d_Twc3' * d_err_d_Twc1);
                H = FillMatrix(H, pose_size, pose_size, pose_3_start_idx, pose_2_start_idx, d_err_d_Twc3' * d_err_d_Twc2);
                
                H = FillMatrix(H, pose_size, pose_size, pose_1_start_idx, pose_2_start_idx, d_err_d_Twc1' * d_err_d_Twc2);
                H = FillMatrix(H, pose_size, pose_size, pose_1_start_idx, pose_3_start_idx, d_err_d_Twc1' * d_err_d_Twc3);
                H = FillMatrix(H, pose_size, pose_size, pose_2_start_idx, pose_3_start_idx, d_err_d_Twc2' * d_err_d_Twc3);
                
                b = FillMatrix(b, pose_size, 1, pose_1_start_idx, 1, -d_err_d_Twc1' * err);
                b = FillMatrix(b, pose_size, 1, pose_2_start_idx, 1, -d_err_d_Twc2' * err);
                b = FillMatrix(b, pose_size, 1, pose_3_start_idx, 1, -d_err_d_Twc3' * err);
                
                error = error + err' * err * 235 * 235;
                
            end
        end
    end
    dx = inv(H(fix_pose_num * pose_size + 1:end, fix_pose_num * pose_size + 1:end)) * b(fix_pose_num * pose_size + 1:end);
    dxMat = reshape(dx, pose_size, []);
    for k = 1 : size(dxMat, 2)
        if is_5dof
            if ~is_left_5dof
                dR = rodrigues(dxMat(1:3,k));
                Twc_stack_noise{k + fix_pose_num, 1}(1:3,1:3) = Twc_stack_noise{k + fix_pose_num, 1}(1:3,1:3) * dR;
                Ag = ProduceOtherOthogonalBasis(Twc_stack_noise{k + fix_pose_num, 1}(1:3,4));
                rot_vec = Ag * dxMat(4:5,k);
                trans_new = rodrigues(rot_vec) * Twc_stack_noise{k + fix_pose_num, 1}(1:3,4);
                len1 = norm(trans_new);
                len2 = norm(Twc_stack_noise{k + fix_pose_num, 1}(1:3,4));
                trans_diff = len1 - len2;
                fprintf(sprintf('########################## trans_diff: %f\n', trans_diff));
                Twc_stack_noise{k + fix_pose_num, 1}(1:3,4) = trans_new;
            else
                if ~is_left_5dof_right_update
                    dr = dxMat(1:3,k);
                    Ag = ProduceOtherOthogonalBasis(Twc_stack_noise{k + fix_pose_num, 1}(1:3,4));
                    rot_vec = Ag * dxMat(4:5,k);
                    trans_new = rodrigues(rot_vec) * Twc_stack_noise{k + fix_pose_num, 1}(1:3,4);
                    dt = trans_new - Twc_stack_noise{k + fix_pose_num, 1}(1:3,4);
                    %                 dt = trans_new - rodrigues(dr) * Twc_stack_noise{k + fix_pose_num, 1}(1:3,4);
                    dT = Exp(dr, dt);
                    %                 dT = [rodrigues(dr) dt;0 0 0 1];
                    if 1
                        len1 = norm(Twc_stack_noise{k + fix_pose_num, 1}(1:3,4));
                        Twc_stack_noise{k + fix_pose_num, 1} = dT * Twc_stack_noise{k + fix_pose_num, 1};
                        len2 = norm(Twc_stack_noise{k + fix_pose_num, 1}(1:3,4));
                        len3 = norm(trans_new);
                        trans_diff = len1 - len2;
                        trans_diff2 = len1 - len3;
                        %                     fprintf(sprintf('########################## trans_diff: %f, trans_diff2: %f\n', trans_diff, trans_diff2));
                        asfdjgdh = 1;
                        %                     Twc_stack_noise{k + fix_pose_num, 1}(1:3,4) = trans_new;
                    else
                        %                     Twc_stack_noise{k + fix_pose_num, 1}(1:3,1:3) = rodrigues(dr) * Twc_stack_noise{k + fix_pose_num, 1}(1:3,1:3);
                        %                     Twc_stack_noise{k + fix_pose_num, 1}(1:3,4) = trans_new;
                    end
                else
                    dR = rodrigues(dxMat(1:3,k));
                    Twc_stack_noise{k + fix_pose_num, 1}(1:3,1:3) = Twc_stack_noise{k + fix_pose_num, 1}(1:3,1:3) * dR;
                    Ag = ProduceOtherOthogonalBasis(Twc_stack_noise{k + fix_pose_num, 1}(1:3,4));
                    rot_vec = Ag * dxMat(4:5,k);
                    trans_new = rodrigues(rot_vec) * Twc_stack_noise{k + fix_pose_num, 1}(1:3,4);
                    Twc_stack_noise{k + fix_pose_num, 1}(1:3,4) = trans_new;
                end
            end
        else
            dT = Exp(dxMat(1:3,k), dxMat(4:6, k));
            if ~use_line
                Twc_stack_noise{k + fix_pose_num, 1} = dT * Twc_stack_noise{k + fix_pose_num, 1};
            else
                Twc_stack_noise{k + fix_pose_num, 1} = Twc_stack_noise{k + fix_pose_num, 1} * dT;
            end
        end
    end
    fprintf(sprintf('iter: %d, err: %f\n', iter, error));
end

% info = [xyz_err err_close_with_gt_pre];
%
% figure,imshow(zeros(480, 640)); hold on;plot(point_trace(:,1), point_trace(:,2),'.g');plot(point_trace(1,1), point_trace(1,2),'or')


end
function [err, X1_predict_normalized, d_err_d_Twc1, d_err_d_Twc2, d_err_d_Twc3] = trifocalTransferLine(Twc1, Twc2, Twc3, X11, X22, X33, is_5dof, use_dist_2_plane)
Rbc1 = Twc1(1:3,1:3);
Rbc2 = Twc2(1:3,1:3);
Rbc3 = Twc3(1:3,1:3);

tbc1 = Twc1(1:3,4);
tbc2 = Twc2(1:3,4);
tbc3 = Twc3(1:3,4);

R21 = Rbc2' * Rbc1;
t21 = Rbc2' * (tbc1 - tbc2);
R31 = Rbc3' * Rbc1;
t31 = Rbc3' * (tbc1 - tbc3);

X1 = X11(1:3);
X2 = X22(1:3);
X3 = X33(1:3);

X1_predict = (R21' * X2) * (t31' * X3) - (R31' * X3) * (t21' * X2);
X1_predict_normalized = X1_predict./norm(X1_predict);
if 0
    err = X1_predict_normalized + X1;
else
    err = cross(X1_predict_normalized, X1);
    if ~use_dist_2_plane
        err = SkewSymMat(X1) * X1_predict_normalized;
    else
        bearing_start1 = X11(4:6);
        bearing_end1 = X11(7:9);
        err = [bearing_start1 bearing_end1]' * X1_predict_normalized;
    end
end

%   predict = (t31 * X1' * R21' - R31 * X1 * t21') * SkewSymMat(X2) * SkewSymMat(t21) * R21 * X1;
%
%   err = predict./norm(predict) + X3;

predict_norm = norm(X1_predict);
predict_norm2 = predict_norm^2;
d_predictN_d_predict = (predict_norm * eye(3) - X1_predict * X1_predict' / predict_norm) / (predict_norm2);
if ~use_dist_2_plane
    d_err_d_predict = SkewSymMat(X1) * d_predictN_d_predict;
else
    d_err_d_predict = [bearing_start1 bearing_end1]' * d_predictN_d_predict;
end
%
% T21
A1 = X2 * t31' * X3;
B1 = R31' * X3;
d_predict_d_R21 = R21' * SkewSymMat(A1) - B1 * t21' * SkewSymMat(X2);
v1_1 = -R31' * X3;
v3_1 = X2;
d_predict_d_t21 = compute_transpose_vec_jac(v1_1, v3_1);
d_predict_d_T21 = [d_predict_d_R21 d_predict_d_t21];

% T31
A2 = R21' * X2 * t31';
B2 = X3 * t21' * X2;
d_predict_d_R31 = A2 * SkewSymMat(X3) - R31' * SkewSymMat(B2);
v1_2 = R21' * X2;
v3_2 = X3;
d_predict_d_t31 = compute_transpose_vec_jac(v1_2, v3_2);
d_predict_d_T31 = [d_predict_d_R31 d_predict_d_t31];
% relative jac, left perturb
d_err_d_T21 = d_err_d_predict * d_predict_d_T21;
d_err_d_T31 = d_err_d_predict * d_predict_d_T31;

T21 = inv(Twc2) * Twc1;
T31 = inv(Twc3) * Twc1;


% jac_Twc1
d_err_d_Tbc1_right_perturb_part1 = d_err_d_T21 * Adj(T21);
d_err_d_Tbc1_right_perturb_part2 = d_err_d_T31 * Adj(T31);
d_err_d_Tbc1_right_perturb_combined = d_err_d_Tbc1_right_perturb_part1 + d_err_d_Tbc1_right_perturb_part2;
if is_5dof
    Jac_trans_2dof = d_err_d_Tbc1_right_perturb_combined(:,4:6) * (-Rbc1' * SkewSymMat(tbc1)) * ProduceOtherOthogonalBasis(tbc1);
    d_err_d_Twc1 =[ d_err_d_Tbc1_right_perturb_combined(:,1:3)  Jac_trans_2dof];
else
    d_err_d_Twc1 = d_err_d_Tbc1_right_perturb_combined;
end
% jac_Twc2
d_err_d_Tbc2_right_perturb = -d_err_d_T21;
if is_5dof
    Jac_trans_2dof = d_err_d_Tbc2_right_perturb(:,4:6) * (-Rbc2' * SkewSymMat(tbc2)) * ProduceOtherOthogonalBasis(tbc2);
    d_err_d_Twc2 = [d_err_d_Tbc2_right_perturb(:,1:3) Jac_trans_2dof];
else
    d_err_d_Twc2 = d_err_d_Tbc2_right_perturb;
end
% jac_Twc3
d_err_d_Tbc3_right_perturb = -d_err_d_T31;
if is_5dof
    Jac_trans_2dof = d_err_d_Tbc3_right_perturb(:,4:6) * (-Rbc3' * SkewSymMat(tbc3)) * ProduceOtherOthogonalBasis(tbc3);
    d_err_d_Twc3 = [d_err_d_Tbc3_right_perturb(:,1:3) Jac_trans_2dof];
else
    d_err_d_Twc3 = d_err_d_Tbc3_right_perturb;
end

end
function Mat2 = GrowMat(Mat1, rows, cols, start_row, start_col)
if size(Mat1,1) < start_row + rows -1
    pad_row = start_row + rows -1 - size(Mat1,1);
    Mat1 = [Mat1; zeros(pad_row, size(Mat1,2))];
%     Mat1 = [Mat1 zeros(size(Mat1,1), pad_row)];
end

if size(Mat1,2) < start_col + cols -1
    pad_col = start_col + cols -1 - size(Mat1,2);
    Mat1 = [Mat1 zeros(size(Mat1,1), pad_col)];
%     Mat1 = [Mat1; zeros(pad_col, size(Mat1,2))];
end
Mat2 = Mat1;

end
function Mat2 = FillMatrix(Mat1, rows, cols, start_row, start_col, data)


% Mat1 = GrowMat(Mat1, rows, cols, start_row, start_col);
Mat1(start_row:start_row+rows-1, start_col:start_col+cols-1) = Mat1(start_row:start_row+rows-1, start_col:start_col+cols-1) + data;

Mat2 = Mat1;
end
function dist_to_line = EpipolarTransfer(Tc1w, Tc2w, K, focal, bearing1, bearing2, use_bearing)
Twc1 = inv(Tc1w);
Twc2 = inv(Tc2w);

Rgc1 = Twc1(1:3,1:3);
Rgc2 = Twc2(1:3,1:3);

pgc1 = Twc1(1:3,4);
pgc2 = Twc2(1:3,4);

T12 = inv(Twc1) * Twc2;

R12 = Rgc1'*Rgc2;
t12 = Rgc1'*(pgc2-pgc1);

R21 = R12';
t21 = -R12'*t12;

invK = inv(K);

m = 1;
if ~use_bearing
    %     uv1_normal = invK*[ bearing1(m,1:2)'; 1 ];
    %     uv2_normal = invK*[ bearing2(m,1:2)'; 1 ];
    uv1_normal = [ bearing1(m,1:2)'; 1 ];
    uv2_normal = [ bearing2(m,1:2)'; 1 ];
else
    uv1_normal =  bearing1(m,1:3)';
    uv2_normal =  bearing2(m,1:3)';
end


F21 = invK' * SkewSymMat(t21) * R21 * invK;
epLine2 = F21*uv1_normal;
epLine2 = epLine2./norm(epLine2(1:2));
dist_to_line = dot(epLine2, uv2_normal);

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
function [err, bearing33, d_err_d_Twc1, d_err_d_Twc2, d_err_d_Twc3] = trifocalTransfer(Tc1w, Tc2w, Tc3w, K, bearing1, bearing2, bearing3, use_bearing, is_5dof, is_left_5dof, is_left_5dof_right_update)

d_err_d_Rwc1 =[];
d_err_d_Rwc2 = [];
d_err_d_Rwc3 = [];
d_err_d_twc1 = [];
d_err_d_twc2 = [];
d_err_d_twc3 = [];

X1 = bearing1';
X2 = bearing2';
X3 = bearing3';




Twc1 = inv(Tc1w);
Twc2 = inv(Tc2w);
Twc3 = inv(Tc3w);

Rbc1 = Twc1(1:3,1:3);
Rbc2 = Twc2(1:3,1:3);
Rbc3 = Twc3(1:3,1:3);

tbc1 = Twc1(1:3,4);
tbc2 = Twc2(1:3,4);
tbc3 = Twc3(1:3,4);


Rgc1 = Twc1(1:3,1:3);
Rgc2 = Twc2(1:3,1:3);
Rgc3 = Twc3(1:3,1:3);

pgc1 = Twc1(1:3,4);
pgc2 = Twc2(1:3,4);
pgc3 = Twc3(1:3,4);

T12 = inv(Twc1) * Twc2;
T13 = inv(Twc1) * Twc3;

R12 = Rgc1'*Rgc2;
t12 = Rgc1'*(pgc2-pgc1);
R13 = Rgc1'*Rgc3;
t13 = Rgc1'*(pgc3-pgc1);

R21 = R12';
t21 = -R12'*t12;

R31 = R13';
t31 = -R13'*t13;



if is_5dof && ~is_left_5dof
    %=================================
    A = (Rbc3' * tbc1 - Rbc3'*tbc3) * X1';
    B = SkewSymMat(X2) * SkewSymMat(Rbc2' * tbc1 - Rbc2' * tbc2);
    C = X1 * (Rbc2' * tbc1 - Rbc2' * tbc2)' * SkewSymMat(X2) * SkewSymMat(Rbc2' * tbc1 - Rbc2' * tbc2);
    d_err_d_Rwc1 = -A * Rbc1' * Rbc2 * B * Rbc2' * Rbc1 * SkewSymMat(X1) + A * SkewSymMat(Rbc1' * Rbc2 * B * Rbc2' * Rbc1 * X1) ...
        + Rbc3' * Rbc1 * C * Rbc2' * Rbc1 * SkewSymMat(X1) + Rbc3' * Rbc1 * SkewSymMat(C * Rbc2' * Rbc1 * X1);
    %==================================
    A = X1' * Rbc1' * Rbc2 * SkewSymMat(X2) * Rbc2';
    B = Rbc3' * Rbc1 * X1;
    C = Rbc2 * SkewSymMat(X2) * Rbc2';
    D = Rbc1 * X1;
    d_err_d_twc1 = Rbc3' * (tbc1 - tbc3) * A * SkewSymMat(D) * SkewSymMat(tbc1) ...
        + Rbc3' * SkewSymMat(tbc1 * A * SkewSymMat(D) * (tbc1 - tbc2)) ...
        - B * (tbc1' - tbc2') * C * SkewSymMat(D) * SkewSymMat(tbc1) ...
        + B * tbc1' * SkewSymMat(C * SkewSymMat(D) * (tbc1 - tbc2));
    d_err_d_twc1 = d_err_d_twc1 * ProduceOtherOthogonalBasis(tbc1);
    %==================================
    A = (Rbc3' * tbc1 - Rbc3' * tbc3) * X1';
    B = Rbc3' * Rbc1 * X1;
    C = SkewSymMat(tbc1 - tbc2) * Rbc1 * X1;
    D = A * Rbc1' - B * (tbc1 - tbc2)';
    d_err_d_Rwc2 = D * SkewSymMat(C) * Rbc2 * SkewSymMat(X2);
    %==================================
    A = (Rbc3' * tbc1 - Rbc3' * tbc3) * X1' * Rbc1' * Rbc2 * SkewSymMat(X2) * Rbc2';
    B = Rbc3' * Rbc1 * X1;
    C = Rbc2 * SkewSymMat(X2) * Rbc2';
    D = Rbc1 * X1;
    d_err_d_twc2 = -A * SkewSymMat(D) * SkewSymMat(tbc2) + B * (tbc1' - tbc2') * C * SkewSymMat(D) * SkewSymMat(tbc2) ...
        - B * tbc2' * SkewSymMat(C * SkewSymMat(D) * (tbc1 - tbc2));
    d_err_d_twc2 = d_err_d_twc2 * ProduceOtherOthogonalBasis(tbc2);
    %=====================================================
    A = X1' * (Rbc2' * Rbc1)' * SkewSymMat(X2) * SkewSymMat(Rbc2' * tbc1 - Rbc2' * tbc2) * Rbc2' * Rbc1 * X1;
    B = Rbc1 * X1 * (Rbc2' * tbc1 - Rbc2' * tbc2)' * SkewSymMat(X2) * SkewSymMat(Rbc2' * tbc1 - Rbc2' * tbc2) * Rbc2' * Rbc1 * X1;
    C = (tbc1 - tbc3) * A - B;
    d_err_d_Rwc3 = SkewSymMat(Rbc3' * C);
    %===================================================
    A = X1' * Rbc1' * Rbc2 * SkewSymMat(X2) * Rbc2' * SkewSymMat(tbc1 - tbc2) * Rbc1 * X1;
    d_err_d_twc3 = Rbc3' * SkewSymMat(tbc3 * A);
    d_err_d_twc3 = d_err_d_twc3 * ProduceOtherOthogonalBasis(tbc3);

else
    %=================================
    d_err_d_Rwc1 =  Rbc3' * (tbc1 - tbc3) * X1' * Rbc1' * SkewSymMat(Rbc2 * X2) * (-1) * (Rbc1 * SkewSymMat(X1) * Rbc1' * (-SkewSymMat(tbc1))    +   Rbc1 * SkewSymMat(X1) * Rbc1' * SkewSymMat(tbc1 - tbc2)   -  SkewSymMat(Rbc1 * SkewSymMat(X1) * Rbc1' * (tbc1 - tbc2)) ) ...
        + Rbc3' * (tbc1 - tbc3) * X1' * Rbc1' * SkewSymMat(SkewSymMat(Rbc2 * X2) * (-1) * Rbc1 * SkewSymMat(X1) * Rbc1' * (tbc1 - tbc2)) ...
        +   (-1) * Rbc3' * SkewSymMat(tbc1 * X1' * Rbc1' * SkewSymMat(Rbc2 * X2) * (-1) * Rbc1 * SkewSymMat(X1) * Rbc1' * (tbc1 - tbc2)) ...
        ...
        -Rbc3' * Rbc1 * X1 * (tbc1' - tbc2') * SkewSymMat(Rbc2 * X2) * (-1) * (Rbc1 * SkewSymMat(X1) * Rbc1' * (-SkewSymMat(tbc1))    +   Rbc1 * SkewSymMat(X1) * Rbc1' * SkewSymMat(tbc1 - tbc2)   -  SkewSymMat(Rbc1 * SkewSymMat(X1) * Rbc1' * (tbc1 - tbc2)) );
    v1 = -Rbc3' * Rbc1 * X1;
    v3 = SkewSymMat(tbc1) * SkewSymMat(Rbc2 * X2) * (-1) *  Rbc1 * SkewSymMat(X1) * Rbc1' * (tbc1 - tbc2);
    jac = compute_transpose_vec_jac(v1, v3);
    d_err_d_Rwc1 = d_err_d_Rwc1 + jac;
    d_err_d_Rwc1 = d_err_d_Rwc1 + (-1) * Rbc3' * (-1) * SkewSymMat(Rbc1 * X1 * (tbc1' - tbc2') * SkewSymMat(Rbc2 * X2) * (-1) * Rbc1 * SkewSymMat(X1) * Rbc1' * (tbc1 - tbc2));
    %=================================
    d_err_d_twc1 = Rbc3' * tbc1 * X1' * Rbc1' * SkewSymMat(Rbc2 * X2) * (-1) * SkewSymMat(Rbc1 * X1) ...
        + Rbc3' *        (X1' * Rbc1' * SkewSymMat(Rbc2 * X2) * SkewSymMat(tbc1) * Rbc1 * X1) ...
        + (-1) * Rbc3' * (X1' * Rbc1' * SkewSymMat(Rbc2 * X2) * SkewSymMat(tbc2) * Rbc1 * X1) ...
        + (-1) * Rbc3' * tbc3 * X1' * Rbc1' * SkewSymMat(Rbc2 * X2) * (-1) * SkewSymMat(Rbc1 * X1) ...
        + (-1) * Rbc3' * Rbc1 * X1 * tbc1' * SkewSymMat(Rbc2 * X2) * (-1) * SkewSymMat(Rbc1 * X1);
    v1 = -Rbc3' * Rbc1 * X1;
    v3 = SkewSymMat(Rbc2 * X2) * SkewSymMat(tbc1) * Rbc1 * X1;
    jac = compute_transpose_vec_jac(v1, v3);
    d_err_d_twc1 = d_err_d_twc1 + jac;
    v1 = Rbc3' * Rbc1 * X1;
    v3 = SkewSymMat(Rbc2 * X2) * SkewSymMat(tbc2) * Rbc1 * X1;
    jac = compute_transpose_vec_jac(v1, v3);
    d_err_d_twc1 = d_err_d_twc1 + jac;
    d_err_d_twc1 = d_err_d_twc1 + Rbc3' * Rbc1 * X1 * tbc2' * SkewSymMat(Rbc2 * X2) * (-1) * SkewSymMat(Rbc1 * X1);
    %=================================
    d_err_d_Rwc2 = (-1) * Rbc3' * (tbc1 - tbc3) * X1' * Rbc1' * Rbc2 * SkewSymMat(X2) * Rbc2' * SkewSymMat(Rbc1 * X1) * SkewSymMat(tbc2) ...
        + Rbc3' * (tbc1 - tbc3) * X1' * Rbc1' * Rbc2 * SkewSymMat(X2) * Rbc2' * (-1) * SkewSymMat(SkewSymMat(Rbc1 * X1) * (tbc1 - tbc2)) ...
        - Rbc3' * (tbc1 - tbc3) * X1' * Rbc1' * (-1) * SkewSymMat(Rbc2 * SkewSymMat(X2) * Rbc2' * SkewSymMat(Rbc1 * X1) * (tbc1 - tbc2)) ...
        + Rbc3' * Rbc1 * X1 * (tbc1' - tbc2') * Rbc2 * SkewSymMat(X2) * Rbc2' * SkewSymMat(Rbc1 * X1) * SkewSymMat(tbc2) ...
        + (-1) * Rbc3' * Rbc1 * X1 * (tbc1' - tbc2') * Rbc2 * SkewSymMat(X2) * Rbc2' * (-1) * SkewSymMat(SkewSymMat(Rbc1 * X1) * (tbc1 - tbc2)) ...
        + Rbc3' * Rbc1 * X1 * (tbc1' - tbc2') * (-1) * SkewSymMat(Rbc2 * SkewSymMat(X2) * Rbc2' * SkewSymMat(Rbc1 * X1) * (tbc1 - tbc2)) ...
        + Rbc3' * Rbc1 * X1 * tbc2'           * (-1) * SkewSymMat(Rbc2 * SkewSymMat(X2) * Rbc2' * SkewSymMat(Rbc1 * X1) * (tbc1 - tbc2));
    %=====================================
    d_err_d_twc2 = (-1) * Rbc3' * (tbc1 - tbc3) * X1' * Rbc1' * SkewSymMat(Rbc2 * X2) * SkewSymMat(Rbc1 * X1) * (-1) ...
        + Rbc3' * Rbc1 * X1 * (tbc1' - tbc2') * SkewSymMat(Rbc2 * X2) * SkewSymMat(Rbc1 * X1) * (-1);
    v1 = Rbc3' * Rbc1 * X1 * (-1);
    v3 = SkewSymMat(Rbc2 * X2) * SkewSymMat(Rbc1 * X1) * (tbc1 - tbc2);
    jac = compute_transpose_vec_jac(v1, v3);
    d_err_d_twc2 = d_err_d_twc2 + jac;
    %=====================================
    d_err_d_Rwc3 = Rbc3' * SkewSymMat(tbc3 * X1' * Rbc1' * SkewSymMat(Rbc2 * X2) * SkewSymMat(tbc1 - tbc2) * Rbc1 * X1) ...
        + (-1) *  Rbc3' * (-1) * SkewSymMat((tbc1 - tbc3) * X1' * Rbc1' * SkewSymMat(Rbc2 * X2) * SkewSymMat(tbc1 - tbc2) * Rbc1 * X1) ...
        + Rbc3' * (-1) * SkewSymMat(Rbc1 * X1 * (tbc1' - tbc2') * SkewSymMat(Rbc2 * X2) * SkewSymMat(tbc1 - tbc2) * Rbc1 * X1);
    %=====================================
    d_err_d_twc3 = -Rbc3' * (X1' * Rbc1' * SkewSymMat(Rbc2 * X2) * SkewSymMat(tbc1 - tbc2) * Rbc1 * X1);
end

if is_left_5dof
    if ~is_left_5dof_right_update
        d_err_d_twc1 = d_err_d_twc1 * (-SkewSymMat(tbc1)) * ProduceOtherOthogonalBasis(tbc1);
        d_err_d_twc2 = d_err_d_twc2 * (-SkewSymMat(tbc2)) * ProduceOtherOthogonalBasis(tbc2);
        d_err_d_twc3 = d_err_d_twc3 * (-SkewSymMat(tbc3)) * ProduceOtherOthogonalBasis(tbc3);
    else
        d_err_d_Twc1 = [d_err_d_Rwc1 d_err_d_twc1] * Adj(Twc1);
        d_err_d_Twc2 = [d_err_d_Rwc2 d_err_d_twc2] * Adj(Twc2);
        d_err_d_Twc3 = [d_err_d_Rwc3 d_err_d_twc3] * Adj(Twc3);
        
        d_err_d_Rwc1 = d_err_d_Twc1(:,1:3);
        d_err_d_twc1 = d_err_d_Twc1(:,4:6);
        
        d_err_d_Rwc2 = d_err_d_Twc2(:,1:3);
        d_err_d_twc2 = d_err_d_Twc2(:,4:6);
        
        d_err_d_Rwc3 = d_err_d_Twc3(:,1:3);
        d_err_d_twc3 = d_err_d_Twc3(:,4:6);
        
        d_err_d_twc11 = d_err_d_twc1 * (-Rbc1' * SkewSymMat(tbc1)) * ProduceOtherOthogonalBasis(tbc1);
        d_err_d_twc22 = d_err_d_twc2 * (-Rbc2' * SkewSymMat(tbc2)) * ProduceOtherOthogonalBasis(tbc2);
        d_err_d_twc33 = d_err_d_twc3 * (-Rbc3' * SkewSymMat(tbc3)) * ProduceOtherOthogonalBasis(tbc3);
        
        
        
        if 0
            %=================================
            A = (Rbc3' * tbc1 - Rbc3'*tbc3) * X1';
            B = SkewSymMat(X2) * SkewSymMat(Rbc2' * tbc1 - Rbc2' * tbc2);
            C = X1 * (Rbc2' * tbc1 - Rbc2' * tbc2)' * SkewSymMat(X2) * SkewSymMat(Rbc2' * tbc1 - Rbc2' * tbc2);
            d_err_d_Rwc1 = -A * Rbc1' * Rbc2 * B * Rbc2' * Rbc1 * SkewSymMat(X1) + A * SkewSymMat(Rbc1' * Rbc2 * B * Rbc2' * Rbc1 * X1) ...
                + Rbc3' * Rbc1 * C * Rbc2' * Rbc1 * SkewSymMat(X1) + Rbc3' * Rbc1 * SkewSymMat(C * Rbc2' * Rbc1 * X1);
            %==================================
            A = X1' * Rbc1' * Rbc2 * SkewSymMat(X2) * Rbc2';
            B = Rbc3' * Rbc1 * X1;
            C = Rbc2 * SkewSymMat(X2) * Rbc2';
            D = Rbc1 * X1;
            d_err_d_twc1 = Rbc3' * (tbc1 - tbc3) * A * SkewSymMat(D) * SkewSymMat(tbc1) ...
                + Rbc3' * SkewSymMat(tbc1 * A * SkewSymMat(D) * (tbc1 - tbc2)) ...
                - B * (tbc1' - tbc2') * C * SkewSymMat(D) * SkewSymMat(tbc1) ...
                + B * tbc1' * SkewSymMat(C * SkewSymMat(D) * (tbc1 - tbc2));
            d_err_d_twc1 = d_err_d_twc1 * ProduceOtherOthogonalBasis(tbc1);
            %==================================
            A = (Rbc3' * tbc1 - Rbc3' * tbc3) * X1';
            B = Rbc3' * Rbc1 * X1;
            C = SkewSymMat(tbc1 - tbc2) * Rbc1 * X1;
            D = A * Rbc1' - B * (tbc1 - tbc2)';
            d_err_d_Rwc2 = D * SkewSymMat(C) * Rbc2 * SkewSymMat(X2);
            %==================================
            A = (Rbc3' * tbc1 - Rbc3' * tbc3) * X1' * Rbc1' * Rbc2 * SkewSymMat(X2) * Rbc2';
            B = Rbc3' * Rbc1 * X1;
            C = Rbc2 * SkewSymMat(X2) * Rbc2';
            D = Rbc1 * X1;
            d_err_d_twc2 = -A * SkewSymMat(D) * SkewSymMat(tbc2) + B * (tbc1' - tbc2') * C * SkewSymMat(D) * SkewSymMat(tbc2) ...
                - B * tbc2' * SkewSymMat(C * SkewSymMat(D) * (tbc1 - tbc2));
            d_err_d_twc2 = d_err_d_twc2 * ProduceOtherOthogonalBasis(tbc2);
            %=====================================================
            A = X1' * (Rbc2' * Rbc1)' * SkewSymMat(X2) * SkewSymMat(Rbc2' * tbc1 - Rbc2' * tbc2) * Rbc2' * Rbc1 * X1;
            B = Rbc1 * X1 * (Rbc2' * tbc1 - Rbc2' * tbc2)' * SkewSymMat(X2) * SkewSymMat(Rbc2' * tbc1 - Rbc2' * tbc2) * Rbc2' * Rbc1 * X1;
            C = (tbc1 - tbc3) * A - B;
            d_err_d_Rwc3 = SkewSymMat(Rbc3' * C);
            %===================================================
            A = X1' * Rbc1' * Rbc2 * SkewSymMat(X2) * Rbc2' * SkewSymMat(tbc1 - tbc2) * Rbc1 * X1;
            d_err_d_twc3 = Rbc3' * SkewSymMat(tbc3 * A);
            d_err_d_twc3 = d_err_d_twc3 * ProduceOtherOthogonalBasis(tbc3);
            
            diff1 = d_err_d_twc11 - d_err_d_twc1;
            diff2 = d_err_d_twc22 - d_err_d_twc2;
            diff3 = d_err_d_twc33 - d_err_d_twc3;
        else
            d_err_d_twc1 = d_err_d_twc11;
            d_err_d_twc2 = d_err_d_twc22;
            d_err_d_twc3 = d_err_d_twc33;
        end
    end
end



PA = [ eye(3), zeros( 3, 1 ) ];
PB = [ R12', -R12'*t12 ]; % T21
PC = [ R13', -R13'*t13 ]; % T31

F12 = R12'*SkewSymMat( t12 );
F12_2 =  SkewSymMat( t12 )' * R12;
invK = inv(K);

z = zeros( 3, size( bearing1, 1 ) );
for m = 1:size( bearing1, 1 )
    if ~use_bearing
        uv1_normal = invK*[ bearing1(m,1:2)'; 1 ]; uv1_normal = uv1_normal/uv1_normal(3);
        uv2_normal = invK*[ bearing2(m,1:2)'; 1 ]; uv2_normal = uv2_normal/uv2_normal(3);
    else
        uv1_normal =  bearing1(m,1:3)';%./norm(bearing1(m,1:3)');
        uv2_normal =  bearing2(m,1:3)';%./norm(bearing2(m,1:3)');
    end
    epLine = F12*uv1_normal;
    F21_check = F12';
    F21 = SkewSymMat(t21) * R21;
    epLine2 = F21*uv1_normal;
    %     epLine = F12_2'*uv1_normal;
    epLine = epLine./norm(epLine);
    if ~use_bearing
        epLineNormal = [epLine(2); -epLine(1); -uv2_normal(1)*epLine(2)+uv2_normal(2)*epLine(1)];
        epLineNormal2 = [epLine2(2); -epLine2(1); -uv2_normal(1)*epLine2(2)+uv2_normal(2)*epLine2(1)];
    else
        epLineNormal = cross(uv2_normal, epLine);
        if 1
            epLineNormal2 = cross(uv2_normal, epLine2);
        else
            epLineNormal2 = epLine2;
        end
    end
    if 1
        epLineNormal = epLineNormal./norm(epLineNormal(1:3));
    end
    if 0
        for k = 1:3
            for i = 1:3
                for j = 1:3
                    Tijk = PB(j,i)*PC(k,4)-PB(j,4)*PC(k,i);
                    z(k,m) = z(k,m)+uv1_normal(i)*epLineNormal(j)*Tijk;
                end
            end
        end
    else
        Tijk = zeros(3,3,3);
        for k = 1:3
            for i = 1:3
                for j = 1:3
                    Tijk(i,j,k) = PB(j,i)*PC(k,4)-PB(j,4)*PC(k,i);
                    %                     z(k,m) = z(k,m)+uv1_normal(i)*epLineNormal(j)*Tijk;
                end
            end
        end
        
        for k = 1:3
            for i = 1:3
                for j = 1:3
                    %                     Tijk = PB(j,i)*PC(k,4)-PB(j,4)*PC(k,i);
                    z(k,m) = z(k,m)+uv1_normal(i)*epLineNormal(j)*Tijk(i,j,k);
                end
            end
        end
        T21 = PB;
        T31 = PC;
        T_check = TFT_from_P([eye(3) zeros(3, 1)],T21,T31);
        for aa = 1 : 3
            T_check2(:,:,aa) = [T21(:,aa) * T31(:,4)' - T21(:,4) * T31(:,aa)'];
            %             T_check_21(:,:,aa) = T21(:,aa) * T31(:,4)';
        end
        err0 = T_check2(:) - T_check(:);
        err1 = SkewSymMat(bearing2(m,:)) * (bearing1(m, 1) * T_check(:,:,1) + bearing1(m, 2) * T_check(:,:,2) + bearing1(m, 3) * T_check(:,:,3)) * SkewSymMat(bearing3(m,:));
        err2 =            epLineNormal2' * (uv1_normal(1)  * T_check(:,:,1) + uv1_normal(2)  * T_check(:,:,2) + uv1_normal(3)  * T_check(:,:,3)) * SkewSymMat(bearing3(m,:));
        
        T_sum = T21(1:3,1:3) * uv1_normal * T31(1:3,4)' - T21(1:3,4) * (T31(1:3,1:3) * uv1_normal)';
        T_sum_transpose = T31(1:3,4) * (T21(1:3,1:3) * uv1_normal)' - (T31(1:3,1:3) * uv1_normal) * T21(1:3,4)';
        T_sum_transpose2 = T31(1:3,4) * (uv1_normal' * T21(1:3,1:3)') - (T31(1:3,1:3) * uv1_normal) * T21(1:3,4)';
        T_sum_transpose3 = T31(1:3,4) * uv1_normal' * T21(1:3,1:3)' - (T31(1:3,1:3) * uv1_normal) * T21(1:3,4)';
        
        %         T_sum2 = T21(1:3,1:3) * ()
        err3 = epLineNormal2' * T_sum * SkewSymMat(bearing3(m,:));
        err4 = T_sum' * epLineNormal2./norm(T_sum' * epLineNormal2) + bearing3(m,:)';
        err5 = T_sum_transpose * epLineNormal2./norm(T_sum_transpose * epLineNormal2) + bearing3(m,:)';
        
        predict = t31 * uv1_normal' * R21' * SkewSymMat(uv2_normal) * SkewSymMat(t21) * R21 * uv1_normal - R31 * uv1_normal * t21' * SkewSymMat(uv2_normal) * SkewSymMat(t21) * R21 * uv1_normal;
        err6 = predict./norm(predict) + bearing3(m,:)';
        err7 = uv1_normal' * SkewSymMat(uv2_normal) - (-(SkewSymMat(uv2_normal) * uv1_normal)');
        %         z(:,m) = (uv1_normal(1) * T(:,:,1) + uv1_normal(2) * T(:,:,2) + uv1_normal(3) * T(:,:,3)) * epLineNormal;
    end
    z(:,m) = K*z(:,m);
    if ~use_bearing
        z(:,m) = z(:,m)/z(3,m);
    else
        z(:,m) = z(:,m)/norm(z(:,m));
        SkewSymMat( z(:,m)) * bearing3(m,:)';
        z(:,m) = -predict'./norm(predict);
    end
end

if 0
    z = z(1:2,:);
    z = z(:);
else
    diff = z' - bearing3;
end
bearing33 = z';
err = predict/norm(predict) + bearing3';

d_bearing_d_pt = compute_d_bearing_d_pt_jac(predict);

d_err_d_Rwc1 = d_bearing_d_pt * d_err_d_Rwc1;
d_err_d_Rwc2 = d_bearing_d_pt * d_err_d_Rwc2;
d_err_d_Rwc3 = d_bearing_d_pt * d_err_d_Rwc3;
d_err_d_twc1 = d_bearing_d_pt * d_err_d_twc1;
d_err_d_twc2 = d_bearing_d_pt * d_err_d_twc2;
d_err_d_twc3 = d_bearing_d_pt * d_err_d_twc3;

d_err_d_Twc1 = [d_err_d_Rwc1 d_err_d_twc1];
d_err_d_Twc2 = [d_err_d_Rwc2 d_err_d_twc2];
d_err_d_Twc3 = [d_err_d_Rwc3 d_err_d_twc3];

end
function T=TFT_from_P(P1,P2,P3)
%TFT_FROM_P Trifocal tensor from the projection matrices.
%
%  General formula to compute the trifocal tensor from any three projection
%  matrices.
%

T=zeros(3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            T(j,k,i)=(-1)^(i+1)*det([P1([1:(i-1) (i+1):3],:);P2(j,:);P3(k,:)]);
        end
    end
end
% T=T/norm(T(:));
end

function A = ProduceOtherOthogonalBasis(n)
N = n;
if (N(1) < 0)
    N(1) = -N(1);
end
if (N(2) < 0)
    N(2) = -N(2);
end
if (N(3) < 0)
    N(3) = -N(3);
end

minIdx = 0;
if (N(1) <= N(2))
    if (N(1) <= N(3))
        minIdx = 0;
    else
        minIdx = 2;
    end
else
    if (N(2) <= N(3))
        minIdx = 1;
    else
        minIdx = 2;
    end
end

A = zeros(3,2);

if minIdx == 0
    A(:,1) = [0, -n(3), n(2)]';
end
if minIdx == 1
    A(:,1) = [n(3), 0, -n(1)]';
end
if minIdx == 2
    A(:,1) = [-n(2), n(1), 0]';
end

A(:,2) = cross(n, A(:,1));
A(:,1) = A(:,1)./norm(A(:,1));
A(:,2) = A(:,2)./norm(A(:,2));

end
function d_bearing_d_pt = compute_d_bearing_d_pt_jac(pt)

d_bearing_d_pt = (norm(pt).*eye(3) - pt * pt'./norm(pt))./(norm(pt)^2);

if 0
    
    ptNorm = norm(pt);
    ptNorm_3_2 = -0.5 * (1 / (ptNorm * ptNorm * ptNorm));
    d_err_d_pt3d = zeros(3,3);
    d_err_d_pt3d(0+1, 0+1) = 1.0 / ptNorm + pt(0+1) * (ptNorm_3_2 * 2 * pt(0+1));
    d_err_d_pt3d(0+1, 1+1) = pt(0+1) * (ptNorm_3_2 * 2 * pt(1+1));
    d_err_d_pt3d(0+1, 2+1) = pt(0+1) * (ptNorm_3_2 * 2 * pt(2+1));
    
    d_err_d_pt3d(1+1, 0+1) = pt(1+1) * (ptNorm_3_2 * 2 * pt(0+1));
    d_err_d_pt3d(1+1, 1+1) = 1.0 / ptNorm + pt(1+1) * (ptNorm_3_2 * 2 * pt(1+1));
    d_err_d_pt3d(1+1, 2+1) = pt(1+1) * (ptNorm_3_2 * 2 * pt(2+1));
    
    d_err_d_pt3d(2+1, 0+1) = pt(2+1) * (ptNorm_3_2 * 2 * pt(0+1));
    d_err_d_pt3d(2+1, 1+1) = pt(2+1) * (ptNorm_3_2 * 2 * pt(1+1));
    d_err_d_pt3d(2+1, 2+1) = 1.0 / ptNorm + pt(2+1) * (ptNorm_3_2 * 2 * pt(2+1));
end

end
function jac = compute_transpose_vec_jac(v1, v3)

% v1 * v2' * v3

v1x = v1(1);
v1y = v1(2);
v1z = v1(3);

v3x = v3(1);
v3y = v3(2);
v3z = v3(3);


jac = zeros(3, 3);
jac(1, 1) = v1x * v3x;
jac(1, 2) = v1x * v3y;
jac(1, 3) = v1x * v3z;

jac(2, 1) = v1y * v3x;
jac(2, 2) = v1y * v3y;
jac(2, 3) = v1y * v3z;

jac(3, 1) = v1z * v3x;
jac(3, 2) = v1z * v3y;
jac(3, 3) = v1z * v3z;

end
function jac = Adj(T)
R = T(1:3,1:3);
jac = zeros(6,6);
jac(1:3,1:3) = R;
jac(4:6, 4:6) = R;
jac(4:6,1:3) = SkewSymMat(T(1:3,4)) * R;
jac(1:3, 4:6) = zeros(3, 3);

end