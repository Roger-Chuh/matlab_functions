function testLineTriangulationIdp()
K = [499 0 321; 0 501 239; 0 0 1];
Rwc1 = rodrigues([0.01 0.02 0.03]);
twc1 = [-100; 0; 100];
% Rwc1 = eye(3);

Rwc{1} = Rwc1;
twc{1} = twc1;
pose_num = 20;


for i = 2 : pose_num
    Rwc{i} = Rwc1*rodrigues(0.01 * (rand(3,1)-0.5));
    twc{i} = twc1 + 1 * (rand(3,1)-0.5);
end



point_num = 20;
split_num = 17;
for pid = 1 : point_num
    pt1 = 100 * rand(3,1);
    pt2 = 100 * rand(3,1);
    pt1 = sign(pt1(3)) * pt1;
    pt2 = sign(pt2(3)) * pt2;
    dir = (pt1 - pt2) / norm(pt1 - pt2);
    pt_stack = [pt1 pt2];
    for j = 1:split_num - 2
        pt = pt1 + 2 * dir;
        pt_stack = [pt_stack pt];
    end
    ptS{pid} = pt_stack;
    
end

host_fid = 1;

for i = 1 : pose_num
    Tcw  = inv([Rwc{i} twc{i}; 0 0 0 1]);
    for pid = 1 : point_num
        obs = [];
        for j = 1 : size(ptS{pid},2)
            point = Tcw(1:3,1:3) * ptS{pid}(:,j) + Tcw(1:3,4);
            obs = [obs; point'./norm(point)];
            if i == host_fid && j == 1
                dir_h1{pid} = point./norm(point);
                uv_start{pid} = point/point(3);
                rou1(pid) = 1./norm(point);
            end
            if i == host_fid && j == 2
                dir_h2{pid} = point./norm(point);
                uv_end{pid} = point/point(3);
                rou2(pid) = 1./norm(point);
            end
        end
        
        
        Obs{i,pid} = obs;
    end
end



for i = 1 : pid
    pii{i} = pi_from_ppp(uv_start{i}, uv_end{i}, [0;0;0]);
    pii{i} = pii{i}./norm(pii{i}(1:3));
    % ni = pii(1:3);
    % ni = ni./norm(ni);
end


for i = 1 : pose_num
    if i == host_fid
        %        continue;
    end
    t1 = twc{i};
    R1 = Rwc{i};
    
    tij = Rwc{host_fid}' * (t1 - twc{host_fid});   % tij
    Rij = Rwc{host_fid}' * R1;          % Rij
    
    for j = 1 : point_num
        p3 = Obs{i,j}(1,:)';
        p3 = p3 / p3(3);
        p4 =  Obs{i,j}(2,:)';
        p4 = p4 / p4(3);
        p3 = Rij * p3 + tij;
        p4 = Rij * p4 + tij;
        pij = pi_from_ppp(p3, p4,tij);
        
        plk = pipi_plk( pii{j}, pij );
        plk = plk./norm(plk(4:6));
        PLK{i}(:,j) = plk;
    end
    
end

err = [];
for fid = 1 : pose_num
    if fid == host_fid
        continue;
    end
    for pid = 1 : point_num
        for sample = 1 : split_num
            [res1, d_err_d_host_pose_1, d_err_d_target_pose_1, d_err_d_rou12_1] = computeResAndJac(Rwc{host_fid}, twc{host_fid}, Rwc{fid},  twc{fid}, PLK{fid}(:,pid), Obs{fid, pid}(sample,:)', dir_h1{pid}, dir_h2{pid}, rou1(pid), rou2(pid));
            err = [err; res1];
        end
    end
    
end

Rwc_gt = Rwc;
twc_gt = twc;

rou1_gt = rou1;
rou2_gt = rou2;


pose_used = pose_num;
for i = 1 : pose_used - 2
    Rwc{i} = Rwc{i} * rodrigues(0.1 * (rand(3,1) - 0.5));
    twc{i} = twc{i} + 2 * (rand(3,1) - 0.5);
end

rou1 = rou1 + 0.1;
rou2 = rou2 + 0.2;

for iter = 1 : 20
    H = zeros(6 * pose_used + 2 * point_num, 6 * pose_used + 2 * point_num);
    b = zeros(6 * pose_used + 2 * point_num, 1);
    err_sum = 0;
    for fid = 1 : pose_num
        if fid == host_fid
            continue;
        end
        for pid = 1 : point_num
            point_start_idx = 2*(pid-1) + 1;
            host_start_idx = 2 * point_num + 1 +(host_fid-1)*6;
            target_start_idx = 2 * point_num + 1 +(fid-1)*6;
            for sample = 1 : split_num
                [res, d_err_d_host_pose, d_err_d_target_pose, d_err_d_rou12] = computeResAndJac(Rwc{host_fid}, twc{host_fid}, Rwc{fid},  twc{fid}, PLK{fid}(:,pid), Obs{fid, pid}(sample,:)', dir_h1{pid}, dir_h2{pid}, rou1(pid), rou2(pid));
                err_sum = err_sum + abs(res);
                H = FillMatrix(H, 2, 2, point_start_idx, point_start_idx, d_err_d_rou12' * d_err_d_rou12);
                if fid <= pose_used
                    H = FillMatrix(H, 6, 6, host_start_idx, host_start_idx, d_err_d_host_pose' * d_err_d_host_pose);
                    H = FillMatrix(H, 6, 6, target_start_idx, target_start_idx, d_err_d_target_pose' * d_err_d_target_pose);
                    
                    H = FillMatrix(H, 6, 2, host_start_idx, point_start_idx, d_err_d_host_pose' * d_err_d_rou12);
                    H = FillMatrix(H, 6, 2, target_start_idx, point_start_idx, d_err_d_target_pose' * d_err_d_rou12);
                    H = FillMatrix(H, 6, 6, target_start_idx, host_start_idx, d_err_d_target_pose' * d_err_d_host_pose);
                    
                    
                    H = FillMatrix(H, 2, 6, point_start_idx, host_start_idx, d_err_d_rou12' * d_err_d_host_pose);
                    H = FillMatrix(H, 2, 6, point_start_idx, target_start_idx, d_err_d_rou12' * d_err_d_target_pose);
                    H = FillMatrix(H, 6, 6, host_start_idx,target_start_idx, d_err_d_host_pose' * d_err_d_target_pose);
                    
                end
                
                b = FillMatrix(b, 2, 1, point_start_idx, 1, -d_err_d_rou12' * res);
                if fid <= pose_used
                    b = FillMatrix(b, 6, 1, host_start_idx, 1, -d_err_d_host_pose' * res);
                    b = FillMatrix(b, 6, 1, target_start_idx, 1, -d_err_d_target_pose' * res);
                end
            end
        end
    end
    fprintf(sprintf('iter: %d, err:%f, rou12 - gt = [%f %f]\n', iter, err_sum, max(abs(rou1 - rou1_gt)), max(abs(rou2 - rou2_gt))));
    for idx = pose_num -1 : pose_num
        H = FillMatrix(H, 6, 6, 2 * point_num + 1 +(idx-1)*6, 2 * point_num + 1 +(idx-1)*6, 1e8 * eye(6));
    end
    dx = inv(H) * b;
    dxx = reshape(dx(1:2 * point_num), 2, [])';
    rou1 = rou1 + dxx(:,1)';
    rou2 = rou2 + dxx(:,2)';
    dx_pose = reshape(dx(2 * point_num+1:end), 6, [])';
    for id = 1 : pose_num
        dx1 = dx_pose(id,:);
        Told = [Rwc{id} twc{id};0 0 0 1];
        dPose = Exp(dx1(1:3)', dx1(4:6)');
        Tnew = dPose * Told;
        Rwc{id} = Tnew(1:3,1:3);
        twc{id} = Tnew(1:3,4);
    end
end


return;
%%



use_host_at_identity = false;
if use_host_at_identity
    plk_h = TransformTth( Rwc{1},  twc{1},  PLK(:,1));
    Rwh = eye(3);
    twh = zeros(3,1);
else
    plk_h = PLK(:,1);
    Rwh = Rwc{1};
    twh = twc{1};
end
err = [];
for i = 1 : 6
    [res1, d_err_d_host_pose_1, d_err_d_target_pose_1, d_err_d_rou12_1] = computeResAndJac(Rwh, twh, Rwc{i},  twc{i}, plk_h, bearings(i,1:3)', dir_h1, dir_h2, rou1, rou2);
    [res2, d_err_d_host_pose_2, d_err_d_target_pose_2, d_err_d_rou12_2] = computeResAndJac(Rwh, twh, Rwc{i},  twc{i}, plk_h, bearings(i,4:6)', dir_h1, dir_h2, rou1, rou2);
    [res3, d_err_d_host_pose_3, d_err_d_target_pose_3, d_err_d_rou12_3] = computeResAndJac(Rwh, twh, Rwc{i},  twc{i}, plk_h, bearings(i,7:9)', dir_h1, dir_h2, rou1, rou2);
    err = [err;[res1;res2; res3]];
end



Rwc_gt = Rwc;
twc_gt = twc;

rou1_gt = rou1;
rou2_gt = rou2;

pose_num = 4;
for i = 1 : pose_num
    Rwc{i} = Rwc{i} * rodrigues(0.01 * (rand(3,1) - 0.5));
    twc{i} = twc{i} + 10 * (rand(3,1) - 0.5);
end

rou1 = rou1 + 0.1;
rou2 = rou2 + 0.2;

host_fid = 1;

for iter = 1 : 10
    H = zeros(2 + 6 * pose_num, 2 + 6 * pose_num);
    b = zeros(2 + 6 * pose_num, 1);
    
    for i = 1 : 6
        if i == host_fid
            continue;
        end
        for j = 1 : 3
            [res, d_err_d_host_pose, d_err_d_target_pose, d_err_d_rou12] = computeResAndJac(Rwc{host_fid}, twc{host_fid}, Rwc{i},  twc{i}, plk_h, bearings(i, 3*j-2 : 3*j)', dir_h1, dir_h2, rou1, rou2);
            H = FillMatrix(H, 2, 2, 1, 1, d_err_d_rou12' * d_err_d_rou12);
            if i <= pose_num
                H = FillMatrix(H, 6, 6, 3+(host_fid-1)*6, 3+(host_fid-1)*6, d_err_d_host_pose' * d_err_d_host_pose);
                H = FillMatrix(H, 6, 6, 3+(i-1)*6, 3+(i-1)*6, d_err_d_target_pose' * d_err_d_target_pose);
                
                H = FillMatrix(H, 6, 2, 3+(host_fid-1)*6, 1, d_err_d_host_pose' * d_err_d_rou12);
                H = FillMatrix(H, 6, 2, 3+(i-1)*6, 1, d_err_d_target_pose' * d_err_d_rou12);
                H = FillMatrix(H, 6, 6, 3+(i-1)*6, 3+(host_fid-1)*6, d_err_d_target_pose' * d_err_d_host_pose);
                
                
                H = FillMatrix(H, 2, 6, 1, 3+(host_fid-1)*6, d_err_d_rou12' * d_err_d_host_pose);
                H = FillMatrix(H, 2, 6, 1, 3+(i-1)*6, d_err_d_rou12' * d_err_d_target_pose);
                H = FillMatrix(H, 6, 6, 3+(host_fid-1)*6,3+(i-1)*6, d_err_d_host_pose' * d_err_d_target_pose);
                
            end
            
            b = FillMatrix(b, 2, 1, 1, 1, -d_err_d_rou12' * res);
            if i <= pose_num
                b = FillMatrix(b, 6, 1, 3+(host_fid-1)*6, 1, -d_err_d_host_pose' * res);
                b = FillMatrix(b, 6, 1, 3+(i-1)*6, 1, -d_err_d_target_pose' * res);
            end
            
            
        end
    end
    
end




khsjaf = 1;
return;





% n v 是第一帧坐标系下的plk坐标。
%【这个例子中，第一帧在世界系下的pose并不是eye(4),12，13，14帧上的线三角化出的plk都是表达在1坐标系下的，并且结果都是一样的】
C0 = [n(:,1) v(:,1)];
kObvNum = 4; 3;
%  Eigen::Matrix<double, 2 * kObvNum, 6> A_temp;
A_temp = zeros(2*kObvNum, 6);
for  i = kObvNum : -1 : 1
    
    temp = zeros(3,6);
    if 1
        temp(:,1:3) = Rwc{i}';
        temp(:,4:6) = SkewSymMat(-Rwc{i}'*twc{i}) * Rwc{i}';
    else
        temp(:,1:3) = Rwc{i};
        temp(:,4:6) = SkewSymMat(twc{i}) * Rwc{i};
    end
    A_temp(2*i-1, : ) = [obs(i,1:2) 1] * temp;
    A_temp(i *2 ,:) = [obs(i,3:4) 1] * temp;
end
A = A_temp;
[U1,S1,V1] = svd(A,0);
para = V1(:,end);
para = para./norm(para(4:6));

% n v 是世界系下【可以理解成第0帧】的plk坐标。
nn = para(1:3);  % V1(1:3,end)./norm(V1(1:3,end));
vv = para(4:6); %V1(4:6,end)./norm(V1(4:6,end));

nn_ = nn;
nn_(1) = -(nn(2)*vv(2) + nn(3)*vv(3))/vv(1);

err = [n(:, 1) nn v(:, 1) vv];
dot(nn,vv)
C = [nn vv];
[U,S,V] = svd(C,0);
% Z的行向量是正交的
Z = S*V';
Z1 = Z(:,1);
Z2 = Z(:,2);
TT = [Z2'; Z1'*[0 -1; 1 0]];
T = [Z(2,1) Z(2,2); Z(1,2) -Z(1,1)];
[U_,S_,V_] = svd(T, 0);
[U_2,S_2,V_2] = svd(TT, 0);
if 1
    V11 = V_(1,end);
    V22 = V_(2,end);
    V11_ = V_2(1,end);
    V22_ = V_2(2,end);
else
    V11 = U_(1,end);
    V22 = U_(2,end);
end

VV = [V11 -V22;V22 V11];
S__ = diag(VV'*S*V');
Z_ = VV*diag(S__);
C_ = U*Z_;

VV_ = [V11_ -V22_;V22_ V11_];
S__1 = diag(VV_'*S*V');
Z_1 = VV_*diag(S__1);
C_1 = U*Z_1;
%         U2*V3*diag()
%
%         S2*VV'
%         V3*

if 0
    line_w = line_to_pose( [C_(:,1); C_(:,2)], Tcw1(1:3,1:3), Tcw1(1:3,4));
    plk_w = plk_w./norm(plk_w(4:6));
    CC = [line_w(1:3)./norm(line_w(1:3)) line_w(4:6)./norm(line_w(4:6))];
else
    %     plk_w = plk_to_pose( Rwc{1}',  -Rwc{1}'*twc{1},  [C0(:,1); C0(:,2)] );
    C00 = 1.3.*C0;
    plk_w = plk_to_pose( Rwc{1},  twc{1},  [C00(:,1); C00(:,2)] );
    %     CC = [plk_w(1:3)./norm(plk_w(1:3)) plk_w(4:6)./norm(plk_w(4:6))];
    plk_w = plk_w./norm(plk_w(4:6));
    CC = [plk_w(1:3) plk_w(4:6)];
end


C0;
[C C_]
norm(C(:) - C_(:))


AA = A;
[a2,b2,c2] = svd(AA,0);
[a1,b1,c1] = svd(AA);
inv(a1)*AA*inv(c1') - b1;
pinv(a2)*AA*inv(c2') - b2;

%figure,line([uvS1(1) uvE1(1) uvS2(1) uvE2(1)], [uvS1(2)  uvE1(2) uvS2(2) uvE2(2)],'Color', [1 0 0])
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
function [res, d_err_d_Thost, d_err_d_Ttarget, d_err_d_rou12] = computeResAndJac(Rwh, twh, Rwt, twt, plk_h_check, obs, dir_h1, dir_h2, rou1, rou2)

p1 = dir_h1./rou1;
p2 = dir_h2./rou2;
ptDiff = p2 - p1;
dir = ptDiff./norm(ptDiff);

d_vh_d_ptDiff = compute_d_bearing_d_pt_jac(ptDiff);

d_ptDiff_d_p1 = -eye(3);
d_ptDiff_d_p2 = eye(3);

d_p1_d_rou1 = compute_d_xyz_d_rou(dir_h1, rou1);
d_p2_d_rou2 = compute_d_xyz_d_rou(dir_h2, rou2);

d_vh_d_p1 = d_vh_d_ptDiff * d_ptDiff_d_p1;
d_vh_d_p2 = d_vh_d_ptDiff * d_ptDiff_d_p2;

d_vh_d_rou1 = d_vh_d_p1 * d_p1_d_rou1;
d_vh_d_rou2 = d_vh_d_p2 * d_p2_d_rou2;

% d_vh_dp12 = d_vh_d_ptDiff * [d_ptDiff_d_p1 d_ptDiff_d_p2];

nh = cross(p1, dir);


d_nh_d_dir = SkewSymMat(p1);
d_dir_d_ptDiff = compute_d_bearing_d_pt_jac(ptDiff);
d_nh_d_p1 = -SkewSymMat(dir) + d_nh_d_dir * d_dir_d_ptDiff * d_ptDiff_d_p1;
d_nh_d_p2 = d_nh_d_dir * d_dir_d_ptDiff * d_ptDiff_d_p2;


d_nh_d_rou1 = d_nh_d_p1 * d_p1_d_rou1;
d_nh_d_rou2 = d_nh_d_p2 * d_p2_d_rou2;

plk_h = [nh; dir];

check = abs(plk_h) - abs(plk_h_check);
% assert(max(abs(check)) < 1e-8);


vh = dir;


Rth = Rwt' * Rwh;
tth = Rwt' * (twh - twt);
d_err_d_host_pose = [];
d_err_d_target_pose = [];
d_err_d_rou12 = [];

plk_t = TransformTth( Rth,  tth,  plk_h);
plane_n = plk_t(1:3) ./ norm(plk_t(1:3));
res = dot(obs,plane_n);


d_err_d_plane_n = obs'; % 1x3

d_plane_n_d_line_t_n = compute_d_bearing_d_pt_jac(plk_t(1:3)); % 3x3

d_line_t_n_d_line_t = [eye(3) zeros(3,3)]; %3x6


nh_check = plk_h_check(1:3);
vh_check = plk_h_check(4:6);

right_perturb = false;

if ~right_perturb % left disturb
    d_nt_d_Rth = -SkewSymMat(Rth * nh) - SkewSymMat(tth) * SkewSymMat(Rth * vh) + SkewSymMat(Rth * vh) * SkewSymMat(tth);
    d_nt_d_tth = -SkewSymMat(Rth * vh);
    d_vt_d_Rth = -SkewSymMat(Rth * vh);
    d_vt_d_tth = zeros(3,3);
else
    d_nt_d_Rth = -Rth * SkewSymMat(nh) - SkewSymMat(tth) * Rth * SkewSymMat(vh);
    d_nt_d_tth = -SkewSymMat(Rth * vh) * Rth;
    d_vt_d_Rth = -Rth * SkewSymMat(vh);
    d_vt_d_tth = zeros(3,3);
end
d_line_t_d_Tth = [d_nt_d_Rth d_nt_d_tth;...
    d_vt_d_Rth d_vt_d_tth];
d_line_t_d_line_h = [Rth SkewSymMat(tth) * Rth; zeros(3,3) Rth];

d_line_h_d_rou12 = [d_nh_d_rou1 d_nh_d_rou2;
    d_vh_d_rou1 d_vh_d_rou2];

d_err_d_Tth = d_err_d_plane_n * d_plane_n_d_line_t_n * d_line_t_n_d_line_t * d_line_t_d_Tth;
d_err_d_rou12 = d_err_d_plane_n * d_plane_n_d_line_t_n * d_line_t_n_d_line_t * d_line_t_d_line_h * d_line_h_d_rou12;

if ~right_perturb
    d_err_d_Thost = d_err_d_Tth * Adj(inv([Rwt twt;0 0 0 1]));
    d_err_d_Ttarget = -d_err_d_Tth * Adj(inv([Rwt twt;0 0 0 1]));
else
    d_err_d_Thost = d_err_d_Tth;
    d_err_d_Ttarget = -d_err_d_Tth * Adj(inv([Rwh twh;0 0 0 1]));
end

end
function d_rou_d_xyz = compute_d_rou_d_xyz(xyz)
d_rou_d_xyz = -0.5 * norm(xyz)^-3 * 2 * [xyz'];
end
function d_xyz_d_rou = compute_d_xyz_d_rou(dir, rou)

d_xyz_d_rou = -dir / rou / rou;
end
function jac = Adj(T)
R = T(1:3,1:3);
jac = zeros(6,6);
jac(1:3,1:3) = R;
jac(4:6, 4:6) = R;
jac(4:6,1:3) = SkewSymMat(T(1:3,4)) * R;
jac(1:3, 4:6) = zeros(3, 3);

end
function d_bearing_d_pt = compute_d_bearing_d_pt_jac(pt)
d_bearing_d_pt = (norm(pt).*eye(3) - pt * pt'./norm(pt))./(norm(pt)^2);
end
function pi = pi_from_ppp(x1, x2, x3)

pi = [cross( x1 - x3 ,( x2 - x3 )); dot(- x3,( cross(x1, x2 )) )]; %// d = - x3.dot( (x1-x3).cross( x2-x3 ) ) = - x3.dot( x1.cross( x2 ) )];


end
function plk = pipi_plk(pi1, pi2)

dp = pi1 * pi2' - pi2 * pi1';

plk = [dp(1,4), dp(2,4), dp(3,4), - dp(2,3), dp(1,3), - dp(1,2)]';
end
function line_c = line_to_pose( line_w, Rcw, tcw)



cp_w = line_w(1:3);
dv_w = line_w(4:6);


cp_c = Rcw * cp_w + tcw;
dv_c = Rcw* dv_w;

line_c(1:3,1) = cp_c;
line_c(4:6,1) = dv_c;


end
function plk_t = TransformTth( Rth,  tth,  plk_h )
nh = plk_h(1:3);
vh = plk_h(4:6);

nt = Rth * nh + SkewSymMat(tth) * Rth * vh;
vt = Rth * vh;


plk_t(1:3,1) = nt;
plk_t(4:6,1) = vt;
end
function plk_c = plk_to_pose( Rcw,  tcw,  plk_w )
nw = plk_w(1:3);
vw = plk_w(4:6);

nc = Rcw * nw + SkewSymMat(tcw) * Rcw * vw;
vc = Rcw * vw;


plk_c(1:3,1) = nc;
plk_c(4:6,1) = vc;



end