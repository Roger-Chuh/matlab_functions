function testPatchNormalOpt()


add_noise = true;

use_idp_res = true;
opt_idp_main = true;
opt_normal = true;
use_fast_livo2_2dof = false;

[xMat, yMat] = meshgrid(-5:5, -5:5);

xy_center = [50 75];
z = 100;

xyz = [xMat(:) yMat(:)] + xy_center;
xyz(:,3) = z;

xyz_center = [xy_center z];

center_id = find(xyz(:,1) == xyz_center(:,1) & xyz(:,2) == xyz_center(:,2) & xyz(:,3) == xyz_center(:,3));

xyz = (rodrigues(0.1 * (rand(3,1)-0.5)) * xyz')' + repmat(100*(rand(1,3)-0.5), size(xyz,1),1);
xyz_center = xyz(center_id,:);


idp_main_gt = 1/norm(xyz_center);
bearing_main = xyz_center' * idp_main_gt;

cam_num = 6;
for i = 1 : cam_num
   if i == 1
       Tbc{i} = eye(4);
   else
       Tbc{i} = [rodrigues(0.1 * (rand(3,1)-0.1)) 10 * (rand(3,1)-0.5);0 0 0 1];
   end
end

xyz_fit = xyz;
xyz_fit(:,4) = 1;
[u,v,w] = svd(xyz_fit);
plane = w(:,4);
plane = plane./norm(plane(1:3));
xyz = xyz;

patch_normal_gt = plane(1:3);

fit_err = xyz_fit * plane;


if add_noise
    if opt_idp_main && opt_normal
        idp_main_opt = idp_main_gt + 0.001;
        patch_normal_opt = rodrigues([0.1 -0.1 0.1]') * patch_normal_gt;
    end
    if opt_idp_main && ~opt_normal
        idp_main_opt = idp_main_gt + 0.001;
        patch_normal_opt = patch_normal_gt;
    end
    if ~opt_idp_main && opt_normal
        idp_main_opt = idp_main_gt;
        patch_normal_opt = rodrigues([0.1 -0.1 0.1]') * patch_normal_gt;
    end
else
    idp_main_opt = idp_main_gt;
    patch_normal_opt = patch_normal_gt;
end

pose_num = 400;

for i = 1 : pose_num
    if i == 1
        Twb{i,1} = eye(4);
    else
        Twb{i,1} = [rodrigues(0.5 * (rand(3,1)-0.5)) 100 * (rand(3,1)-0.5);0 0 0 1];
    end
end
Twch = Twb{1,1} * Tbc{1};
for fid = 1 : pose_num
    for cid = 1 : cam_num
        for pid = 1 : size(xyz,1)
            Twct = Twb{fid,1} * Tbc{cid};
            Tth = inv(Twct) * Twch;
            xyz_ccs = Tth(1:3, 1:3) * xyz(pid,:)' + Tth(1:3, 4);
            idp_each_obs = 1/norm(xyz_ccs);
            bearing_each = xyz_ccs * idp_each_obs;
            idp_each_obs_stack{fid,pid}{cid} = idp_each_obs;
            bearing_each_stack{fid,pid}{cid} = bearing_each;
        end
    end
end

Err = [];
for iter = 1 : 20
    
    H = zeros(3,3);
    b = zeros(3,1);
    err_sum = 0;
    
    for fid = 2 : pose_num
        Twch = Twb{1,1} * Tbc{1};
        bearing_main_host = bearing_each_stack{1, center_id}{1};
        for cid = 1 : cam_num
            for pid = 1 : length(size(xyz, 1))
                if 0
                    xyz_ccs = Tcw{fid,1}(1:3, 1:3) * xyz(pid,:)' + Tcw{fid,1}(1:3, 4);
                    idp_each_obs = 1/norm(xyz_ccs);
                    bearing_each = xyz_ccs * idp_each_obs;
                else
                    idp_each_obs = idp_each_obs_stack{fid, pid}{cid};
                    bearing_each = bearing_each_stack{fid, pid}{cid};
                end
                bearing_each_host = bearing_each_stack{1, pid}{1};
                bearing_each_obs = bearing_each_stack{fid, pid}{cid};
                Twct = Twb{fid,1} * Tbc{cid};
                Tth = inv(Twct) * Twch;
                [err, d_err_each_d_idp_main, d_err_each_d_normal_2dof] = computeJac(Tth,idp_main_opt, patch_normal_opt, idp_each_obs, bearing_each_host, bearing_main_host, bearing_each_obs, use_fast_livo2_2dof, use_idp_res);
                err_sum = err_sum + norm(err);
                H(1:2, 1:2) = H(1:2, 1:2) + d_err_each_d_normal_2dof' * d_err_each_d_normal_2dof;
                H(3, 3) = H(3, 3) + d_err_each_d_idp_main' * d_err_each_d_idp_main;
                H(1:2,3) = H(1:2,3) + d_err_each_d_normal_2dof' * d_err_each_d_idp_main;
                H(3, 1:2) = H(3,1:2) + d_err_each_d_idp_main' * d_err_each_d_normal_2dof;
                b(1:2) = b(1:2) - d_err_each_d_normal_2dof' * err;
                b(3) = b(3) - d_err_each_d_idp_main' * err;
            end
        end
    end
    Err = [Err; err_sum];
    fprintf(sprintf('iter: %d, err_sum: %f\n', iter, err_sum));
    for k = 1: 2
       H(k,k) = H(k,k) + 1000 * H(k,k);
    end
    if opt_idp_main && opt_normal
        dx = inv(H) * b;
        if ~use_fast_livo2_2dof
            rot_vec = ProduceOtherOthogonalBasis(patch_normal_opt) * dx(1:2);
            patch_normal_opt = rodrigues(rot_vec) * patch_normal_opt;
        else
            [B_old, m_old, b_old, M_check_old] = GetB(patch_normal_opt, bearing_main_host, idp_main_opt);
            m_new = m_old + dx(1:2);
            M_check_new = B_old * m_new + b_old;
        end
        idp_main_opt = idp_main_opt + dx(3);
    end
   
    if opt_idp_main && ~opt_normal
        dx = inv(H(3,3)) * b(3);
        idp_main_opt = idp_main_opt + dx;
    end
    if ~opt_idp_main && opt_normal
        dx = inv(H(1:2,1:2)) * b(1:2);
        rot_vec = ProduceOtherOthogonalBasis(patch_normal_opt) * dx(1:2);
        patch_normal_opt = rodrigues(rot_vec) * patch_normal_opt;
    end
end


Err'
end
function [B, m, b, M_check] = GetB(patch_normal_opt, bearing_main_host, idp_main_opt)
xyz_main_host = bearing_main_host./idp_main_opt;

M = patch_normal_opt / (patch_normal_opt' * xyz_main_host);
Mx = M(1);
My = M(2);
m = [Mx;My];
B = [1 0; 0 1;-xyz_main_host(1)/xyz_main_host(3) -xyz_main_host(2)/xyz_main_host(3)];
b = [0;0;1/xyz_main_host(3)];
M_check = B * m + b;
end
function [err, d_err_each_d_idp_main, d_err_each_d_normal_2dof] = computeJac(Tth,idp_main_opt, patch_normal_opt, idp_each_obs, bearing_each_host, bearing_main_host, bearing_each_obs, use_fast_livo2_2dof, use_idp_res)

res_amplify_ratio = 1;

if 0
    xyz_main_host = bearing_main_host./idp_main_opt;
    
    M = patch_normal_opt / (patch_normal_opt' * xyz_main_host);
    Mx = M(1);
    My = M(2);
    m = [Mx;My];
    B = [1 0; 0 1;-xyz_main_host(1)/xyz_main_host(3) -xyz_main_host(2)/xyz_main_host(3)];
    b = [0;0;1/xyz_main_host(3)];
    M_check = B * m + b;
else
    [B, m, b, M_check] = GetB(patch_normal_opt, bearing_main_host, idp_main_opt);
end
if ~use_fast_livo2_2dof
    idp_each_check = bearing_each_host' * M_check;
    idp_each = idp_main_opt * ((bearing_each_host' * patch_normal_opt) /(bearing_main_host' * patch_normal_opt));
else
    idp_each = bearing_each_host' * M_check;
    idp_each_check = idp_main_opt * ((bearing_each_host' * patch_normal_opt) /(bearing_main_host' * patch_normal_opt));
end
xyz_each_host = bearing_each_host / idp_each;
xyz_each_target = Tth(1:3,1:3) *xyz_each_host + Tth(1:3,4);
if 0
    err = 1 ./ norm(xyz_each_target) - idp_each_obs;
else
    err = xyz_each_target./norm(xyz_each_target) - bearing_each_obs;
end

d_idp_each_d_idp_main = ((bearing_each_host' * patch_normal_opt) / (bearing_main_host' * patch_normal_opt));

d_err_d_xyz_target = compute_d_bearing_d_pt_jac(xyz_each_target);
d_xyz_target_d_xyz_host = Tth(1:3,1:3);
d_xyz_host_d_idp_each = compute_d_xyz_d_rou(bearing_each_host, idp_each);
d_err_each_d_idp_main = d_err_d_xyz_target * d_xyz_target_d_xyz_host * d_xyz_host_d_idp_each * d_idp_each_d_idp_main;


d_idp_each_d_normal = idp_main_opt * bearing_each_host' / (bearing_main_host' * patch_normal_opt) +...
    idp_main_opt * bearing_each_host' * patch_normal_opt * (-(bearing_main_host' * patch_normal_opt)^-2) * bearing_main_host';

d_idp_each_d_normal_2dof = d_idp_each_d_normal * (-SkewSymMat(patch_normal_opt)) * ProduceOtherOthogonalBasis(patch_normal_opt);
d_err_each_d_normal_2dof = d_err_d_xyz_target * d_xyz_target_d_xyz_host * d_xyz_host_d_idp_each * d_idp_each_d_normal_2dof;
d_idp_each_d_normal_fast_livo_2dof = bearing_each_host' * B;
if use_fast_livo2_2dof
   d_err_each_d_normal_2dof =  d_err_d_xyz_target * d_xyz_target_d_xyz_host * d_xyz_host_d_idp_each * d_idp_each_d_normal_fast_livo_2dof;
end



if use_idp_res
    err = 1 ./ norm(xyz_each_target) - idp_each_obs;
    d_err_d_xyz_target = compute_d_rou_d_xyz(xyz_each_target);
    d_err_each_d_idp_main = d_err_d_xyz_target * d_xyz_target_d_xyz_host * d_xyz_host_d_idp_each * d_idp_each_d_idp_main;
    d_err_each_d_normal_2dof = d_err_d_xyz_target * d_xyz_target_d_xyz_host * d_xyz_host_d_idp_each * d_idp_each_d_normal_2dof;
end
err = res_amplify_ratio * err;
d_err_each_d_idp_main = res_amplify_ratio * d_err_each_d_idp_main;
d_err_each_d_normal_2dof = res_amplify_ratio * d_err_each_d_normal_2dof;
end
function d_rou_d_xyz = compute_d_rou_d_xyz(xyz)
d_rou_d_xyz = -0.5 * norm(xyz)^-3 * 2 * [xyz'];
end
function d_xyz_d_rou = compute_d_xyz_d_rou(dir, rou)

d_xyz_d_rou = -dir / rou / rou;
end
function d_bearing_d_pt = compute_d_bearing_d_pt_jac(pt)
d_bearing_d_pt = (norm(pt).*eye(3) - pt * pt'./norm(pt))./(norm(pt)^2);
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