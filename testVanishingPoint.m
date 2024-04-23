function testVanishingPoint(varargin)
close all

if nargin < 1
    use_vec3_err = false;
elseif nargin == 1
     use_vec3_err = varargin{1};
end


change_err_order = false;true;
% use_vec3_err = true;

colorStack = {'k';'m';'c'};
colorStack = repmat(colorStack,100,1);

L2R = [rodrigues([0.02; 0.01; -0.003]) [-50.1 0.01 -0.1]';0 0 0 1];

markerSize = 800; 150; 100; 200;
markerRow = 5;3; 15; 4; 7; 15;20; 7;
markerCol = 5;3; 15; 4; 7; 15;20; 7;

[xMat0, yMat0] = meshgrid(markerSize : markerSize : (markerCol - 0) * markerSize, markerSize : markerSize : (markerRow - 0) * markerSize);
xMat0 = xMat0';
yMat0 = yMat0';
marker1 = [xMat0(:) yMat0(:) zeros(markerRow*markerCol,1)];

Cam2MarkerRot = rodrigues([-0.01 1.0 -0.02] + 0.05*(rand(1,3)-0.5));
Cam2MarkerRot = rodrigues([-0.01 1.0 -0.02] + 0.00*(rand(1,3)-0.5));
Cam2MarkerTrans = [-100 6500 -8000]';
Cam2MarkerTrans = [0 0 -8000]';
Cam2Marker = [Cam2MarkerRot Cam2MarkerTrans; 0 0 0 1];
Marker2CamL = inv(Cam2Marker);
Marker2CamR = L2R*Marker2CamL;


width = 640; height = 480;
K = eye(3);
K(1,1) = 236.8820281-0;
K(2,2) = 236.7437341-0;
K(1,3) = 325.0065396;
K(2,3) = 241.9557923;
dist = [0.2124214875, -0.1719834148, 0.06211290885, -0.009899091649]';
% dist = [0.06314965524 -0.07006363688 0.01760113423 0.0001556135281]';
% dist = [0.46314965524 -0.47006363688 0.1760113423 0.001556135281]';

param = [K(1,1) K(2,2) K(1,3) K(2,3) dist']';


param = [234.9461165 235.0992169 331.2555939 233.9358652 0.2195841015 -0.1960459177 0.1109940203 -0.048751368 0.01353063659 -0.001607160663 -0.007680295854 0.004871400778 0.008343518773 -0.002910213568 -0.006609047968 0.0008683913441 0.001013424359 -0.0003548134221 -0.01310959013 -0.02949166679 ]';



T00 = eye(4);
T01 = [ -7.4656856642149472e-01, -7.0503228692542969e-02, -6.6156229515842346e-01, -3.5247847596384795e-02;
    -6.2170843696788058e-01, -2.8011666159265608e-01,7.3144601666913212e-01, 4.1311291893071443e-02;
    -2.3688392734482649e-01, 9.5737346455921379e-01,1.6529384236437439e-01, -3.4150037619299412e-02;
    0., 0., 0., 1. ];

T02 = [ -4.3634354777464124e-01, 4.7520621068004289e-02, -8.9852440083102414e-01, -6.8945219190976056e-02;
    8.9764389729039329e-01, 9.1761591914208740e-02, -4.3106292337277669e-01, -1.1696462645325466e-01;
    6.1965651555951469e-02, -9.9464647027106146e-01, -8.2696174062438135e-02, -6.6681057254321219e-02;
    0., 0., 0., 1. ];

T03 = [ 9.2370408189253994e-01, 3.7671281461209621e-01, -6.9700964140337329e-02, -2.1413263321725769e-02;
    -3.7601452313871109e-01, 8.5662781853225878e-01, -3.5327306563822758e-01, -1.0905148943810981e-01;
    -7.3374706022088163e-02, 3.5252834754625934e-01, 9.3292010198755815e-01, -2.2677364423121626e-02;
    0., 0., 0., 1. ];

T01(1:3,4) = 1000.*T01(1:3,4);
T02(1:3,4) = 1000.*T02(1:3,4);
T03(1:3,4) = 1000.*T03(1:3,4);

[xMat, yMat] = meshgrid(1 : width, 1 : height);
pix = [xMat(:) yMat(:)];

RotVec = [];
uMat_stack = {};
vMat_stack = {};
Marker2CamL0 = Marker2CamL;
for frame_id = 1 : 20
    Marker2CamL = [rodrigues(0.5 * (rand(3,1)-0.5)) 100 * (rand(3,1)-0.5);0 0 0 1] * Marker2CamL0;
    [pt2dL00, pt3dL] = TransformAndProject(marker1, K, Marker2CamL(1:3,1:3), Marker2CamL(1:3,4));
    
    for i = 1 : size(marker1,1)
        [pt2d_show(i,:)] = ProjectKB20(pt3dL(i,1), pt3dL(i,2), pt3dL(i,3),param(1),param(2),param(3),param(4),param(5),param(6),param(7),param(8),param(9),param(10),param(11),param(12),param(13),param(14),param(15),param(16),param(17),param(18),param(19),param(20));
        [pt3d_show(i,:)] = UnProjectKB20(pt2d_show(i,1), pt2d_show(i,2),param(1),param(2),param(3),param(4),param(5),param(6),param(7),param(8),param(9),param(10),param(11),param(12),param(13),param(14),param(15),param(16),param(17),param(18),param(19),param(20));
        [pt2d_show_check(i,:)] = ProjectKB20(pt3d_show(i,1), pt3d_show(i,2), pt3d_show(i,3),param(1),param(2),param(3),param(4),param(5),param(6),param(7),param(8),param(9),param(10),param(11),param(12),param(13),param(14),param(15),param(16),param(17),param(18),param(19),param(20));
    end
    figure(20),imshow(zeros(480, 640));hold on; plot(pt2d_show(:,1), pt2d_show(:,2), 'or','MarkerSize',5,'LineWidth',5);
    plot(pt2d_show(1,1), pt2d_show(1,2), 'xg','MarkerSize',10,'LineWidth',10);
    Rcw = Marker2CamL(1:3,1:3);
    rot_vec = rodrigues(Rcw);
    RotVec = [RotVec;rot_vec];
    uMat = reshape(pt2d_show(:,1), markerRow, markerRow);
    vMat = reshape(pt2d_show(:,2), markerRow, markerRow);
    uMat_stack{frame_id,1} = uMat;
    vMat_stack{frame_id,1} = vMat;
end



computeError([param; RotVec], uMat_stack, vMat_stack, use_vec3_err, change_err_order);





end
function computeError(param, uMat_stack, vMat_stack, use_vec3_err, change_err_order)
param0 = param;
offset = 0;2;
fix_f = false;
fix_center = false;
fix_rot = false;
add_noise = true;

frame_num = length(uMat_stack);

if add_noise
    param(offset+1:20) = param(offset+1:20) + 0.001;
end
if fix_f
    param(1:2) = param0(1:2);
    if add_noise
%         param(3:4) = param0(3:4) + [10 -10]';
    end
else
    if add_noise
        param(1:2) = param0(1:2) + [5 -5]';
    end
end
if fix_center
    param(3:4) = param0(3:4);
    if add_noise
%         param(3:4) = param0(3:4) + [10 -10]';
    end
else
    if add_noise
        param(3:4) = param0(3:4) + [5 -5]';
    end
end
if add_noise && ~fix_rot
    for frame_id = 1 : frame_num
        param(21 + 3*(frame_id-1):23 + 3*(frame_id-1)) = param(21 + 3*(frame_id-1):23 + 3*(frame_id-1)) + (0.05 * (rand(3,1)-0.5)); [-0.01 0.02 -0.03]';
    end
end
% rotMat = rodrigues(param(end-2:end));
% intr = param(1:end-3);
% combinations = combnk(1:size(uMat,1), 2);
comb = nchoosek(1:size(uMat_stack{1,1},1), 2);


loss = [];
for iter = 1 : 10
    err = 0;
    
    intr = param(1:20);
    
    H = zeros(length(param), length(param));
    b = zeros(length(param),1);
    rotMat_stack = {};
    for frame_id = 1 : frame_num
        uMat = uMat_stack{frame_id, 1};
        vMat = vMat_stack{frame_id, 1};
        rotMat = rodrigues(param(21 + 3*(frame_id-1):23 + 3*(frame_id-1)));
        X1 = rotMat(:,1)';
        X2 = rotMat(:,2)';
        rotMat_stack{frame_id,1} = rotMat;
        for i = 1 : size(uMat, 1)
            uv_list_hori = [uMat(:,i) vMat(:,i)];
            uv_list_vert = [uMat(i,:); vMat(i,:)]';
            for j = 1 : size(comb,1)
                if 1
                    uv1_hori = uv_list_hori(comb(j,1),:);
                    uv2_hori = uv_list_hori(comb(j,2),:);
                    [hori_bearing1, ~, d_pt3d_d_param_hori1] = UnProjectKB20(uv1_hori(1), uv1_hori(2),intr(1),intr(2),intr(3),intr(4),intr(5),intr(6),intr(7),intr(8),intr(9),intr(10),intr(11),intr(12),intr(13),intr(14),intr(15),intr(16),intr(17),intr(18),intr(19),intr(20));
                    [hori_bearing2, ~, d_pt3d_d_param_hori2] = UnProjectKB20(uv2_hori(1), uv2_hori(2),intr(1),intr(2),intr(3),intr(4),intr(5),intr(6),intr(7),intr(8),intr(9),intr(10),intr(11),intr(12),intr(13),intr(14),intr(15),intr(16),intr(17),intr(18),intr(19),intr(20));
                    
                    d_unproject_d_bearing_hori1 = compute_d_bearing_d_pt_jac(hori_bearing1');
                    d_unproject_d_bearing_hori2 = compute_d_bearing_d_pt_jac(hori_bearing2');
                    
                    if fix_f
                        d_pt3d_d_param_hori1(:,1:2) = 0;
                        d_pt3d_d_param_hori2(:,1:2) = 0;
                    end
                    if fix_center
                        d_pt3d_d_param_hori1(:,3:4) = 0;
                        d_pt3d_d_param_hori2(:,3:4) = 0;
                    end
                    d_pt3d_d_param_hori1 = d_pt3d_d_param_hori1(:,1:end-offset);
                    d_pt3d_d_param_hori2 = d_pt3d_d_param_hori2(:,1:end-offset);
                    if 0
                        d_pt3d_d_param_hori1 = d_unproject_d_bearing_hori1 * d_pt3d_d_param_hori1;
                        d_pt3d_d_param_hori2 = d_unproject_d_bearing_hori2 * d_pt3d_d_param_hori2;
                    end
                    hori_plane_norm = cross(hori_bearing1, hori_bearing2);
                    hori_norm = norm(hori_plane_norm);
                    hori_plane_norm_change_order = cross(X1, hori_bearing2);
                    hori_norm_change_order = norm(hori_plane_norm_change_order);
                    %         hori_plane_norm = hori_plane_norm./hori_norm; 这里不能归一化
                    hori_plane_norm_check = cross(hori_bearing1, X1);
                    hori_plane_norm_check = hori_plane_norm_check./norm(hori_plane_norm_check);
                    if 0
                        assert(norm(hori_plane_norm_check - hori_plane_norm./hori_norm) < 0.00001);
                    end
                    err_hori = dot(hori_plane_norm, X1) / hori_norm; %真正的残差要归一化% <----
                    res_hori = dot(hori_plane_norm, X1); %这是中间结果，不要归一化
                    err_hori_change_order = dot(hori_plane_norm_change_order, hori_bearing1) / hori_norm_change_order; %真正的残差要归一化% <----
                    res_hori_change_order = dot(hori_plane_norm_change_order, hori_bearing1); %这是中间结果，不要归一化
                    d_hori_plane_d_param = [-SkewSymMat(hori_bearing2) * d_pt3d_d_param_hori1 + SkewSymMat(hori_bearing1) * d_pt3d_d_param_hori2];
                    
                    if ~use_vec3_err
                        if 0
                            d_hori_err_d_hori_plane1 = [X1(1)/hori_norm - res_hori * hori_plane_norm(1) / hori_norm^3;
                                X1(2)/hori_norm - res_hori * hori_plane_norm(2) / hori_norm^3;
                                X1(3)/hori_norm - res_hori * hori_plane_norm(3) / hori_norm^3]';
                        else
                            d_hori_err_d_hori_plane = compute_d_dist_d_plane(res_hori, X1, hori_plane_norm);
                            d_hori_err_d_hori_plane_change_order = compute_d_dist_d_plane(res_hori_change_order, hori_bearing1, hori_plane_norm_change_order);
                        end
                        d_hori_err_d_param = d_hori_err_d_hori_plane * d_hori_plane_d_param;% <----
                        d_hori_err_d_X1 = hori_plane_norm./hori_norm;
                        d_X1_d_R = -SkewSymMat(X1);
                        d_hori_err_d_R = d_hori_err_d_X1 * d_X1_d_R; % <---
                     %%
                        d_hori_err_d_bearing2 = d_hori_err_d_hori_plane_change_order * SkewSymMat(X1) * d_pt3d_d_param_hori2;
                        d_hori_err_d_bearing1 = hori_plane_norm_change_order./hori_norm_change_order * d_pt3d_d_param_hori1;
                        d_hori_err_d_param_change_order = d_hori_err_d_bearing1 + d_hori_err_d_bearing2;
                        d_X1_d_R_change_order = -SkewSymMat(X1);
                        d_hori_err_d_R_change_order = d_hori_err_d_hori_plane_change_order * (-SkewSymMat(hori_bearing2)) * (d_X1_d_R_change_order);
                        if change_err_order
                            err_hori = err_hori_change_order;
                            d_hori_err_d_param = d_hori_err_d_param_change_order;
                            d_hori_err_d_R = d_hori_err_d_R_change_order;
                        end
                    else
                        normalized_plane_hori = hori_plane_norm'./hori_norm;
                        pt_in_plane_hori = X1' - (X1 * normalized_plane_hori) .* normalized_plane_hori;
                        pt_in_plane_norm_hori = pt_in_plane_hori./norm(pt_in_plane_hori);
                        err_hori = pt_in_plane_norm_hori - X1';
                        d_err_d_pt_hori = compute_d_bearing_d_pt_jac(pt_in_plane_hori);
                        d_pt_d_plane_norm_hori = compute_d_pt_d_plane_norm_jac(X1', normalized_plane_hori); 
                        d_plane_norm_d_plane_hori = compute_d_bearing_d_pt_jac(hori_plane_norm');
                        d_hori_err_d_hori_plane = d_err_d_pt_hori * d_pt_d_plane_norm_hori * d_plane_norm_d_plane_hori;
                        d_hori_err_d_param = d_hori_err_d_hori_plane * d_hori_plane_d_param;% <----
                        d_pt_in_plane_d_X1 = [1-normalized_plane_hori(1)^2 -normalized_plane_hori(1)*normalized_plane_hori(2) -normalized_plane_hori(1)*normalized_plane_hori(3);
                            -normalized_plane_hori(1)*normalized_plane_hori(2) 1-normalized_plane_hori(2)^2 -normalized_plane_hori(2)*normalized_plane_hori(3);
                            -normalized_plane_hori(1)*normalized_plane_hori(3) -normalized_plane_hori(2)*normalized_plane_hori(3) 1-normalized_plane_hori(3)^2];
                        d_hori_err_d_X1 = -eye(3) + d_err_d_pt_hori * d_pt_in_plane_d_X1;
                        d_X1_d_R = -SkewSymMat(X1);
                        d_hori_err_d_R = d_hori_err_d_X1 * d_X1_d_R; % <---
                        %                 d_hori_err_d_R = zeros(size(d_hori_err_d_R));
                     %%
                        normalized_plane_hori_change_order = hori_plane_norm_change_order'./hori_norm_change_order;
                        pt_in_plane_hori_change_order = hori_bearing1' - (hori_bearing1 * normalized_plane_hori_change_order) .* normalized_plane_hori_change_order;
                        pt_in_plane_norm_hori_change_order = pt_in_plane_hori_change_order./norm(pt_in_plane_hori_change_order);
                        err_hori_change_order = pt_in_plane_norm_hori_change_order - hori_bearing1';
                        d_err_d_pt_hori_change_order = compute_d_bearing_d_pt_jac(pt_in_plane_hori_change_order);
                        d_pt_d_plane_norm_hori_change_order = compute_d_pt_d_plane_norm_jac(hori_bearing1', normalized_plane_hori_change_order); 
                        d_plane_norm_d_plane_hori_change_order = compute_d_bearing_d_pt_jac(hori_plane_norm_change_order');
                        
                        d_X1_d_R_chnage_order = -SkewSymMat(X1);
                        d_pt_in_plane_d_hori_bearing1 = [1-normalized_plane_hori_change_order(1)^2 -normalized_plane_hori_change_order(1)*normalized_plane_hori_change_order(2) -normalized_plane_hori_change_order(1)*normalized_plane_hori_change_order(3);
                            -normalized_plane_hori_change_order(1)*normalized_plane_hori_change_order(2) 1-normalized_plane_hori_change_order(2)^2 -normalized_plane_hori_change_order(2)*normalized_plane_hori_change_order(3);
                            -normalized_plane_hori_change_order(1)*normalized_plane_hori_change_order(3) -normalized_plane_hori_change_order(2)*normalized_plane_hori_change_order(3) 1-normalized_plane_hori_change_order(3)^2];
                        if 0
                            d_pt_in_plane_d_hori_bearing1_check = eye(3) - normalized_plane_hori_change_order * normalized_plane_hori_change_order';
                        end
                        
                        d_hori_err_d_bearing1 = (-eye(3) + d_err_d_pt_hori_change_order * d_pt_in_plane_d_hori_bearing1) * d_pt3d_d_param_hori1;
                        d_hori_err_d_X1_change_order = d_err_d_pt_hori_change_order * d_pt_d_plane_norm_hori_change_order * d_plane_norm_d_plane_hori_change_order * (-SkewSymMat(hori_bearing2));
                        d_hori_err_d_R_change_order = d_hori_err_d_X1_change_order * d_X1_d_R_chnage_order;
                        d_hori_err_d_bearing2 = d_err_d_pt_hori_change_order * d_pt_d_plane_norm_hori_change_order * d_plane_norm_d_plane_hori_change_order * (SkewSymMat(X1) * d_pt3d_d_param_hori2);
                        d_hori_err_d_param_change_order = d_hori_err_d_bearing1 + d_hori_err_d_bearing2;
                        if change_err_order
                            err_hori = err_hori_change_order;
                            d_hori_err_d_param = d_hori_err_d_param_change_order;
                            d_hori_err_d_R = d_hori_err_d_R_change_order;
                        end
                    end
                    if fix_rot
                        d_hori_err_d_R = zeros(size(d_hori_err_d_R));
                    end
                    H(1:20-offset, 1:20-offset) = H(1:20-offset, 1:20-offset) + d_hori_err_d_param' * d_hori_err_d_param;
                    H(21 + 3*(frame_id-1):23 + 3*(frame_id-1),21 + 3*(frame_id-1):23 + 3*(frame_id-1)) = H(21 + 3*(frame_id-1):23 + 3*(frame_id-1),21 + 3*(frame_id-1):23 + 3*(frame_id-1)) + d_hori_err_d_R' * d_hori_err_d_R;
                    H(1:20-offset, 21 + 3*(frame_id-1):23 + 3*(frame_id-1)) = H(1:20-offset, 21 + 3*(frame_id-1):23 + 3*(frame_id-1)) + d_hori_err_d_param' * d_hori_err_d_R;
                    H(21 + 3*(frame_id-1):23 + 3*(frame_id-1), 1:20-offset) = H(21 + 3*(frame_id-1):23 + 3*(frame_id-1), 1:20-offset) + d_hori_err_d_R' * d_hori_err_d_param;
                    b(1:20-offset) = b(1:20-offset) - d_hori_err_d_param' * err_hori;
                    b(21 + 3*(frame_id-1):23 + 3*(frame_id-1)) = b(21 + 3*(frame_id-1):23 + 3*(frame_id-1)) - d_hori_err_d_R' * err_hori;
                    err = err + sum(err_hori.^2);
                end
                if 1
                    uv1_vert = uv_list_vert(comb(j,1),:);
                    uv2_vert = uv_list_vert(comb(j,2),:);
                    [vert_bearing1, ~, d_pt3d_d_param_vert1] = UnProjectKB20(uv1_vert(1), uv1_vert(2),intr(1),intr(2),intr(3),intr(4),intr(5),intr(6),intr(7),intr(8),intr(9),intr(10),intr(11),intr(12),intr(13),intr(14),intr(15),intr(16),intr(17),intr(18),intr(19),intr(20));
                    [vert_bearing2, ~, d_pt3d_d_param_vert2] = UnProjectKB20(uv2_vert(1), uv2_vert(2),intr(1),intr(2),intr(3),intr(4),intr(5),intr(6),intr(7),intr(8),intr(9),intr(10),intr(11),intr(12),intr(13),intr(14),intr(15),intr(16),intr(17),intr(18),intr(19),intr(20));
                    d_unproject_d_bearing_vert1 = compute_d_bearing_d_pt_jac(vert_bearing1');
                    d_unproject_d_bearing_vert2 = compute_d_bearing_d_pt_jac(vert_bearing2');
                    if fix_f
                        d_pt3d_d_param_vert1(:,1:2) = 0;
                        d_pt3d_d_param_vert2(:,1:2) = 0;
                    end
                    if fix_center
                        d_pt3d_d_param_vert1(:,3:4) = 0;
                        d_pt3d_d_param_vert2(:,3:4) = 0;
                    end
                    
                    d_pt3d_d_param_vert1 = d_pt3d_d_param_vert1(:,1:end-offset);
                    d_pt3d_d_param_vert2 = d_pt3d_d_param_vert2(:,1:end-offset);
                    if 0
                        d_pt3d_d_param_vert1 = d_unproject_d_bearing_vert1 * d_pt3d_d_param_vert1;
                        d_pt3d_d_param_vert2 = d_unproject_d_bearing_vert2 * d_pt3d_d_param_vert2;
                    end
                    vert_plane_norm = cross(vert_bearing1, vert_bearing2);
                    vert_norm = norm(vert_plane_norm);
                    vert_plane_norm_change_order = cross(X2, vert_bearing2);
                    vert_norm_change_order = norm(vert_plane_norm_change_order);
                    %         vert_plane_norm = vert_plane_norm./vert_norm; 这里不能归一化
                    vert_plane_norm_check = cross(vert_bearing1, X2);
                    vert_plane_norm_check = vert_plane_norm_check./norm(vert_plane_norm_check);
                    if 0
                        assert(norm(vert_plane_norm_check - vert_plane_norm./vert_norm) < 0.00001);
                    end
                    err_vert = dot(vert_plane_norm, X2) / vert_norm;%真正的残差要归一化 % <----
                    res_vert = dot(vert_plane_norm, X2); %这是中间结果，不要归一化
                    err_vert_change_order = dot(vert_plane_norm_change_order, vert_bearing1) / vert_norm_change_order; %真正的残差要归一化% <----
                    res_vert_change_order = dot(vert_plane_norm_change_order, vert_bearing1); %这是中间结果，不要归一化
                    d_vert_plane_d_param = [-SkewSymMat(vert_bearing2) * d_pt3d_d_param_vert1 + SkewSymMat(vert_bearing1) * d_pt3d_d_param_vert2];
                    
                    if ~use_vec3_err
                        if 0
                            d_vert_err_d_vert_plane1 = [X2(1)/vert_norm - res_vert * vert_plane_norm(1) / vert_norm^3;
                                X2(2)/vert_norm - res_vert * vert_plane_norm(2) / vert_norm^3;
                                X2(3)/vert_norm - res_vert * vert_plane_norm(3) / vert_norm^3]';
                        else
                            d_vert_err_d_vert_plane = compute_d_dist_d_plane(res_vert, X2, vert_plane_norm);
                            d_vert_err_d_vert_plane_change_order = compute_d_dist_d_plane(res_vert_change_order, vert_bearing1, vert_plane_norm_change_order);
                        end
                        d_vert_err_d_param = d_vert_err_d_vert_plane * d_vert_plane_d_param;% <----
                        d_vert_err_d_X2 = vert_plane_norm./vert_norm;
                        d_X2_d_R = -SkewSymMat(X2);
                        d_vert_err_d_R = d_vert_err_d_X2 * d_X2_d_R;% <----
                        %%
                        d_vert_err_d_bearing2 = d_vert_err_d_vert_plane_change_order * SkewSymMat(X2) * d_pt3d_d_param_vert2;
                        d_vert_err_d_bearing1 = vert_plane_norm_change_order./vert_norm_change_order * d_pt3d_d_param_vert1;
                        d_vert_err_d_param_change_order = d_vert_err_d_bearing1 + d_vert_err_d_bearing2;
                        d_X2_d_R_change_order = -SkewSymMat(X2);
                        d_vert_err_d_R_change_order = d_vert_err_d_vert_plane_change_order * (-SkewSymMat(vert_bearing2)) * (d_X2_d_R_change_order);
                        if change_err_order
                            err_vert = err_vert_change_order;
                            d_vert_err_d_param = d_vert_err_d_param_change_order;
                            d_vert_err_d_R = d_vert_err_d_R_change_order;
                        end
                    else
                        normalized_plane_vert = vert_plane_norm'./vert_norm;
                        pt_in_plane_vert = X2' - (X2 * normalized_plane_vert) .*(normalized_plane_vert);
                        err_vert = (pt_in_plane_vert./norm(pt_in_plane_vert) - X2');
                        d_err_d_pt_vert = compute_d_bearing_d_pt_jac(pt_in_plane_vert);
                        d_pt_d_plane_norm_vert = compute_d_pt_d_plane_norm_jac(X2', normalized_plane_vert);
                        d_plane_norm_d_plane_vert = compute_d_bearing_d_pt_jac(vert_plane_norm');
                        d_vert_err_d_vert_plane = d_err_d_pt_vert * d_pt_d_plane_norm_vert * d_plane_norm_d_plane_vert;
                        d_vert_err_d_vert_plane1 = d_err_d_pt_vert * d_pt_d_plane_norm_vert * d_plane_norm_d_plane_vert * [-SkewSymMat(vert_bearing2) * d_pt3d_d_param_vert1];
                        d_vert_err_d_vert_plane2 = d_err_d_pt_vert * d_pt_d_plane_norm_vert * d_plane_norm_d_plane_vert * [SkewSymMat(vert_bearing1) * d_pt3d_d_param_vert2];
                        d_vert_err_d_param = d_vert_err_d_vert_plane * d_vert_plane_d_param;% <----
                        d_vert_err_d_param_comb = d_vert_err_d_vert_plane1 + d_vert_err_d_vert_plane2;
                        if 1
                            d_pt_in_plane_d_X2 = [1-normalized_plane_vert(1)^2 -normalized_plane_vert(1)*normalized_plane_vert(2) -normalized_plane_vert(1)*normalized_plane_vert(3);
                                -normalized_plane_vert(1)*normalized_plane_vert(2) 1-normalized_plane_vert(2)^2 -normalized_plane_vert(2)*normalized_plane_vert(3);
                                -normalized_plane_vert(1)*normalized_plane_vert(3) -normalized_plane_vert(2)*normalized_plane_vert(3) 1-normalized_plane_vert(3)^2];
                            d_pt_in_plane_d_X2_check = eye(3) - normalized_plane_vert * normalized_plane_vert';
                        end
                        d_vert_err_d_X2 = -eye(3) + d_err_d_pt_vert * d_pt_in_plane_d_X2;
                        d_X2_d_R = -SkewSymMat(X2);
                        d_vert_err_d_R = d_vert_err_d_X2 * d_X2_d_R; % <---
                        %                 d_vert_err_d_R = zeros(size(d_vert_err_d_R));
                        %%
                        normalized_plane_vert_change_order = vert_plane_norm_change_order'./vert_norm_change_order;
                        pt_in_plane_vert_change_order = vert_bearing1' - (vert_bearing1 * normalized_plane_vert_change_order) .* normalized_plane_vert_change_order;
                        pt_in_plane_norm_vert_change_order = pt_in_plane_vert_change_order./norm(pt_in_plane_vert_change_order);
                        err_vert_change_order = pt_in_plane_norm_vert_change_order - vert_bearing1';
                        d_err_d_pt_vert_change_order = compute_d_bearing_d_pt_jac(pt_in_plane_vert_change_order);
                        d_pt_d_plane_norm_vert_change_order = compute_d_pt_d_plane_norm_jac(vert_bearing1', normalized_plane_vert_change_order); 
                        d_plane_norm_d_plane_vert_change_order = compute_d_bearing_d_pt_jac(vert_plane_norm_change_order');
                        
                        d_X2_d_R_chnage_order = -SkewSymMat(X2);
                        d_pt_in_plane_d_vert_bearing1 = [1-normalized_plane_vert_change_order(1)^2 -normalized_plane_vert_change_order(1)*normalized_plane_vert_change_order(2) -normalized_plane_vert_change_order(1)*normalized_plane_vert_change_order(3);
                            -normalized_plane_vert_change_order(1)*normalized_plane_vert_change_order(2) 1-normalized_plane_vert_change_order(2)^2 -normalized_plane_vert_change_order(2)*normalized_plane_vert_change_order(3);
                            -normalized_plane_vert_change_order(1)*normalized_plane_vert_change_order(3) -normalized_plane_vert_change_order(2)*normalized_plane_vert_change_order(3) 1-normalized_plane_vert_change_order(3)^2];
                        d_vert_err_d_bearing1 = (-eye(3) + d_err_d_pt_vert_change_order * d_pt_in_plane_d_vert_bearing1) * d_pt3d_d_param_vert1;
                        d_vert_err_d_X2_change_order = d_err_d_pt_vert_change_order * d_pt_d_plane_norm_vert_change_order * d_plane_norm_d_plane_vert_change_order * (-SkewSymMat(vert_bearing2));
                        d_vert_err_d_R_change_order = d_vert_err_d_X2_change_order * d_X2_d_R_chnage_order;
                        d_vert_err_d_bearing2 = d_err_d_pt_vert_change_order * d_pt_d_plane_norm_vert_change_order * d_plane_norm_d_plane_vert_change_order * (SkewSymMat(X2) * d_pt3d_d_param_vert2);
                        d_vert_err_d_param_change_order = d_vert_err_d_bearing1 + d_vert_err_d_bearing2;
                        if change_err_order
                            err_vert = err_vert_change_order;
                            d_vert_err_d_param = d_vert_err_d_param_change_order;
                            d_vert_err_d_R = d_vert_err_d_R_change_order;
                        end
                    end
                    if fix_rot
                        d_vert_err_d_R = zeros(size(d_vert_err_d_R));
                    end
                    H(1:20-offset, 1:20-offset) = H(1:20-offset, 1:20-offset) + d_vert_err_d_param' * d_vert_err_d_param;
                    H(21 + 3*(frame_id-1):23 + 3*(frame_id-1),21 + 3*(frame_id-1):23 + 3*(frame_id-1)) = H(21 + 3*(frame_id-1):23 + 3*(frame_id-1),21 + 3*(frame_id-1):23 + 3*(frame_id-1)) + d_vert_err_d_R' * d_vert_err_d_R;
                    H(1:20-offset, 21 + 3*(frame_id-1):23 + 3*(frame_id-1)) = H(1:20-offset, 21 + 3*(frame_id-1):23 + 3*(frame_id-1)) + d_vert_err_d_param' * d_vert_err_d_R;
                    H(21 + 3*(frame_id-1):23 + 3*(frame_id-1), 1:20-offset) = H(21 + 3*(frame_id-1):23 + 3*(frame_id-1), 1:20-offset) + d_vert_err_d_R' * d_vert_err_d_param;
                    b(1:20-offset) = b(1:20-offset) - d_vert_err_d_param' * err_vert;
                    b(21 + 3*(frame_id-1):23 + 3*(frame_id-1)) = b(21 + 3*(frame_id-1):23 + 3*(frame_id-1)) - d_vert_err_d_R' * err_vert;
                    err = err + sum(err_vert.^2);
                end
            end
        end
    end
    param_diff = param - param0;
    if fix_rot
       last_dim = 20;
    else
        last_dim = length(param);
    end
    if ~fix_f && ~fix_center
        dx = inv(H(1:last_dim-offset,1:last_dim-offset)) * b(1:last_dim-offset);
        param(offset+1:20) = param(offset+1:20) + dx(1:20);
    end
    if fix_f && ~fix_center
        dx = inv(H(3:last_dim,3:last_dim)) * b(3:last_dim);
        dx = [0;0;dx];
        param(1:20) = param(1:20) + dx(1:20);
    end
    if ~fix_f && fix_center
        dx = inv(H([1:2 5:last_dim],[1:2 5:last_dim])) * b([1:2 5:last_dim]);
        dx = [dx(1:2);0;0;dx(3:end)];
        param(1:20) = param(1:20) + dx(1:20);
    end
    if fix_f && fix_center
        dx = inv(H([5:last_dim],[5:last_dim])) * b([5:last_dim]);
        dx = [0;0;0;0;dx];
        param(1:20) = param(1:20) + dx(1:20);
    end
    if ~fix_rot
        for frame_id = 1 : frame_num
            rotMatNew = rodrigues(dx(21 + 3*(frame_id-1):23 + 3*(frame_id-1))) * rotMat_stack{frame_id,1};
            param(21 + 3*(frame_id-1):23 + 3*(frame_id-1)) = rodrigues(rotMatNew);
        end
    end
    loss = [loss; err];
    fprintf(sprintf('iter: %d, err: %f\n', iter, err));
end

figure,subplot(1,2,1);plot(loss);title('loss');subplot(1,2,2);plot(diff(loss));title('loss diff');


end
function d_bearing_d_pt = compute_d_pt_d_plane_norm_jac(target_bearing, plane_norm)

d_bearing_d_pt = zeros(3,3);

d_bearing_d_pt(1,1) = -2*target_bearing(1)*plane_norm(1) - target_bearing(2)*plane_norm(2) - target_bearing(3)*plane_norm(3);
d_bearing_d_pt(1,2) = -target_bearing(2)*plane_norm(1);
d_bearing_d_pt(1,3) = -target_bearing(3)*plane_norm(1);

d_bearing_d_pt(2,1) = -target_bearing(1)*plane_norm(2);
d_bearing_d_pt(2,2) = -target_bearing(1)*plane_norm(1) - 2*target_bearing(2)*plane_norm(2) - target_bearing(3)*plane_norm(3);
d_bearing_d_pt(2,3) = -target_bearing(3)*plane_norm(2);

d_bearing_d_pt(3,1) = -target_bearing(1)*plane_norm(3);
d_bearing_d_pt(3,2) = -target_bearing(2)*plane_norm(3);
d_bearing_d_pt(3,3) = -target_bearing(1)*plane_norm(1) - target_bearing(2)*plane_norm(2) - 2*target_bearing(3)*plane_norm(3);



J = -(dot(target_bearing, plane_norm) *eye(3) + plane_norm * target_bearing');

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
function J_err_plane = compute_d_dist_d_plane(res, p3d_target, plane)
plane_norm = norm(plane);
J_err_plane(0+1) = p3d_target(0+1)/plane_norm - res * plane(0+1) / plane_norm^3;
J_err_plane(1+1) = p3d_target(1+1)/plane_norm - res * plane(1+1) / plane_norm^3;
J_err_plane(2+1) = p3d_target(2+1)/plane_norm - res * plane(2+1) / plane_norm^3;
end
function plane_norm_norm = Check(do_col, markerRow, MarkerId, pt3d_show, pt3dL, colorStack)
for i = 1 :1: markerRow
    
    %     markerId = (i*markerRow - markerRow+1): i*markerRow;
    if do_col
        markerId = MarkerId(:,i)';
    else
        markerId = MarkerId(i,:);
    end
    pt3d_show_list = pt3d_show(markerId,:);
    
    line_stack = [];
    for j = 1 : 1: size(pt3d_show_list,1)
        line_stack = [line_stack; [0 0 0];pt3d_show_list(j,:) ];
    end
    line_stack = [line_stack; [0 0 0]];
    figure(10), plot3(pt3d_show_list(:,1), pt3d_show_list(:,2), pt3d_show_list(:,3), '-sg','MarkerSize',5,'LineWidth',5);
    plot3(pt3dL(markerId,1)./5000, pt3dL(markerId,2)./5000, pt3dL(markerId,3)./5000, '-sb','MarkerSize',5,'LineWidth',5);
    line(line_stack(:,1),line_stack(:,2),line_stack(:,3),'Color',colorStack{i},'LineWidth', 2,'MarkerSize', 2);
    
    [~,~,plane_norm] = svd(pextend(pt3d_show_list(:,1:3)')',0);
    plane_norm = plane_norm(:,4);
    plane_norm = plane_norm./norm(plane_norm(1:3));
    plane_norm = sign(plane_norm(2)).*plane_norm;
    
    plotQuiver(plane_norm(1:3)',colorStack{i});
    %     plotQuiver(plane_norm(1:3)');
    hori_van_dir_err(i,1) = sum(abs(dot(repmat(plane_norm,1,size(pt3d_show_list,1)),pextend(pt3d_show_list(:,1:3)'))));
    plane_norms(i,:) = plane_norm(1:3)';
end


[~,~,plane_norm_norm] = svd(pextend(plane_norms(:,1:3)')',0);
plane_norm_norm = plane_norm_norm(:,4);
plane_norm_norm = plane_norm_norm./norm(plane_norm_norm(1:3));
plane_norms_dir_err = abs(dot(repmat(plane_norm_norm,1,size(plane_norms,1)),pextend(plane_norms(:,1:3)')));

figure(10),hold on; plotQuiver(plane_norm_norm(1:3)','r');




end


function [pt2d, inImageFlag] = projectKB8(pt3d, K, dis, R, t)
pt3d = ([R t; 0 0 0 1] * pt3d')';
pt3d = pt3d(:,1:3);

% inImageFlag = pt3d(:,3) > 0.005;

fx = K(1,1);
fy = K(2,2);
cx = K(1,3);
cy = K(2,3);
k1 = dis(1);
k2 = dis(2);
k3 = dis(3);
k4 = dis(4);

x = pt3d(:,1);
y = pt3d(:,2);
z = pt3d(:,3);

r2 = x .* x + y .* y;
r = sqrt(r2);

theta = atan2(r, z);
theta2 = theta .* theta;
theta4 = theta2 .* theta2;
theta6 = theta4 .* theta2;
theta8 = theta6 .* theta2;

r_theta = theta .* (1 + k1 .* theta2 + k2 .* theta4 + k3 .* theta6 + k4 .* theta8);
if r > 1e-8
    norm_inv = 1./r;
else
    norm_inv = ones(length(r), 1);
end
% norm_inv = r > 1e-8 ? double(1.0) / r : 1;

mx = r_theta .* x .* norm_inv;
my = r_theta .* y .* norm_inv;

pt2d(:,1) = fx .* mx + cx;
pt2d(:,2) = fy .* my + cy;



inImageFlag = find(pt2d(:,1) > 1 & pt2d(:,1) < 640-1 & pt2d(:,2) > 1 & pt2d(:,2) < 480-1 & pt3d(:,3)> 10);

end

function [theta,d_func_d_theta1] = solveTheta( k1,  k2,  k3,  k4,   r_theta ,d_func_d_theta, ITER)


theta = r_theta;
for  i = ITER:-1:1
    theta2 = theta .* theta;
    
    func = k4 .* theta2;
    func = func + k3;
    func = func .* theta2;
    func = func + k2;
    func = func .* theta2;
    func = func + k1;
    func = func .* theta2;
    func = func + 1.0;
    func = func .* theta;
    
    d_func_d_theta = (9) .* k4 .* theta2;
    d_func_d_theta = d_func_d_theta + (7) .* k3;
    d_func_d_theta = d_func_d_theta .* theta2;
    d_func_d_theta = d_func_d_theta + (5) .* k2;
    d_func_d_theta = d_func_d_theta .* theta2;
    d_func_d_theta = d_func_d_theta + (3) .* k1;
    d_func_d_theta = d_func_d_theta .* theta2;
    d_func_d_theta = d_func_d_theta + (1);
    
    
    %theta = theta + (r_theta - func) / d_func_d_theta;
    theta_fix = (r_theta - func) ./ d_func_d_theta;
    theta = theta+ theta_fix;
    if 0 % (abs(theta_fix) < 1e-8)
        break;
    end
end
d_func_d_theta1 = d_func_d_theta;
end
function p3d = unprojectKB8( fx,  fy,  cx,  cy,  k1,  k2,  k3,  k4, proj)


mx = (proj(:,1) - cx) ./ fx;
my = (proj(:,2) - cy) ./ fy;




theta = 0;
sin_theta = 0;
cos_theta = 1;
thetad = sqrt(mx .* mx + my .* my);

thetad = min(max(thetad,-3.141592653/2),3.141592653/2);
scaling = 1.0;
d_func_d_theta = 0;
if (thetad > 1e-8)
    [theta,d_func_d_theta1] = solveTheta(k1, k2, k3, k4, thetad, d_func_d_theta,3);
    d_func_d_theta = d_func_d_theta1;
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    scaling = sin_theta ./ thetad;
end

p3d = [mx .* scaling my .* scaling cos_theta];


end
function [BW10, BW20, BW30, matAll] = CalcRegion(distance, K, dist, pt3d, validId,T01,T02,T03, idVec)
% function [mat1, mat2, mat3] = CalcRegion(K, dist, pt3d, validId,T01,T02,T03, idVec)
% 10 12 13
validId = find(validId);
[xMat, yMat] = meshgrid(1:640, 1:480);
pix = [xMat(:) yMat(:)];
mat1 = -1.*ones(480, 640);
mat2 = mat1;
mat3 = mat1;
xyz0 = distance.*pt3d(validId,:);
xyz1 = ((T01)*pextend(xyz0'))'; xyz1 = xyz1(:,1:3);
xyz2 = ((T02)*pextend(xyz0'))'; xyz2 = xyz2(:,1:3);
xyz3 = ((T03)*pextend(xyz0'))'; xyz3 = xyz3(:,1:3);

[pt2d_1_in_0, validFlag10] = projectKB8(pextend(xyz1')', K, dist, eye(3), [0 0 0]');
BW10 = CalcValidPix(pt2d_1_in_0, validFlag10);
mat1(validId(validFlag10)) = idVec(1);
BW10(BW10 > 0) = idVec(1);
% figure(51),imshow(zeros(480,640));hold on;plot(pt2d_1_in_0(validFlag10,1), pt2d_1_in_0(validFlag10,2),'.r');title('1 in 0');

[pt2d_2_in_0, validFlag20] = projectKB8(pextend(xyz2')', K, dist, eye(3), [0 0 0]');
BW20 = CalcValidPix(pt2d_2_in_0, validFlag20);
mat2(validId(validFlag20)) = idVec(2);
BW20(BW20 > 0) = idVec(2);
% figure(52),imshow(zeros(480,640));hold on;plot(pt2d_2_in_0(validFlag20,1), pt2d_2_in_0(validFlag20,2),'.r');title('2 in 0');

[pt2d_3_in_0, validFlag30] = projectKB8(pextend(xyz3')', K, dist, eye(3), [0 0 0]');
BW30 = CalcValidPix(pt2d_3_in_0, validFlag30);
mat3(validId(validFlag30)) = idVec(3);
BW30(BW30 > 0) = idVec(3);
% figure(50),imshow(zeros(480,640));hold on;plot(pt2d_3_in_0(validFlag30,1), pt2d_3_in_0(validFlag30,2),'.r');title('3 in 0');

matAll = BW10 + BW20 + BW30;
end
function BW2 = CalcValidPix(pt2d_2_in_0, validFlag20)
pixels_hit = unique(round(pt2d_2_in_0(validFlag20,:)),'rows');
ind = sub2ind([480 640], pixels_hit(:,2),pixels_hit(:,1));
validMask = zeros(480, 640);
validMask(ind) = 1;
se = strel('square',3);
BW2 = imdilate(validMask,se);
BW2(BW2==0) = -1;
end
function [pt2d, d_uv_d_pt3d, d_uv_d_param] = project(X,Y,Z, fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2)
global Param
d_uv_d_pt3d = []; d_uv_d_param = [];

% X = pt3d(:,1);
% Y = pt3d(:,2);
% Z = pt3d(:,3);




param = [fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2];
% fx = param(1);
% fy = param(2);
% cx = param(3);
% cy = param(4);
% k1 = param(5);
% k2 = param(6);
% k3 = param(7);
% k4 = param(8);
% k5 = param(9);
% k6 = param(10);
% p1 = param(11);
% p2 = param(12);
% s1 = param(13);
% s2 = param(14);
% s3 = param(15);
% s4 = param(16);

if 1
    a = X./Z;
    b = Y./Z;
    r = sqrt(a.^2 + b.^2);
    th = atan(r);
else
    a = X;
    b = Y;
    r = sqrt(X.^2 + Y.^2);
    th = atan2(r, Z);
end


th2 = th.^2;
th4 = th2.^2;
th6 = th2.*th4;
th8 = th4.*th4;
th10 = th6.*th4;
th12 = th6.*th6;

thd = th.*(1 + k1.*th2 + k2.*th4 + k3.*th6 + k4.*th8 + k5.*th10 + k6.*th12);
x_r = a./r.*thd;
y_r = b./r.*thd;
r_d = sqrt(x_r.^2 + y_r.^2);
uvDistorted = [x_r      +      p1.*(2.*x_r.^2 + r_d.^2) + 2.*x_r.*y_r.*p2     +     s1.*r_d.^2 + s2.*r_d.^4   + s5.*r_d.^6 ...
    y_r      +      p2.*(2.*y_r.^2 + r_d.^2) + 2.*x_r.*y_r.*p1     +     s3.*r_d.^2 + s4.*r_d.^4   + s6.*r_d.^6];


[matTilt, dMatTiltdTauX, dMatTiltdTauY, invMatTilt, dInvMatTiltdTauX, dInvMatTiltdTauY] = computeTiltProjectionMatrix(t1, t2);
uvTilted =    (matTilt * [uvDistorted ones(size(uvDistorted,1),1)]')';
uvTilted = uvTilted./repmat(uvTilted(:,3),1,3);


% pt2d = [fx.*uvDistorted(:,1) + cx ...
%     fy.*uvDistorted(:,2) + cy];
pt2d = [fx.*uvTilted(:,1) + cx ...
    fy.*uvTilted(:,2) + cy];




duvDistorted_dxryr = compute_duvDistorted_dxryr([x_r y_r], x_r.^2 + y_r.^2, param);


d_thd_d_th = 1 + 3*k1.*th2 + 5*k2.*th4 + 7*k3.*th6 + 9*k4.*th8 + 11*k5.*th10 + 13*k6.*th12;

d_x_r_d_a = b.^2./r.^3.*thd + a.^2./(r.^2 + r.^4)*d_thd_d_th;
d_x_r_d_b = -a.*b./r.^3.*thd + a.*b./(r.^2 + r.^4)*d_thd_d_th;
d_y_r_d_a = -a.*b./r.^3.*thd + a.*b./(r.^2 + r.^4)*d_thd_d_th;
d_y_r_d_b = a.^2./r.^3.*thd + b.^2./(r.^2 + r.^4)*d_thd_d_th;
d_xr_yr_d_ab = [d_x_r_d_a d_x_r_d_b; d_y_r_d_a d_y_r_d_b];

d_ab_d_xyz = [ 1./Z    0     -X./Z.^2;
    0    1./Z   -Y./Z.^2];




stx = sin(t1);
ctx = cos(t1);
sty = sin(t2);
cty = cos(t2);

tt1 = ctx;
tt4 = -stx*sty;
tt5 = cty;
tt7 = sty;
tt8 = -stx*cty;
tt9 = ctx*cty;

u_d = uvDistorted(:,1);
v_d = uvDistorted(:,2);

d_ut_d_ud = tt1/(tt7*u_d + tt8*v_d + tt9) - tt1*u_d*tt7/(tt7*u_d + tt8*v_d + tt9)^2;

d_ut_d_vd = -tt1*u_d*tt8/(tt7*u_d + tt8*v_d + tt9)^2;

d_vt_d_ud = tt4/(tt7*u_d + tt8*v_d + tt9) - (tt4*u_d+tt5*v_d)*tt7/(tt7*u_d + tt8*v_d + tt9)^2;

d_vt_d_vd = tt5/(tt7*u_d + tt8*v_d + tt9) - (tt4*u_d+tt5*v_d)*tt8/(tt7*u_d + tt8*v_d + tt9)^2;

d_uvTilted_d_uvDistorted = [d_ut_d_ud d_ut_d_vd; d_vt_d_ud d_vt_d_vd];
d_uv_d_uvTilted = [fx 0;0 fy];

d_uv_d_pt3d = d_uv_d_uvTilted * d_uvTilted_d_uvDistorted * duvDistorted_dxryr * d_xr_yr_d_ab * d_ab_d_xyz;




%%%%
d_uv_d_fxfycxcy = [uvTilted(:,1)     0            1   0;
    0         uvTilted(:,2)  0   1];

d_xr_yr_d_thd = [a./r; b./r];

d_thd_d_k1k2k3k4k5k6 = [th.^3   th.^5   th.^7   th.^9   th.^11   th.^13];

d_uv_d_k1k2k3k4k5k6 = d_uv_d_uvTilted * d_uvTilted_d_uvDistorted * duvDistorted_dxryr * d_xr_yr_d_thd * d_thd_d_k1k2k3k4k5k6;

duvDistorted_d_p1p2 = [(2.*x_r.^2 + r_d.^2)        2.*x_r.*y_r;
    2.*x_r.*y_r         (2.*y_r.^2 + r_d.^2)];

duvDistorted_d_s1s2s3s4s5s6 = [r_d.^2   r_d.^4    0        0     r_d.^6     0 ;
    0        0     r_d.^2      r_d.^4     0      r_d.^6];

d_uv_d_p1p2 = d_uv_d_uvTilted * d_uvTilted_d_uvDistorted  * duvDistorted_d_p1p2;

d_uv_d_s1s2s3s4s5s6 = d_uv_d_uvTilted * d_uvTilted_d_uvDistorted  * duvDistorted_d_s1s2s3s4s5s6;


temp = sty*u_d-stx*cty*v_d+ctx*cty;

d_ut_d_tx = -stx*u_d/temp - ctx*u_d*(-ctx*cty*v_d-stx*cty)/temp^2;

d_ut_d_ty = -ctx*u_d*(cty*u_d+stx*sty*v_d-ctx*sty)/temp^2;

d_vt_d_tx = (-ctx*sty*u_d+cty*v_d)/temp - (cty*v_d-stx*sty*u_d)*(-ctx*cty*v_d-stx*cty)/temp^2;

d_vt_d_ty = (-stx*cty*u_d-sty*v_d)/temp-(cty*v_d-stx*sty*u_d)*(cty*u_d+stx*sty*v_d-ctx*sty)/temp^2;

d_uvTilted_d_t1t2 = [d_ut_d_tx d_ut_d_ty; d_vt_d_tx d_vt_d_ty];

d_uv_d_t1t2 = d_uv_d_uvTilted * d_uvTilted_d_t1t2;

d_uv_d_param = [d_uv_d_fxfycxcy d_uv_d_k1k2k3k4k5k6 d_uv_d_p1p2 d_uv_d_s1s2s3s4s5s6 d_uv_d_t1t2];



end
function [pt3d, d_pt3d_d_uv, d_pt3d_d_param] = unproject(u,v, fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2)
global Param data

d_pt3d_d_uv = []; d_pt3d_d_param = [];
% u = pt2d(:,1);
% v = pt2d(:,2);
% Z = pt3d(:,3);

param = [fx,fy,cx,cy,k1,k2,k3,k4,k5,k6,p1,p2,s1,s2,s3,s4,s5,s6,t1,t2];
% fx = param(1);
% fy = param(2);
% cx = param(3);
% cy = param(4);
% k1 = param(5);
% k2 = param(6);
% k3 = param(7);
% k4 = param(8);
% k5 = param(9);
% k6 = param(10);
% p1 = param(11);
% p2 = param(12);
% s1 = param(13);
% s2 = param(14);
% s3 = param(15);
% s4 = param(16);

uvTilted = [(u - cx)./fx (v - cy)./fy];

[matTilt, dMatTiltdTauX, dMatTiltdTauY, invMatTilt, dInvMatTiltdTauX, dInvMatTiltdTauY] = computeTiltProjectionMatrix(t1, t2);

uvDistorted =    (invMatTilt * [uvTilted ones(size(uvTilted,1),1)]')';
uvDistorted = uvDistorted./repmat(uvDistorted(:,3),1,3);

uvDistorted = uvDistorted(:,1:2);

if 1 %isempty(data.xr_yr)
    [xr_yr, duvDistorted_dxryr] = compute_xr_yr_from_uvDistorted(uvDistorted, param);
    data.xr_yr = xr_yr;
    data.duvDistorted_dxryr = duvDistorted_dxryr;
else
    xr_yr = data.xr_yr;
    duvDistorted_dxryr = data.duvDistorted_dxryr;
end
xr_yrNorm = norm(xr_yr);
if (xr_yrNorm == 0.0)
    pt3d = [0 0 1];
else
    [theta, dthD_dth] = getThetaFromNorm_xr_yr(xr_yrNorm, param);
end

if 0
    pt3d = [tan(theta) / xr_yrNorm * xr_yr 1];
else
    thd = xr_yrNorm;
    scaling = sin(theta)./thd;
    x = xr_yr(:,1)*scaling;
    y = xr_yr(:,2)*scaling;
    z = cos(theta);
    pt3d = [x y z];
    
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    mx = xr_yr(:,1);
    my = xr_yr(:,2);
    
    x_r = mx;
    y_r = my;
    r_d = sqrt(x_r.^2 + y_r.^2);
    
    d_thetad_d_mx = mx / thd;
    d_thetad_d_my = my / thd;
    
end


theta2 = theta * theta;
d_scaling_d_thetad = (thd * cos_theta / dthD_dth - sin_theta) / (thd * thd);
d_cos_d_thetad = -sin_theta / dthD_dth;
d_scaling_d_k1 = -cos_theta * theta * theta2 / (dthD_dth * thd);
d_cos_d_k1 =  -d_cos_d_thetad * theta * theta2;


d_pt3d_d_k1 = [mx * d_scaling_d_k1;
    my * d_scaling_d_k1;
    d_cos_d_k1];
d_pt3d_d_k2 = d_pt3d_d_k1 .* theta2;
d_pt3d_d_k3 = d_pt3d_d_k2 .* theta2;
d_pt3d_d_k4 = d_pt3d_d_k3 .* theta2;
d_pt3d_d_k5 = d_pt3d_d_k4 .* theta2;
d_pt3d_d_k6 = d_pt3d_d_k5 .* theta2;

d_X_d_mx = scaling + mx * d_scaling_d_thetad * d_thetad_d_mx;
d_X_d_my = mx * d_scaling_d_thetad * d_thetad_d_my;

d_Y_d_mx = my * d_scaling_d_thetad * d_thetad_d_mx;
d_Y_d_my = scaling + my * d_scaling_d_thetad * d_thetad_d_my;

d_Z_d_mx = d_cos_d_thetad * d_thetad_d_mx;
d_Z_d_my = d_cos_d_thetad * d_thetad_d_my;


d_pt3d_d_mxmy = [d_X_d_mx d_X_d_my;
    d_Y_d_mx d_Y_d_my;
    d_Z_d_mx d_Z_d_my];

d_pt3d_d_xryr = d_pt3d_d_mxmy;

d_xryr_duvDistorted = inv(duvDistorted_dxryr);
if 0
    d_tilt_d_fxfycxcy = [-uvDistorted(:,1)./fx           0             -1./fx     0;
        0          -uvDistorted(:,2)./fy    0      -1./fy];
else
    if 0
        d_tilt_d_fxfycxcy = [-(uvDistorted(:,1)-cx)./fx./fx           0             -1./fx     0;
            0          -(uvDistorted(:,2)-cy)./fy./fy    0      -1./fy];
        
    else
        d_tilt_d_fxfycxcy = [-(u-cx)./fx./fx           0             -1./fx     0;
            0          -(v-cy)./fy./fy    0      -1./fy];
    end
end
d_tilt_d_uv = [1./fx      0;
    0        1./fy];





stx = sin(t1);
ctx = cos(t1);
sty = sin(t2);
cty = cos(t2);
ttx = tan(t1);
tty = tan(t2);

tt1 = 1/ctx;
tt4 = ttx*tty;
tt5 = 1/cty;
tt7 = -tty;
tt8 = ttx/cty;
tt9 = 1/ctx/cty;

u_t = uvTilted(:,1);
v_t = uvTilted(:,2);

d_ud_d_ut = tt1/(tt7*u_t + tt8*v_t + tt9) - tt1*u_t*tt7/(tt7*u_t + tt8*v_t + tt9)^2;

d_ud_d_vt = -tt1*u_t*tt8/(tt7*u_t + tt8*v_t + tt9)^2;

d_vd_d_ut = tt4/(tt7*u_t + tt8*v_t + tt9) - (tt4*u_t+tt5*v_t)*tt7/(tt7*u_t + tt8*v_t + tt9)^2;

d_vd_d_vt = tt5/(tt7*u_t + tt8*v_t + tt9) - (tt4*u_t+tt5*v_t)*tt8/(tt7*u_t + tt8*v_t + tt9)^2;

d_uvDistorted_d_uvTilted = [d_ud_d_ut d_ud_d_vt; d_vd_d_ut d_vd_d_vt];





d_udvd_d_tilt = d_uvDistorted_d_uvTilted;




temp = -tty*u_t + ttx/cty*v_t+1/ctx/cty;


d_ud_d_tx = (stx/ctx^2*u_t/temp)-(1/ctx*u_t*(v_t/ctx^2/cty+cty*stx/(ctx*cty)^2)/temp^2);

d_ud_d_ty = -u_t/ctx*(-u_t/cty^2+ttx*sty*v_t/cty^2+ctx*sty/(ctx*cty)^2)/temp^2;

d_vd_d_tx = (1/ctx^2*tty*u_t/temp) - ((ttx*tty*u_t+1/cty*v_t)*(1/cty/ctx^2*v_t+cty*stx/(ctx*cty)^2)/temp^2);

d_vd_d_ty = ((ttx*u_t/cty^2 + sty*v_t/cty^2)/temp) - ((ttx*tty*u_t+1/cty*v_t)*(-u_t/cty^2+ttx*sty*v_t/cty^2+ctx*sty/(ctx*cty)^2)/temp^2);

d_udvd_d_t1t2 = [d_ud_d_tx d_ud_d_ty; d_vd_d_tx d_vd_d_ty];
%%
duvDistorted_d_p1p2 = -[(2.*x_r.^2 + r_d.^2)        2.*x_r.*y_r;
    2.*x_r.*y_r         (2.*y_r.^2 + r_d.^2)];
%%
duvDistorted_d_s1s2s3s4s5s6 = -[r_d.^2   r_d.^4    0        0     r_d.^6    0;
    0        0     r_d.^2   r_d.^4      0     r_d.^6];


d_xryr_d_uv = d_xryr_duvDistorted * d_udvd_d_tilt * d_tilt_d_uv;
d_xryr_d_fxfycxcy = d_xryr_duvDistorted * d_udvd_d_tilt* d_tilt_d_fxfycxcy;
% d_xryr_d_t1t2 = d_xryr_duvDistorted * d_udvd_d_tilt* d_tilt_d_t1t2;
d_xryr_d_t1t2 = d_xryr_duvDistorted * d_udvd_d_t1t2;
d_xryr_d_p1p2 = d_xryr_duvDistorted * duvDistorted_d_p1p2;
d_xryr_d_s1s2s3s4s5s6 = d_xryr_duvDistorted * duvDistorted_d_s1s2s3s4s5s6;


d_pt3d_d_param = [d_pt3d_d_mxmy * d_xryr_d_fxfycxcy ...
    [d_pt3d_d_k1 d_pt3d_d_k2 d_pt3d_d_k3 d_pt3d_d_k4 d_pt3d_d_k5 d_pt3d_d_k6] ...
    d_pt3d_d_mxmy * [d_xryr_d_p1p2 d_xryr_d_s1s2s3s4s5s6 d_xryr_d_t1t2]];

d_pt3d_d_uv = d_pt3d_d_mxmy * d_xryr_d_uv;

end
function [xr_yr, duvDistorted_dxryr] = compute_xr_yr_from_uvDistorted(uvDistorted, param)
global test_Jac
fx = param(1);
fy = param(2);
cx = param(3);
cy = param(4);
k1 = param(5);
k2 = param(6);
k3 = param(7);
k4 = param(8);
k5 = param(9);
k6 = param(10);
p1 = param(11);
p2 = param(12);
s1 = param(13);
s2 = param(14);
s3 = param(15);
s4 = param(16);
s5 = param(17);
s6 = param(18);
t1 = param(19);
t2 = param(20);
%initial guess:
xr_yr = uvDistorted;
max_iter = 20;
for i = 1 : max_iter
    uvDistorted_est = xr_yr;
    xr_yr_squaredNorm = xr_yr(:,1).^2 + xr_yr(:,2).^2;
    x_r = xr_yr(:,1);
    y_r = xr_yr(:,2);
    r_d = sqrt(x_r.^2 + y_r.^2);
    uvDistorted_est = uvDistorted_est + [ p1.*(2.*x_r.^2 + r_d.^2) + 2.*x_r.*y_r.*p2     +     s1.*r_d.^2 + s2.*r_d.^4 + s5.*r_d.^6 ...
        p2.*(2.*y_r.^2 + r_d.^2) + 2.*x_r.*y_r.*p1     +     s3.*r_d.^2 + s4.*r_d.^4 + s6.*r_d.^6];
    duvDistorted_dxryr = compute_duvDistorted_dxryr(xr_yr, xr_yr_squaredNorm, param);
    correction = (inv(duvDistorted_dxryr) * (uvDistorted' - uvDistorted_est'))';
    xr_yr  = xr_yr + correction;
    if(~test_Jac)
        if(norm(correction) < 1e-20)
            break;
        end
    end
end



end
function duvDistorted_dxryr = compute_duvDistorted_dxryr(xr_yr, xr_yr_squaredNorm, param)
fx = param(1);
fy = param(2);
cx = param(3);
cy = param(4);
k1 = param(5);
k2 = param(6);
k3 = param(7);
k4 = param(8);
k5 = param(9);
k6 = param(10);
p1 = param(11);
p2 = param(12);
s1 = param(13);
s2 = param(14);
s3 = param(15);
s4 = param(16);
s5 = param(17);
s6 = param(18);
t1 = param(19);
t2 = param(20);
duvDistorted_dxryr(1,1,:) = 1 + 6 .* xr_yr(:,1).*p1 + 2.*xr_yr(:,2).*p2   +   2.*(s1 + 2.*s2.*xr_yr_squaredNorm + 3.*s5.*xr_yr_squaredNorm.^2 ).*xr_yr(:,1);
duvDistorted_dxryr(1,2,:) = 2.*xr_yr(:,2).*p1 + 2.*xr_yr(:,1).*p2         +   2.*(s1 + 2.*s2.*xr_yr_squaredNorm + 3.*s5.*xr_yr_squaredNorm.^2 ).*xr_yr(:,2);
duvDistorted_dxryr(2,1,:) = 2.*xr_yr(:,2).*p1 + 2.*xr_yr(:,1).*p2         +   2.*(s3 + 2.*s4.*xr_yr_squaredNorm + 3.*s6.*xr_yr_squaredNorm.^2 ).*xr_yr(:,1);
duvDistorted_dxryr(2,2,:) = 1 + 6 .* xr_yr(:,2).*p2 + 2.*xr_yr(:,1).*p1   +   2.*(s3 + 2.*s4.*xr_yr_squaredNorm + 3.*s6.*xr_yr_squaredNorm.^2 ).*xr_yr(:,2);



end
function [th, dthD_dth] = getThetaFromNorm_xr_yr(th_radialDesired, param)
global test_Jac Param


fx = Param(1);
fy = Param(2);
cx = Param(3);
cy = Param(4);
k1 = Param(5);
k2 = Param(6);
k3 = Param(7);
k4 = Param(8);
k5 = Param(9);
k6 = Param(10);
p1 = Param(11);
p2 = Param(12);
s1 = Param(13);
s2 = Param(14);
s3 = Param(15);
s4 = Param(16);
t1 = param(17);
t2 = param(18);

param = Param;


%  initial guess
th = th_radialDesired;



startK = 5;
max_iter = 20;
for i = 1 : max_iter
    
    if 1
        thetaSq = th * th;
        th_radial = 1;
        dthD_dth = 1;
        theta2is = thetaSq;
        for j = 0:5
            th_radial = th_radial + theta2is * param(startK + j);
            dthD_dth = dthD_dth + (2 * j + 3) * param(startK + j) * theta2is;
            theta2is = theta2is * thetaSq;
        end
    else
        th2 = th.^2;
        th4 = th2.^2;
        th6 = th2.*th4;
        th8 = th4.*th4;
        th10 = th6.*th4;
        th12 = th6.*th6;
        dthD_dth = 1 + 3*k1.*th2 + 5*k2.*th4 + 7*k3.*th6 + 9*k4.*th8 + 11*k5.*th10 + 13*k6.*th12;
    end
    
    th_radial = th_radial * th;
    %     if(~test_Jac)
    if (abs(dthD_dth) > 1e-20)
        step = (th_radialDesired - th_radial) / dthD_dth;
    else
        
        if (th_radialDesired - th_radial) * dthD_dth > 0.0
            step =  1e-19;
        else
            step = -1e-19;
        end
    end
    
    th = th + step;
    if(norm(step) < 1e-20)
        break;
    end
    
    if (abs(th) >=3.1415926 / 2.0)
        
        th = (0.999) * 3.1415926 / 2.0;
    end
    %     else
    %          step = (th_radialDesired - th_radial) / dthD_dth;
    %           th = th + step;
    %     end
end

end
function [matTilt, dMatTiltdTauX, dMatTiltdTauY, invMatTilt, dInvMatTiltdTauX, dInvMatTiltdTauY] =  computeTiltProjectionMatrix(tauX, tauY)

cTauX = cos(tauX);
sTauX = sin(tauX);
cTauY = cos(tauY);
sTauY = sin(tauY);
matRotX = [1,0,0;0,cTauX,sTauX;0,-sTauX,cTauX];
matRotY = [cTauY,0,-sTauY;0,1,0;sTauY,0,cTauY];
matRotXY = matRotY * matRotX;
matProjZ = [matRotXY(3,3),0,-matRotXY(1,3);0,matRotXY(3,3),-matRotXY(2,3);0,0,1];

% Matrix for trapezoidal distortion of tilted image sensor
matTilt = matProjZ * matRotXY;

% Derivative with respect to tauX
dMatRotXYdTauX = matRotY * [0,0,0;0,-sTauX,cTauX;0,-cTauX,-sTauX];
dMatProjZdTauX = [dMatRotXYdTauX(3,3),0,-dMatRotXYdTauX(1,3);0,dMatRotXYdTauX(3,3),-dMatRotXYdTauX(2,3);0,0,0];
dMatTiltdTauX = (matProjZ * dMatRotXYdTauX) + (dMatProjZdTauX * matRotXY);
dMatTiltdTauX_check = [-sin(tauX) 0 0; -cos(tauX)*sin(tauY) 0  0; 0 -cos(tauX)*cos(tauY) -sin(tauX)*cos(tauY)];
err1 = dMatTiltdTauX_check - dMatTiltdTauX;

% Derivative with respect to tauY
dMatRotXYdTauY = [-sTauY,0,-cTauY;0,0,0;cTauY,0,-sTauY] * matRotX;
dMatProjZdTauY = [dMatRotXYdTauY(3,3),0,-dMatRotXYdTauY(1,3);0,dMatRotXYdTauY(3,3),-dMatRotXYdTauY(2,3);0,0,0];
dMatTiltdTauY = (matProjZ * dMatRotXYdTauY) + (dMatProjZdTauY * matRotXY);
dMatTiltdTauY_check = [0 0 0; -sin(tauX)*cos(tauY) -sin(tauY) 0; cos(tauY) sin(tauX)*sin(tauY) -cos(tauX)*sin(tauY)];
err2 = dMatTiltdTauY_check - dMatTiltdTauY;

invZ = 1./matRotXY(3,3);
invMatProjZ = [invZ,0,invZ*matRotXY(1,3);0,invZ,invZ*matRotXY(2,3);0,0,1];
invMatTilt = matRotXY'*invMatProjZ;
invMatTilt_check = [1/cos(tauX) 0 0; tan(tauX)*tan(tauY) 1/cos(tauY) 0; -tan(tauY) tan(tauX)/cos(tauY) 1/(cos(tauX)*cos(tauY))];
err3 = invMatTilt_check - invMatTilt;

dInvMatTiltdTauX = [sin(tauX)/(cos(tauX) * cos(tauX))               0                                     0;
    tan(tauY)/(cos(tauX) * cos(tauX))                0                                     0;
    0                   1/cos(tauX)/cos(tauX)/cos(tauY)  sin(tauX)/cos(tauX)/cos(tauX)/cos(tauY)];

dInvMatTiltdTauY = [            0                                   0                                        0;
    tan(tauX)/cos(tauY)/cos(tauY)       sin(tauY)/cos(tauY)/cos(tauY)                          0;
    -1/cos(tauY)/cos(tauY)     tan(tauX)*sin(tauY)/cos(tauY)/cos(tauY)    sin(tauY)/cos(tauX)/cos(tauY)/cos(tauY)];
end