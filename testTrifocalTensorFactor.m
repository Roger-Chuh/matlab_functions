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

pose_num = 7;

T_wc_host = [reshape(poseMat(1, 1:9), 3, 3) poseMat(1, 10:12)';0 0 0 1];

use_bearing = true;false;
is_5dof = false;true;
point_trace = pextend(host_uv')';
% point_trace_gt = [host_uv 1];

rho_diff = zeros(pose_num, 1);
rho_noise = 0.0 * (rand(pose_num, 1) - 0.5);
rho_diff(3:end) = rho_noise(3:end);
T_cw_stack{1, 1} = inv(T_wc_host);
% xyz_err = 1./(host_rho + rho_diff) - 1/host_rho;
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
Twc_stack_gt = Twc_stack;

for id = 1 : pose_num
    T_wc_cur = Twc_stack{id, 1};
    T_cw_stack{id, 1} = inv(T_wc_cur);
    if id >= 3
        target_braring_predict = zeros(size(target_braring, 1), 3);
        for pid = 1 : size(target_braring, 1)
            [err, target_braring_predict(pid,:)] = trifocalTransfer(T_cw_stack{id-2, 1}, T_cw_stack{id-1, 1}, T_cw_stack{id, 1}, eye(3), uv_stack{id-2}(pid,:), uv_stack{id-1}(pid,:), uv_stack{id}(pid,:), use_bearing, is_5dof);
        end
    end
end


% info = [xyz_err err_close_with_gt_pre];
% 
% figure,imshow(zeros(480, 640)); hold on;plot(point_trace(:,1), point_trace(:,2),'.g');plot(point_trace(1,1), point_trace(1,2),'or')


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
function [err, bearing33, d_err_d_Rwc1, d_err_d_Rwc2, d_err_d_Rwc3, d_err_d_twc1, d_err_d_twc2, d_err_d_twc3] = trifocalTransfer(Tc1w, Tc2w, Tc3w, K, bearing1, bearing2, bearing3, use_bearing, is_5dof)

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



if is_5dof
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