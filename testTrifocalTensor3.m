function testTrifocalTensor3()

close all;

data = load('G:\matlab\data\direct\gt\D2_011\4\output.txt');

data = data(:,2:8);

poseMat = [];
for i = 1 : size(data,1)
    data1 = data(i,:);
    xyzw = data1(4:7);
    trans = data1(1:3);
    R = quat2rotm(xyzw([4 1 2 3]));
    poseMat = [poseMat; [reshape(R,1,9), trans]];
end


figure,plotPath(poseMat(1:10:100,:));


intrMat = [500 0 320; 0 500 240; 0 0 1];

host_uv = [100 10];
host_rho = 0.5;

pose_num = 20;

T_wc_host = [reshape(poseMat(1, 1:9), 3, 3) poseMat(1, 10:12)';0 0 0 1];

use_bearing = true;false;

point_trace = [host_uv 1];
point_trace_gt = [host_uv 1];

rho_diff = zeros(pose_num, 1);
rho_noise = 0.1 * (rand(pose_num, 1) - 0.5);
rho_diff(3:end) = rho_noise(3:end);
T_cw_stack{1, 1} = inv(T_wc_host);
xyz_err = 1./(host_rho + rho_diff) - 1/host_rho;
for id = 2 : pose_num
    T_wc_cur = [reshape(poseMat(id, 1:9), 3, 3) poseMat(id, 10:12)';0 0 0 1];
    T_cw_stack{id, 1} = inv(T_wc_cur);
    T_th = inv(T_wc_cur) * T_wc_host;
    rho = host_rho + rho_diff(id);
    rho_gt = host_rho;
    xyz_in_host = inv(intrMat) * [host_uv 1]'./rho;
    xyz_in_host_gt = inv(intrMat) * [host_uv 1]'./rho_gt;
    
    xyz_in_cur = T_th(1:3,1:3) * xyz_in_host + T_th(1:3,4);
    xyz_in_cur_gt = T_th(1:3,1:3) * xyz_in_host_gt + T_th(1:3,4);
    
    if xyz_in_cur_gt(3) < 0.5
        continue;
    end
    
    cur_cur = pflat(intrMat * xyz_in_cur);
    cur_cur_gt = pflat(intrMat * xyz_in_cur_gt);
    
    point_trace = [point_trace; cur_cur(1:2)' 1];
    point_trace_gt = [point_trace_gt; cur_cur_gt(1:2)' 1];
    if id >= 3
        dist_to_line_close = EpipolarTransfer(T_cw_stack{id-1, 1}, T_cw_stack{id, 1}, intrMat, intrMat(1, 1),point_trace(id-1,:), point_trace(id,:), use_bearing);
        dist_to_line_far = EpipolarTransfer(T_cw_stack{1, 1}, T_cw_stack{id, 1}, intrMat, intrMat(1, 1),point_trace(1,:), point_trace(1,:), use_bearing);
        dist_to_line_close_gt = EpipolarTransfer(T_cw_stack{id-1, 1}, T_cw_stack{id, 1}, intrMat, intrMat(1, 1),point_trace_gt(id-1,:), point_trace_gt(id,:), use_bearing);
        
        dist_to_line_close_with_gt_pre = EpipolarTransfer(T_cw_stack{id-1, 1}, T_cw_stack{id, 1}, intrMat, intrMat(1, 1),point_trace_gt(id-1,:), point_trace(id,:), use_bearing);
        
        
        bearing33_close = trifocalTransfer(T_cw_stack{id-2, 1}, T_cw_stack{id-1, 1}, T_cw_stack{id, 1}, intrMat, point_trace(id-2,:), point_trace(id-1,:), point_trace(id,:), use_bearing);
        point_trans_err_close = norm(bearing33_close -  cur_cur');
        if 0
            bearing33_close_with_gt_pre = trifocalTransfer(T_cw_stack{id-2, 1}, T_cw_stack{id-1, 1}, T_cw_stack{id, 1}, intrMat, point_trace_gt(id-2,:), point_trace_gt(id-1,:), point_trace(id,:), use_bearing);
        else
            bearing33_close_with_gt_pre = trifocalTransfer(T_cw_stack{id-2, 1}, T_cw_stack{id-1, 1}, T_cw_stack{id, 1}, intrMat, point_trace_gt(id-2,:), point_trace(id-1,:), point_trace(id,:), use_bearing);
        end
        point_trans_err_close_with_gt_pre = norm(bearing33_close_with_gt_pre -  cur_cur');
        bearing33_far = trifocalTransfer(T_cw_stack{1, 1}, T_cw_stack{2, 1}, T_cw_stack{id, 1}, intrMat, point_trace(1,:), point_trace(2,:), point_trace(id,:), use_bearing);
        point_trans_err_far = norm(bearing33_far -  cur_cur');
        
        bearing33_gt = trifocalTransfer(T_cw_stack{id-2, 1}, T_cw_stack{id-1, 1}, T_cw_stack{id, 1}, intrMat, point_trace_gt(id-2,:), point_trace_gt(id-1,:), point_trace_gt(id,:), use_bearing);
        point_trans_err_gt = norm(bearing33_gt -  cur_cur_gt');
        reproj_err = norm(cur_cur_gt - cur_cur);
        err_close(id,:) = [dist_to_line_close point_trans_err_close reproj_err];
        err_far(id,:) = [dist_to_line_far point_trans_err_far reproj_err];
        err_close_with_gt_pre(id,:) = [dist_to_line_close_with_gt_pre point_trans_err_close_with_gt_pre reproj_err];
    end
end


info = [xyz_err err_close_with_gt_pre];

figure,imshow(zeros(480, 640)); hold on;plot(point_trace(:,1), point_trace(:,2),'.g');plot(point_trace(1,1), point_trace(1,2),'or')


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
function bearing33 = trifocalTransfer(Tc1w, Tc2w, Tc3w, K, bearing1, bearing2, bearing3, use_bearing)
Twc1 = inv(Tc1w);
Twc2 = inv(Tc2w);
Twc3 = inv(Tc3w);

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
    end
end

if 0
    z = z(1:2,:);
    z = z(:);
else
    diff = z' - bearing3;
end
bearing33 = z';
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