function testMargPoseToIdp()


pose_dim = 6;
point_dim = 1;
point_num = 10;
all_dim = pose_dim + point_num;

H = zeros(all_dim, all_dim);
b = zeros(all_dim, 1);


H11_each = cell(point_num,1);
b1_each = cell(point_num,1);
for pid = 1 : point_num
    pose_jac = rand(3, 6);
    point_jac = rand(3, 1);
    err = rand(3, 1);
    
    H = FillMatrix(H, pose_dim, pose_dim, 1, 1, pose_jac' * pose_jac);
    b = FillMatrix(b, pose_dim, 1, 1, 1, -pose_jac' * err);
    H11_each{pid,1} = pose_jac' * pose_jac;
    b1_each{pid,1} = -pose_jac' * err;
    H = FillMatrix(H, pose_dim, point_dim, point_dim, pose_dim + pid, pose_jac' * point_jac);
    H = FillMatrix(H, point_dim, pose_dim, pose_dim + pid, point_dim, point_jac' * pose_jac);
    
    H = FillMatrix(H, point_dim, point_dim, pose_dim + pid, pose_dim + pid, point_jac' * point_jac);
    b = FillMatrix(b, point_dim, 1, pose_dim + pid, 1, -point_jac' * err);
end
dx = inv(H) * b;
H11 = H(1:pose_dim, 1:pose_dim);
H12 = H(1:pose_dim, pose_dim+1:all_dim);
H22 = diag(H(pose_dim + 1:all_dim, pose_dim + 1:all_dim));
b1 = b(1:pose_dim);
b2 = b(pose_dim + 1 : all_dim);

H11_orig = H11;
b1_orig = b1;

H11 = H11 - H12 * diag(1./H22) * H12';
b1 = b1 - H12 * diag(1./H22) * b2;
dx1 = inv(H11) * b1;
dx2 = diag(1./H22) * (b2 - H12' * dx1);
[dx - [dx1;dx2]]'


cov_full = inv(H);
cov_pose = inv(H11);
cov_point = diag(1./H22) + diag(1./H22) * H12' * cov_pose * H12 * diag(1./H22);
cov_point_diff = cov_point - cov_full(pose_dim+1:end, pose_dim+1:end)



H22_temp = zeros(point_num, 1);
b2_temp = zeros(point_num, 1);
dx2_temp = zeros(point_num, 1);
dx1_temp = zeros(pose_dim, 1);
H12_temp = zeros(pose_dim, point_num);

for pid = 1 : point_num   
    H22_temp(pid) = H22(pid);
    H12_temp(:, pid) = H12(:,pid);
    b2_temp(pid) = b2(pid);
    
%     H22_temp(pid) = H22_temp(pid) - H12_temp(:,pid)' * inv(H11_each{pid,1}) * H12_temp(:,pid);
%     b2_temp(pid) = b2_temp(pid) - H12_temp(:,pid)' * inv(H11_each{pid,1}) * b1_each{pid,1};
    H22_temp(pid) = H22_temp(pid) - H12_temp(:,pid)' * inv(H11_orig) * H12_temp(:,pid);
    b2_temp(pid) = b2_temp(pid) - H12_temp(:,pid)' * inv(H11_orig) * b1_orig;
    dx2_temp(pid) = 1./H22_temp(pid) * b2_temp(pid);
    dx1_temp_temp = inv(H11_orig) * (b1_orig - H12_temp(:,pid) * dx2_temp(pid)); 
    dx1_temp = dx1_temp + dx1_temp_temp;
    
end


dim_marg = pose_dim;
HH = zeros(size(H));
HH(1 : (all_dim - dim_marg),1 : (all_dim - dim_marg)) = H(dim_marg + 1: end, dim_marg + 1: end);
HH((all_dim - dim_marg)+1:end,(all_dim - dim_marg)+1:end) = H(1: dim_marg, 1: dim_marg);
HH((all_dim - dim_marg)+1:all_dim,1 : (all_dim - dim_marg)) = H(1: dim_marg, dim_marg + 1: end);
HH(1 : (all_dim - dim_marg), (all_dim - dim_marg)+1:end) = H(dim_marg + 1: end, 1: dim_marg);

HH21 = HH((all_dim - dim_marg)+1:all_dim,1 : (all_dim - dim_marg));
HH11 = HH(1 : (all_dim - dim_marg),1 : (all_dim - dim_marg));
HH22 = HH((all_dim - dim_marg)+1:end,(all_dim - dim_marg)+1:end);
bb1 = b(dim_marg + 1: end);

bb2 = b(1: dim_marg);

HH11 = HH11 - HH21' * inv(HH22) * HH21;
bb11 = bb1 - HH21' *inv(HH22) * bb2;
dx11 = inv(HH11) * bb11;
dx22 = inv(HH22) * (bb2 - HH21 * dx11);

%把pose marg到idp会导致idp之间存在关联，hessian变稠密，HH11变稠密了
dx_new = [dx22; dx11];
[dx_new - dx]'

end
function Mat2 = FillMatrix(Mat1, rows, cols, start_row, start_col, data)


Mat1 = GrowMat(Mat1, rows, cols, start_row, start_col);
Mat1(start_row:start_row+rows-1, start_col:start_col+cols-1) = Mat1(start_row:start_row+rows-1, start_col:start_col+cols-1) + data;

Mat2 = Mat1;
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