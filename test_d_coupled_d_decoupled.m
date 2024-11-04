function test_d_coupled_d_decoupled()

pose_num = 10;
Tbc = [rodrigues(0.1 * (rand(3, 1) - 0.5)) 1 * (rand(3, 1) - 0.5); 0 0 0 1];

for i = 1 : pose_num
    Twbs{i, 1} = [rodrigues(0.1 * (rand(3, 1) - 0.5)) 1 * (rand(3, 1) - 0.5); 0 0 0 1];
    Twcs{i, 1} = Twbs{i, 1} * Tbc;
end


% FROM BASALT-VIO
% TEST(SophusUtilsCase, RelPoseTestDecoupledLeftIncSE3) {
%   Sophus::SE3d T_w_i = Sophus::SE3d::exp(Sophus::Vector6d::Random());
%   Sophus::SE3d T_w_j = Sophus::SE3d::exp(Sophus::Vector6d::Random());
% 
%   for (const double meas_error : {1e0, 1e-1, 1e-2, 1e-4}) {
%     Sophus::SE3d T_ij_meas =
%         T_w_i.inverse() * T_w_j *
%         Sophus::SE3d::exp(Sophus::Vector6d::Random() * meas_error);
% 
%     Sophus::SE3d T_j_i = T_w_j.inverse() * T_w_i;
%     Sophus::Vector6d res = Sophus::se3_logd(T_ij_meas * T_j_i);
% 
%     Sophus::Matrix6d J_T_w_i, d_left_perturb_d_so3_trans_i;
%     Sophus::Matrix6d J_T_w_j, d_left_perturb_d_so3_trans_j, J_T_w_j_comp,
%         d_left_perturb_d_so3_trans_j2;
%     Sophus::Matrix6d rr_i, rr_i_comp;
%     Sophus::Matrix6d rr_j, rr_j_comp, rr_j_comp2;
% 
%     Sophus::Matrix6d d_err_d_Twi_left, d_err_d_Twj_left,
%         d_err_d_Twi_left_decoupled, d_err_d_Twj_left_decoupled, Jrinv;
% 
%     /// left perturb SE3
%     Sophus::rightJacobianInvSE3Decoupled(res, J_T_w_i);
%     /// left perturb SE3
%     J_T_w_j = -J_T_w_i * T_j_i.inverse().Adj();
%     J_T_w_j_comp = -J_T_w_i * T_w_i.inverse().Adj();
% 
%     rr_i.setZero();
%     rr_i.topLeftCorner<3, 3>() = rr_i.bottomRightCorner<3, 3>() =
%         T_w_i.so3().inverse().matrix();
% 
%     d_left_perturb_d_so3_trans_i.setIdentity();
%     d_left_perturb_d_so3_trans_i.topRightCorner<3, 3>() =
%         Sophus::SO3d::hat(T_w_i.translation());
%     rr_i_comp = T_w_i.inverse().Adj() * d_left_perturb_d_so3_trans_i;
%     std::cout << "rr_i:\n" << rr_i << std::endl;
%     std::cout << "rr_i_comp:\n" << rr_i_comp << std::endl;
%     std::cout << "rr_i_diff:\n" << rr_i - rr_i_comp << std::endl;
% 
%     rr_j.setZero();
%     rr_j.topLeftCorner<3, 3>() = rr_j.bottomRightCorner<3, 3>() =
%         T_w_j.so3().inverse().matrix();
% 
%     d_left_perturb_d_so3_trans_j.setIdentity();
%     if (false) {
%       d_left_perturb_d_so3_trans_j.topRightCorner<3, 3>() =
%           Sophus::SO3d::hat(T_j_i.inverse().translation());
%     } else {
%       d_left_perturb_d_so3_trans_j.topRightCorner<3, 3>() =
%           Sophus::SO3d::hat(T_w_j.translation());
%     }
%     rr_j_comp = T_w_j.inverse().Adj() * d_left_perturb_d_so3_trans_j;
%     std::cout << "rr_j:\n" << rr_j << std::endl;
%     std::cout << "rr_j_comp:\n" << rr_j_comp << std::endl;
%     std::cout << "rr_j_diff:\n" << rr_j - rr_j_comp << std::endl;
% 
%     Sophus::Vector6d x0;
%     x0.setZero();
%     printf("aaaaa\n");
% 
%     {
%       Sophus::rightJacobianInvSE3Decoupled(res, Jrinv);
%       d_err_d_Twi_left = Jrinv * T_w_i.inverse().Adj();
%       d_err_d_Twj_left = -Jrinv * T_w_i.inverse().Adj();
%       Sophus::Matrix6d d_left_se3_d_left_so3_trans_i,
%           d_left_se3_d_left_so3_trans_j;
%       d_left_se3_d_left_so3_trans_i.setIdentity();
%       d_left_se3_d_left_so3_trans_j.setIdentity();
%       d_left_se3_d_left_so3_trans_i.topRightCorner<3, 3>() =
%           Sophus::SO3d::hat(T_w_i.translation());
%       d_left_se3_d_left_so3_trans_j.topRightCorner<3, 3>() =
%           Sophus::SO3d::hat(T_w_j.translation());
%       d_err_d_Twi_left_decoupled =
%           d_err_d_Twi_left * d_left_se3_d_left_so3_trans_i;
%       d_err_d_Twj_left_decoupled =
%           d_err_d_Twj_left * d_left_se3_d_left_so3_trans_j;
% 
%       Sophus::Matrix6d J_w_i_basalt = J_T_w_i * rr_i;
%       Sophus::Matrix6d J_w_j_basalt = J_T_w_j * rr_j;
% 
%       std::cout << "J_w_i_basalt:\n"
%                 << J_w_i_basalt << "\nJ_w_i_roger:\n"
%                 << d_err_d_Twi_left_decoupled << "\ndiff_i:\n"
%                 << J_w_i_basalt - d_err_d_Twi_left_decoupled << std::endl;
%       std::cout << "J_w_j_basalt:\n"
%                 << J_w_j_basalt << "\nJ_w_j_roger:\n"
%                 << d_err_d_Twj_left_decoupled << "\ndiff_j:\n"
%                 << J_w_j_basalt - d_err_d_Twj_left_decoupled << std::endl;
%     }
% 
%     test_jacobian(
%         "d_res_d_T_w_i", J_T_w_i * rr_i,
%         [&](const Sophus::Vector6d &x) {
%           Sophus::SE3d T_w_i_new;
%           T_w_i_new.so3() = Sophus::SO3d::exp(x.tail<3>()) * T_w_i.so3();
%           T_w_i_new.translation() = T_w_i.translation() + x.head<3>();
% 
%           return Sophus::se3_logd(T_ij_meas * T_w_j.inverse() * T_w_i_new);
%         },
%         x0);
% 
%     test_jacobian(
%         "d_res_d_T_w_j", J_T_w_j * rr_j,
%         [&](const Sophus::Vector6d &x) {
%           Sophus::SE3d T_w_j_new;
%           T_w_j_new.so3() = Sophus::SO3d::exp(x.tail<3>()) * T_w_j.so3();
%           T_w_j_new.translation() = T_w_j.translation() + x.head<3>();
% 
%           return Sophus::se3_logd(T_ij_meas * T_w_j_new.inverse() * T_w_i);
%         },
%         x0);
%   }
% }



end

function jac = Adj(T)
R = T(1:3,1:3);
jac = zeros(6,6);
jac(1:3,1:3) = R;
jac(4:6, 4:6) = R;
jac(4:6,1:3) = SkewSymMat(T(1:3,4)) * R;
jac(1:3, 4:6) = zeros(3, 3);

end
function J = rightJacobianInvSE3Decoupled(phi)
J = zeros(6, 6);
J(1:3,1:3) = JrInv(phi(1:3));
J(4:6,4:6) = rodrigues(phi(1:3));
end
function J = JrInv(phi)
EPSILON = 1e-6;
EPSILONSQRT = sqrt(EPSILON);

J = eye(3);

phi_norm2 = sum(phi.^2);
phi_hat = SkewSymMat(phi);
phi_hat2 = phi_hat * phi_hat;

J = J + phi_hat / 2;
if (phi_norm2 > EPSILON)
    phi_norm = sqrt(phi_norm2);
    
    
    
    if (phi_norm < 3.141592653 - EPSILONSQRT)
        
        J = J + phi_hat2 * (1 / phi_norm2 - (1 + cos(phi_norm)) / (2 * phi_norm * sin(phi_norm)));
    else
        
        J = J + phi_hat2 / (3.141592653 * 3.141592653);
    end
    
else
    
    J = J + phi_hat2 / 12;
    
end

end
function drdp = LogSE3(T)
drdp = zeros(6,1);
R = T(1:3,1:3);
omega = rodrigues(R);
theta = norm(omega);
drdp(1:3) = omega;
Omega = SkewSymMat(omega);
if (theta < 1e-10)
    V_inv = eye(3) - 0.5 * Omega + (1. / 12.) * (Omega * Omega);
    drdp(4:6) = V_inv * T(1:3,4);
else
    half_theta = (0.5) * theta;
    V_inv = (eye(3) - (0.5) * Omega + ((1) - theta * cos(half_theta) / ((2) * sin(half_theta))) / (theta * theta) * (Omega * Omega));
    drdp(4:6) = V_inv * T(1:3,4);
end
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

function [res, J_T_w_i_left_perturb, J_T_w_j_left_perturb, J_T_w_i_left_perturb_decoupled, J_T_w_j_left_perturb_decoupled] = compute_d_left_perturb_d_so3_trans_jac(T_ij_meas, T_w_i, T_w_j)
if 0
    T_w_i = [rodrigues(0.1 * (rand(3, 1) - 0.5)) 1 * (rand(3, 1) - 0.5); 0 0 0 1];
    T_w_j = [rodrigues(0.1 * (rand(3, 1) - 0.5)) 1 * (rand(3, 1) - 0.5); 0 0 0 1]; 
    T_ij_meas = inv(T_w_i) * T_w_j * [rodrigues(0.001 * (rand(3, 1) - 0.5)) 0.001 * (rand(3, 1) - 0.5); 0 0 0 1];
end

T_j_i = inv(T_w_j) * T_w_i;
res = LogSE3(T_ij_meas * T_j_i);

Jrinv = rightJacobianInvSE3Decoupled(res);

J_T_w_i_left_perturb = Jrinv * Adj(inv(T_w_i));
J_T_w_j_left_perturb = Jrinv * (-Adj(inv(T_w_i)));

d_left_se3_d_left_so3_trans_i = compute_d_left_se3_d_left_so3_trans(T_w_i);
d_left_se3_d_left_so3_trans_j = compute_d_left_se3_d_left_so3_trans(T_w_j);

J_T_w_i_left_perturb_decoupled = J_T_w_i_left_perturb * d_left_se3_d_left_so3_trans_i;
J_T_w_j_left_perturb_decoupled = J_T_w_j_left_perturb * d_left_se3_d_left_so3_trans_j;
end

function d_left_se3_d_left_so3_trans = compute_d_left_se3_d_left_so3_trans(T_w_i)
d_left_se3_d_left_so3_trans = eye(6);
d_left_se3_d_left_so3_trans(4:6,1:3) = SkewSymMat(T_w_i(1:3,4));

end