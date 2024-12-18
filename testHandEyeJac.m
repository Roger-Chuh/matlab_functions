function testHandEyeJac()

if 0
    test_d_left_perturb_d_so3_trans();
end


r = [0.1; 0.2; -0.3];

R1 = JrInv(r);
R2 = Jr(r);
R2L = Jr(-r);
% SkewSymMat(r) * R2

dr = [0.1;0.2;-0.1];
dt = [0.2;0.1;-0.2];
dT1 = Exp(dr, dt);
dT2 = [rodrigues(dr) Jr(-dr) * dt;0 0 0 1];

R = rodrigues(rand(3, 1));
R * rodrigues(dr) - rodrigues(dr) * R

use_decoupled_left_SE3 = false;


pose_num = 6;
Tbc = [rodrigues(0.1 * (rand(3, 1) - 0.5)) 1 * (rand(3, 1) - 0.5); 0 0 0 1];

for i = 1 : pose_num
    Twbs{i, 1} = [rodrigues(0.1 * (rand(3, 1) - 0.5)) 1 * (rand(3, 1) - 0.5); 0 0 0 1];
    Twcs{i, 1} = Twbs{i, 1} * Tbc;
end

max_iter = 20;


Tbc_est = Tbc * [rodrigues(0.5 * (rand(3, 1) - 0.5)) 10 * (rand(3, 1) - 0.5); 0 0 0 1];
for iter = 1 : max_iter
    
    H = zeros(6, 6);
    b = zeros(6, 1);
    err_sum = 0;
    for i = 1 : pose_num
        Twb1 = Twbs{i, 1};
        Twc1 = Twcs{i, 1};
        for j = i + 1 : pose_num
            Twb2 = Twbs{j, 1};
            Twc2 = Twcs{j, 1};
            if 0
                T_delta_b = inv(Twb2) * Twb1;
                T_delta_c = inv(Twc2) * Twc1;
            else
                T_delta_b = inv(Twb1) * Twb2;
                T_delta_c = inv(Twc1) * Twc2;
            end
            [err, d_err_d_Tbc] = computeJac(T_delta_b, T_delta_c, Tbc_est, use_decoupled_left_SE3);
            H = H + d_err_d_Tbc' * d_err_d_Tbc;
            b = b - d_err_d_Tbc' * err;
            err_sum = err_sum + norm(err);
        end
    end
    fprintf(sprintf('iter: %d, err_sum: %f\n', iter, err_sum));
    dx = inv(H) * b;
    dT = Exp(dx(1:3), dx(4:6));
    if ~use_decoupled_left_SE3
        Tbc_est = Tbc_est * dT;
    else
        if 1
            Tbc_est(1:3,1:3) = rodrigues(dx(1:3)) * Tbc_est(1:3,1:3);
            Tbc_est(1:3,4) = Tbc_est(1:3,4) + dx(4:6);
        else
            Tbc_est = dT * Tbc_est;
        end
    end
    
end


T1 = [rodrigues(rand(3, 1)) rand(3, 1);0 0 0 1];
T2 = [rodrigues(rand(3, 1)) rand(3, 1);0 0 0 1];
T3 = [rodrigues(rand(3, 1)) rand(3, 1);0 0 0 1];
adj_diff = Adj(T1 * T2 * T3) - Adj(T1 * T2) * Adj(T3);


end

function [err, d_err_d_Tbc] = computeJac(T_delta_b, T_delta_c, Tbc, use_decoupled_left_SE3)
d_err_d_Tbc = [];
T_err = inv(Tbc) * T_delta_b * Tbc * inv(T_delta_c);
err = LogSE3(T_err);

d_err_d_Tbc1 = -rightJacobianInvSE3Decoupled(-err);
d_err_d_Tbc2 = rightJacobianInvSE3Decoupled(err) * Adj(T_delta_c);

d_err_d_Tbc = d_err_d_Tbc1 + d_err_d_Tbc2;

if use_decoupled_left_SE3
    d_err_d_Tbc_right = d_err_d_Tbc;
    d_err_d_Tbc_left = d_err_d_Tbc_right * Adj(inv(Tbc));
    d_left_se3_d_left_so3_trans = compute_d_left_se3_d_left_so3_trans(Tbc);
    d_err_d_Tbc = d_err_d_Tbc_left * d_left_se3_d_left_so3_trans;
end

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
function J = Jr(phi)
EPSILON = 1e-10;

J = eye(3);

phi_norm2 = sum(phi.^2);
phi_hat = SkewSymMat(phi);
phi_hat2 = phi_hat * phi_hat;

if (phi_norm2 > EPSILON)
    phi_norm = sqrt(phi_norm2);
    phi_norm3 = phi_norm2 * phi_norm;
    
    J = J - phi_hat * (1 - cos(phi_norm)) / phi_norm2;
    J = J + phi_hat2 * (phi_norm - sin(phi_norm)) / phi_norm3;
else
    % sin and cos Taylor expansion around 0
    J = J - phi_hat / 2;
    J = J + phi_hat2 / 6;
end
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
