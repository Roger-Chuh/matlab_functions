function testRotMatMinusJac()

R0 = rodrigues(rand(3, 1));

R_noise = R0 * rodrigues(0.01 * (rand(3, 1) - 0.5));






end

function [err, d_err_d_R] = ccomputeJac(R_obs, R_est)



end