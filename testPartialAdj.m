function testPartialAdj()

J1 = rand(8, 6);
J1_bak = J1;
J2 = rand(8,2);
res = rand(8,1);
Adj = rand(6, 6);

H = zeros(8,8);
b = zeros(8,1);

H(1:6, 1:6) = H(1:6, 1:6) + J1' * J1;
H(7:8, 7:8) = H(7:8, 7:8) + J2' * J2;
b(1:6) = b(1:6) + J1' * res;
b(7:8) = b(7:8) + J2' * res;
H(7:8, 1:6) = H(7:8, 1:6) + J2' * J1;
H(1:6, 7:8) = H(1:6, 7:8) + J1' * J2;

H_bak = H;
b_bak = b;
% % % % % % % % % % % % % % % 
J1 = J1 * Adj;

H = zeros(8,8);
b = zeros(8,1);

H(1:6, 1:6) = H(1:6, 1:6) + J1' * J1;
H(7:8, 7:8) = H(7:8, 7:8) + J2' * J2;
b(1:6) = b(1:6) + J1' * res;
b(7:8) = b(7:8) + J2' * res;
H(7:8, 1:6) = H(7:8, 1:6) + J2' * J1;
H(1:6, 7:8) = H(1:6, 7:8) + J1' * J2;

H_new = H;
b_new = b;

H_comp = zeros(8,8);
b_comp = zeros(8,1);

H_comp(7:8,7:8) = H_bak(7:8,7:8);

H_comp(1:6, 1:6) = Adj' * H_bak(1:6, 1:6) * Adj;
H_comp(7:8, 1:6) = H_bak(7:8, 1:6) * Adj;
H_comp(1:6, 7:8) = Adj' * H_bak(1:6, 7:8);

H_new - H_comp

b_comp(1:6) = Adj' * b_bak(1:6);
b_comp(7:8) = b_bak(7:8);

b_new' - b_comp'


T1 = [rodrigues(rand(3, 1)) rand(3, 1);0 0 0 1];
adj1 = Adjoint(T1);
adj2 = Adjoint(inv(T1));
inv(adj1) - adj2

end

function jac = Adjoint(T)
R = T(1:3,1:3);
jac = zeros(6,6);
jac(1:3,1:3) = R;
jac(4:6, 4:6) = R;
jac(4:6,1:3) = SkewSymMat(T(1:3,4)) * R;
jac(1:3, 4:6) = zeros(3, 3);

end
