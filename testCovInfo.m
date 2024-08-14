function testCovInfo()
clc;
clear;
R = rodrigues([1 2 3]);
J = 100*rand(8,3);
H = J' * J;

% inv((J*R') * (J*R')') - R * inv(H) * R'

inv((J*R')' * (J*R')) - R * inv(H) * R'
J_other = 1000*rand(3,3);
inv(J_other * H * inv(J_other)) - J_other * inv(H) * inv(J_other)

inv(inv(J_other) * H * (J_other)) - inv(J_other) * inv(H) * (J_other)

% (J * eye(3) * J')
(inv(J_other') * H * inv(J_other)) - inv((J_other) * inv(H) * J_other')

point = rand(3,1);
point = point ./norm(point);
point = sign(point(3)) .* point;
v1 = [0;0;1];
R = RotationBetweenPoints(v1, point);

R * v1 - point

R' * point

end
function R = RotationBetweenPoints(v1, v2)
v_hat = SkewSymMat(cross(v1, v2));
c = dot(v1, v2);

R = eye(3) + v_hat + ((v_hat * v_hat) / (1 + c));
end