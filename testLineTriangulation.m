function testLineTriangulation()


K = [499 0 321; 0 501 239; 0 0 1];

fix_first_at_identity = false;

Rwc1 = rodrigues([0.01 0.02 0.03]);
if fix_first_at_identity
    Rwc1 = eye(3);
end
Rwc2 = Rwc1*rodrigues([0.01 0.02 0.01]);
Rwc3 = Rwc2*rodrigues([-0.03 -0.01 0.01]);
Rwc4 = Rwc3*rodrigues([0.04 0.04 0.01]);


twc1 = [-100; 0; 100];
if fix_first_at_identity
    twc1 = [-0; 0; 0];
end
twc2 = [-0; 0; 200];
twc3 = [100; 0; 100];
twc4 = [-200; 0; 300];


ptS = [50; 20; 1000];
ptE = [100; -50; 800];

Tcw1 = inv([Rwc1 twc1; 0 0 0 1]);
Tcw2 = inv([Rwc2 twc2; 0 0 0 1]);
Tcw3 = inv([Rwc3 twc3; 0 0 0 1]);
Tcw4 = inv([Rwc4 twc4; 0 0 0 1]);

uvS1 = pflat(Tcw1(1:3,1:3)*ptS + Tcw1(1:3,4));
uvE1 = pflat(Tcw1(1:3,1:3)*ptE + Tcw1(1:3,4));
uvS2 = pflat(Tcw2(1:3,1:3)*ptS + Tcw2(1:3,4));
uvE2 = pflat(Tcw2(1:3,1:3)*ptE + Tcw2(1:3,4));
uvS3 = pflat(Tcw3(1:3,1:3)*ptS + Tcw3(1:3,4));
uvE3 = pflat(Tcw3(1:3,1:3)*ptE + Tcw3(1:3,4));
uvS4 = pflat(Tcw4(1:3,1:3)*ptS + Tcw4(1:3,4));
uvE4 = pflat(Tcw4(1:3,1:3)*ptE + Tcw4(1:3,4));


pii = pi_from_ppp(uvS1, uvE1, [0;0;0]);
pii = pii./norm(pii(1:3));
ni = pii(1:3);
ni = ni./norm(ni);

twc{1} = twc1;
twc{2} = twc2;
twc{3} = twc3;
twc{4} = twc4;

Rwc{1} = Rwc1;
Rwc{2} = Rwc2;
Rwc{3} = Rwc3;
Rwc{4} = Rwc4;

t0 =twc{1};
R0 = Rwc{1};

obs = [uvS1(1:2)' uvE1(1:2)'; uvS2(1:2)' uvE2(1:2)';uvS3(1:2)' uvE3(1:2)';uvS4(1:2)' uvE4(1:2)'];
obs = obs;
% obs = obs + 0.001*(rand(size(obs))-0.5);

for i = 2 : 4
    t1 = twc{i};
    R1 = Rwc{i};
    
    tij = R0' * (t1 - t0);   % tij
    Rij = R0' * R1;          % Rij
    
    
    p3 = [obs(i,1:2)';1];%(obsj_tmp(0), obsj_tmp(1), 1);
    p4 = [obs(i,3:4)';1];  %Vector3d p4(obsj_tmp(2), obsj_tmp(3), 1);
    p3 = p3 ./norm(p3);
    p4 = p4 ./norm(p4);
    p3 = Rij * p3 + tij;
    p4 = Rij * p4 + tij;
    pij = pi_from_ppp(p3, p4,tij);
    pij = pij./norm(pij(1:3));
    plk = pipi_plk( pii, pij );
    plk = plk./norm(plk(4:6));
    PLK(:,i-1) = plk;
    
    n(:,i-1) = plk(1:3); %plk(1:3)./norm(plk(1:3));
    v(:,i-1) = plk(4:6); %plk(4:6)./norm(plk(4:6));
    
end


% n v 是第一帧坐标系下的plk坐标。
%【这个例子中，第一帧在世界系下的pose并不是eye(4),12，13，14帧上的线三角化出的plk都是表达在1坐标系下的，并且结果都是一样的】


C0 = [n(:,1) v(:,1)];
kObvNum = 4; 3;
%  Eigen::Matrix<double, 2 * kObvNum, 6> A_temp;
A_temp = zeros(2*kObvNum, 6);
for  i = kObvNum : -1 : 1
    
    temp = zeros(3,6);
    if 1
        temp(:,1:3) = Rwc{i}';
        temp(:,4:6) = SkewSymMat(-Rwc{i}'*twc{i}) * Rwc{i}';
    else
        temp(:,1:3) = Rwc{i};
        temp(:,4:6) = SkewSymMat(twc{i}) * Rwc{i};
    end
    A_temp(2*i-1, : ) = [obs(i,1:2) 1] * temp;
    A_temp(i *2 ,:) = [obs(i,3:4) 1] * temp;
end
A = A_temp;
[U1,S1,V1] = svd(A,0);
para = V1(:,end);
para = para./norm(para(4:6));
nn = para(1:3);  % V1(1:3,end)./norm(V1(1:3,end));
vv = para(4:6); %V1(4:6,end)./norm(V1(4:6,end));

nn_ = nn;
nn_(1) = -(nn(2)*vv(2) + nn(3)*vv(3))/vv(1);

err = [n(:, 1) nn v(:, 1) vv];
dot(nn,vv)
C = [nn vv];
[U,S,V] = svd(C,0);
% Z的行向量是正交的
Z = S*V';
Z1 = Z(:,1);
Z2 = Z(:,2);
TT = [Z2'; Z1'*[0 -1; 1 0]];
T = [Z(2,1) Z(2,2); Z(1,2) -Z(1,1)];
[U_,S_,V_] = svd(T, 0);
[U_2,S_2,V_2] = svd(TT, 0);
if 1
    V11 = V_(1,end);
    V22 = V_(2,end);
    V11_ = V_2(1,end);
    V22_ = V_2(2,end);
else
    V11 = U_(1,end);
    V22 = U_(2,end);
end

VV = [V11 -V22;V22 V11];
S__ = diag(VV'*S*V');
Z_ = VV*diag(S__);
C_ = U*Z_;

VV_ = [V11_ -V22_;V22_ V11_];
S__1 = diag(VV_'*S*V');
Z_1 = VV_*diag(S__1);
C_1 = U*Z_1;
%         U2*V3*diag()
%
%         S2*VV'
%         V3*

if 0
    line_w = line_to_pose( [C_(:,1); C_(:,2)], Tcw1(1:3,1:3), Tcw1(1:3,4));
    plk_w = plk_w./norm(plk_w(4:6));
    CC = [line_w(1:3)./norm(line_w(1:3)) line_w(4:6)./norm(line_w(4:6))];
else
%     plk_w = plk_to_pose( Rwc{1}',  -Rwc{1}'*twc{1},  [C0(:,1); C0(:,2)] );
C00 = 1.3.*C0;
    plk_w = plk_to_pose( Rwc{1},  twc{1},  [C00(:,1); C00(:,2)] );
    %     CC = [plk_w(1:3)./norm(plk_w(1:3)) plk_w(4:6)./norm(plk_w(4:6))];
    plk_w = plk_w./norm(plk_w(4:6));
    CC = [plk_w(1:3) plk_w(4:6)];
end


dist_pose1_origin_to_line = CC(:,1) - cross(twc1, CC(:,2));
dist_pose1_origin_to_line_check = C0(:,1) - cross(zeros(3,1), C0(:,2));

dist_diff = norm(dist_pose1_origin_to_line_check) - norm(dist_pose1_origin_to_line);

C0;
[C C_]
norm(C(:) - C_(:))


AA = A;
[a2,b2,c2] = svd(AA,0);
[a1,b1,c1] = svd(AA);
inv(a1)*AA*inv(c1') - b1;
pinv(a2)*AA*inv(c2') - b2;

%figure,line([uvS1(1) uvE1(1) uvS2(1) uvE2(1)], [uvS1(2)  uvE1(2) uvS2(2) uvE2(2)],'Color', [1 0 0])
end
function pi = pi_from_ppp(x1, x2, x3)

pi = [cross( x1 - x3 ,( x2 - x3 )); dot(- x3,( cross(x1, x2 )) )]; %// d = - x3.dot( (x1-x3).cross( x2-x3 ) ) = - x3.dot( x1.cross( x2 ) )];


end
function plk = pipi_plk(pi1, pi2)

dp = pi1 * pi2' - pi2 * pi1';

plk = [dp(1,4), dp(2,4), dp(3,4), - dp(2,3), dp(1,3), - dp(1,2)]';
end
function line_c = line_to_pose( line_w, Rcw, tcw)



cp_w = line_w(1:3);
dv_w = line_w(4:6);


cp_c = Rcw * cp_w + tcw;
dv_c = Rcw* dv_w;

line_c(1:3,1) = cp_c;
line_c(4:6,1) = dv_c;


end
function plk_c = plk_to_pose( Rcw,  tcw,  plk_w )
nw = plk_w(1:3);
vw = plk_w(4:6);

nc = Rcw * nw + SkewSymMat(tcw) * Rcw * vw;
vc = Rcw * vw;


plk_c(1:3,1) = nc;
plk_c(4:6,1) = vc;



end