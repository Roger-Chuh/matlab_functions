function testOnlineCalibrationIdp()

Tbci = [rodrigues([0.01 0.02 0.03]) [0.01 0.02 0.03]';0 0 0 1];
Tbcj = [rodrigues([0.01 -0.02 0.03]) [-0.01 0.02 0.03]';0 0 0 1];

Twbi = [rodrigues([0.01 0.02 -0.03]) [0.01 0.02 -0.03]';0 0 0 1];
Twbj = [rodrigues([0.01 -0.02 -0.03]) [0.01 -0.02 -0.03]';0 0 0 1];


K = [400 0 320; 0 400 240; 0 0 1];

uv = [350 150];

bearing = inv(K) * [uv 1]';
bearing = bearing./norm(bearing);

rho = 0.2;

T_th = inv(Twbj * Tbcj) * (Twbi * Tbci);

X2 = T_th(1:3,1:3) * bearing + rho * T_th(1:3,4);

Rbci = Tbci(1:3,1:3);
Rbcj = Tbcj(1:3,1:3);
tbci = Tbci(1:3,4);
tbcj = Tbcj(1:3,4);

Rwbi = Twbi(1:3,1:3);
Rwbj = Twbj(1:3,1:3);
twbi = Twbi(1:3,4);
twbj = Twbj(1:3,4);
Tbjbi = inv(Twbj) * Twbi;
Rbjbi = Tbjbi(1:3,1:3);



X2_check = Rbcj' * Rbjbi * Rbci * bearing + rho * (Rbcj' * Rbjbi * tbci + Rbcj'* Rwbj' * twbi - Rbcj' * tbcj - Rbcj' * Rwbj' * twbj);
% X2_check2 =  Rbcj' * Rbjbi * Rbci * bearing + rho * Rbcj' * (Rbjbi * tbci + Rwbj' * twbi - tbcj - Rwbj' * twbj);
X2_check22 =  Rbcj' * (Rbjbi * Rbci * bearing + rho * (Rbjbi * tbci + Rwbj' * twbi - tbcj - Rwbj' * twbj));

X2_check - X2

end