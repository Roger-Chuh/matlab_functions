function testSE3SIM3()
Tc0ci = [rodrigues(rand(3, 1)) rand(3, 1); 0 0 0 1];
Tp0c0 = [rodrigues(rand(3, 1)) rand(3, 1); 0 0 0 1];
Tp0pj = [rodrigues(rand(3, 1)) rand(3, 1); 0 0 0 1];
s = 1.2;

X1 = rand(4,1);

X2 = inv(Tc0ci) * (1/s) * (inv(Tp0c0) * (s * Tp0pj * X1));
X22 = inv(Tc0ci) * ((1) * (inv(Tp0c0) * (1 * Tp0pj * X1)));


Xp0 = s * (Tp0pj(1:3,1:3) * X1(1:3) +  Tp0pj(1:3,4));
Xc0 = inv(Tp0c0) * [Xp0;1];
Xc0 = 1/s * Xc0(1:3);
Xci = inv(Tc0ci) * [Xc0;1]


Xci__ = Tc0ci(1:3,1:3)' * (Tp0c0(1:3,1:3)' * (Tp0pj(1:3,1:3) * X1(1:3) + Tp0pj(1:3,4)) - 1/s * (Tp0c0(1:3,1:3)' * Tp0c0(1:3,4))) - Tc0ci(1:3,1:3)' * Tc0ci(1:3,4);
Xci(1:3) - Xci__

return

s = 1;
Xp0_ = s * (Tp0pj(1:3,1:3) * X1(1:3) +  Tp0pj(1:3,4));
Xc0_ = inv(Tp0c0) * [Xp0_;1];
Xc0_ = 1/s * Xc0_(1:3);
Xci_ = inv(Tc0ci) * [Xc0_;1]

end