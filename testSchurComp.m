function testSchurComp()


H = rand(10, 10);
H = 0.5 * (H + H');
dx = rand(10,1);
b = H * dx;

dim_marg = 3;
dim_all = size(H, 1);



HH = zeros(size(H));

HH(1 : (dim_all - dim_marg),1 : (dim_all - dim_marg)) = H(dim_marg + 1: end, dim_marg + 1: end);
HH((dim_all - dim_marg)+1:end,(dim_all - dim_marg)+1:end) = H(1: dim_marg, 1: dim_marg);
HH((dim_all - dim_marg)+1:dim_all,1 : (dim_all - dim_marg)) = H(1: dim_marg, dim_marg + 1: end);
HH(1 : (dim_all - dim_marg), (dim_all - dim_marg)+1:end) = H(dim_marg + 1: end, 1: dim_marg);

H21 = HH((dim_all - dim_marg)+1:dim_all,1 : (dim_all - dim_marg));
H11 = HH(1 : (dim_all - dim_marg),1 : (dim_all - dim_marg));
H22 = HH((dim_all - dim_marg)+1:end,(dim_all - dim_marg)+1:end);
b1 = b(dim_marg + 1: end);

b2 = b(1: dim_marg);

H11 = H11 - H21' * inv(H22) * H21;
b11 = b1 - H21' *inv(H22) * b2;
dx1 = inv(H11) * b11;
dx2 = inv(H22) * (b2 - H21 * dx1);

dx_new = [dx2; dx1];
dx_new - dx

end