function testTrifocalTensor2()

use_bearing = true;


intrMat = [500 0 320; 0 500 240;0 0 1];
scale = 10;
xyz = [10 20 1000; 10 22 1100; 15 20 1020; -10 30 900; 17 17 800; 20 -12 950];
xyz = scale.* xyz;

Tc1w = [rodrigues([0.01 0.02 -0.03]) scale * [10 20 -30]';0 0 0 1];
Tc2w = [rodrigues([0.01 -0.02 0.03]) scale * [-10 20 -30]';0 0 0 1];
Tc3w = [rodrigues([-0.01 0.02 -0.03]) scale * [10 -20 30]';0 0 0 1];


pix1 = TransformAndProject(xyz, intrMat,Tc1w(1:3,1:3), Tc1w(1:3,4));
pix2 = TransformAndProject(xyz, intrMat,Tc2w(1:3,1:3), Tc2w(1:3,4));
pix3 = TransformAndProject(xyz, intrMat,Tc3w(1:3,1:3), Tc3w(1:3,4));

bearing1 = (inv(intrMat) * pextend(pix1'))';
% [bearing1, ~] = NormalizeVector(bearing1);
bearing2 = (inv(intrMat) * pextend(pix2'))';
% [bearing2, ~] = NormalizeVector(bearing2);
bearing3 = (inv(intrMat) * pextend(pix3'))';
% [bearing3, ~] = NormalizeVector(bearing3);

if use_bearing
    [bearing1, ~] = NormalizeVector(bearing1);
    [bearing2, ~] = NormalizeVector(bearing2);
    [bearing3, ~] = NormalizeVector(bearing3);
end

T = TFT_from_P(Tc1w(1:3,:),Tc2w(1:3,:),Tc3w(1:3,:));


for j = 1 : size(bearing1,1)
    err1(j, 1) = sum(sum(abs(SkewSymMat(bearing2(j,:)) * (bearing1(j, 1) * T(:,:,1) + bearing1(j, 2) * T(:,:,2) + bearing1(j, 3) * T(:,:,3)) * SkewSymMat(bearing3(j,:)))));
    [ point3 ] = pointTransfer( T, bearing1(j,:), bearing2(j,:) );
    point3 = point3./point3(3);
end


bearing33 = trifocalTransfer(Tc1w, Tc2w, Tc3w, eye(3), bearing1, bearing2, bearing3, use_bearing);

end
function T=TFT_from_P(P1,P2,P3)
%TFT_FROM_P Trifocal tensor from the projection matrices.
%
%  General formula to compute the trifocal tensor from any three projection
%  matrices.
%

T=zeros(3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            T(j,k,i)=(-1)^(i+1)*det([P1([1:(i-1) (i+1):3],:);P2(j,:);P3(k,:)]);
        end
    end
end
% T=T/norm(T(:));
end
function [ point3 ] = pointTransfer( Tri, point1, point2 )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
point3 = zeros(3,1);
i = 1;
j = 2;

for l=1:3,
    v1 = 0;
    v2 = 0;
    for k=1:3,
        v1 = v1+point1(k)*Tri(k,j,l);
        v2 = v2+point1(k)*Tri(k,i,l);
    end;
    
    point3(l) = point2(i)*v1 - point2(j)*v2;
end;
end
function bearing33 = trifocalTransfer(Tc1w, Tc2w, Tc3w, K, bearing1, bearing2, bearing3, use_bearing)
Twc1 = inv(Tc1w);
Twc2 = inv(Tc2w);
Twc3 = inv(Tc3w);

Rgc1 = Twc1(1:3,1:3);
Rgc2 = Twc2(1:3,1:3);
Rgc3 = Twc3(1:3,1:3);

pgc1 = Twc1(1:3,4);
pgc2 = Twc2(1:3,4);
pgc3 = Twc3(1:3,4);


R12 = Rgc1'*Rgc2;
t12 = Rgc1'*(pgc2-pgc1);
R13 = Rgc1'*Rgc3;
t13 = Rgc1'*(pgc3-pgc1);

PA = [ eye(3), zeros( 3, 1 ) ];
PB = [ R12', -R12'*t12 ];
PC = [ R13', -R13'*t13 ];

F12 = R12'*SkewSymMat( t12 );
F12_2 =  SkewSymMat( t12 )' * R12;
invK = inv(K);

z = zeros( 3, size( bearing1, 1 ) );
for m = 1:size( bearing1, 1 )
    if ~use_bearing
        uv1_normal = invK*[ bearing1(m,1:2)'; 1 ]; uv1_normal = uv1_normal/uv1_normal(3);
        uv2_normal = invK*[ bearing2(m,1:2)'; 1 ]; uv2_normal = uv2_normal/uv2_normal(3);
    else
        uv1_normal =  bearing1(m,1:3)';%./norm(bearing1(m,1:3)');
        uv2_normal =  bearing2(m,1:3)';%./norm(bearing2(m,1:3)');
    end
    epLine = F12*uv1_normal;
    %     epLine = F12_2'*uv1_normal;
    epLine = epLine./norm(epLine);
    if ~use_bearing
        epLineNormal = [epLine(2); -epLine(1); -uv2_normal(1)*epLine(2)+uv2_normal(2)*epLine(1)];
    else
        epLineNormal = cross(uv2_normal, epLine);
    end
    if 1
        epLineNormal = epLineNormal./norm(epLineNormal(1:3));
    end
    for k = 1:3
        for i = 1:3
            for j = 1:3
                Tijk = PB(j,i)*PC(k,4)-PB(j,4)*PC(k,i);
                z(k,m) = z(k,m)+uv1_normal(i)*epLineNormal(j)*Tijk;
            end
        end
    end
    
    z(:,m) = K*z(:,m);
    if ~use_bearing
        z(:,m) = z(:,m)/z(3,m);
    else
        z(:,m) = z(:,m)/norm(z(:,m));
    end
end

if 0
    z = z(1:2,:);
    z = z(:);
else
    diff = z' - bearing3;
end
bearing33 = z';
end