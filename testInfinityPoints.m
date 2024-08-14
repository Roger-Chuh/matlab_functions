function testInfinityPoints()

close all

T = [rodrigues([-0.9 0.2 0.3]) [-100 -100 200]';0 0 0 1];
% T = [rodrigues([0 0 0]) [-100 -100 200]';0 0 0 1];
K = [300 0 320; 0 300 240; 0 0 1];

point_num = 20;
x_ = randperm(640);
y_ = randperm(480);
pix = [x_(1:point_num)' y_(1:point_num)'];

pix = [[314,340;214,419;118,349;443,194;263,107;68,253;550,299;376,241;613,212;450,329;315,129;587,435;610,392;252,478;217,364;298,80;64,402;291,246;143,378;629,79]];

bearing = (inv(K) * pextend(pix'))';


figure,imshow(ones(480, 640));hold on;
depth_near = 1000;100;-20;100;
use_idp = true;

if ~use_idp
    step = 10;
    list = step : step : 10000000;
    list = list(1:5);
else
    
    idepth_near = 1 / depth_near;
    step1 = 0.0001;
    step2 = -step1;
    list = step1 : step1 : 10000;
    list = list(1:17);
    list = fliplr(list);
end
list = list-mean(list);
cnt = 1;
bearings = zeros(point_num, 10);
list = [list(1) 0  list(end)];
for k = list %1 : 100
    if ~use_idp
        depth_far = depth_near + k; %step * k; 1000;
    else
        depth_far = 1 / (idepth_near + k);
    end
    xyz_near = bearing.*depth_near;
    xyz_far = bearing.*depth_far;
    
    for i = 1 : point_num
        [pt2d_near(i,:), tgtPt3d_host(i,:)] = TransformAndProject(xyz_near(i,:), K, T(1:3, 1:3), T(1:3, 4));
        [pt2d_far(i,:), tgtPt3d_far] = TransformAndProject(xyz_far(i,:), K, T(1:3, 1:3), T(1:3, 4));
        [pt2d_near_rot(i,:), tgtPt3d] = TransformAndProject(xyz_near(i,:), K, T(1:3, 1:3), zeros(3,1));
        [pt2d_far_rot(i,:), tgtPt3d_list(i,:)] = TransformAndProject(xyz_far(i,:), K, T(1:3, 1:3), zeros(3,1));
    end
    [tgtPt3d_list, ~] = NormalizeVector(tgtPt3d_list);
    plot(pt2d_near(:,1), pt2d_near(:,2), 'or');
    plot(pt2d_far(:,1), pt2d_far(:,2), 'og');
    plot(pt2d_near_rot(:,1), pt2d_near_rot(:,2), 'sc');
    plot(pt2d_far_rot(:,1), pt2d_far_rot(:,2), 'sy');
    if(cnt == 1)
        plot(pt2d_far(:,1), pt2d_far(:,2), 'oc');
        bearing_far = (inv(K) * pextend(pt2d_far'))';
        [bearing_far, ~] = NormalizeVector(bearing_far);
        bearing_host = (inv(K) * pextend(pt2d_near'))';
        for j = 1 : point_num
%            ang(j,1) = CalcDegree(bearing_far(j,:), bearing(j,:)./norm(bearing(j,:))); 
%            ang(j,1) = acosd(dot(bearing_far(j,:), bearing(j,:)./norm(bearing(j,:)))); 
           ang(j,1) = acosd(dot(bearing_far(j,:), tgtPt3d_host(j,:)./norm(tgtPt3d_host(j,:)))); 
%            tgtPt3d_host(i,:)
        end
        [~, epi_len(:,1)] = NormalizeVector(pt2d_far  - pt2d_near);
        bearings(:,1:3) = bearing_far;
        bearings(:,4:5) = pt2d_far;
    end
    if(cnt == length(list))
        plot(pt2d_far(:,1), pt2d_far(:,2), 'om');
        bearing_far = (inv(K) * pextend(pt2d_far'))';
        [bearing_far, ~] = NormalizeVector(bearing_far);
        for j = 1 : point_num
%            ang(j,2) = CalcDegree(bearing_far(j,:), bearing(j,:)./norm(bearing(j,:))); 
%            ang(j,2) = acosd(dot(bearing_far(j,:), bearing(j,:)./norm(bearing(j,:)))); 
           ang(j,2) = acosd(dot(bearing_far(j,:), tgtPt3d_host(j,:)./norm(tgtPt3d_host(j,:)))); 
        end
        [~, epi_len(:,2)] = NormalizeVector(pt2d_far  - pt2d_near);
        bearings(:,6:8) = bearing_far;
        bearings(:,9:10) = pt2d_far;
    end
    drawnow;
    cnt = cnt + 1;
end
uint16([ang epi_len])
[ang(:,1)./ang(:,2) epi_len(:,1)./epi_len(:,2)]
pole = pflat(K * T(1:3,4));

Lines = [];
for i = 1 : point_num
    pt1 = pt2d_near_rot(i,:);
    pt2 = pt2d_near(i,:);
    pt3 = pt2d_far(i,:);
    
    line_para1 = GetLine(pt1, pt2);
    err1_1 = dot(line_para1, [pt1 1]);
    err1_2 = dot(line_para1, [pt2 1]);
    err1_3 = dot(line_para1, pole);
    
    line_para2 = GetLine(pt1, pt3);
    err2_1 = dot(line_para2, [pt1 1]);
    err2_2 = dot(line_para2, [pt3 1]);
    err2_3 = dot(line_para2, pole);
    
    line_para1 - line_para2;
    Lines = [Lines; line_para1];
end

points = lineToBorderPoints(Lines,[1000 1000]);
for k = 1 : size(points,1)
    line(points(k,[1,3])',points(k,[2,4])','Color', [1 0 1]);plot(pole(1), pole(2),'*b');
    text(pt2d_near(k,1)+1,pt2d_near(k,2)-1, num2str(k));
end

check_id = 4;
pt0 = pt2d_near(check_id,:);
pt1 = bearings(check_id,4:5);
pt2 = bearings(check_id,9:10);
vec0 = inv(K)*[pt0 1]';
vec1 = inv(K)*[pt1 1]';
vec2 = inv(K)*[pt2 1]';
len1 = norm(pt0 - pt1);
len2 = norm(pt0 - pt2);
ang1 = CalcDegree(vec0',vec1');
ang2 = CalcDegree(vec0',vec2');

[[ang1 / ang2] [len1 / len2]]
[ang(check_id,1)./ang(check_id,2) epi_len(check_id,1)./epi_len(check_id,2)]
end
function line_para1 = GetLine(pt1, pt2)
slope1 = (pt1(2) - pt2(2)) / (pt1(1) - pt2(1));
line_para1 = [slope1, -1, -slope1 * pt1(1) + pt1(2)];
line_para1 = line_para1./norm(line_para1(1:2));
end