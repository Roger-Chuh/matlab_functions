function testMultiViewTriangulate()


Twc{1,1} = [rodrigues([0.1 0.02 0.03]) [10 20 30]';0 0 0 1];
Twc{2,1} = [rodrigues([0.1 -0.02 0.03]) [10 20 -30]';0 0 0 1];
Twc{3,1} = [rodrigues([-0.1 0.02 0.03]) [-10 20 30]';0 0 0 1];
Twc{4,1} = [rodrigues([0.1 -0.02 -0.03]) [-10 -20 30]';0 0 0 1];



XYZ = [10 20 1300];


K = [400 0 320; 0 400 240;0 0 1];

Tcw = {};
for i = 1 : length(Twc)
    [ptIcs(i,:), tgtPt3d] = TransformAndProject(XYZ, K, Twc{i,1}(1:3, 1:3)', -Twc{i,1}(1:3, 1:3)'*Twc{i,1}(1:3, 4));
    
    bearing = inv(K) * [ptIcs(i,:) 1]';
    bearings(i,:) = bearing'./norm(bearing);
    Tcw{i,1} = inv(Twc{i,1});
end
 
pt3d = triangulate(Tcw, bearings);


end
function point_3d = triangulate(poses, points)

design_matrix = zeros(length(poses)*2, 4);
for (i = 1 : length(poses))
    p0x = points(i,1);
    p0y = points(i,2);
    p0z = points(i,3);
    design_matrix(i*2-1,:) = p0x * poses{i,1}(3,:) - p0z*poses{i,1}(1,:);
    design_matrix(i*2,:) = p0y * poses{i,1}(3,:) - p0z*poses{i,1}(2,:);
end


[A, B, C] = svd(design_matrix,0);
[x, y, z] = svd(design_matrix' * design_matrix,0);
triangulated_point = C (:, end);
point_3d = triangulated_point(1:3) ./ triangulated_point(4);


errs = design_matrix*[point_3d; 1];
error =  norm(errs)/ length(errs);

end