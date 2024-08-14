function testSE2Warp()

host = rand(2, 8);

s = 0.9;
theta = deg2rad(3);
t = [10;20];
% t = [0;0];
target = zeros(size(host));
for i = 1 : size(host,2)
    target(:, i) = s *[cos(theta) -sin(theta);sin(theta) cos(theta)] * host(:,i) + t;
    
end
T = [s *[cos(theta) -sin(theta);sin(theta) cos(theta)] t; 0 0 1];
t1 = T * pextend(host);
t2 = inv(T) * t1;
pextend(host) - t2
inv(T) - [[cos(theta) -sin(theta);sin(theta) cos(theta)]'/s [0;0]; 0 0 1]
theta0 = theta;
s0 = s;
target0 = target;


theta = theta + 2;
s = s - 2;
for iter = 1 : 10
   H = zeros(2, 2);
   b = zeros(2, 1);
   err_sum = 0;
   for i = 1 : size(host, 2)
       target_ = s *[cos(theta) -sin(theta);sin(theta) cos(theta)] * host(:,i) + t;
       err = target_ - target(:,i);
       d_err_d_theta = [t(2) - target_(2);target_(1) - t(1)];
       d_err_d_s = [(target_(1) - t(1))/s; (target_(2) - t(2))/s];
       d_err_d_state = [d_err_d_theta d_err_d_s];
       H = H + d_err_d_state' * d_err_d_state;
       b = b - d_err_d_state' * err;
       err_sum = err_sum + norm(err);
   end
   fprintf(sprintf('iter: %d, err: %f\n', iter, err_sum)); 
   
   dx = inv(H) * b;
   theta  = theta + dx(1);
   s  = s + dx(2);
end

end