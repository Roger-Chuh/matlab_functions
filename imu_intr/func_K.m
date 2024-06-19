function [y] = func_K(x)


%y = zeros(3,6);

y = [x(1) x(2) x(3) 0 0 0;0 0 0 x(1) x(2) 0;0 0 0 0 0 x(3)];


end

