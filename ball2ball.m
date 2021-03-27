function [forces_x,forces_y,forces_z] = ball2ball(loc,q)
n = size(loc,1);
if size(q,1) ~= n
    forces_x= NaN;
    forces_y = NaN;
    forces_z= NaN;
    return
end

qq = q*q';
k = 9e9;
Dists = dist(loc');
% the matrix describe the reletive distance in x y z direction
Dxs = loc(:,1)'-loc(:,1);
Dys = loc(:,2)'-loc(:,2);
Dzs = loc(:,3)'-loc(:,3);

% The force in x y z direction
Forces_x = k.*qq.*Dxs./(abs(Dists).^3);
Forces_y = k.*qq.*Dys./(abs(Dists).^3);
Forces_z = k.*qq.*Dzs./(abs(Dists).^3);

% make the diagonal elements equal to 0
Forces_x(isnan(Forces_x)) = 0;
Forces_y(isnan(Forces_y)) = 0;
Forces_z(isnan(Forces_z)) = 0;

forces_x = sum(Forces_x);
forces_y = sum(Forces_y);
forces_z = sum(Forces_z);

end