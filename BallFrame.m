% For the spherical frame work

% initial set up for the cluster
clear
T = 50;
dt = 0.01;
N = floor(T/dt);
R = 0.03;
a = 1.5*R;
q0 = 1e-7;
n = 300;
mass = 1e2;
ratio = 100;
qc = 0;
col = 1;

% location of the balls 
x = zeros(1,n);
y = zeros(1,n);
z = zeros(1,n);
r = zeros(1,n);
E = zeros(1,N);
q = q0*ones(n,1);
Dr = 0.1;

for i=1:n
    r(i)=2*R;
    
    % put the stars in the sphere with radius R
    while r(i) > R
        x(i) = R*2*(rand-0.5);
        y(i) = R*2*(rand-0.5);
        z(i) = R*2*(rand-0.5);
        r(i) = sqrt((x(i))^2 + (y(i))^2 + (z(i))^2);
    end
end


% velocity of the balls
vx = zeros(1,n);
vy = zeros(1,n);
vz = zeros(1,n);

% plot the initial set up
zscaled = abs(q);                                                 
%cn = ceil(max(abs(zscaled)));                                        
cn = ceil(max(q));     
cm = colormap(winter(cn));   
fig = scatter3(x,y,z, [], cm(ceil(zscaled),:), 'filled');
axis equal
axis([-a a -a a -a a]);

for i = 1:N
    [Fx, Fy, Fz] = ball2ball([x',y',z'],q);
    vx = vx+dt.*Fx./mass-Dr.*vx.^2.*sign(vx);
    vy = vy+dt.*Fy./mass-Dr.*vy.^2.*sign(vy);
    vz = vz+dt.*Fz./mass-Dr.*vz.^2.*sign(vz);
    x = x+dt.*vx;
    y = y+dt.*vy;
    z = z+dt.*vz;
    rr = sqrt(x.^2+y.^2+z.^2);
    
    % if touch the boundary, do the reflection
    inx = find(sqrt(x.^2+y.^2+z.^2)>R);
    for j = inx
        r = sqrt(x(j)^2+y(j)^2+z(j)^2);
        x(j) = x(j)*R/r; y(j) = y(j)*R/r; z(j) = z(j)*R/r;
        m = [vx(j);vy(j);vz(j)]; n = [x(j);y(j);z(j)];
        nm = m-(2*(m'*n)./(n'*n))*n;
        vx(j) = col*nm(1); vy(j) = col*nm(2); vz(j) = col*nm(3);
        qtotal = qc+q(j);
        qc = qtotal*(ratio/(ratio+1));
        q(j) = qtotal/(ratio+1);
    end
    inx = find(x<0&y<0&z<0);
%     fig.XData = x(inx);
%     fig.YData = y(inx);
%     fig.ZData = z(inx);

%     fig.XData = x;
%     fig.YData = y;
%     fig.ZData = z; 
%     drawnow

    scatter3(x,y,z,[],[zeros(size(x,2), 2),rr'./R],'filled');
    axis equal
    axis([-a a -a a -a a]);
    
    drawnow    
    E(i) = sum(vx.^2+vy.^2+vz.^2);
end