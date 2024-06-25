function laplaceNBC(N)
%laplaceNBC(N) VERSION 12-9-2023
% solves the 2D Laplace equation
%       u_xx + u_yy = 0
% on the unit square with (N+1) dx by N dy, with dx = dy = h.
% Grid points are labeled i = 1,...,N+1,N+2 and j = 1,...,N+1.
% In u(i,j), i labels x and j labels y.
% BCs are Dirichlet with u = 0 on boundary except for a Neumann BC
%     u_x(x=1,y) = -1.
% Uses SOR iteration.
% Try: laplaceNBC(100)

tic;

fprintf('SOR iteration\n');
fprintf('N = %g\n',N);

h = 1/N; % dx = dy = h
k = (1:N+1)';
x = (k-1)/N;
y = (k-1)/N;
u = zeros(N+2,N+1); % includes ghost points at right boundary
u(N+2,k) = u(N,k) - 2*h; % Neumann BC u_x(x=1,y) = -1 using ghost point 

% calculate optimal omega for SOR
mu = cos(pi/(N+1)); % Jacobi spectral radius
omega = 2*(1-sqrt(1-mu^2))/mu^2;
fprintf('SOR omega = %g\n',omega);

SUM = 0;
for i = 2:N+1
    for j = 2:N
        residual = -4*u(i,j) + u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1);
        SUM = SUM + abs(residual);  
    end
end
normresidual0 = SUM/(N+1)^2;
normresidual = normresidual0;
EPSILON = 10^-5*h^2;
fprintf('EPSILON for convergence = %g\n',EPSILON);
iter = 0;
while normresidual > EPSILON*normresidual0
    iter = iter + 1;
    for i = 2:N+1
        for j = 2:N
            residual = -4*u(i,j) + u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1);
            u(i,j) = u(i,j) + omega*residual/4;
        end
    end
    SUM = 0;
    for i = 2:N+1
        for j = 2:N
            residual = -4*u(i,j) + u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1);
            SUM = SUM + abs(residual);
        end
    end
    normresidual = SUM/(N+1)^2;
    u(N+2,k) = u(N,k) - 2*h; % Neumann BC u_x(x=1,y) = -1 
end
fprintf('||r_k||/||r_0|| = %g\n',normresidual/normresidual0);
fprintf('SOR iterations = %g\n',iter);

toc;

u = u(1:N+1,:); % omit ghost points for graphics
u = u'; % transpose to plot

figure;
surf(x,y,u);
colormap('jet');
shading interp;
set(gca,'fontsize',24); 
xlabel('x'); ylabel('y'); zlabel('u');

figure;
hh = pcolor(x,y,u);
set(hh,'edgecolor','none','facecolor','interp');
axis equal;
axis off;
set(gca,'fontsize',24);
title('Potential');
colormap hsv;
shading interp;
colorbar;

end









