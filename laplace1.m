function laplace1(N)
%laplace1(N) VERSION 12-9-2023
% solves the 2D Laplace equation 
%        u_xx + u_yy = 0 
% on the unit square with N dx by N dy, with dx = dy = h. 
% Grid points are labeled i = 1,...,N+1 and j = 1,...,N+1.
% In u(i,j), i labels x and j labels y.
% BCs are Dirichlet with u = 0 on boundary except u(x=1,y) = 4 y (1-y).
% Uses Jacobi iteration.
% Try: laplace1(100)

tic;

fprintf('Jacobi iteration\n');
fprintf('N = %g\n',N);

h = 1/N; % dx = dy = h
k = (1:N+1)';
x = (k-1)/N;
y = (k-1)/N;
u = zeros(N+1,N+1);
u(N+1,k) = 4*y(k).*(1 - y(k)); % nonzero BC u(x=1,y)

central = u(2:N,2:N);
north = u(2:N,3:N+1);
south = u(2:N,1:N-1);
east = u(3:N+1,2:N);
west = u(1:N-1,2:N);
residual = -4*central + north + south + east + west;
normresidual0 = norm(residual,1)/(N-1)^2;
normresidual = normresidual0;
EPSILON = 10^-5*h^2;
fprintf('EPSILON for convergence = %g\n',EPSILON);

iter = 0;
while normresidual > EPSILON*normresidual0 && iter < 10^6
    iter = iter + 1;
    u(2:N,2:N) = u(2:N,2:N) + residual/4;
    central = u(2:N,2:N);
    north = u(2:N,3:N+1);
    south = u(2:N,1:N-1);
    east = u(3:N+1,2:N);
    west = u(1:N-1,2:N);
    residual = -4*central + north + south + east + west;
    normresidual = norm(residual,1)/(N-1)^2;
end
fprintf('||r_k||/||r_0|| = %g\n',normresidual/normresidual0);
fprintf('Jacobi iterations = %g\n',iter);

toc;

u = u'; % transpose to plot

figure;
surf(x,y,u);
colormap jet;
shading interp;
set(gca,'fontsize',24); 
xlabel('x'); ylabel('y'); zlabel('Potential');

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









