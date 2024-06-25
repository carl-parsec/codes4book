function laplace6(N,method)
%laplace6(N,method) VERSION 8-23-2023
% solves the 2D Laplace equation
%       u_xx + u_yy = 0
% on the unit square with N dx by N dy, with dx = dy = h.
% Method (built-in MATLAB functions):
%     1 CG
%     2 PCG with ICHOL
%     3 GMRES with ICHOL
% Grid points are labeled i = 1,...,N+1 and j = 1,...,N+1.
% In u(i,j), i labels x and j labels y.
% BCs are Dirichlet with u = 0 on boundary except u(x=1,y) = 4 y (1-y).
% Uses CG, PCG, or GMRES
% Try: laplace6(200,2)

tic;

% h = 1/N; % dx = dy = h
k = (1:N+1)';
x = (k-1)/N;
y = (k-1)/N;
ufull = zeros(N+1,N+1);
ufull(N+1,k) = 4*y(k).*(1 - y(k)); % nonzero BC u(x=1,y)

b = zeros((N-1)^2,1);
for k = 1:N-1
    n = (N-1)*k;
    b(n) = 4*y(k)*(1 - y(k)); % solving -del^2 u = Au = -b
end

A = delsq(numgrid('S',N+1)); % Note this is -laplacian = -del^2
EPSILON = 10^-6;

if method == 1 % CG
    fprintf('CG method\n');
    u = pcg(A,b,EPSILON,10000);
elseif method == 2 % PCG
    fprintf('PCG method with ICHOL\n');
    % L = ichol(A);
    L = ichol(A,struct('michol','on')); % much faster
    u = pcg(A,b,EPSILON,1000,L,L');
elseif method == 3
    fprintf('GMRES method with ICHOL\n');
    % L = ichol(A);
    L = ichol(A,struct('michol','on')); % much faster
    u = gmres(A,b,50,EPSILON,1000,L,L');
else
    error('incorrect method type = %g',method);
end
fprintf('N = %g\n',N);
fprintf('EPSILON for convergence = %g\n',EPSILON);

toc;

u = reshape(u,N-1,N-1);
ufull(2:N,2:N) = u;
ufull = ufull'; % transpose to plot

figure;
surf(x,y,ufull);
colormap('jet');
shading flat;
set(gca,'fontsize',24); 
xlabel('x'); ylabel('y'); zlabel('Potential');

figure;
hh = pcolor(x,y,ufull);
set(hh,'edgecolor','none','facecolor','interp');
axis equal;
axis off;
set(gca,'fontsize',24);
title('Potential');
colormap hsv;
colorbar;

end









