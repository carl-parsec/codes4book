function laplace5(N)
%laplace5(N) VERSION 12-9-2023
% solves the 2D Laplace equation
%       u_xx + u_yy = 0
% on the unit square with N dx by N dy, with dx = dy = h.
% Grid points are labeled i = 1,...,N+1 and j = 1,...,N+1.
% In u(i,j), i labels x and j labels y.
% BCs are Dirichlet with u = 0 on boundary except u(x=1,y) = 4 y (1-y).
% Uses Preconditioned Conjugate Gradient Method with
%     M^-1 = Minv = (I - L*Dinv)*(I - Dinv*L')
% due to Ament, Knittel, Weiskopf, & Strasser (2010).
% Try: laplace5(100)

tic;

fprintf("PCG method with Minv = (I - L*Dinv)*(I - Dinv*L')\n");
fprintf('N = %g\n',N);

% h = 1/N; % dx = dy = h
k = (1:N+1)';
x = (k-1)/N;
y = (k-1)/N;
ufull = zeros(N+1,N+1);
ufull(N+1,k) = 4*y(k).*(1 - y(k)); % nonzero BC u(x=1,y)
u = ufull(2:N,2:N); 
u = reshape(u,(N-1)^2,1);

b = zeros((N-1)^2,1);
for k = 1:N-1
    n = (N-1)*k;
    b(n) = -4*y(k)*(1 - y(k));
end

% number grid points in (N-1)*(N-1) interior
G = numgrid('S',N+1);
% calculate Laplacian
A = -delsq(G);
% Matrix_A = full(A) % Set N = 4
% Vector_b = b

L = tril(A,-1);
D = diag(diag(A));
% Dinv = inv(D);
% Minv = (speye((N-1)^2) - L*Dinv)*(speye((N-1)^2) - Dinv*L');
% In general:
% replace inv(A)*b with A\b 
% replace b*inv(A) with b/A 
Minv = (speye((N-1)^2) - L/D)*(speye((N-1)^2) - D\L');

residual = b - A*u;
z = Minv*residual;
d = z;
zr_old = z'*residual;
normresidual0 = norm(residual,1)/(N-1)^2;
normresidual = normresidual0;
% EPSILON = 10^-5*h^2;
EPSILON = 10^-6;
fprintf('EPSILON for convergence = %g\n',EPSILON);
iter = 0;
while iter < length(b) && normresidual > EPSILON*normresidual0
    iter = iter + 1;
    Ad = A*d;
    alpha = zr_old/(d'*Ad);
    u = u + alpha*d;
    residual = residual - alpha*Ad;
    z = Minv*residual;
    zr_new = z'*residual;
    d = z + (zr_new/zr_old)*d;
    zr_old = zr_new;
    normresidual = norm(residual,1)/(N-1)^2;
end
fprintf('||r_k||/||r_0|| = %g\n',normresidual/normresidual0);
fprintf('PCG iterations = %g\n',iter);

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









