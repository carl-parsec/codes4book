function nonlin_laplace(N)
%nonlin_laplace(N) VERSION 12-9-2023
% solves the 2D nonlinear Laplace equation
%     div (epsilon(u) grad u) = 0, epsilon(u) = 1 + a u^2
% on the unit square with N dx by N dy, with dx = dy = h.
% Grid points are labeled i = 1,...,N+1 and j = 1,...,N+1.
% In u(i,j), i labels x and j labels y.
% BCs are Dirichlet with u = 0 on boundary except u(x=1,y) = 4 y (1-y).
% Uses Newton iteration and a sparse direct solve.
% Try: nonlin_laplace(100)

tic;

EPSILON = 10^-12;

fprintf('Newton method for 2D nonlinear Laplace equation\n');
fprintf('N = %g\n',N);
fprintf("EPSILON for Newton convergence = %g\n",EPSILON);

% h = 1/N; % dx = dy = h
k = (1:N+1)';
x = (k-1)/N;
y = (k-1)/N;
u = zeros(N+1,1);
u(N+1,k) = 4*y(k).*(1 - y(k)); % nonzero BC u(x=1,y)

F = zeros((N-1)^2,1);
delta_u = zeros(N+1,N+1);
% J = zeros((N-1)^2,(N-1)^2);
B = zeros((N-1)^2,5);

% calculate initial residual F and normresidual0
for i = 2:N
    for j = 2:N
        m = index(N,i,j);
        F(m) = residual(u,i,j);
    end
end
norm_residual0 = norm(F,1)/(N-1)^2;
norm_residual = norm_residual0;
fprintf("\t||r_0|| = %g\n",norm_residual0);

% Newton iteration
iter = 0;
while norm_residual > EPSILON*norm_residual0 && iter < 10
    iter = iter + 1;
    % calculate Jacobian J = dF/du
    for i = 2:N
        for j = 2:N
            m = index(N,i,j);
            B(m,3) = 0.5*depsilon(u(i,j))*(u(i+1,j) + u(i-1,j) + ...
                u(i,j+1) + u(i,j-1) - 4*u(i,j)) - ...
                0.5*(epsilon(u(i,j)) + epsilon(u(i+1,j))) - ...
                0.5*(epsilon(u(i-1,j)) + epsilon(u(i,j))) - ...
                0.5*(epsilon(u(i,j)) + epsilon(u(i,j+1))) - ...
                0.5*(epsilon(u(i,j-1)) + epsilon(u(i,j)));

            if i < N
                B(m,2) = 0.5*depsilon(u(i+1,j))*(u(i+1,j) - u(i,j)) + ...
                    0.5*(epsilon(u(i,j)) + epsilon(u(i+1,j)));
            end

            if i > 2
                B(m,4) = -0.5*depsilon(u(i-1,j))*(u(i,j) - u(i-1,j)) + ...
                    0.5*(epsilon(u(i-1,j)) + epsilon(u(i,j)));
            end

            if j < N
                B(m,1) = 0.5*depsilon(u(i,j+1))*(u(i,j+1) - u(i,j)) + ...
                    0.5*(epsilon(u(i,j)) + epsilon(u(i,j+1)));
            end

            if j > 2
                B(m,5) = -0.5*depsilon(u(i,j-1))*(u(i,j) - u(i,j-1)) + ...
                    0.5*(epsilon(u(i,j-1)) + epsilon(u(i,j)));
            end
        end
    end
    if N <= 4
        fprintf('\nB = \n');
        disp(B);
        fprintf('\n');
    end

    % A = spdiags(B,d,m,n) creates an m-by-n sparse matrix from the
    %     columns of B and places them along the diagonals specified by d.
    % D2 = spdiags([e -2*e e],[-1 0 1],N+1,N+1)/h^2;
    J = spdiags(B,[-(N-1) -1 0 1 (N-1)],(N-1)^2,(N-1)^2);

    % solve for change in states
    du = -J'\F; % transpose J for the linear solve
    % EPSILON_iter = 10^-6;
    % [L,U] = ilu(J',struct('type','ilutp','droptol',1e-6));
    % du = gmres(J',-F,50,EPSILON_iter,1000,L,U);
    for i = 2:N
        for j = 2:N
            m = index(N,i,j);
            delta_u(i,j) = du(m);
        end
    end
    u = u + delta_u;

    % calculate residual F and normresidual
    for i = 2:N
        for j = 2:N
            m = index(N,i,j);
            F(m) = residual(u,i,j);
        end
    end
    norm_residual = norm(F,1)/(N-1)^2;
    fprintf('Newton iteration = %d\n',iter);
    fprintf("\t||r|| = %g\n",norm_residual);
    fprintf("\t||r||/||r0|| = %g\n",norm_residual/norm_residual0);
end

if N <= 4
    fprintf('\nJ^T = \n');
    disp(full(J'));
    fprintf('\n');
end

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

function m = index(N,i,j)
m = (j-2)*(N-1) + (i-1);  % = (N-2)*(N-1) + (N-1) = (N-1)^2 at i,j = N
% m = (i-2)*(N-1) + (j-1); % linear convergence
end

function r = residual(u,i,j)
r = -(0.5*(epsilon(u(i-1,j)) + epsilon(u(i,j))) + ...
    0.5*(epsilon(u(i,j)) + epsilon(u(i+1,j))) + ...
    0.5*(epsilon(u(i,j-1)) + epsilon(u(i,j))) + ...
    0.5*(epsilon(u(i,j)) + epsilon(u(i,j+1))))*u(i,j) + ...
    0.5*(epsilon(u(i,j)) + epsilon(u(i+1,j)))*u(i+1,j) + ...
    0.5*(epsilon(u(i-1,j)) + epsilon(u(i,j)))*u(i-1,j) + ...
    0.5*(epsilon(u(i,j)) + epsilon(u(i,j+1)))*u(i,j+1) + ...
    0.5*(epsilon(u(i,j-1)) + epsilon(u(i,j)))*u(i,j-1);
end

function e = epsilon(u)
a = 2.5;
e = 1 + a*u^2;
end

function eprime = depsilon(u)
a = 2.5;
eprime = 2*a*u;
end









