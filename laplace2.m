function laplace2(N)
%laplace2(N) VERSION 12-9-2023
% solves the 2D Laplace equation 
%       u_xx + u_yy = 0 
% on the unit square with N dx by N dy, with dx = dy = h. 
% Grid points are labeled i = 1,...,N+1 and j = 1,...,N+1.
% In u(i,j), i labels x and j labels y.
% BCs are Dirichlet with u = 0 on boundary except u(x=1,y) = 4 y (1-y).
% Uses Gauss-Seidel iteration.
% Try: laplace2(100)

tic;

fprintf('Gauss-Seidel iteration\n');
fprintf('N = %g\n',N);

h = 1/N; % dx = dy = h
k = (1:N+1)';
x = (k-1)/N;
y = (k-1)/N;
u = zeros(N+1,N+1);
u(N+1,k) = 4*y(k).*(1 - y(k)); % nonzero BC u(x=1,y)

SUM = 0;
for i = 2:N
    for j = 2:N
        residual = -4*u(i,j) + u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1); 
        SUM = SUM + abs(residual);
    end
end
normresidual0 = SUM/(N-1)^2;
normresidual = normresidual0;
EPSILON = 10^-5*h^2;
fprintf('EPSILON for convergence = %g\n',EPSILON);

iter = 0;
while normresidual > EPSILON*normresidual0 && iter < 10^6
    iter = iter + 1;
    SUM = 0;
    for j = 2:N
        for i = 2:N
            residual = -4*u(i,j) + u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1);
            u(i,j) = u(i,j) + residual/4;
            SUM = SUM + abs(residual);
        end
    end

    % To converge on new residual rather than current residual, uncomment:
    SUM = 0;
    for i = 2:N
        for j = 2:N
            residual = -4*u(i,j) + u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1);
            SUM = SUM + abs(residual);
        end
    end

    normresidual = SUM/(N-1)^2;
end
fprintf('||r_k||/||r_0|| = %g\n',normresidual/normresidual0);
fprintf('Gauss-Seidel iterations = %g\n',iter);

toc;

u = u'; % transpose to plot

figure;
surf(x,y,u);
colormap jet;
shading interp;
set(gca,'fontsize',24); 
xlabel('x'); ylabel('y'); zlabel('u');

figure;
hh = pcolor(x,y,u);
set(hh,'edgecolor','none','facecolor','interp');
axis equal;
axis off;
set(gca,'fontsize',24);
title('u(x,y)');
colormap hsv;
shading interp;
colorbar;

u = u';
uexact = fourier_sol(x,y,N);
numerr = sum(sum(abs(uexact - u)))/(N+1)^2;
fprintf('|err|_1 = %g\n',numerr);

end

function ue = fourier_sol(x,y,N)
Nfourier = 100;
f = zeros(Nfourier,1);
fsum = zeros(N+1,N+1);
% m = (1:Nfourier)';
for m = 1:Nfourier
    f(m) = 32/(pi*(2*m-1))^3; % 4/(pi*(2*m-1));
end
k = (1:N+1)';
fsum(N+1,k) = 4*y(k).*(1-y(k));
for i = 2:N
    for j = 2:N
        fsum(i,j) = 0;
        for m = 1:Nfourier
            fsum(i,j) = fsum(i,j) + f(m)*sin((2*m-1)*pi*y(j))*...
                sinh((2*m-1)*pi*x(i))/sinh((2*m-1)*pi);
        end
    end
end
ue = fsum;
end









