function laplace0(N)
%laplace0(N) VERSION 8-23-2023
% solves the 2D Laplace equation
%       u_xx + u_yy = 0
% on the unit square with N dx by N dy, with dx = dy = h.
% Grid points are labeled i = 1,...,N+1 and j = 1,...,N+1.
% In u(i,j), i labels x and j labels y.
% BCs are Dirichlet with u = 0 on boundary except u(x=1,y) = 4 y (1-y).
% Uses banded direct solve.
% Try: laplace0(100)

tic;

fprintf("banded direct solve\n");
fprintf('N = %g\n',N);

% h = 1/N; % dx = dy = h
k = (1:N+1)';
x = (k-1)/N;
y = (k-1)/N;
ufull = zeros(N+1,N+1);
ufull(N+1,k) = 4*y(k).*(1 - y(k)); % nonzero BC u(x=1,y)

b = zeros((N-1)^2,1);
for k = 1:N-1
    n = (N-1)*k;
    b(n) = -4*y(k)*(1 - y(k));
end

% number grid points in (N-1)*(N-1) interior
G = numgrid('S',N+1); % S = square
% calculate Laplacian
A = -delsq(G); % -delsq is the Laplacian
% Matrix_A = full(A) % Set N = 4
% Vector_b = b

u = A\b;

toc;

u = reshape(u,N-1,N-1);
ufull(2:N,2:N) = u;
ufull = ufull'; % transpose to plot

figure;
surf(x,y,ufull);
colormap('jet');
shading flat;
set(gca,'fontsize',24); 
xlabel('x'); ylabel('y'); zlabel('u');

figure;
hh = pcolor(x,y,ufull);
set(hh,'edgecolor','none','facecolor','interp');
axis equal;
axis off;
set(gca,'fontsize',24);
title('Potential');
colormap hsv;
colorbar;

k = (1:N-1)';
x = k/N;
y = k/N;

uexact = fourier_sol(x,y,N);
numerr = sum(sum(abs(uexact - u)))/(N+1)^2;
fprintf('|err|_1 = %g\n',numerr);

    function ue = fourier_sol(x,y,N)
        fsum = zeros(N-1,N-1);
        Nfourier = 64;
        m = (1:Nfourier)';
        f = 4./(pi*(2*m-1));
        for i = 1:N-1
            for j = 1:N-1
                fsum(i,j) = 0;
                for m = 1:Nfourier
                    fsum(i,j) = fsum(i,j) + f(m)*sin((2*m-1)*pi*y(j))*...
                        sinh((2*m-1)*pi*x(i))/sinh((2*m-1)*pi);
                end
            end
            ue = fsum;
        end
    end

end









