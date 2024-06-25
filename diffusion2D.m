function diffusion2D(N,tf)
%diffusion2D(N,tf) VERSION 8-22-2023
% solves the 2D diffusion (heat) equation 
%     u_t = u_xx + u_yy
% on the unit square with N dx x N dy ((N+1)*(N+1) grid points) to time tf
% with homogeneous Dirichlet BCs and a Gaussian IC.
% Uses TRBDF2 method and fixed dt for movie.
% Try: diffusion2D(100,0.01)

tic;

fprintf('TRBDF2 method\n');
fprintf('N = %g, tf = %g\n',N,tf);

GAMMA = 2 - sqrt(2);
h = 1/N;
k = (1:N+1)';
x = (k-1)*h; y = (k-1)*h;  
% dt = h^2/2; % forward Euler stability limit
dt = sqrt(N)*h^2/2; % fixed dt for movie
steps = round(tf/dt);

u = zeros(N-1,N-1); % interior solution values
U = zeros(N+1,N+1); % full region
i = 1:N-1; j = 1:N-1;
u(i,j) = exp(-40*((y(i+1)-0.5).^2))*exp(-40*((x(j+1)'-0.5).^2)); % ICs

% number grid points in (N-1)*(N-1) interior
G = numgrid('S',N+1);
% calculate Laplacian D2
D2 = -delsq(G)/h^2; % Note the MINUS sign!
% full(h^2*D2) % set N = 4

CONST = GAMMA/2;
CONST1 = (1 - GAMMA)/(2 - GAMMA);
CONST2 = 1/(GAMMA*(2 - GAMMA));
CONST3 = (1 - GAMMA)^2/(GAMMA*(2 - GAMMA));
for n = 1:steps % timestep loop
    u = reshape(u,(N-1)^2,1);   
    % u = (speye((N-1)^2) - dt*D2)\u; % backward Euler
    umid = (speye((N-1)^2) - CONST*dt*D2)\((speye((N-1)^2) + ...
        CONST*dt*D2)*u); % TR
    u = (speye((N-1)^2) - CONST1*dt*D2)\(CONST2*umid - CONST3*u); % BDF2   
    u = reshape(u,N-1,N-1);
    
    % for movie:
    U(i+1,j+1) = u(i,j);
    surf(x,y,U,'FaceColor','interp','EdgeColor','none','FaceLighting',...
        'gouraud');
    view(-50,30);
    camlight left;
    set(gca,'fontsize',24);
    xlabel('x'); ylabel('y'); zlabel('u');
    axis([0 1 0 1 0 1.1]);
    getframe;
end

toc;

end









