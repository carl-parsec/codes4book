function diffusion0(N,steps)
%diffusion0(N,steps) VERSION 8-22-2023
% solves the diffusion (heat) equation u_t = u_xx 
% with N dx (N+1 grid points 1,...,N+1) to time t = steps*dt
% with Dirichlet BCs u(-1,t) = 0 = u(1,t) and 
% ICs u(x,0) = 1 - x^2 or sin(pi*x).
%   For heat equation, u(x,t) = T(x,t) - T_amb.
%   For diffusion equation, u(x,t) = rho(x,t) - rho_amb.
% Uses backward Euler method.
% Try: diffusion0(100,50)

tic;

fprintf('N = %g, steps = %g\n',N,steps);
h = 2/N; % h = dx
% dt = h^2/2; % forward Euler stability limit dt <= h^2/2
dt = sqrt(N)*h^2/2; % fixed timestep
fprintf('tf = %g\n',steps*dt);
M = struct('cdata',cell(1,steps),'colormap',cell(1,steps)); % for movie

% Note there are N-1 interior points
e = ones(N-1,1);
D2 = spdiags([e -2*e e],[-1 0 1],N-1,N-1)/h^2; % 2nd derivative matrix
if N <= 16
    display(h^2*D2);
    display(full(h^2*D2));
end

j = (1:N-1)'; % ' denotes transpose
x = -1 + h*j;
% ICs and allocates memory for u
u = 1 - x.^2; 
% u = sin(pi*x); 

figure; % for movie
for n = 1:steps % timestep loop     
    % u = u + dt*D2*u; % forward Euler   
    u = (speye(N-1) - dt*D2)\u; % backward Euler
    
    plot([-1; x; 1],[0; u; 0],'b-','LineWidth',2); % prepend/append BCs
    set(gca,'fontsize',24);
    xlabel('x'); ylabel('u');
    axis([-1 1 0 1.1]); % for 1 - x^2 ICs
    % axis([-1 1 -1.1 1.1]); % for sin(pi*x) ICs
    M(n) = getframe; % movie
end

toc;

movie(M,1,5); % replay movie

figure;
plot([-1; x; 1],[0; u; 0],'b-','LineWidth',2); % prepend/append BCs
set(gca,'fontsize',24);
xlabel('x'); ylabel('u');
axis([-1 1 0 1.1]); % for 1 - x^2 ICs

end









