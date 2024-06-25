function diffusion2(N,tf)
%diffusion2(N,tf) VERSION 8-22-2023
% solves the diffusion (heat) equation u_t = u_xx 
% with N dx (N+1 grid points 1,...,N+1) to time tf
% with Dirichlet BCs u(-1,t) = 0 = u(1,t) and
% ICs u(x,0) = 1 - x^2 or sin(pi*x).
% Uses TRBDF2 method with dynamic timestep.
% Try: diffusion2(100,0.25)

tic;

fprintf('TRBDF2 method\n');
fprintf('N = %g, tf = %g\n',N,tf);

t0 = 0;
dt_max = (tf - t0)/10; % require at least 10 timesteps
dt_min = 16*eps*(tf - t0); 
nmax = 100000; % max number of timesteps
GAMMA = 2 - sqrt(2);
EPSILON_REL = 10^-6; EPSILON_ABS = 10^-12;

h = 2/N;
% dt = h^2/2; % forward Euler stability limit
dt = sqrt(N)*h^2/2;

% Note there are N-1 interior points
e = ones(N-1,1);
D2 = spdiags([e -2*e e],[-1 0 1],N-1,N-1)/h^2;

j = (1:N-1)';
x = -1 + h*j;
% ICs and allocates memory for u
% u = 1 - x.^2;
u = sin(pi*x);

CONST = GAMMA/2;
CONST1 = (1 - GAMMA)/(2 - GAMMA);
CONST2 = 1/(GAMMA*(2 - GAMMA));
CONST3 = (1 - GAMMA)^2/(GAMMA*(2 - GAMMA));
k = (-3*GAMMA^2 + 4*GAMMA - 2)/(12*(2 - GAMMA));
C = 2*k;

t = t0;
n = 0;
while t < tf && n < nmax % timestep loop
    n = n + 1;
    r = Inf;
    while r > 2
        if 1.1*dt >= tf - t
            dt = tf - t;
        end
        tnew = t + dt;
        
        umid = (speye(N-1) - CONST*dt*D2)\((speye(N-1) + CONST*dt*D2)*u); % TR
        unew = (speye(N-1) - CONST1*dt*D2)\(CONST2*umid - CONST3*u); % BDF2
        
        norm_solution = norm(u,1)/(N-1);
        e_l = norm(C*dt*D2*(u/GAMMA - umid/(GAMMA*(1-GAMMA)) + ...
            unew/(1-GAMMA)),1)/(N-1); % divided-difference formula
        r = e_l/(EPSILON_REL*norm_solution + EPSILON_ABS);
        dt = min(min(2*dt,dt/r^(1/3)),dt_max);
        if dt <= dt_min
            warning('MATLAB:dt_min','dt = %e < dt_min at t = %e.\n',dt,t);
            n = nmax; r = 0; % to exit from while loops
        end 
    end
    t = tnew;
    u = unew;
    
    plot([-1; x; 1],[0; u; 0],'b-','LineWidth',4);
    set(gca,'fontsize',24);
    xlabel('x'); ylabel('u');
    axis([-1 1 -1.1 1.1]);
    getframe;
end

toc;

end









