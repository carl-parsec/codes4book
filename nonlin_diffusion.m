function nonlin_diffusion(N,tf)
%nonlin_diffusion(N,tf) VERSION 12-8-2023
% solves the nonlinear diffusion equation
%     u_t = (d(u) u_x)_x
% where d(u) = a u + b is the diffusion coefficient,
% with N dx from time t0 = 0 to tf with Neumann BCs
%     u_x(0,t) = 0 = u_x(xmax,t)
% and ICs u(x,0) = 1.5*exp(-2*(x-x0).^2) (Gaussian implant).
%                                            x-x-- --x-x
% Grid pts for Neumann BCs with ghost pts: 2 1 2 ... N N+1 N
% Uses TRBDF2 method with dynamic timestep and Newton's method.
% Try: nonlin_diffusion(100,0.5)

% Computational Units
%   Fundamental units
%       L_scale = 0.1 micron = 10^-5 cm
%       t_scale = 1000 sec
%       u_scale = 10^20 cm^-3
%   Derived units
%       d_scale = 10^-13 cm^2/sec = L_scale^2/t_scale

tic;

VERBOSE = 1; % 1 (verbose mode) or 0 (quiet mode)

fprintf('nonlinear diffusion using TRBDF2 method\n');
fprintf('N = %g\n',N);

xmax = 10; % = 1 micron
t0 = 0;
dt_max = (tf - t0)/10; % require at least 10 timesteps
dt_min = 16*eps*(tf - t0); 
nmax = 100000; % max number of timesteps
GAMMA = 2 - sqrt(2);
EPSILON_REL = 10^-6; EPSILON_ABS = 10^-15;

h = xmax/N;
dt = sqrt(N)*h^2/2;

j = (1:N+1)';
x = h*(j-1);
% ICs and allocates memory for u
x0 = xmax/10;
u = 1.5*exp(-2*(x-x0).^2); % Gaussian implant
Q0 = total_number(u,h);
fprintf('total number of particles at t = 0: Q0 = %g\n',Q0);

% Note: There are N+1 interior points for the Neumann problem
e = ones(N+1,1);
Dplus = spdiags([-e e],[0 1],N+1,N+1)/h;
Dplus(N+1,N) = 1/h; % for Neumann BC at right
Dminus = spdiags([-e e],[-1 0],N+1,N+1)/h;
Dminus(1,2) = -1/h;  % for Neumann BC at left

D2 = spdiags([e -2*e e],[-1 0 1],N+1,N+1)/h^2;
D2(1,2) = 2/h^2; % for Neumann BC at left
D2(N+1,N) = 2/h^2; % for Neumann BC at right

CONST = GAMMA/2;
CONST1 = (1 - GAMMA)/(2 - GAMMA);
CONST2 = 1/(GAMMA*(2 - GAMMA));
CONST3 = (1 - GAMMA)^2/(GAMMA*(2 - GAMMA));
k = (-3*GAMMA^2 + 4*GAMMA - 2)/(12*(2 - GAMMA));

t = t0;
n = 0;
tstep = 0;
figure; % for movie
while t < tf && n < nmax % timestep loop
    n = n + 1;
    tstep = tstep + 1;
    if VERBOSE == 1
        fprintf('tstep = %g, t = %g, initial dt = %g\n',tstep,t,dt);
    end
    r = Inf;
    REDO = 0;
    while r > 2
        if REDO > 0 && VERBOSE == 1
           fprintf('\tREDO = %g\n',REDO); 
        end
        if 1.1*dt >= tf - t
            dt = tf - t;
        end
        tnew = t + dt;
        if VERBOSE == 1
            fprintf('    t = %g, r = %g, dt = %g\n',tnew,r,dt);
        end

        umid = TRstep(dt,h,N,u,Dplus,Dminus,CONST,VERBOSE);
        unew = BDF2step(dt,h,N,u,umid,Dplus,Dminus,CONST1,CONST2,CONST3,VERBOSE);

        norm_solution = norm(u,1)/(N+1);
        e_l = norm(2*k*dt*d(u).*D2*(u/GAMMA - umid/(GAMMA*(1-GAMMA)) + ...
            unew/(1-GAMMA)),1)/(N+1); % linearized divided-difference formula
        r = e_l/(EPSILON_REL*norm_solution + EPSILON_ABS);
        dt = min(min(2*dt,dt/r^(1/3)),dt_max);
        if dt <= dt_min
            warning('MATLAB:dt_min','dt = %e < dt_min at t = %e.\n',dt,t);
            n = nmax; r = 0; % to exit from while loops
        end 
        REDO = REDO + 1;
    end
    t = tnew;
    u = unew;
    if VERBOSE == 1
        fprintf('    new t = %g, new r = %g, next dt = %g\n',t,r,dt);
    end

    plot(x,u,'b-','LineWidth',2);
    set(gca,'fontsize',24);
    xlabel('x'); ylabel('u');
    axis([0 xmax 0 1.5]);
    getframe;
end
Qf = total_number(u,h);
fprintf('total number of particles at tf = %g: Qf = %g\n',t,Qf);
fprintf('    (Qf - Q0)/Q0 = %g\n',(Qf - Q0)/Q0);

toc;

end

function umid = TRstep(dt,h,N,u,Dplus,Dminus,CONST,VERBOSE)
EPSILON = 10^-9;

umid = u;
iter = 0;
R = umid - u - CONST*dt*((dplus(u).*Dplus - dminus(u).*Dminus)*u + ...
    (dplus(umid).*Dplus - dminus(umid).*Dminus)*umid)/h;
normresidual0 = norm(R,1)/(N+1);
if VERBOSE == 1
    fprintf('\t\tTR step\n');
    fprintf('\t||residual_0|| = %g\n',normresidual0);
end
norm_residual = Inf;    
while norm_residual > EPSILON*normresidual0 && iter < 8
    iter = iter + 1;
    J = speye(N+1) - ...
        CONST*dt*((dplus(umid).*Dplus - dminus(umid).*Dminus) + ...
        d_dplus(umid).*(Dplus*umid) - d_dminus(umid).*(Dminus*umid))/h;
    umid = umid - J\R;
    % umid = (speye(N+1) - CONST*dt*d_D2)\((speye(N+1) + CONST*dt*d_D2)*u);
    R = umid - ...
        (u + CONST*dt*((dplus(u).*Dplus - dminus(u).*Dminus)*u + ...
        (dplus(umid).*Dplus - dminus(umid).*Dminus)*umid)/h);
    norm_residual = norm(R,1)/(N+1);
    if VERBOSE == 1
        fprintf('\tNewton iteration %d: ',iter);
        fprintf('||residual||/||residual0|| = %g\n',...
            norm_residual/normresidual0);
    end
end
end
 
function unew = BDF2step(dt,h,N,u,umid,Dplus,Dminus,CONST1,CONST2,CONST3,VERBOSE)
EPSILON = 10^-9;

unew = umid;
iter = 0;
R = unew - (-CONST3*u + CONST2*umid  + ...
    CONST1*dt*((dplus(unew).*Dplus - dminus(unew).*Dminus)*unew)/h);
normresidual0 = norm(R,1)/(N+1);
if VERBOSE == 1
    fprintf('\t\tBDF2 step\n');
    fprintf('\t||residual_0|| = %g\n',normresidual0);
end
norm_residual = Inf;
while norm_residual > EPSILON*normresidual0 && iter < 8
    iter = iter + 1;
    J = speye(N+1) - CONST1*dt*((dplus(unew).*Dplus - ...
        dminus(unew).*Dminus) + d_dplus(unew).*(Dplus*unew) - ...
        d_dminus(unew).*(Dminus*unew))/h;
    unew = unew - J\R;
    % unew = (speye(N+1) - CONST1*dt*d_D2)\(CONST2*umid - CONST3*u);
    R = unew - (-CONST3*u + CONST2*umid  + ...
        CONST1*dt*((dplus(unew).*Dplus - dminus(unew).*Dminus)*unew)/h);
    norm_residual = norm(R,1)/(N+1);
    if VERBOSE == 1
        fprintf('\tNewton iteration %d: ',iter);
        fprintf('||residual||/|residual0|| = %g\n',...
            norm_residual/normresidual0);
    end
end
end

% diffusion coefficient d(u) = a u + b
% a = 1.5; b = 0.02;
% a = 0; b = 1.5; % for linear diffusion

function d = d(u)
a = 1.5; b = 0.02;
d = a*u + b;
end

function dplus = dplus(u)
N = size(u,1) - 1;
a = 1.5; b = 0.02;
e = ones(N+1,1);
aplus = a*spdiags([e e],[0 1],N+1,N+1)/2;
aplus(N+1,N) = a/2; % for Neumann BC at right
dplus = aplus*u + b;
end

function d_dplus = d_dplus(u)
% d dplus/du
N = size(u,1) - 1;
a = 1.5;
e = ones(N+1,1);
aplus = a*spdiags([e e],[0 1],N+1,N+1)/2;
aplus(N+1,N) = a/2; % for Neumann BC at right
d_dplus = aplus;
end

function dminus = dminus(u)
N = size(u,1) - 1;
a = 1.5; b = 0.02;
e = ones(N+1,1);
aminus = a*spdiags([e e],[-1 0],N+1,N+1)/2;
aminus(1,2) = a/2; % for Neumann BC at left
dminus = aminus*u + b;
end

function d_dminus = d_dminus(u)
% d dminus/du
N = size(u,1) - 1;
a = 1.5;
e = ones(N+1,1);
aminus = a*spdiags([e e],[-1 0],N+1,N+1)/2;
aminus(1,2) = a/2; % for Neumann BC at left
d_dminus = aminus;
end

function Q = total_number(u,h)
N = size(u,1) - 1;
Q = 0;
for i = 1:N
    Q = Q + (u(i) + u(i+1))/2;
end
Q = Q*h;
end







