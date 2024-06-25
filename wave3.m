function wave3(N,steps,CFL)
%wave3(N,steps,CFL,method) VERSION 8-24-2023
% solves the first-order wave equation
%     u_t + f_x = 0, f(u) = u
% with N+1 grid points for 0 <= x <= 1, h = 1/N, dt = CFL*h/umax, and
%     t_f = steps*dt.
% Uses WENO3-LF and through-flow BCs.  The exact solution is superimposed.
% Try: wave3(200,120,1)

% WENO3 solver developed by Jeremiah Jones and Carl Gardner, Arizona State
% University.
%     Email comments to carl.gardner@asu.edu
% WENO implementation based on a Fortran code of Guan-Shan Jiang and
% Chi-Wang Shu.

% The inviscid Burgers equation is written as a conservation law
%     q_t + f_x = 0
% which is then approximated by 
%     (dq/dt)(j) = -(F(j+1/2)-F(j-1/2))/dx
% q is an (N+1)x1 column vector

tic;

fprintf('WENO3 method\n');
fprintf('N = %g, steps = %g, CFL = %g\n',N,steps,CFL);

% if N is odd, reset N to N+1
if mod(N,2) 
    N = N + 1;
    fprintf('N reset = %g to be even\n',N);
end

dx = 1/N;
x = linspace(0,1,N+1); t = 0;
q = zeros(N+1,1); % N+1 dimensional column vector for u(x,t)
% j = 0:2*N/10;
j = 0:N/10;
q(j+1) = (sin(10*j*pi/N)).^4; % ICs

figure; % for movie
for n = 1:steps % timestep loop
    alpha = 1; % max characteristic speed
    dt = CFL*dx/alpha;
    t = t + dt;

    % RK3 timestep for dq/dt = rhs = -df/dx
    q1 = q + dt*rhs(dx,q,alpha);
    alpha1 = max(abs(q1));
    q2 = 0.75*q + 0.25*(q1 + dt*rhs(dx,q1,alpha1));
    alpha2 = max(abs(q2));
    q = (q + 2*(q2 + dt*rhs(dx,q2,alpha2)))/3;

    % movie
    plot(x,qexact(t),'r.',x,q,'b-','LineWidth',2,'MarkerSize',16);
    % plot(x,q,'b-','LineWidth',2);
    set(gca,'fontsize',24);
    axis([0 1 -0.5 1.1]);
    xlabel('x'); ylabel('u');
    getframe;
end
fprintf('tf = %g\n',t);

toc;

numerr = sum(abs(qexact(t) - q))/(N+1);
fprintf('|err|_1 = %g\n',numerr);

    function qe = qexact(t)
        qe = zeros(N+1,1);
        s = 1;
        Nleft = round(s*t*N);
        Nright = round(s*t*N + N/10);
        for i = 1:N
            if i < Nleft
                qe(i) = 0;
            elseif i > Nright
                qe(i) = 0;
            else
                qe(i+1) = sin(10*pi*(i - Nleft)/N).^4;
            end
        end
    end
end

function f = flux(q)
% exact flux
f = q;
end

function RHS = rhs(dx,q,alpha)
N = size(q,1) - 1;
F = WENO_flux(q,alpha);
RHS = -(F(2:N+2) - F(1:N+1))/dx;  
end

function F = WENO_flux(q,alpha)
% WENO method using WENO3()
N = size(q,1) - 1;
% through-flow BCs with ghost points
qghost = [q(1,:);q(1,:);q;q(N+1,:);q(N+1,:)];

f = flux(qghost);

% Lax-Friedrichs flux splitting
% there are N+5 elements of fplus and fminus
fplus = 0.5*(f + alpha*qghost); 
fminus = 0.5*(f - alpha*qghost);

% compute positive WENO flux 
                        % F(j+1/2) and F(j-1/2) stencils
fjm1 = fplus(1:N+2,:);  % f(j-1)       f(j-2)
fj = fplus(2:N+3,:);    % f(j)         f(j-1)
fjp1 = fplus(3:N+4,:);  % f(j+1)       f(j)

Fplus = WENO3(fjm1,fj,fjp1);

% compute negative WENO flux
                        % F(j+1/2) and F(j-1/2) stencils
fj = fminus(2:N+3,:);   % f(j)         f(j-1)
fjp1 = fminus(3:N+4,:); % f(j+1)       f(j)
fjp2 = fminus(4:N+5,:); % f(j+2)       f(j+1)

Fminus = WENO3(fjp2,fjp1,fj);

% compute WENO flux
F = Fplus + Fminus;
end

function F = WENO3(fjm1,fj,fjp1)
epsilon = 10^-6;
c1 = 1/3; 
c2 = 2/3; 

% compute WENO flux 
beta1 = (fj - fjm1).^2;        % left smoothness indicator
beta2 = (fj - fjp1).^2;        % right smoothness indicator
w1 = c1./(epsilon + beta1).^2; % non-scaled left weight
w2 = c2./(epsilon + beta2).^2; % non-scaled right weight 
sum = w1 + w2;                 % sum of non-scaled weights
w1 = w1./sum;                  % scaled left weight
w2 = w2./sum;                  % scaled right weight
f1 = 0.5*(3*fj - fjm1);        % numerical flux on left stencil 
f2 = 0.5*(fj + fjp1);          % numerical flux on right stencil 
F = w1.*f1 + w2.*f2;           % WENO flux    
end









