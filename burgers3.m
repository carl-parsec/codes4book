function burgers3(N,steps,CFL,prob)
%burgers3(N,steps,CFL,method,prob) VERSION 8-24-2023
% solves the inviscid Burgers equation 
%     u_t + f_x = 0, f(u) = u^2/2
% with N+1 grid points for 0 <= x <= 1, h = 1/N, dt = CFL*h/umax, and
%     t_f = steps*dt.
% prob = 1 (shock), 2 (rarefaction), 3 (merging shocks), or 4 (N wave) 
% Uses WENO3-LF (about 47 lines vs 13 lines for conservative upwind) 
% and through-flow BCs.  The exact solution is superimposed.
% Try: burgers3(200,120,0.9,1) for Riemann problem shock
%      burgers3(200,100,0.9,2) for Riemann problem rarefaction
%      burgers3(200,180,0.9,3) for merging shocks
%      burgers3(200,90,0.9,4)  for N wave
%      burgers3(200,120,1,0)   for line

% WENO3 solver developed by Jeremiah Jones and Carl Gardner, Arizona State 
% University.
%     Email comments to carl.gardner@asu.edu
% WENO implementation based on a Fortran code of Guan-Shan Jiang and 
% Chi-Wang Shu.

% The inviscid Burgers equation is written as a conservation law
%     q_t + f_x = 0
% which is then approximated by 
%     (dq/dt)(j) = -(F(j+1/2)-F(j-1/2))/dx
% q is an (N+1)x1 column vector.

tic;

fprintf('WENO3 method\n');

% if N is odd, reset N to N+1
if mod(N,2) 
    N = N + 1;
    fprintf('N reset = %g to be even\n',N);
end

fprintf('N = %g, steps = %g, CFL = %g\n',N,steps,CFL);

dx = 1/N;
x = linspace(0,1,N+1); t = 0;
q = zeros(N+1,1); % N+1 dimensional column vector for u(x,t)

% ICs: 
if prob == 0
    fprintf('line\n');
    j = 1:N+1; q(j) = x(j) - 0.5;
elseif prob == 1
    fprintf('Riemann problem shock wave\n');
    j = 1:N/2; q(j) = 2; 
    j = N/2+1:N+1; q(j) = 1;
elseif prob == 2
    fprintf('Riemann problem rarefaction wave\n');
    j = 1:N/2+1; q(j) = 1;  
	j = N/2+2:N+1; q(j) = 2;    
elseif prob == 3
    fprintf('merging shocks\n');
    j = 1:N/4; q(j) = 3;
    j = N/4+1:N/2; q(j) = 2;
    j = N/2+1:N+1; q(j) = 1;
elseif prob == 4
    fprintf('N wave\n');
    j = N/4+1:3*N/4+1; q(j) = 4*(4*(j-N/4-1)/N - 1);
else
    error('incorrect problem type: %g',prob);
end

figure; % for movie
for n = 1:steps % timestep loop
    alpha = max(abs(q)); % max characteristic speed
    dt = CFL*dx/alpha;
    t = t + dt;

    % RK3 timestep for dq/dt = rhs = -df/dx
    q1 = q + dt*rhs(dx,q,alpha);
    alpha1 = max(abs(q1));
    q2 = 0.75*q + 0.25*(q1 + dt*rhs(dx,q1,alpha1));
    alpha2 = max(abs(q2));
    q = (q + 2*(q2 + dt*rhs(dx,q2,alpha2)))/3;

    % movie
    if prob == 0
        plot(x,qexact(t),'r.',x,q,'b-','LineWidth',2,'MarkerSize',12);
        axis([0 1 -0.5 0.5]);
    elseif prob == 3 % merging shocks
        plot(x,qexact(t),'r-',x,q,'b-','LineWidth',2);
        axis([0 1 0 3.2]);
    elseif prob == 4 % N wave
        plot(x,qexact(t),'r-',x,q,'b-','LineWidth',2);
        axis([0 1 -4 4]);
    else % RP
        plot(x,qexact(t),'r-',x,q,'b-','LineWidth',2);
        axis([0 1 0 2.2]); % RP
    end
    % legend('exact','WENO3','Location','South');
    set(gca,'fontsize',24);
    xlabel('x'); ylabel('u');
    getframe;
end
fprintf('tf = %g\n',t);

toc;

numerr = sum(abs(qexact(t) - q))/(N+1);
fprintf('|err|_1 = %g\n',numerr);

    function qe = qexact(t)
        qe = zeros(N+1,1);
        if prob == 0
            % line
            for i = 1:N+1
                qe(i) = (x(i)-0.5)/(1+t);
            end
        elseif prob == 1
            % shock wave
            s = 1.5;
            Ns = round(N/2 + s*t*N);
            Ns = min(N+1,Ns);
            i = 1:Ns; qe(i) = 2;
            i = Ns+1:N+1; qe(i) = 1;
        elseif prob == 2
            % rarefaction wave
            xleft = 0.5 + t;
            xright  = 0.5 + 2*t;
            for i = 1:N+1
                if (x(i) <= xleft)
                    qe(i) = 1;
                elseif (x(i) >= xright)
                    qe(i) = 2;
                elseif t > 0
                    qe(i) = (x(i)-0.5)/t;
                end
            end
        elseif prob == 3
            % ICs above were:
            %     j = 1:N/4; u(j) = 3;
            %     j = N/4+1:N/2; u(j) = 2;
            %     j = N/2+1:N+1; u(j) = 1;
            % merging shocks
            s1 = 2.5;
            Ns1 = round(N/4 + s1*t*N);
            Ns1 = min(N+1,Ns1);
            s2 = 1.5;
            Ns2 = round(N/2 + s2*t*N);
            Ns2 = min(N+1,Ns2);
            if Ns1 < Ns2
                i = 1:Ns1; qe(i) = 3;
                i = Ns1+1:Ns2; qe(i) = 2;
                i = Ns2+1:N+1; qe(i) = 1;
            else
                tm = 1/(4*(s1-s2));
                Nm = round(N/4 + s1*tm*N);
                s3 = 2;
                Ns3 = Nm + round(s3*(t-tm)*N);
                Ns3 = min(N+1,Ns3);
                i = 1:Ns3; qe(i) = 3;
                i = Ns3+1:N+1; qe(i) = 1;
            end
        elseif prob == 4
            % N wave
            t = t + 1/16;
            xleft = 0.5 - sqrt(t);
            xright  = 0.5 + sqrt(t);
            for i = 1:N+1
                if (x(i) < xleft)
                    qe(i) = 0;
                elseif (x(i) > xright)
                    qe(i) = 0;
                else
                    qe(i) = (x(i)-0.5)/t;
                end
            end
        end
    end

end

function f = flux(q)
% exact flux
f = q.^2/2;
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









