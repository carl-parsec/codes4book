function weno3(rhol,ul,pl,rhor,ur,pr,gamma,tf,N,cfl,method)
%weno3(rhol,ul,pl,rhor,ur,pr,gamma,tf,N,cfl,method) VERSION 8-20-2023
% solves 1D gas dynamical Riemann Problems using 3rd-order WENO (method = 0)
% or positivity-preserving (PP) 3rd-order WENO (method = 1) with 
% Lax-Friedrichs flux splitting. 
% Try weno3(5,30,0.4127,0.05,0,0.4127,5/3,0.01,500,0.5,0 or 1)
% to see the PP method in action!
% BCs are through-flow.  The exact RP solution is superimposed.
%     rhol,ul,pl =  left RP parameters: density, velocity, pressure
%     rhor,ur,pr = right RP parameters: density, velocity, pressure
%     gamma = polytropic gas constant
%     tf = final time
%     N = number of dx (with N+1 grid points)
%     cfl = CFL factor (cfl <= 0.5 for positivity-preserving method)
%     method = 0 WENO3
%            = 1 positivity-preserving WENO3
% Examples:
%     1. weno3(5,30,0.4127,0.5,0,0.4127,5/3,0.015,500,0.5,0) Jet RP
%          shock/contact/shock
%          Ha & Gardner, J Sci Comp 34 (2008) 247
%        weno3(5,30,0.4127,0.05,0,0.4127,5/3,0.01,500,0.5,1) for PP
%     2. weno3(1,0,1,0.125,0,0.1,1.4,0.2,200,0.9,0) Sod RP
%          rarefaction/contact/shock
% From Liska & Wendroff, SIAM J Sci Comp 25 (2003) 995:
%     3. weno3(1,0.75,1,0.125,0,0.1,1.4,0.2,200,0.9,0) Toro/Sod RP
%          rarefaction/contact/shock
%     4. weno3(1,-19.59745,1000,1,-19.59745,0.01,1.4,0.01,200,0.9,0)
%          Stationary Contact RP (Toro)  
%          rarefaction/contact/shock
%     5. weno3(1,-2,0.4,1,2,0.4,1.4,0.1,400,0.5,0) Near-Vacuum RP (Toro)
%          rarefaction/contact(0 strength)/rarefaction
%     6. weno3(1,1,10^-6,1,-1,10^-6,5/3,1,200,0.75,0) Noh RP
%          shock/contact(0 strength)/shock
%     7. weno3(0.1261192,8.9047029,782.92899,6.591493,...
%          2.2654207,3.1544874,1.4,0.0055,800,0.9,0) Peak RP (Kucharik)
%          rarefaction/contact/shock
%     8. weno3(5.99924,19.5975,460.894,5.99242,-6.19633,46.095,1.4,...
%          0.02,400,0.9,0) Two Strong Shocks RP (Toro)
%          shock/contact/shock

% Developed by Jeremiah Jones and Carl Gardner, Arizona State University.
%     Email comments to carl.gardner@asu.edu
% WENO implementation based on a Fortran code of Guang-Shan Jiang and 
% Chi-Wang Shu.

% NOTES:
%     0. To speed the code up significantly, comment out the movie!
%     1. The MATLAB code is optimized for readability instead of speed or
%        memory usage.
%     2. 5th-order WENO gives better resolution of the RPs (good student
%        project!).
%     3. If transformation to characteristic variables is not used with
%        WENO fluxes, the Jet, Stationary Contact, and Noh RP examples
%        above fail, while the Peak RP example gives p and u overshoots
%        and oscillations (good HW problem!).  Search for comments
%        1/5-5/5 for "code: char variables".
%     4. The computed internal energy ~ p/rho is not as accurate as 
%        rho, u, p, especially when p and rho are both near 0.

% Reverse RP Examples for Testing Code:
%     weno3(0.125,0,0.1,1,0,1,1.4,0.2,200,0.9,0) reverse Sod RP
%     weno3(0.125,0,0.1,1,-0.75,1,1.4,0.2,200,0.9,0) reverse Toro/Sod RP
%     weno3(0.5,0,0.4127,5,-30,0.4127,5/3,0.015,500,0.5,0) reverse Jet RP
% Single Shock RP for Testing Code:
%     (There are some small start-up oscillations.)
%     weno3(2,1,3,1,0,1,1.5,0.2,400,0.9,0)
%     weno3(3.86,0,10.33,1,-2.63,1,1.4,0.2,200,0.9,0) % LeVeque p. 117 of
%       Nonlinear Conservation Laws and Finite Volume Methods for
%       Astrophysical Fluid Flow

% The Euler equations of 1D gas dynamics are written as conservation laws
%     q_t + f_x = 0
% which are then approximated by 
%     (dq/dt)(j) = -(F(j+1/2)-F(j-1/2))/dx
% Conserved quantities solution array
%     q(j,eq#) 
% is an (N+1)x3 matrix:
%     rho = q(:,1) = mass density
%      mx = q(:,2) = x momentum density
%       E = q(:,3) = energy density
% Pressure is denoted by p.

tic;

if method == 0
    fprintf('WENO3 method\n');
elseif method == 1
    fprintf('positivity-preserving WENO3 method\n');
    if cfl > 0.5 % cfl <= 1/2 for positivity-preserving method
        cfl = 0.5;
        fprintf('cfl reset to 0.5 for positivity-preserving method\n');
    end
else
    error('incorrect method = %g',method);
end

fprintf('(rhol,ul,pl) = (%g, %g, %g)\n',rhol,ul,pl);
fprintf('(rhor,ur,pr) = (%g, %g, %g)\n',rhor,ur,pr);
fprintf('gamma = %g\n',gamma);
fprintf('tf = %g, N = %g, CFL factor = %g\n',tf,N,cfl);

% if N is odd, reset N to N+1
if mod(N,2) 
    N = N + 1;
    fprintf('N reset = %g to be even\n',N);
end

dx = 1/N;
x = (-0.5:dx:0.5)'; % grid points

% allocate memory for solution array
q = zeros(N+1,3); % 3 conserved quantities for 1D gas dynamics

% left RP parameters
ml = rhol*ul;
El = pl/(gamma - 1) + 0.5*rhol*ul^2;

% right RP parameters
mr = rhor*ur;
Er = pr/(gamma - 1) + 0.5*rhor*ur^2;

% set up initial conditions
q(1:N/2,1) = rhol;
q(1:N/2,2) = ml;
q(1:N/2,3) = El;

% q(N/2+1,1) = rhol;
% q(N/2+1,2) = ml;
% q(N/2+1,3) = El;
% or
q(N/2+1,1) = (rhol + rhor)/2;
q(N/2+1,2) = (ml + mr)/2;
q(N/2+1,3) = (El + Er)/2;

q(N/2+2:N+1,1) = rhor;
q(N/2+2:N+1,2) = mr;
q(N/2+2:N+1,3) = Er;

% timestep loop
figure; % for movie
n = 0;
t = 0;
while t < tf
    n = n + 1;
    p = pressure(q,gamma);
    % compute dt 
    alpha = max_char_speed(q,p,gamma);
    dt = cfl*dx/alpha;
    if 1.01*dt >= tf - t
        dt = tf - t;
    end
    t = t + dt;
    lambda = dt/dx;
    
    % RK3 timestep for dq/dt = rhs = -df/dx
    q1 = q + dt*rhs(dx,q,p,gamma,alpha,lambda,method);
    p1 = pressure(q1,gamma);
    alpha1 = max_char_speed(q1,p1,gamma);
    q2 = 0.75*q + 0.25*(q1 + dt*rhs(dx,q1,p1,gamma,alpha1,lambda,method));
    p2 = pressure(q2,gamma);
    alpha2 = max_char_speed(q2,p2,gamma);
    q = (q + 2*(q2 + dt*rhs(dx,q2,p2,gamma,alpha2,lambda,method)))/3;
    
    % movie
    plot(x,q(:,1),'b-','LineWidth',2);
    xlabel('x'); ylabel('Density'); 
    set(gca,'fontsize',24);
    xticks([-0.5 -0.25 0 0.25 0.5]);
    getframe;
end
fprintf('number of timesteps = %g\n',n);

toc;

% exact RP solution [rho,v,p]:
qexact = exact_RP_sol(x,rhol,ul,pl,rhor,ur,pr,gamma,tf);

figure;
plot(x,qexact(:,1),'r-',x,q(:,1),'b-','LineWidth',2);
set(gca,'fontsize',24);
xticks([-0.5 -0.25 0 0.25 0.5]);
xlabel('x'); ylabel('Density'); 
% legend('exact','WENO','Location','North');

figure;
plot(x,qexact(:,2),'r-',x,q(:,2)./q(:,1),'b-','LineWidth',2);
set(gca,'fontsize',24);
xticks([-0.5 -0.25 0 0.25 0.5]);
xlabel('x'); ylabel('Velocity');  
% legend('exact','WENO','Location','South'); 

figure;
p = pressure(q,gamma);
plot(x,qexact(:,3)./((gamma - 1)*qexact(:,1)),'r-',...
    x,p./((gamma - 1)*q(:,1)),'b-','LineWidth',2);
set(gca,'fontsize',24);
xticks([-0.5 -0.25 0 0.25 0.5]);
xlabel('x'); ylabel('Internal Energy'); 
% legend('exact','WENO','Location','South');

figure;
p = pressure(q,gamma);
plot(x,qexact(:,3),'r-',x,p,'b-','LineWidth',2);
set(gca,'fontsize',24);
xticks([-0.5 -0.25 0 0.25 0.5]);
xlabel('x'); ylabel('Pressure'); 
% legend('exact','WENO','Location','SouthWest');

numerr = sum(abs(qexact(:,2) - q(:,2)./q(:,1)))/(N+1);
fprintf('|err|_1 in u = %g\n',numerr);

end

function f = flux(q,p)
% exact flux for 1D gas dynamics
u = q(:,2)./q(:,1);
f(:,1) = q(:,2);
f(:,2) = q(:,2).^2./q(:,1) + p; % or q(:,1).*u.^2 + p;
f(:,3) = u.*(q(:,3) + p);
end

function p = pressure(q,gamma)
p = (gamma - 1)*(q(:,3) - 0.5*q(:,2).^2./q(:,1));
end

function alpha = max_char_speed(q,p,gamma)
% computes maximum characteristic speed
if numel(find(p < 0)) > 0
    fprintf('try positivity-preserving method or decrease cfl and dx\n');
    error('pressure < 0 in max_char_speed');
end 
u = q(:,2)./q(:,1);
c = sqrt(gamma*p./q(:,1));
alpha = max(abs(u) + c);
end

function RHS = rhs(dx,q,p,gamma,alpha,lambda,method)
% alpha = max(abs(u) + c), lambda = dt/dx
N = size(q,1) - 1;
F = WENO_flux(q,p,gamma,alpha);
% F = LF_flux(q,p,alpha); % to implement LF simulations
if method == 1 % PP WENO3
    F_LF = LF_flux(q,p,alpha);
    Fstar = Flimited_flux(q,p,F,F_LF,gamma,lambda,1);
    F = Flimited_flux(q,p,Fstar,F_LF,gamma,lambda,2);
end
RHS = -(F(2:N+2,:) - F(1:N+1,:))/dx;  
end

function F_LF = LF_flux(q,p,alpha)
% LF method
N = size(q,1) - 1;
% through-flow BCs with ghost points
qghost = [q(1,:);q;q(N+1,:)];
pghost = [p(1);p;p(N+1)];

f = flux(qghost,pghost);
F_LF = 0.5*(f(1:N+2,:) + f(2:N+3,:) - alpha*(qghost(2:N+3,:) - ...
    qghost(1:N+2,:)));
end

function F = WENO_flux(q,p,gamma,alpha)
% WENO method using WENO3()
N = size(q,1) - 1;
% through-flow BCs with ghost points
qghost = [q(1,:);q(1,:);q;q(N+1,:);q(N+1,:)];
pghost = [p(1);p(1);p;p(N+1);p(N+1)];

f = flux(qghost,pghost);

% Lax-Friedrichs flux splitting
% there are N+5 elements of fplus and fminus
fplus = 0.5*(f + alpha*qghost); 
fminus = 0.5*(f - alpha*qghost);

% 1/5 compute left and right Roe matrices (code: char variables)
L = zeros(3,3,N+2);
R = zeros(3,3,N+2);
for j = 1:N+2
    [LL,RR] = roe_matrix(qghost(j+1:j+2,:),pghost(j+1:j+2),gamma);
    L(:,:,j) = LL;
    R(:,:,j) = RR;
end

% compute positive WENO flux 
                        % F(j+1/2) and F(j-1/2) stencils
fjm1 = fplus(1:N+2,:);  % f(j-1)       f(j-2)
fj = fplus(2:N+3,:);    % f(j)         f(j-1)
fjp1 = fplus(3:N+4,:);  % f(j+1)       f(j)

% 2/5 transform to characteristic variables for Fplus (code: char variables)
for j = 1:N+2
   fj(j,:) = fj(j,:)*L(:,:,j)';
   fjp1(j,:) = fjp1(j,:)*L(:,:,j)';
   fjm1(j,:) = fjm1(j,:)*L(:,:,j)';
end

Fplus = WENO3(fjm1,fj,fjp1);

% compute negative WENO flux
                        % F(j+1/2) and F(j-1/2) stencils
fj = fminus(2:N+3,:);   % f(j)         f(j-1)
fjp1 = fminus(3:N+4,:); % f(j+1)       f(j)
fjp2 = fminus(4:N+5,:); % f(j+2)       f(j+1)

% 3/5 transform to characteristic variables for Fminus (code: char variables)
for j = 1:N+2
   fjp2(j,:) = fjp2(j,:)*L(:,:,j)';
   fjp1(j,:) = fjp1(j,:)*L(:,:,j)';
   fj(j,:) = fj(j,:)*L(:,:,j)';
end

Fminus = WENO3(fjp2,fjp1,fj);

% compute WENO flux in characteristic variables
F = Fplus + Fminus;

% 4/5 transform WENO flux back to physical variables (code: char variables)
for j = 1:N+2
   F(j,:) = F(j,:)*R(:,:,j)';
end
end

% 5/5 (code: char variables)
function [L,R] = roe_matrix(q,p,gamma)
% See P. L. Roe, Journal of Computational Physics 135, 250-258 (1997)
% Originally published 1981
% q and p have only 2 rows j and j+1: 2x3 matrix
% computes left and right Roe matrices at x(j+1/2) based on Roe average of
% two states q(j) and q(j+1)                  

u = q(:,2)./q(:,1);
c = sqrt(gamma*p./q(:,1));
h = (q(:,3) + p)./q(:,1); % enthalpy
w = sqrt(q(:,1));

% mean values
t = w(1)/(w(1) + w(2));
um = t*u(1) + (1 - t)*u(2);
cm = t*c(1) + (1 - t)*c(2);
hm = t*h(1) + (1 - t)*h(2);

% define Roe matrices with mean values
R = [1, 1, 1; um - cm, um, um + cm; hm - cm*um, 0.5*um^2, hm + cm*um];
% L = R\eye(3);
b = (gamma - 1)/cm^2;
L = [0.5*(b*0.5*um^2 + um/cm), -0.5*(b*um + 1/cm), 0.5*b;
     1 - b*0.5*um^2, b*um, -b;
     0.5*(b*0.5*um^2 - um/cm), -0.5*(b*um - 1/cm), 0.5*b];
end

function F = WENO3(fjm1,fj,fjp1)
epsilon = 10^-6;
c1 = 1/3; 
c2 = 2/3; 

% compute WENO flux 
beta1 = (fj - fjm1).^2;        % left smoothness indicator
beta2 = (fj - fjp1).^2;        % right smoothness indicator
w1 = c1./(epsilon + beta1).^2; % unscaled left weight
w2 = c2./(epsilon + beta2).^2; % unscaled right weight 
sum = w1 + w2;                 % sum of unscaled weights
w1 = w1./sum;                  % scaled left weight
w2 = w2./sum;                  % scaled right weight
f1 = 0.5*(3*fj - fjm1);        % numerical flux on left stencil 
f2 = 0.5*(fj + fjp1);          % numerical flux on right stencil 
F = w1.*f1 + w2.*f2;           % WENO flux    
end

function Flimited = Flimited_flux(q,p,F,F_LF,gamma,lambda,type)
% Positivity-preserving method
% See X. Y. Hu, N. A. Adams, & C.-W. Shu, Journal of Computational Physics 
% 242 (2013) 169-180
N = size(q,1) - 1;
qplus = q - 2*lambda*F(2:N+2,:);
qminus = q + 2*lambda*F(1:N+1,:);
qplus_LF = q - 2*lambda*F_LF(2:N+2,:);
qminus_LF = q + 2*lambda*F_LF(1:N+1,:);
if type == 1
    epsilon = min(10^-13,min(q(:,1)));
    gplus = qplus(:,1);
    gminus = qminus(:,1);
    gplus_LF = qplus_LF(:,1);
    gminus_LF = qminus_LF(:,1);
elseif type == 2
    epsilon = min(10^-13,min(p));
    gplus = pressure(qplus,gamma);
    gminus = pressure(qminus,gamma);
    gplus_LF = pressure(qplus_LF,gamma);
    gminus_LF = pressure(qminus_LF,gamma);
else
    error('incorrect type = %g in Flimited_flux',type);
end
theta = positivity_limiter(gplus,gminus,gplus_LF,gminus_LF,epsilon);
Flimited = (1-theta).*F_LF + theta.*F;

% Ntheta = numel(find(theta < 1 - 10^-6));
% if Ntheta > 0
%     fprintf('theta < 1 - 10^-6 at %g points\n',Ntheta);
% end
end

function theta = positivity_limiter(gplus,gminus,gplus_LF,gminus_LF,epsilon)
N = size(gplus,1) - 1;
theta = ones(N+2,1);
gplus = [gplus(1);gplus;gplus(N+1)];
gminus = [gminus(1);gminus;gminus(N+1)];
gplus_LF = [gplus_LF(1);gplus_LF;gplus_LF(N+1)];
gminus_LF = [gminus_LF(1);gminus_LF;gminus_LF(N+1)];
for i=1:N+2
    thetaplus = 1; thetaminus = 1;
    if gplus(i) < epsilon
        thetaplus = (epsilon - gplus_LF(i))/(gplus(i) - gplus_LF(i));
    end
    if gminus(i+1) < epsilon
        thetaminus = (epsilon - gminus_LF(i+1))/...
            (gminus(i+1) - gminus_LF(i+1));
    end
    theta(i) = min(thetaplus,thetaminus);
    if theta(i) > 1
        fprintf('WARNING: theta = %g > 1, resetting to 1\n',theta(i));
        % error('theta > 1 computed in positivity_limiter');
        theta(i) = min(1,theta(i));
    elseif theta(i) < 0
        fprintf('WARNING: theta = %g < 0, resetting to 0\n',theta(i));
        % error('theta < 0 computed in positivity_limiter: decrease cfl');
        theta(i) = max(0,theta(i));
    end
end 
end

function q = exact_RP_sol(x,rhol,ul,pl,rhor,ur,pr,gamma,t)
% q = [rho,u,p](j)
SHOCK = 1;
RAREFACTION = 2;
EPSILON_WAVE = 10^-6;
N = size(x,1) - 1;
q = zeros(N+1,3);

[pstar,ustar,Ml,Mr,lwave,rwave] = ...
    exact_RP(rhol,ul,pl,rhor,ur,pr,gamma,SHOCK,RAREFACTION);
fprintf('pstar = %g, ustar = %g\n',pstar,ustar);

figure;
if lwave == SHOCK
    us = ul - Ml/rhol;
    rhostarII = Ml/(ustar - us);
    X = us*t;
    fprintf('left wave = SHOCK, shock speed = %g\n',us);
    if abs((pl - pstar)/(pl + pstar)) > EPSILON_WAVE
        line([0 X],[0 t],'color','r','LineWidth',2); % from (0,0) to (X,t)
    else
        fprintf('L SHOCK appears to be of 0 strength\n');
    end
elseif lwave == RAREFACTION
    cl = sqrt(gamma*pl/rhol);
    cstar = cl + 0.5*(ul - ustar)*(gamma - 1);
    rhostarII = gamma*pstar/cstar^2;
    % x1 = (ul - cl)*t; % left boundary of rarefaction
    % x2 = (ustar - cstar)*t; % right boundary of rarefaction
    s1 = ul - cl;
    s2 = ustar - cstar;
    fprintf('left wave = RAREFACTION, %g (left) = boundary speed = %g (right)\n',...
        s1,s2);
    if abs((pl - pstar)/(pl + pstar)) > EPSILON_WAVE
        XL1 = s1*t;
        XL2 = s2*t;
        hold on;
        xx = [0 XL1 XL2];
        yy = [0  t   t];
        fill(xx,yy,'c');
    else
        fprintf('L RAREFACTION appears to be of 0 strength\n');
    end
end
fprintf('CONTACT speed = %g\n',ustar);
if rwave == SHOCK
    us = ur + Mr/rhor;
    rhostarI = Mr/(us - ustar);
    X = us*t;
    fprintf('right wave = SHOCK, shock speed = %g\n',us);
    if abs((pr - pstar)/(pr + pstar)) > EPSILON_WAVE
        line([0 X],[0 t],'color','r','LineWidth',2); % from (0,0) to (X,t)
    else
        fprintf('R SHOCK appears to be of 0 strength\n');
    end
elseif rwave == RAREFACTION
    cr = sqrt(gamma*pr/rhor);
    cstar = cr + 0.5*(ustar - ur)*(gamma - 1);
    rhostarI = gamma*pstar/cstar^2;
    % x1 = (ustar + cstar)*t; % left boundary of rarefaction
    % x2 = (ur + cr)*t; % right boundary of rarefaction
    s1 = ur + cr;
    s2 = ustar + cstar;
    fprintf('right wave = RAREFACTION, %g (left) = boundary speed = %g (right)\n',...
        s2,s1);
    if abs((pr - pstar)/(pr + pstar)) > EPSILON_WAVE
        XR1 = s1*t;
        XR2 = s2*t;
        hold on;
        xx = [0 XR1 XR2];
        yy = [0  t   t];
        fill(xx,yy,'c');
    else
        fprintf('R RAREFACTION appears to be of 0 strength\n');
    end
end
X = ustar*t;
if abs((rhostarI - rhostarII)/(rhostarI + rhostarII)) > EPSILON_WAVE
    line([0 X],[0 t],'color','b','LineWidth',2); % from (0,0) to (X,t)
else
    fprintf('CONTACT appears to be of 0 strength\n');
end
fprintf('rhoI = %g, rhoII = %g\n',rhostarI,rhostarII);
set(gca,'fontsize',24);
axis([-0.5 0.5 0 t]);
xticks([-0.5 -0.25 0 0.25 0.5]);
xlabel('x'); ylabel('t'); 
title('Exact RP Solution');

xc = ustar*t; % location of contact
for j = 1:N+1
    if x(j) >= xc
        if rwave == SHOCK
            us = ur + Mr/rhor;
            if x(j) > us*t
                q(j,:) = [rhor,ur,pr];
            else
                rhostar = Mr/(us - ustar);
                q(j,:) = [rhostar,ustar,pstar];
            end
        else % rwave == RAREFACTION
            cr = sqrt(gamma*pr/rhor);
            cstar = cr + 0.5*(ustar - ur)*(gamma - 1);
            x1 = (ustar + cstar)*t; % left boundary of rarefaction
            x2 = (ur + cr)*t; % right boundary of rarefaction
            if x(j) < x1
                rhostar = gamma*pstar/cstar^2;
                q(j,:) = [rhostar,ustar,pstar];
            elseif x(j) > x2
                q(j,:) = [rhor,ur,pr];
            else % inside rarefaction
                uj = x(j)/t;
                cbar = (2*cr + (gamma - 1)*(uj - ur))/(1 + gamma);
                ubar = uj - cbar;
                rhobar = (cbar^2*rhor^gamma/(gamma*pr))^(1/(gamma - 1));
                pbar = cbar^2*rhobar/gamma;
                q(j,:) = [rhobar,ubar,pbar];
            end
        end
    else % x(j) < xc
        if lwave == SHOCK
            us = ul - Ml/rhol;
            if x(j) < us*t
                q(j,:) = [rhol,ul,pl];
            else
                rhostar = Ml/(ustar - us);
                q(j,:) = [rhostar,ustar,pstar];
            end
        else % lwave == RAREFACTION
            cl = sqrt(gamma*pl/rhol);
            cstar = cl + 0.5*(ul - ustar)*(gamma - 1);
            x1 = (ul - cl)*t; % left boundary of rarefaction
            x2 = (ustar - cstar)*t; % right boundary of rarefaction 
            if x(j) < x1
                 q(j,:) = [rhol,ul,pl];
            elseif x(j) > x2
                rhostar = gamma*pstar/cstar^2;
                q(j,:) = [rhostar,ustar,pstar];
            else % inside rarefaction
                uj = x(j)/t;
                cbar = (2*cl + (gamma - 1)*(ul - uj))/(1 + gamma);
                ubar = uj + cbar;
                rhobar = (cbar^2*rhol^gamma/(gamma*pl))^(1/(gamma - 1));
                pbar = cbar^2*rhobar/gamma;
                q(j,:) = [rhobar,ubar,pbar];
            end            
        end
    end
end
end

function [pstar,ustar,Ml,Mr,lwave,rwave] = ...
    exact_RP(rhol,ul,pl,rhor,ur,pr,gamma,SHOCK,RAREFACTION)
% based on Chorin, J Comp Phys 22 (1976) 517, modification of Godunov
MAX_ITER_alpha = 50;
MAX_ITER = 50;
EPS = 8*eps;
MIN_PRESSURE = realmin;

% initial guess
pstar = 0.5*(pl + pr);
coefl = sqrt(pl*rhol);
coefr = sqrt(pr*rhor);
Ml = coefl*phi(pstar/pl,gamma);
Mr = coefr*phi(pstar/pr,gamma);

% iterate at least twice to avoid spurious convergence when pl = pr    
pstar = ((ul - ur)*Ml*Mr + pr*Ml + pl*Mr)/(Ml + Mr);
if pstar < MIN_PRESSURE 
    pstar = MIN_PRESSURE;
end
Ml = coefl*phi(pstar/pl,gamma);
Mr = coefr*phi(pstar/pr,gamma);

% iterative loop using Godunov's method, as adapted by Chorin
n = 1;
alpha = 2;
for n1 = 1:MAX_ITER_alpha
alpha = alpha/2;
for n2 = 1:MAX_ITER
    n = n + 1;
    Ml_old = Ml;
    Mr_old = Mr;  
    ptilde = ((ul - ur)*Ml*Mr + pr*Ml + pl*Mr)/(Ml + Mr);
    if ptilde <= MIN_PRESSURE 
        ptilde = MIN_PRESSURE;
    end
    pstar = alpha*ptilde + (1 - alpha)*(pstar);
    Ml = coefl*phi(pstar/pl,gamma);
    Mr = coefr*phi(pstar/pr,gamma);
    if abs(Ml - Ml_old) <= EPS*Ml && abs(Mr - Mr_old) <= EPS*Mr
        ustar = (pl - pr + Mr*ur + Ml*ul)/(Ml + Mr);
        if pstar >= pl      
            lwave = SHOCK;
        else
            lwave = RAREFACTION;
        end
        if pstar >= pr       
            rwave = SHOCK;
        else
            rwave = RAREFACTION;
        end
        fprintf('exact RP solution converged in %g iterations\n',n);
        if pstar == MIN_PRESSURE
            fprintf('pstar = MIN_PRESSURE: likely error in exact RP\n');
        end
        return;
    end
end
end
error('exact RP solution failed to converge in %g iterations',n);
end

function PHI = phi(r,gamma)
EPS = 8*eps;
if abs(1 - r) <= EPS
    PHI = sqrt(gamma);
elseif r < 1
    beta = 0.5*(gamma - 1)/gamma;
    PHI = 0.5*(gamma - 1)*(1 - r)/(sqrt(gamma)*(1 - r^beta));
else % r > 1
    PHI = sqrt(0.5*(gamma + 1)*r + 0.5*(gamma - 1));
end
end









