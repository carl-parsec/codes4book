function jets(rhojet,ujet,vjet,Tjet,rhoamb,uamb,vamb,Tamb,gamma, ...
    xmax,ymax,jetwidth,Nx,Ny,tf,cfl)
%jets(rhojet,ujet,vjet,Tjet,rhoamb,uamb,vamb,Tamb,gamma,...
%   xmax,ymax,jetwidth,Nx,Ny,tf,cfl) VERSION 10-31-2024 
% Copyright by Carl L. Gardner 2023-2024
% simulates 2D gas dynamical jets using 3rd-order WENO with Lax-Friedrichs
% flux splitting.
%       rho in H/cm^3, velocity = (u,v) in km/s, T = temperature in K
%           Note: 2.652 x 10^19 molecules/cm^3 in air at STP
%       rhojet,ujet,vjet,TjeT = jet parameters 
%       rhoamb,uamb,vamb,Tamb = ambient parameters
%       gamma = polytropic gas constant
%       xmax, ymax, jetwidth = geometry of jet grid in 10^11 km  
%       Nx = number of dx (with Nx+1 grid points)
%       Ny = number of dy (with Ny+1 grid points)
%       tf = final time in yr
%       cfl = CFL factor
%
% Try: jets(500,300,0,10^3,50,0,0,10^4,5/3,1,1,0.05,200,200,11,0.25)
%           heavy Mach 80 astrophysical jet
%        jets(500,10,0,200,50,0,0,10^4,5/3,4,2,0.1,400,200,1750,0.9)
%           heavy Mach 6 astrophysical jet with jet p < amb p (shock
%           diamonds)
%           
% Developed by Jeremiah Jones and Carl Gardner, Arizona State University
%     Email comments to carl.gardner@asu.edu
% WENO implementation based on a Fortran code of Guan-Shan Jiang and
% Chi-Wang Shu

% NOTES:
%     0. To speed the code up significantly, comment out the movie!
%     1. The MATLAB code is optimized for readability instead of speed or
%        memory usage.
%     2. 5th-order WENO gives better resolution (good student project!).
%     3. Does NOT use positivity-preserving scheme (another good student 
%        project!).
%     4. If characteristic variables are not used, maximum Mach number for
%        Mach 80 jet is about 7

% The Euler equations of 2D gas dynamics are written as conservation laws
%     q_t + f_x + g_y = 0
% which are then approximated by (i labels x and j labels y)
%     (dq/dt)(i,j) = -(F(i+1/2,j)-F(i-1/2,j))/dx - (G(i,j+1/2)-G(i,j-1/2))/dy
% Conserved quantities solution array
%     q(i,j,eq#) 
% is an (Nx+1)x(Ny+1)x4 matrix:
%     rho = q(:,:,1) = mass density
%      mx = q(:,:,2) = x momentum density = rho*u
%      my = q(:,:,3) = y momentum density = rho*v
%       E = q(:,:,4) = energy density

tic;

fprintf('2D WENO3 method\n');

% UNITS for astrophysical stellar jets
lbar = 10^11; % km
tbar = 10^10; % s (1 yr ~ pi*10^7 s)
rhobar = 100; % H/cm^3
nbar = 100;   % H/cm^3 since mH = 1 in our units
% Derived Units
ubar = lbar/tbar; % 10 km/s
% mH = 1 in computational units
Ebar = 104.454; % eV/cm^3
K_per_eV = 11604.519; % 11,604.519 K = 1 eV (NIST)
Tbar_eV = Ebar/nbar; % 1.04454 eV
Tbar_K = Tbar_eV*K_per_eV; % 12,121.384 K

X = 1; % X SWEEP
Y = 2; % Y SWEEP

fprintf('jet: rho = %g H/cm^3, (u,v) = (%g,%g) km/s, T = %g K\n',...
    rhojet,ujet,vjet,Tjet);
fprintf('ambient: rho = %g H/cm^3, (u,v) = (%g,%g) km/s, T = %g K\n',...
    rhoamb,uamb,vamb,Tamb);
fprintf('gamma = %g\n',gamma);
fprintf('x = [0,%g] x 10^11 km, y = [0,%g] x 10^11 km\n',xmax,ymax);
fprintf('jetwidth = %g x 10^11 km\n',jetwidth);
fprintf('Nx = %g, Ny = %g\n',Nx,Ny);
fprintf('tf = %g yr\n',tf);
fprintf('CFL factor = %g\n',cfl);

% SWITCH TO COMPUTATIONAL UNITS:
rhojet = rhojet/rhobar;
ujet = ujet/ubar;
vjet = vjet/ubar;
Tjet = Tjet/Tbar_K;
rhoamb = rhoamb/rhobar;
uamb = uamb/ubar;
vamb = vamb/ubar;
Tamb = Tamb/Tbar_K;
tf = tf*(3600*24*365.242199)/tbar;
fprintf('\tCOMPUTATIONAL UNITS\n');
fprintf('\tjet: rho = %g, (u,v) = (%g,%g), T = %g\n',rhojet,ujet,vjet,Tjet);
fprintf('\tambient: rho = %g, (u,v) = (%g,%g), T = %g\n',rhoamb,uamb,vamb,Tamb);
fprintf('\ttf = %g\n',tf);

dx = xmax/Nx;
dy = ymax/Ny;
fprintf('\tdx = %g, dy = %g\n',dx,dy);
x = (0:dx:xmax)'; y = (0:dy:ymax)'; % grid points

% allocate memory for solution array
q = zeros(Nx+1,Ny+1,4); % 4 conserved quantities for 2D gas dynamics

% ambient parameters
mxamb = rhoamb*uamb;
myamb = rhoamb*vamb;
pamb = rhoamb*Tamb;
Eamb =  pamb/(gamma - 1) + 0.5*rhoamb*(uamb^2 + vamb^2);

% jet parameters
fprintf('\tjetwidth = %g\n',jetwidth);
yjet1 = ymax/2 - jetwidth/2;
yjet2 = ymax/2 + jetwidth/2;
fprintf('\t[yjet1,yjet2] = [%g,%g]\n',yjet1,yjet2);
j1 = ceil(yjet1/dy) + 1;
j2 = floor(yjet2/dy) + 1;
fprintf('\tj1 = %g, j2 = %g\n',j1,j2);
if abs(y(j1) - yjet1) > 10^-15 || abs(y(j2) - yjet2) > 10^-15
    fprintf('\tadjusted [yjet1,yjet2] = [%g,%g] to hit grid points\n',...
        y(j1),y(j2));
end
mxjet = rhojet*ujet;
myjet = rhojet*vjet;
pjet = rhojet*Tjet;
Ejet = pjet/(gamma - 1) + 0.5*rhojet*(ujet^2 + vjet^2);

fprintf('\tP_jet = %g, P_amb = %g\n',pjet,pamb);

% set up initial conditions
q(1:Nx+1,1:Ny+1,1) = rhoamb;
q(1:Nx+1,1:Ny+1,2) = mxamb;
q(1:Nx+1,1:Ny+1,3) = myamb;
q(1:Nx+1,1:Ny+1,4) = Eamb;

q = enforce_jet_BCs(j1,j2,q,rhojet,mxjet,myjet,Ejet);
p = pressure(q,gamma);
camb = sqrt(gamma*pamb/rhoamb);
cjet = sqrt(gamma*pjet/rhojet);
fprintf('c in jet = %g km/s, c in ambient = %g km/s\n',cjet*ubar,camb*ubar);
% fprintf('max(c) = %g, min(c) = %g\n',cmax,cmin);
fprintf('M_jet wrt jet = %g, M_jet wrt ambient = %g\n',ujet/cjet,ujet/camb);

% timestep loop
figure; % for movie
n = 0;
t = 0;
while t < tf
    n = n + 1;
    % compute dt
    [alphax,alphay] = max_char_speed(q,p,gamma);
    dt = cfl/(alphax/dx + alphay/dy);
    if 1.01*dt >= tf - t
        dt = tf - t;
    end
    t = t + dt;

    % RK3 timestep for dq/dt = -(df/dx + dg/dy)
    q1 = q + dt*rhs(X,Nx,Ny,dx,q,p,gamma,alphax) + ...
        + dt*rhs(Y,Nx,Ny,dy,q,p,gamma,alphay);
    q1 = enforce_jet_BCs(j1,j2,q1,rhojet,mxjet,myjet,Ejet);
    p1 = pressure(q1,gamma);
    q2 = 0.75*q + 0.25*(q1 + dt*rhs(X,Nx,Ny,dx,q1,p1,gamma,alphax) + ...
        dt*rhs(Y,Nx,Ny,dy,q1,p1,gamma,alphay));
    q2 = enforce_jet_BCs(j1,j2,q2,rhojet,mxjet,myjet,Ejet);
    p2 = pressure(q2,gamma);
    q = (q + 2*(q2 + dt*rhs(X,Nx,Ny,dx,q2,p2,gamma,alphax) + ...
        + dt*rhs(Y,Nx,Ny,dy,q2,p2,gamma,alphay)))/3;
    q = enforce_jet_BCs(j1,j2,q,rhojet,mxjet,myjet,Ejet);
    p = pressure(q,gamma);

    % movie
    hh = pcolor(x,y,log10(q(:,:,1)*rhobar)');
    set(hh,'edgecolor','none','facecolor','interp');
    axis equal;
    axis tight;
    set(gca,'fontsize',24);
    title('Density'); % Log_{10}(n/(H/cm^3))
    colorbar;
    colormap(jet);
    getframe;

    tyr = t*tbar/(3600*24*365.242199);
    if mod(n,100) == 0
        fprintf('t = %g yr, number of timesteps = %g\n',tyr,n);
    end
end
fprintf('tf = %g yr, number of timesteps = %g\n',tyr,n);

toc;

figure;
hh = pcolor(x,y,log10(q(:,:,1)*rhobar)');
set(hh,'edgecolor','none','facecolor','interp');
axis equal;
axis tight;
set(gca,'fontsize',24);
title('Density'); % Log_{10}(n/(H/cm^3))
colorbar;
colormap(jet);

figure;
hh = pcolor(x,y,log10(p)');
set(hh,'edgecolor','none','facecolor','interp');
axis equal;
axis tight;
set(gca,'fontsize',24);
title('Pressure');
colorbar;
colormap(jet);

T = p./q(:,:,1);
figure;
hh = pcolor(x,y,log10(T*Tbar_K)'); 
set(hh,'edgecolor','none','facecolor','interp');
axis equal;
axis tight;
set(gca,'fontsize',24);
title('Temperature'); % Log_{10}(T/K)
colorbar;
colormap(jet);

end

function q = enforce_jet_BCs(j1,j2,q,rhojet,mxjet,myjet,Ejet)
% enforce jet inflow BCs
j = j1:j2;
q(1,j,1) = rhojet;
q(1,j,2) = mxjet;
q(1,j,3) = myjet;
q(1,j,4) = Ejet;
end

function f = f_flux_1d(q,p)
% exact flux for X sweep for 2D gas dynamics
u = q(:,2)./q(:,1);
v = q(:,3)./q(:,1);
f(:,1) = q(:,2);
f(:,2) = q(:,1).*u.^2 + p;
f(:,3) = q(:,1).*u.*v;
f(:,4) = u.*(q(:,4) + p);
end

function g = g_flux_1d(q,p)
% exact flux for Y sweep for 2D gas dynamics
u = q(:,2)./q(:,1);
v = q(:,3)./q(:,1);
g(:,1) = q(:,3);
g(:,2) = q(:,1).*u.*v;
g(:,3) = q(:,1).*v.^2 + p;
g(:,4) = v.*(q(:,4) + p);
end

function p = pressure(q,gamma)
p = (gamma - 1)*(q(:,:,4) - 0.5*(q(:,:,2).^2 + q(:,:,3).^2)./q(:,:,1));
end

function [alphax,alphay] = max_char_speed(q,p,gamma)
% computes maximum characteristic speed in x and y directions
if (numel(find(p < 0)) > 0)
    fprintf('try positivity-preserving method or decrease cfl and dx\n');
    error('pressure < 0 in max_char_speed');
end 
u = q(:,:,2)./q(:,:,1);
v = q(:,:,3)./q(:,:,1);
c = sqrt(gamma*p./q(:,:,1));
alphax = max(max(abs(u) + c));
alphay = max(max(abs(v) + c));
end

function RHS = rhs(SWEEP,Nx,Ny,h,q,p,gamma,alpha)
X = 1; % X SWEEP
Y = 2; % Y SWEEP

RHS = zeros(Nx+1,Ny+1,4);
if SWEEP == X
    q1d = zeros(Nx+1,4);
    p1d = zeros(Nx+1,1);
    for j = 1:Ny+1
        for i = 1:Nx+1
            q1d(i,:) = q(i,j,:);
            p1d(i) = p(i,j);
        end
        F = WENO_flux_1d(X,Nx,q1d,p1d,gamma,alpha);
        RHS(1:Nx+1,j,:) = -(F(2:Nx+2,:) - F(1:Nx+1,:))/h;
    end
elseif SWEEP == Y
    q1d = zeros(Ny+1,4);
    p1d = zeros(Ny+1,1);
    for i = 1:Nx+1
        for j = 1:Ny+1
            q1d(j,:) = q(i,j,:);
            p1d(j) = p(i,j);
        end
        G = WENO_flux_1d(Y,Ny,q1d,p1d,gamma,alpha);
        RHS(i,1:Ny+1,:) = -(G(2:Ny+2,:) - G(1:Ny+1,:))/h;
    end
end
end

function F = WENO_flux_1d(SWEEP,N,q,p,gamma,alpha)
% F = WENO_flux(SWEEP,N,q,p,~,alpha) % gamma
% WENO method using WENO3()
% through-flow BCs with ghost points
X = 1; % X SWEEP
Y = 2; % Y SWEEP

qghost = [q(1,:);q(1,:);q;q(N+1,:);q(N+1,:)];
pghost = [p(1);p(1);p;p(N+1);p(N+1)];
if SWEEP == X
    f = f_flux_1d(qghost,pghost);
elseif SWEEP == Y
    f = g_flux_1d(qghost,pghost);
end

% Lax-Friedrichs flux splitting
% there are N+5 elements of fplus and fminus
fplus = 0.5*(f + alpha*qghost); 
fminus = 0.5*(f - alpha*qghost);

% 1/5 compute left and right Roe matrices (code: char variables)
L = zeros(4,4,N+2);
R = zeros(4,4,N+2);
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

% compute WENO flux
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
% q and p have only 2 rows j and j+1: q is a 2x4 matrix
% computes left and right Roe matrices at x(j+1/2) based on Roe average of
% two states q(j) and q(j+1)                  

u = q(:,2)./q(:,1);
v = q(:,3)./q(:,1);
c = sqrt(gamma*p./q(:,1));
h = (q(:,4) + p)./q(:,1); % enthalpy
w = sqrt(q(:,1));

% mean values
t = w(1)/(w(1) + w(2));
um = t*u(1) + (1 - t)*u(2);
vm = t*v(1) + (1 - t)*v(2);
cm = t*c(1) + (1 - t)*c(2);
hm = t*h(1) + (1 - t)*h(2);

% define Roe matrices with mean values
R = [1, 0, 1, 1; 
    um - cm, 0, um, um + cm;
    vm, 1, vm, vm; 
    hm - cm*um, vm, 0.5*(um^2+vm^2), hm + cm*um];
% L = R\eye(4);
b = (gamma - 1)/cm^2;
L = [0.5*(b*0.5*(um^2+vm^2) + um/cm), -0.5*(b*um + 1/cm), -0.5*b*vm, 0.5*b;
    -vm, 0, 1, 0;
    1 - b*0.5*(um^2+vm^2), b*um, b*vm, -b;
    0.5*(b*0.5*(um^2+vm^2) - um/cm), -0.5*(b*um - 1/cm), -0.5*b*vm, 0.5*b];
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









