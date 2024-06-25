function chd(prob,V,N,cfl)
%chd(prob,V,N,cfl) VERSION 8-20-2023 Copyright by Carl L. Gardner 2023
% solves the 1D CHD model using 3rd-order WENO with Lax-Friedrichs flux
% splitting (using characteristic variables for fluxes), with TRBDF2 for 
% heat conduction and a tridiagonal solve for Poisson's equation.
%     prob = problem type
%     V = applied bias in volts
%     N = number of dx
%     cfl = CFL factor
%
% Problem types:  
%                channel     type       models  article
%   1  Si   300K 0.4  micron (subsonic) BW      Electron Shock Wave
%   2  Si   300K 0.1  micron (shock)    BW      Electron Shock Wave
%   3  Si    77K 1.0  micron (shock)    BW      Electron Shock Wave
%   4  Si    77K 1.0  micron (shock)    BWS     Monte Carlo
%   5  GaAs 300K 0.25 micron (shock)    CONST   NTK vs CLAWPACK
%
%                                   EPSILON = 10^-6
%                                   (2021 MacBook M1 Pro)
%                          timesteps  approx wall time (without movie)
% Try: chd(1,1,240,1)      1850         15 s
%      chd(2,1,240,1)      1686         15 s
%      chd(3,1,960,1)      5991        190 s
%      chd(4,1,960,1)      5823        180 s
%      chd(5,1,300,1)      3133         30 s

% Developed by Carl Gardner, SoMSS, Arizona State University.
%     Email comments to carl.gardner@asu.edu
% WENO3-LF code developed by Jeremiah Jones and Carl Gardner, SoMSS, ASU,
% based on a Fortran code of Guan-Shan Jiang and Chi-Wang Shu.

% The CHD equations are written as conservation laws
%     q_t + f_x = S
% which are then approximated by
%     (dq/dt)(j) = -(F(j+1/2)-F(j-1/2))/dx + S(j)
% Conserved quantities solution array
%     q(j,eq#) 
% is an (N+1)x3 matrix:
%     rho = q(:,1) = mass density = mn
%      mx = q(:,2) = momentum density = mnu
%       W = q(:,3) = energy density = 3nT/2 + mnu^2/2
% Pressure is denoted by p.
%      VP = VP(:) = -e phi = Poisson potential energy, e > 0
%       E = E(:) = dVP/dx = electric field

%                        Computational Units
%      Fundamental units
%              L = 0.1 micron = 10^-5 cm
%              energy = 1 eV
%              t = 10^-13 sec
%              N0 = 10^18 cm^-3
%      Derived units
%              u = 10^8 cm/sec
%              L^-3 = 10^15 cm^-3
%              mass = 10^-16 eV/(cm/sec)^2
%
%      m_e = 0.51099906 MeV/c_light^2 = 5.68563 10^-16 eV/(cm/sec)^2
%        where the speed of light c_light = 2.99792458 10^10 cm/sec
%      e = 4.803250 10^-10 esu (unrationalized)
%      eV = 1.60217733 10^-12 erg = 1.60217733 10^-19 J
%
%      epsilon_0 = 8.854 10^-12 F/m
%      e^2/epsilon_0 = 0.180955 10^-5 eV cm = 2.8816 10^-27 J m
%      epsilon = 11.7 for Si, 12.9 for GaAs

tic;

fprintf('WENO3/TRBDF2 method\n');

EPSILON = 10^-6; % 10^-6 to determine steady state
MAX_STEPS = 10^5; % or 10^6
EPSILON_SOR = 10^-12; % for SOR convergence if implemented
fprintf('EPSILON to determine steady state = %g\n',EPSILON);
% fprintf('EPSILON_SOR = %g\n',EPSILON_SOR);

gamma = 5/3;
e2 = 180.955; % (q_e^2/epsilon0)*(density scale factor of 1000)
m_e = 5.68563; % e- mass in 10^-16 eV/(cm/sec)^2

if prob == 1
    fprintf('\t Si 300K 0.4 micron channel (subsonic)\n');
    % Nd = 5*10^17 | 2*10^15 | 5*10^17 /cm^3
    models = 'BW';
    T0 = 0.025852; % eV
    Lchan = 4; L1 = 1; L2 = 1;
    xleft = 0; xright = Lchan; delta = 0.5;
    Nleft = 0.5; Nmid = 0.002; Nright = 0.5;
    m = 0.26*m_e;
    v_s = 0.1035;
    kappa0 = 0.4; % 0.4;
    epsilon = 11.7;
    mu_min = 0.080;
    delta_mu = 1.350;
    Nref = 0.112;
    aa = 0.72;
elseif prob == 2
    fprintf('\t Si 300K 0.1 micron channel (shock)\n');
    % Nd = 5*10^17 | 2*10^15 | 5*10^17 /cm^3
    models = 'BW';
    T0 = 0.025852; % eV
    Lchan = 1; L1 = 1; L2 = 1;
    % Lchan = 1; L1 = 2; L2 = 2;
    xleft = 0; xright = Lchan; delta = 0.125;
    Nleft = 0.5; Nmid = 0.002; Nright = 0.5;
    m = 0.26*m_e;
    v_s = 0.1035;
    kappa0 = 0.4;
    epsilon = 11.7;
    mu_min = 0.080;
    delta_mu = 1.350;
    Nref = 0.112;
    aa = 0.72;
elseif prob == 3
    fprintf('\t Si 77K 1 micron channel (shock)\n');
    % Nd = 10^18 | 10^15 | 10^18 /cm^3
    models = 'BW';
    T0 = 0.00663535; % eV
    Lchan = 10; L1 = 1; L2 = 1;
    xleft = 0; xright = Lchan; delta = 0.5;
    Nleft = 1; Nmid = 0.001; Nright = 1;
    m = 0.24*m_e;
    v_s = 0.1257;
    kappa0 = 0.4;
    epsilon = 11.7;
    mu_min = 0.080*(300/200)^0.45*(200/77)^0.15;
    delta_mu = 1.430*(300/77)^2 - 0.080;
    Nref = 0.112*(77/300)^3.2;
    aa = 0.72*(77/300)^0.065;
elseif prob == 4
    fprintf('\t Si 77K 1 micron channel (shock) MC\n');
    % Nd = 10^18 | 10^15 | 10^18 /cm^3
    models = 'BWS';
    T0 = 0.00663535; % eV
    Lchan = 10; L1 = 1; L2 = 1;
    xleft = 0; xright = Lchan; delta = 0.5; % 0.2 in MC
    Nleft = 1; Nmid = 0.001; Nright = 1;
    m = 0.24*m_e;
    v_s = 0.12;
    kappa0 = 0.05;
    epsilon = 11.7;
    tau_p0 = 16.7;
    Nref = 0.00114;
    aa = 0.659;
elseif prob == 5
    fprintf('\t GaAs 300K 0.25 micron channel (shock)\n');
    % Nd = 5*10^17 | 2*10^15 | 5*10^17 /cm^3
    models = 'CONST';
    T0 = 0.025852; % eV
    Lchan = 2.5; L1 = 2.5; L2 = 2.5;
    xleft = 0; xright = Lchan; delta = 0.5;
    Nleft = 0.5; Nmid = 0.002; Nright = 0.5;
    m = 0.063*m_e;
    v_s = 0.2;
    kappa0 = 0.4; % 0.4; 0 in NTK vs CLAWPACK
    epsilon = 12.9;
    tau_p0 = 2;
else
    error('incorrect problem type: %g',prob);
end

Ldevice = L1 + Lchan + L2;
dx = Ldevice/N;
x = (-L1:dx:(Lchan+L2))'; % grid points

% chd(prob,V,N,cfl)
fprintf('V = %g, T0 = %g eV, kappa0 = %g, N = %g, dx = %g, CFL = %g\n',...
    V,T0,kappa0,N,dx,cfl);
fprintf('source | channel | drain = %g | %g | %g = %g in 0.1 micron\n',...
    L1,Lchan,L2,Ldevice);
fprintf('Nd = %g | %g | %g in 10^18/cm^3\n',Nleft,Nmid,Nright);

e = ones(N+1,1);
Dplus = spdiags([-e e],[0 1],N+1,N+1);
Dminus = spdiags([-e e],[-1 0],N+1,N+1);

% doping Nd
Nd = zeros(N+1,1); 
for j = 1:N+1
    if x(j) <= xleft - delta
        Nd(j) = Nleft;
    elseif x(j) > xleft - delta && x(j) < xleft + delta
        Nd(j) = 0.5*(Nleft - Nmid)*tanh(-4*(x(j) - xleft)/delta) + ... 
            0.5*(Nleft + Nmid);
    elseif x(j) >= xleft + delta && x(j) <= xright - delta
        Nd(j) = Nmid;
    elseif x(j) > xright - delta && x(j) < xright + delta
        Nd(j) = 0.5*(Nmid - Nright)*tanh(-4*(x(j) - xright)/delta) + ...
            0.5*(Nmid + Nright);
    elseif x(j) >= xright + delta
        Nd(j) = Nright;
    end
end

figure;
plot(x,Nd,'b-','LineWidth',2,'MarkerSize',12);
xlim([-L1 (Lchan+L2)]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('N_D'); 

% tau and kappa models 
if strcmp(models,'BW')
    fprintf('BW models: ');
    mu_n0 = mu_min + delta_mu./(1 + (Nd/Nref).^aa);
    tau_p0 = m*mu_n0; % tau_p0 = 0.625 | 2 | 0.625 for 300K Si
    tau_w0 = tau_p0/2;
    kappa = kappa0*mu_n0*T0;
elseif strcmp(models,'BWS')
    fprintf('BWS models: ');
    kappa = ones(N+1,1)*kappa0*tau_p0*T0/m;
    tau_p0 = tau_p0./(1 + (Nd/Nref).^aa);
    tau_w0 = tau_p0/2;
elseif strcmp(models,'modBW')
    fprintf('modified BW models: ');
    tau_p0 = ones(N+1,1)*tau_p0;
    tau_w0 = tau_p0/2;
    kappa = kappa0*tau_p0*T0/m; 
elseif strcmp(models,'CONST')
    fprintf('CONST models: ');
    tau_p0 = ones(N+1,1)*tau_p0;
    tau_w0 = tau_p0;
    kappa = kappa0*tau_p0*T0/m; 
else
    error('incorrect models type'); 
end
% kappa will be multiplied by density n (incorporated into kappaD2 below)
fprintf('tau_p0 = %g | %g | %g in 0.1 picosec\n',...
    tau_p0(1),tau_p0(round(N/2)),tau_p0(N+1));

figure;
plot(x,tau_p0/m,'r-','LineWidth',2,'MarkerSize',12);
xlim([-L1 (Lchan+L2)]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('Mobility'); 

% allocate memory for solution array
q = zeros(N+1,3); % 3 conserved quantities for 1D gas dynamics

% set up initial conditions
q(:,1) = m*Nd;
q(:,2) = m*Nd(round(N/2))*v_s;
T = T0*ones(N+1,1);
q(:,3) = energy_density(q,T,m);
VP = -V*(x+L1)/Ldevice;

% initial Poisson solver
VP = Poisson(VP,Nd,q(:,1)/m,e2,epsilon,dx,EPSILON_SOR);
VPleft = [0;VP(1:N)];
VPright = [VP(2:N+1);-V];
E = (VPright - VPleft)/(2*dx);
E(1) = E(2); % dE/dx = 0 BCs
E(N+1) = E(N);

% timestep loop
% figure; % for movie
n = 0;
t = 0;
dqdt = Inf; % ||dq/dt||_1
while dqdt > EPSILON && n < MAX_STEPS
    qold = q;
    n = n + 1;
    p = pressure(q,gamma);
    % compute dt 
    alpha = max_char_speed(q,p,gamma);
    dt = cfl*dx/alpha;
    t = t + dt;
    
    % RK3 timestep for gas dynamics dq/dt = rhs = -df/dx + S
    q1 = q + dt*rhs(dx,q,p,gamma,alpha) + ...
        dt*source(q,tau_p0,tau_w0,T0,T,E,m,v_s,models);
    % BCs
    q1(1,1) = m*Nd(1); q1(N+1,1) = m*Nd(N+1);

    p = pressure(q1,gamma);
    alpha = max_char_speed(q1,p,gamma);
    T = temperature(q1,m);
    q2 = 0.75*q + 0.25*(q1 + dt*rhs(dx,q1,p,gamma,alpha) + ...
        dt*source(q1,tau_p0,tau_w0,T0,T,E,m,v_s,models));
    % BCs
    q2(1,1) = m*Nd(1); q2(N+1,1) = m*Nd(N+1);
    
    p = pressure(q2,gamma);
    alpha = max_char_speed(q2,p,gamma);
    T = temperature(q2,m);
    q = (q + 2*(q2 + dt*rhs(dx,q2,p,gamma,alpha) + ...
        dt*source(q2,tau_p0,tau_w0,T0,T,E,m,v_s,models)))/3;
    % BCs
    q(1,1) = m*Nd(1); q(N+1,1) = m*Nd(N+1);
    T = temperature(q,m);

    % heat conduction term
    % + heat conduction term: kappa = (kappa0*tau_p0*T0/m)*n
    % with Neumann BCs dT/dx = 0
    kappa_mid = 0.5*(kappa(1:N).*q(1:N,1)+kappa(2:N+1).*q(2:N+1,1));
    kappa_plus = [kappa_mid;kappa_mid(N)];
    kappa_minus = [kappa_mid(1);kappa_mid];
    kappaD2 = 2*(kappa_plus.*Dplus - kappa_minus.*Dminus)./...
        (3*q(:,1)*dx^2);
    % for dT/dx = 0 BCs
    kappaD2(1,2) = 2*kappaD2(1,2);
    kappaD2(N+1,N) = 2*kappaD2(N+1,N);
    T = heat_eq(T,dt,kappaD2);
    q(:,3) = energy_density(q,T,m);

    % Poisson solver
    VP = Poisson(VP,Nd,q(:,1)/m,e2,epsilon,dx,EPSILON_SOR);
    VPleft = [0;VP(1:N)];
    VPright = [VP(2:N+1);-V];
    E = (VPright - VPleft)/(2*dx);
    E(1) = E(2); % dE/dx = 0 BCs
    E(N+1) = E(N);
        
    dqdt = sum(sum(abs((q - qold)/dt)))/(3*N);

%     % movie
%     if prob == 3
%         plot(x,T/T0,'r-','LineWidth',2); title('T/T_0');
%         xlim([-L1 (Lchan+L2)]); ylim([0 20]);
%     elseif prob == 1 || prob == 2
%         plot(x,T/T0,'r-','LineWidth',2); title('T/T_0');
%         xlim([-L1 (Lchan+L2)]); ylim([0 10]);
%     else
%         plot(x,T/T0,'r-','LineWidth',2); title('T/T_0');
%         xlim([-L1 (Lchan+L2)]); ylim([0 12]);
%     end
%     %   OR
%     % plot(x,q(:,1)/m,'b-','LineWidth',2); title('Density');
%     % xlim([-L1 (Lchan+L2)]); ylim([0 1]);
%     %   OR
%     % plot(x,q(:,2)./q(:,1),'c-','LineWidth',2); title('Velocity');    
%     % xlim([-L1 (Lchan+L2)]); ylim([-0.5 2.5]);
%     
%     % PLUS
%     set(gca,'fontsize',24);
%     getframe;    
end

fprintf('final ||dq/dt|| = %g\n',dqdt);
fprintf('number of timesteps = %g\n',n);
fprintf('tf = %g\n',t);

toc;

figure;
plot(x,q(:,1)/m - Nd,'b-','LineWidth',2);
xlim([-L1 (Lchan+L2)]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('n - N_D'); 

figure;
plot(x,q(:,1)/m,'b-','LineWidth',2);
xlim([-L1 (Lchan+L2)]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('Density');

figure;
plot(x,log10(q(:,1)/m),'b-','LineWidth',2);
xlim([-L1 (Lchan+L2)]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('Log_{10}(Density)');

MAXq2 = max(q(:,2));
figure;
plot(x,q(:,2),'r-','LineWidth',2);
xlim([-L1 (Lchan+L2)]); ylim([0 1.2*MAXq2]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('Momentum Density');

u = q(:,2)./q(:,1);
figure;
plot(x,q(:,3),'b-','LineWidth',2);
xlim([-L1 (Lchan+L2)]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('Energy Density');

figure;
plot(x,u,'c-',x,u,'b.','LineWidth',2,'MarkerSize',12);
xlim([-L1 (Lchan+L2)]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('Velocity');

Mach = u./sqrt(T/m);
figure;
plot(x,Mach,'c-','LineWidth',2);
xlim([-L1 (Lchan+L2)]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('Mach Number');

figure;
plot(x,T/T0,'b-','LineWidth',2);
xlim([-L1 (Lchan+L2)]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('Temperature');

figure;
plot(x,VP,'b-','LineWidth',2);
xlim([-L1 (Lchan+L2)]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('Potential Energy');

figure;
plot(x,E,'c-','LineWidth',2);
xlim([-L1 (Lchan+L2)]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('Electric Field');

if dqdt > EPSILON
    error('steady state not reached: ||dq/dt|| = %g > EPSILON = %g\n',...
        dqdt,EPSILON);
end

end

function f = flux(q,p) 
% exact flux for 1D gas dynamics
u = q(:,2)./q(:,1);
f(:,1) = q(:,2);
f(:,2) = q(:,2).^2./q(:,1) + p; % or q(:,1).*v.^2 + p;
f(:,3) = u.*(q(:,3) + p);  
end

function p = pressure(q,gamma)
p = (gamma - 1)*(q(:,3) - 0.5*q(:,2).^2./q(:,1)); 
end

function T = temperature(q,m)
T = (2/3)*(q(:,3) - 0.5*q(:,2).^2./q(:,1))./(q(:,1)/m);
end

function W = energy_density(q,T,m)
W = 1.5*(q(:,1)/m).*T + 0.5*q(:,2).^2./q(:,1);
end

function alpha = max_char_speed(q,p,gamma)
% computes maximum characteristic speed
if (numel(find(p < 0)) > 0)
    fprintf('try positivity-preserving method or decrease cfl and dx\n');
    error('pressure < 0 in max_char_speed');
end 
u = q(:,2)./q(:,1);
c = sqrt(gamma*p./q(:,1));
alpha = max(abs(u) + c);
end

function RHS = rhs(dx,q,p,gamma,alpha)
N = size(q,1) - 1;
F = WENO_flux(q,p,gamma,alpha);
RHS = -(F(2:N+2,:) - F(1:N+1,:))/dx;  
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

% compute left and right Roe matrices
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

% transform to characteristic variables
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

% transform to characteristic variables
for j = 1:N+2
   fjp2(j,:) = fjp2(j,:)*L(:,:,j)';
   fjp1(j,:) = fjp1(j,:)*L(:,:,j)';
   fj(j,:) = fj(j,:)*L(:,:,j)';
end

Fminus = WENO3(fjp2,fjp1,fj);

% compute WENO flux in characteristic variables
F = Fplus + Fminus;

% transform WENO flux back to physical variables
for j = 1:N+2
   F(j,:) = F(j,:)*R(:,:,j)';
end
end

function [L,R] = roe_matrix(q,p,gamma)
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
w1 = c1./(epsilon + beta1).^2; % non-scaled left weight
w2 = c2./(epsilon + beta2).^2; % non-scaled right weight 
sum = w1 + w2;                 % sum of non-scaled weights
w1 = w1./sum;                  % scaled left weight
w2 = w2./sum;                  % scaled right weight
f1 = 0.5*(3*fj - fjm1);        % numerical flux on left stencil 
f2 = 0.5*(fj + fjp1);          % numerical flux on right stencil 
F = w1.*f1 + w2.*f2;           % WENO3 flux    
end

function S = source(q,tau_p0,tau_w0,T0,T,E,m,v_s,models)
N = size(q,1) - 1;

W0 = 1.5*(q(:,1)/m)*T0;
if strcmp(models,'BW')
    tau_p = tau_p0*T0./T;
    tau_w = (tau_w0*T0./T).*(1 + 3*T.^2./(m*v_s^2*(T + T0)));
elseif strcmp(models,'BWS') || strcmp(models,'modBW')
    tau_p = tau_p0*T0./T;
    tau_w = (tau_w0*T0./T).*(1 + 3*T./(m*v_s^2));
elseif strcmp(models,'CONST')
    tau_p = tau_p0;
    tau_w = tau_w0;
end

S(:,1) = zeros(N+1,1);
S(:,2) = -q(:,2)./tau_p - q(:,1).*E/m;
S(:,3) = -(q(:,3) - W0)./tau_w - q(:,2).*E/m; 
end

function VP = Poisson(VP,Nd,n,e2,epsilon,dx,~) % EPSILON_SOR
% VP = -e phi, where e > 0
N = size(VP,1) - 1;

% tridiagonal solve
b = zeros(N-1,1);
b(N-1) = -VP(N+1)/dx^2;
e = ones(N-1,1);
Dplus = spdiags([-e e],[0 1],N-1,N-1);
Dminus = spdiags([-e e],[-1 0],N-1,N-1);
D2 = (Dplus - Dminus)/dx^2;
VP(2:N) = D2\(e2*(Nd(2:N) - n(2:N))/epsilon + b);

% % SOR solver
% % calculate optimal omega for SOR
% mu_J = cos(pi/(N+1)); % Jacobi spectral radius
% omega = 2*(1-sqrt(1-mu_J^2))/mu_J^2;
% 
% normresidual = Inf;
% iter = 0;
% while normresidual > EPSILON_SOR && iter < 5000
%     iter = iter + 1;
%     sum = 0;
%     % interior points:
%     for j = 2:N
%         residual = -2*VP(j) + VP(j+1) + VP(j-1) - ...
%             e2*(Nd(j) - n(j))*dx^2/epsilon;
%         sum = sum + abs(residual);
%         VP(j) = VP(j) + omega*residual/2;
%     end
%     normresidual = sum/N;
% end
% if normresidual > EPSILON_SOR
%     fprintf('SOR failed to converge in %g iterations\n',iter);
%     error('SOR failed to converge: ||r|| = %g > EPSILON_SOR = %g\n',...
%         normresidual,EPSILON_SOR);  
% end
end

function T = heat_eq(T,dt,kappaD2)
N = size(T,1) - 1;
GAMMA = 2 - sqrt(2);
C = GAMMA/2;
C1 = (1-GAMMA)/(2-GAMMA);
C2 = 1/(GAMMA*(2-GAMMA));
C3 = (1-GAMMA)^2/(GAMMA*(2-GAMMA));

% with Neumann BCs dT/dx = 0 in kappaD2
% TRBDF2
Tmid = (speye(N+1)-C*dt*kappaD2)\((speye(N+1)+C*dt*kappaD2)*T); % TR
T = (speye(N+1)-C1*dt*kappaD2)\(C2*Tmid-C3*T); % BDF2
end









