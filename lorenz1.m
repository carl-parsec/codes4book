function lorenz1(tf)
% lorenz1(tf) VERSION 8-21-2023
% solves the Lorenz equations from t = 0 to tf with 
% x(0) = x0, y(0) = y0, z(0) = z0 given below,
% using fourth/fifth-order Runge-Kutta or TRBDF2.
% Here w(t) = [x(t), y(t), z(t)] and the Lorenz equations are 
%     dw/dt = F(w)
% In the plots, the red dots = the equilibria and the cyan dot = the ICs.
% The Jacobian for TRBDF2 is supplied in J.
% Try: lorenz1(40)

tic;

fprintf('tf = %g\n',tf);
tspan = [0 tf];

sigma = 10;
b = 8/3;
r = 28; % standard chaotic attractor, D = 2.06
% r = 24.737;% Moler problem
% r = 10;    % -> fixed point
% r = 99.65; % complicated periodic
% r = 100.5; % complicated periodic
% r = 160;   % complicated periodic
% r = 350;   % simplest periodic orbit

% equilibria
eta = sqrt(b*(r-1)); % or = -sqrt(b*(r-1))
xeq = eta;
yeq = eta;
zeq = r - 1;

w0 = [0; 1; 0]; % standard ICs
% w0 = [xeq; yeq + eps; zeq]; % near unstable equilibrium inside attractor
% w0 = [0; eps; 0]; % near very unstable equilibrium at origin

% Runge-Kutta 4/5
options = odeset('RelTol',10^-12,'AbsTol',10^-15);
[t,w] = ode45(@F,tspan,w0,options);

% TRBDF2 with analytical Jacobian
% options = odeset('RelTol',10^-6,'AbsTol',10^-9,'Jacobian',@jacobian);
% [t,w] = ode23tb(@F,tspan,w0,options);

% TRBDF2 with MATLAB calculated finite-difference Jacobian
% options = odeset('RelTol',10^-6,'AbsTol',10^-9);
% [t,w] = ode23tb(@F,tspan,w0,options);

toc;

N = size(t,1);
M = round(0.5*N); % to discard transients in attractor for first M steps

figure;
plot(t,w(:,2),'r-','LineWidth',2);
set(gca,'fontsize',24);
xlabel('t'); ylabel('y');
title('Lorenz equations');

figure;
plot(t,w(:,1),'b-','LineWidth',2);
set(gca,'fontsize',24);
xlabel('t'); ylabel('x');
title('Lorenz equations');

figure;
plot(t,w(:,3),'c-','LineWidth',2);
set(gca,'fontsize',24);
xlabel('t'); ylabel('z');
title('Lorenz equations');

figure;
plot(w(M:N,1),w(M:N,3),'b-',xeq,zeq,'r.',-xeq,zeq,'r.',...
    w(1,1),w(1,3),'c.','MarkerSize',24,'LineWidth',2);
set(gca,'fontsize',24);
xlabel('x'); ylabel('z');
title('Lorenz attractor');

figure;
plot3(w(1,1),w(1,3),w(1,2),'c.',...
    xeq,zeq,yeq,'r.',-xeq,zeq,-yeq,'r.',0,0,0,'r.',...
    w(M:N,1),w(M:N,3),w(M:N,2),'b-','MarkerSize',24,'LineWidth',2);
set(gca,'fontsize',24);
xlabel('x'); ylabel('z'); zlabel('y');
title('3D Lorenz attractor');

figure;
comet3(w(:,1),w(:,3),w(:,2));

    function wprime = F(~,w)
        wprime = [sigma*(w(2) - w(1))
                  r*w(1) - w(2) - w(1)*w(3)
                  w(1)*w(2) - b*w(3)];
    end

    function J = jacobian(~,w) % dF/dw
        J = [-sigma   sigma   0
             r-w(3)   -1      -w(1)
             w(2)     w(1)    -b   ];
    end

end









