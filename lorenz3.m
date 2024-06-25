function lorenz3(tf)
%lorenz3(tf) VERSION 8-21-2023
% solves the Lorenz equations from t = 0 to tf with 3 different methods
% or ICs labeled 1, 2, 3.
%     w(t) = [x(t); y(t); z(t)] 
% and the Lorenz equations are
%     dw/dt = F(w)
% Try: lorenz3(40)

tic;

fprintf('tf = %g\n',tf);
tspan = [0 tf];
w0 = [0; 1; 0];

sigma = 10;
b = 8/3;
r = 28; % r = 28 standard chaotic attractor
% equilibria
eta = sqrt(b*(r-1)); % or = -sqrt(b*(r-1))
xeq = eta;
yeq = eta;
zeq = r - 1;

% Runge-Kutta 4/5
options = odeset('RelTol',10^-6,'AbsTol',10^-9);
[t1,w1] = ode45(@F,tspan,w0,options);
% TRBDF2 with analytical Jacobian
options = odeset('RelTol',10^-6,'AbsTol',10^-9,'Jacobian',@jacobian);
[t2,w2] = ode23tb(@F,tspan,w0,options);
% TRBDF2 with MATLAB calculated "finite-difference" Jacobian
options = odeset('RelTol',10^-6,'AbsTol',10^-9);
[t3,w3] = ode23tb(@F,tspan,w0,options);

% % Runge-Kutta 4/5 with different RelTol
% options = odeset('RelTol',10^-6,'AbsTol',10^-15);
% [t1,w1] = ode45(@F,tspan,w0,options);
% options = odeset('RelTol',10^-9,'AbsTol',10^-15);
% [t2,w2] = ode45(@F,tspan,w0,options);
% options = odeset('RelTol',10^-12,'AbsTol',10^-15);
% [t3,w3] = ode45(@F,tspan,w0,options);

% % Runge-Kutta 4/5 with different ICs
% w0_1 = [0; 1; 0];
% w0_2 = [1; 4; -1];
% w0_3 = [-1; -4; 1];
% options = odeset('RelTol',10^-9,'AbsTol',10^-12);
% [t1,w1] = ode45(@F,tspan,w0_1,options);
% options = odeset('RelTol',10^-9,'AbsTol',10^-12);
% [t2,w2] = ode45(@F,tspan,w0_2,options);
% options = odeset('RelTol',10^-9,'AbsTol',10^-12);
% [t3,w3] = ode45(@F,tspan,w0_3,options);

toc;

figure;
plot(t1,w1(:,2),'r-',t2,w2(:,2),'b-',t3,w3(:,2),'c-','LineWidth',2);
xlim([0 tf]);
set(gca,'fontsize',24);
xlabel('t'); ylabel('y');
legend('1','2','3','Location','NorthEast');
title('Lorenz equations');

N1 = size(t1,1);
M1 = round(0.5*N1); % to discard transients in attractor for first M steps

figure;
plot3(w1(M1:N1,1),w1(M1:N1,3),w1(M1:N1,2),'r-',...
    xeq,zeq,yeq,'m.',-xeq,zeq,-yeq,'m.','MarkerSize',24,'LineWidth',2);
grid on;
set(gca,'fontsize',24);
xlabel('x'); ylabel('z'); zlabel('y');
title('Lorenz attractor 1');

N2 = size(t2,1);
M2 = round(0.5*N2);

figure;
plot3(w2(M2:N2,1),w2(M2:N2,3),w2(M2:N2,2),'b-',...
    xeq,zeq,yeq,'m.',-xeq,zeq,-yeq,'m.','MarkerSize',24,'LineWidth',2);
grid on;
set(gca,'fontsize',24);
xlabel('x'); ylabel('z'); zlabel('y');
title('Lorenz attractor 2');

N3 = size(t3,1);
M3 = round(0.5*N3);

figure;
plot3(w3(M3:N3,1),w3(M3:N3,3),w3(M3:N3,2),'c-',...
    xeq,zeq,yeq,'m.',-xeq,zeq,-yeq,'m.','MarkerSize',24,'LineWidth',2);
grid on;
set(gca,'fontsize',24);
xlabel('x'); ylabel('z'); zlabel('y');
title('Lorenz attractor 3');

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









