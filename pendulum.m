function pendulum(gamma)
%pendulum(gamma) VERSION 8-20-2023 
% plots the phase diagram/direction field for the nonlinear pendulum  
%       u'' + gamma u' + sin(u) = 0, (m = 1, omega0 = 1)
% where u = angle theta, along with 3 solution curves.
% Set y(1) = u and y(2) = u'.
% Try: pendulum(0.05)

tic;

fprintf('damped nonlinear pendulum: gamma = %g\n',gamma);
ya0 = [pi/3; 1];    % 1st IC
yb0 = [-8*pi/5; 2]; % 2nd IC
yc0 = [8*pi/5; -2]; % 3rd IC
t1 = 0; t2 = 100;

tspan = [t1,t2];
[ta,ya] = ode45(@f,tspan,ya0);
[tb,yb] = ode45(@f,tspan,yb0);
[tc,yc] = ode45(@f,tspan,yc0);

toc;

figure;
plot(tc/(2*pi),yc(:,1)/(2*pi),'m-',tb/(2*pi),yb(:,1)/(2*pi),'b-',...
    ta/(2*pi),ya(:,1)/(2*pi),'r-','MarkerSize',24,'LineWidth',2);
set(gca,'fontsize',24);
xlabel('t/(2\pi)'); ylabel('\theta/(2\pi)');
% xlabel('t'); ylabel('u');
xlim([t1 t2]/(2*pi));

[y1,y2] = meshgrid(-10:0.5:10,-3:0.5:3);
dy1 = y2; dy2 = -sin(y1) - gamma*y2;

figure;
quiver(y1/(2*pi),y2/(2*pi),dy1/(2*pi),dy2/(2*pi));
hold on;
plot(ya(:,1)/(2*pi),ya(:,2)/(2*pi),'r-',yb(:,1)/(2*pi),yb(:,2)/(2*pi),'b-',...
    yc(:,1)/(2*pi),yc(:,2)/(2*pi),'m-','MarkerSize',24,'LineWidth',1);
set(gca,'fontsize',24);
xlabel('\theta/(2\pi)'); ylabel('d\theta/(2\pidt)');
% xlabel('u'); ylabel('du/dt');
axis equal;
xlim([-10,10]/(2*pi));
ylim([-3,3]/(2*pi));

    function yprime = f(~,y)
        yprime = [y(2); -sin(y(1)) - gamma*y(2)];
    end

end









