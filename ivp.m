function ivp(steps,tf,y0)
%ivp1(steps,tf,y0) VERSION 8-20-2023 
% solves the initial value problem dy/dt = -y
% with y(0) = y0 using forward Euler, backward Euler, and TR.
% Try: ivp(10,2,1)

tic; 

fprintf('steps = %g, tf = %g, y(0) = %g\n',steps,tf,y0);
dt = tf/steps;
fprintf('dt = %g\n',dt);
y1 = zeros(steps+1,1); % set y's to be steps+1 dimensional column vectors
y2 = zeros(steps+1,1);
y3 = zeros(steps+1,1);
t_pts = linspace(0,tf,steps+1);

y1(1) = y0; y2(1) = y0; y3(1) = y0;
for n = 1:steps % timestep loop
    y1(n+1) = (1 - dt)*y1(n); % forward Euler
    y2(n+1) = y2(n)/(1 + dt); % backward Euler
    y3(n+1) = (1 - dt/2)*y3(n)/(1 + dt/2); % trapezoidal rule (TR)
end

toc;

figure;
plot(t_pts,y1,'r-',t_pts,y2,'b-',t_pts,y3,'c-',t_pts,exp(-t_pts),'k.',...
    'MarkerSize',24,'LineWidth',2);
set(gca,'fontsize',24);
xlim([0 tf]);
xlabel('t'); ylabel('y');
legend('FE','BE','TR','exact','Location','NorthEast');

figure;
plot(t_pts,y2,'b-',t_pts,y3,'c-',t_pts,exp(-t_pts),'k.',...
    'MarkerSize',24,'LineWidth',2);
set(gca,'fontsize',24);
xlim([0 tf]);
xlabel('t'); ylabel('y');
legend('BE','TR','exact','Location','NorthEast');

numerr = sum(abs(exp(-t_pts)' - y1))/(steps+1);
fprintf('|FE_err|_1 = %g\n',numerr);

numerr = sum(abs(exp(-t_pts)' - y2))/(steps+1);
fprintf('|BE_err|_1 = %g\n',numerr);

numerr = sum(abs(exp(-t_pts)' - y3))/(steps+1);
fprintf('|TR_err|_1 = %g\n',numerr);

end









