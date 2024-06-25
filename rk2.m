function rk2(w0,steps,tf,prob)
%rk2(w0,steps,tf,type) VERSION 8-20-2023
% solves dw/dt = f(w) from t = 0 to tf = steps*dt
% with the column vector of ICs w(0) = w0
% using second-order Runge-Kutta with a fixed timestep.
% The function f(w) is implemented below for problem types:
%   1 damped harmonic oscillator y'' + 2 p y' + omega2 y = 0
%      rk2([1; -0.05],200,8*pi,1)
%   2 van der Pol oscillator
%      rk2([1; 0],2000,50,2)

tic;

omega2 = 1; % for harmonic oscillator
p = 0.05; % with friction
% p = 0; % no friction

if prob == 1
    fprintf('damped harmonic oscillator\n');
    fprintf('omega^2 = %g, p = %g\n',omega2,p);
elseif prob == 2
    fprintf('van der Pol oscillator\n');
else
    error('incorrect problem type: %g',prob);
end

fprintf('w0 = \n');
disp(w0);
fprintf('steps = %g, tf = %g\n',steps,tf);
t0 = 0;
dt = (tf - t0)/steps;
fprintf('dt = %g\n',dt);
wsol = zeros(steps+1,2); % store solution in rows
tsol = linspace(t0,tf,steps+1);

w = w0;
wsol(1,:) = transpose(w);
for n = 1:steps % timestep loop
    w = w + dt*f(w + dt*f(w)/2); % RK2
    wsol(n+1,:) = transpose(w);
end

toc;

if prob == 1
    figure;
    plot(tsol,exp(-p*tsol).*cos(sqrt((-p^2+omega2))*tsol),'c-',...
        tsol,wsol(:,1),'b.','MarkerSize',24,'LineWidth',2);
    legend('exact','RK2','Location','NorthEast');
    set(gca,'fontsize',24);
    xlim([t0 tf]);
    xlabel('t'); ylabel('y');

    numerr = sum(abs((exp(-p*tsol).*cos(sqrt((-p^2+omega2))*tsol))' - ...
        wsol(:,1)))/(steps+1);
    fprintf('|err|_1 = %g\n',numerr);
else
    figure;
    plot(tsol,wsol(:,1),'r-','MarkerSize',12,'LineWidth',2);
    set(gca,'fontsize',24);
    xlim([t0 tf]);
    xlabel('t'); ylabel('y');
end

figure;
plot(wsol(:,1),wsol(:,2),'b-',wsol(1,1),wsol(1,2),'r.',...
    'MarkerSize',24,'LineWidth',2);
axis square;
set(gca,'fontsize',24);
xlabel('y'); ylabel('dy/dt');
title('attractor');

    function wprime = f(w)
        if prob == 1 % damped harmonic oscillator
            wprime = [w(2); -omega2*w(1) - 2*p*w(2)];
        elseif prob == 2 % van der Pol oscillator
            wprime = [w(2); -w(1) + (4 - w(1)^2)*w(2)];
        else
            error('incorrect type = %g',prob);
        end
    end

end









