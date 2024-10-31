function lorenz2(delta,tf)
%lorenz2(delta,tf) VERSION 8-21-2023
% solves the Lorenz equations from t = 0 to tf with
% y(0) = 1 vs y(0) = 1 + delta and x(0) = 0 = z(0) 
% using fourth/fifth-order Runge-Kutta.
% Here w(t) = [x(t); y(t); z(t)] and the Lorenz equations are
%     dw/dt = F(w)
% Try: lorenz2(0.001,40)

tic;

fprintf('tf = %g\n',tf);
fprintf('y(0) = 1 + delta, delta = %g\n',delta);
tspan = [0 tf];
w0 = [0; 1; 0];

opts = odeset('RelTol',10^-12,'AbsTol',10^-15);
[t1,w1] = ode45(@F,tspan,w0,opts);
[t2,w2] = ode45(@F,tspan,w0+[0; delta; 0],opts);

toc;

figure;
plot(t2,w2(:,2),'r-',t1,w1(:,2),'b-','LineWidth',2);
set(gca,'fontsize',24);
xlabel('t'); ylabel('y');
legend('y_0 = 1+\delta','y_0 = 1','Location','NorthWest');
title('Lorenz equations');

    function wprime = F(~,w)
        wprime = [10*(w(2) - w(1))
                  28*w(1) - w(2) - w(1)*w(3)
                  w(1)*w(2) - 8*w(3)/3];
    end

end










