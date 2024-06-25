function shaw(u0,tf)
%shaw(u0,tf) VERSION 8-21-2023 
% solves the Shaw equations from t = 0 to tf using fourth/fifth-order 
% Runge-Kutta.
% Try: shaw(-0.73,100)

tic;

tspan = [0 tf];
w0 = [u0; 0; 0];

fprintf('w0 = \n');
disp(w0);
fprintf('tf = %g\n',tf);

opts = odeset('RelTol',10^-6,'AbsTol',10^-9);
% opts = odeset('RelTol',10^-12,'AbsTol',10^-15);
[t,w] = ode45(@F,tspan,w0,opts);

toc;

figure;
plot(t,w(:,1),'r-','LineWidth',2);
set(gca,'fontsize',24);
xlabel('t'); ylabel('u');
title('Shaw oscillator');

figure;
plot(w(:,1),w(:,2),'b-',u0,0,'c.','MarkerSize',24,'LineWidth',2);
axis square;
set(gca,'fontsize',24);
xlabel('u'); ylabel('v'); 
title('Shaw attractor with transients');

N = size(t,1);
M = round(0.5*N); % to discard transients in attractor for first M steps

figure;
plot(w(M:N,1),w(M:N,2),'b-','MarkerSize',24,'LineWidth',2);
axis square;
set(gca,'fontsize',24);
xlabel('u'); ylabel('v');
title('Shaw attractor');

figure;
plot3(w(M:N,1),w(M:N,3),w(M:N,2),'b-','MarkerSize',24,'LineWidth',2);
grid on;
set(gca,'fontsize',24);
xlabel('u'); ylabel('t'); zlabel('v');
title('3D Shaw attractor');

figure;
comet3(w(:,1),w(:,3),w(:,2));

    function wprime = F(~,w)
        wprime = [0.7*w(2) + 10*w(1)*(0.1 - w(2)^2)
                  -w(1) + 0.25*sin(1.57*w(3))
                  1];
    end

end









