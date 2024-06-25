function layer(epsilon,N)
%layer(epsilon,N) VERSION 8-21-2023 
% solves the nonlinear BVP
%     epsilon*y'' + 2*y' + exp(y) = 0, y(0) = 0, y(1) = 0 
% using central differences with N dx by Newton iteration. 
% The solution has a boundary layer at the left with width 
%       -epsilon ln(epsilon)/2
% Try: layer(0.1,100)

tic;

fprintf('epsilon = %g, N = %g\n',epsilon,N);
EPSILON = 10^-10;
h = 1/N; % h = dx
x = linspace(0,1,N+1);
y = zeros(N-1,1);

e = ones(N-1,1);
D2 = spdiags([e -2*e e],[-1 0 1],N-1,N-1)/h^2;
D1 = spdiags([-e e],[-1 1],N-1,N-1)/(2*h);

F = epsilon*D2*y + 2*D1*y + exp(y);
normresidual0 = norm(F,1)/(N-1);
fprintf('\t||residual_0|| = %g\n',normresidual0);
normresidual = Inf;
iter = 0;
while normresidual > EPSILON*normresidual0 && iter < 20
    iter = iter + 1;
    J = epsilon*D2 + 2*D1 + spdiags(exp(y),0,N-1,N-1);
    y = y - J\F;
    F = epsilon*D2*y + 2*D1*y + exp(y);
    normresidual = norm(F,1)/(N-1);
    fprintf('Newton iteration = %d\n',iter);
    fprintf('\t||residual|| = %g\n',normresidual);
end
y = [0; y; 0]; % add in BCs to y

toc;

% To superimpose multiple plots uncomment:
% hold on;
% and comment:
figure;
% 
plot(x,log(2./(x+1)) - log(2)*exp(-2*x/epsilon),'b-',x,y,'c.',...
    'MarkerSize',12,'LineWidth',2);
xlim([0 1]); ylim([0 0.7]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('y');
legend('computed','approx analytical','Location','NorthEast');
title('Nonlinear Boundary Layer Solution');

% pts = linspace(10^-6,0.2,201);
% epsilon_pts = [0.01;0.05;0.1;0.2];
% widths = [0.025;0.085;0.14;0.22];
% figure;
% plot(pts,-pts.*(log(pts))/2,'r-',epsilon_pts,widths,'b.',...
%     'MarkerSize',24,'LineWidth',2);
% xlim([0 0.2]); ylim([0 0.25]);
% set(gca,'fontsize',24);
% xlabel('\epsilon'); ylabel('${\cal B}$','Interpreter','latex');

end









