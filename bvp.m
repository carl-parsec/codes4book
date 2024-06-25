function bvp(N)
%bvp(N) VERSION 8-21-2023
% solves the linear boundary value problem 
%       y'' + 2y' + y = 0, y(0) = 1, y(1) = 0 
% using central differences with N dx.
% The exact solution is y(x) = (1 - x) exp(-x).
% Try: bvp(50)

tic;

fprintf('N = %g\n',N);
h = 1/N; % h = dx
xpts = linspace(0,1,N+1);

e = ones(N-1,1);
D2 = spdiags([e -2*e e],[-1 0 1],N-1,N-1)/h^2;
D1 = spdiags([-e e],[-1 1],N-1,N-1)/(2*h);
% discretized ODE BVP is A y + b = 0
A = D2 + 2*D1 + speye(N-1);
b = [1/h^2 - 1/h; zeros(N-3,1); 0]; % incorporates BCs
y = -A\b;

toc;

% For displaying full matrices, uncomment full_... lines and try bvp(8)
% full_h2_D2 = full(h^2*D2)
% full_2h_D1 = full(2*h*D1)
% full_A = full(A)

y = [1; y; 0]; % add in BCs to y
figure;
plot(xpts,(1-xpts).*exp(-xpts),'r-',xpts,y,'b.','MarkerSize',12,...
    'LineWidth',2);
xlim([0 1]);
set(gca,'fontsize',24);
xlabel('x'); ylabel('y');
legend('exact','computed','Location','NorthEast');
title('Linear BVP Solution');

yexact = (1-xpts).*exp(-xpts);
numerr = sum(abs(yexact' - y))/(N+1);
fprintf('|err|_1 = %g\n',numerr);

end
