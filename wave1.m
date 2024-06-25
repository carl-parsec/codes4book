function wave1(N,steps,CFL)
%wave1(N,steps,CFL) VERSION 8-24-2023
% solves the wave equation u_t + u_x = 0
% with N dx for 0 <= x <= 1, h = 1/N, dt = CFL*h, and
%     t_f = steps*dt.
% For the ICs, N >= 20 should be an integer multiple of 10.
% Uses upwind method and through-flow BCs.
% The exact solution (red dots) is compared with the computed solution 
% (blue line).
% Try: wave1(200,150,1)
%      wave1(200,150,0.9)

tic;

fprintf('upwind method\n');
fprintf('N = %g, steps = %g, CFL = %g\n',N,steps,CFL);

h = 1/N;
dt = CFL*h;
u = zeros(N+1,1); % set u to be an N+1 dimensional column vector
j = 0:2*N/10;
u(j+1) = sin(10*j*pi/N); % ICs
x = linspace(0,1,N+1);

t = 0;
figure; % for movie
for n = 1:steps % timestep loop
    t = n*dt;

    % with periodic BCs
    % grid pts: ... N+1 | 1 2 ... N N+1 | 1 ...
    % upwind = [u(N+1);u(1:N)]; % with periodic BCs

    % with through-flow BCs
    % grid pts: ... 1 | 1 2 ... N N+1 | N+1 ...
    upwind = [u(1);u(1:N)]; % with through-flow BCs
    u = u - dt*(u - upwind)/h; % upwind method
    
    ue = uexact(t);
    plot(x,u,'b-',x,ue,'r.','LineWidth',2,'MarkerSize',16);
    set(gca,'fontsize',24);
    xlabel('x'); ylabel('u');
    axis([0 1 -1 1]);
    getframe;
end
fprintf('tf = %g\n',t);

toc;

ue = uexact(t);
numerr = sum(abs(ue - u))/(N+1);
fprintf('|err|_1 = %g\n',numerr);

    function ue = uexact(t)
        ue = zeros(N+1,1);
        s = 1;
        Nleft = round(s*t*N);
        Nright = round(s*t*N + 2*N/10);

        for i = 1:N
            if i < Nleft
                ue(i) = 0;
            elseif i > Nright
                ue(i) = 0;
            else
                ue(i+1) = sin(10*pi*(i - Nleft)/N);
            end
        end
    end
end









