function wave2(N,steps,CFL,method)
%wave2(N,steps,CFL,method) VERSION 8-24-2023
% solves the wave equation w_tt - w_xx = 0 as a 
%   first-order system (u = w_x, v = w_t):
%      u_t = v_x 
%      v_t = u_x
%   with N+1 grid points for 0 <= x <= 1, h = 1/N, dt = CFL*h, and
%       t_f = steps*dt.
%   Methods: 1 LF
%            2 LW
%   For the ICs, N >= 20 should be an integer multiple of 10.
%   Uses periodic BCs.
%   Try: wave2(200,280,1,2)

tic;

if method == 1 % LF
    fprintf('LF method\n');
elseif method == 2 % LW
    fprintf('LW method\n');
else
    error('incorrect method type: %g',method);
end

fprintf('N = %g, steps = %g, CFL = %g\n',N,steps,CFL);
h = 1/N;
dt = CFL*h;
u = zeros(N+1,1); v = zeros(N+1,1); % N+1 dimensional column vectors

j = N/2-N/20:N/2+N/20;
u(j+1) = -cos(10*j*pi/N); % ICs
x = linspace(0,1,N+1);

t = 0;
% figure; % for movie
for n = 1:steps % timestep loop
    t = t + dt;
    % with periodic BCs
    % grid pts: ... N+1 | 1 2 ... N N+1 | 1 ...
    uleft = [u(N+1);u(1:N)];
    uright = [u(2:N+1);u(1)];
    vleft = [v(N+1);v(1:N)];
    vright = [v(2:N+1);v(1)];

    if method == 1 % LF
        u = (uleft + uright)/2 + dt*(vright - vleft)/(2*h);
        v = (vleft + vright)/2 + dt*(uright - uleft)/(2*h);
    elseif method == 2 % LW
        % LF partial step: timestep = dt/2
        umidright = (u + uright)/2 + dt*(vright - v)/(2*h);
        vmidright = (v + vright)/2 + dt*(uright - u)/(2*h);
        umidleft = (u + uleft)/2 + dt*(v - vleft)/(2*h);
        vmidleft = (v + vleft)/2 + dt*(u - uleft)/(2*h);
        % leapfrog partial step: timestep = dt
        u = u + dt*(vmidright - vmidleft)/h;
        v = v + dt*(umidright - umidleft)/h;
    end

    plot(x,u,'r-','LineWidth',2);
    set(gca,'fontsize',24);
    xlabel('x'); ylabel('u');
    axis([0 1 0 1]);
    getframe;
end
fprintf('tf = %g\n',t);

toc;

end









