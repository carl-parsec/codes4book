function rk2dyn(w0,tf,prob)
%rk2dyn(w0,tf,type) VERSION 12-8-2023
% solves dw/dt = f(w) from t = 0 to tf
% with the column vector of ICs w(0) = w0 
% using second-order Runge-Kutta with a dynamic timestep.
%
% The function f(w) is implemented below for problem types:
%     1 damped harmonic oscillator y'' + 2 p y' + omega2 y = 0
%         rk2dyn([1; -0.05],8*pi,1)
%     2 van der Pol oscillator
%         rk2dyn([1; 0],50,2)
%     3 singularity for y' = y^2, y(0) = 1
%         rk2dyn(1,2,3)
%     4 Lorenz Equations
%         rk2dyn([0; 1; 0],40,4)

tic;

VERBOSE = 1; % 1 (verbose mode) or 0 (quiet mode)

omega2 = 1; % for harmonic oscillator
p = 0.05; % with friction
% p = 0; % no friction

if prob == 1
    fprintf('damped harmonic oscillator\n');
    fprintf('omega^2 = %g, p = %g\n',omega2,p);
elseif prob == 2
    fprintf('van der Pol oscillator\n');
elseif prob == 3
    fprintf('singularity\n');
elseif prob == 4
    fprintf('Lorenz Equations\n');
else
    error('incorrect problem type: %g',prob);
end

fprintf('w0 = \n');
disp(w0);

t0 = 0;
t = t0;
dt_max = (tf - t0)/10;
dt_min = 16*eps*(tf - t0); 
Nmax = 1000000; % max number of timesteps
EPSILON_REL = 10^-3; % 10^-6 for more accuracy
EPSILON_ABS = 10^-6; % 10^-9

w = w0;
wsol(1,:) = transpose(w);
tsol(1) = t;

% calculate initial dt
norm_solution = norm(w,1);
norm_f = norm(f(w),1);
tbar = (norm_solution + EPSILON_ABS/EPSILON_REL)/...
    (norm_f + EPSILON_ABS/EPSILON_REL);
dt = EPSILON_REL^(1/3)*tbar; % since e_l ~ dt^3
dt = min(dt,dt_max);

n = 0;
total_redo = 0;
SINGULARITY = 0;
while t < tf && n < Nmax && ~SINGULARITY % timestep loop
    n = n + 1;
    r = Inf;
    redo = -1;
    while r > 2
        redo = redo + 1;
        if redo > 0
            if VERBOSE == 1
                fprintf('REDOING TIMESTEP %d\n',n);
            end
            total_redo = total_redo + 1;
        end  
        if 1.1*dt >= tf - t
            dt = tf - t;
        end
        tnew = t + dt;
        w1 = w + dt*f(w + dt*f(w)/2); % RK2
        wmid = w + dt*f(w + dt*f(w)/4)/2; % first dt/2 for 2*RK2
        wnew = wmid + dt*f(wmid + dt*f(wmid)/4)/2; % second dt/2 for 2*RK2
        norm_solution = norm(w,1);
        e_l = (4/3)*norm(wnew - w1,1);
        r = e_l/(EPSILON_REL*norm_solution + EPSILON_ABS);
        if VERBOSE == 1
            fprintf('n = %d, t = %g, dt = %g, r = %g,\n',n,t,dt,r); 
        end
        dt = min(min(2*dt,dt/r^(1/3)),dt_max);
        if dt < dt_min
            warning('MATLAB:dt_min',...
                'possible singularity: dt = %e < dt_min at t = %g',dt,t); 
            SINGULARITY = 1; r = 0; % to exit from while loops
        end        
    end
    t = tnew;
    w = wnew;
    tsol(end+1) = t;
    wsol(end+1,:) = transpose(w);
end
fprintf('number of timesteps = %d\n',n);
fprintf('REDO timestep = %d times\n',total_redo);
fprintf('total number of steps = %d\n',n + total_redo);
fprintf('tf = %g\n',tf); 

toc;

if prob == 1
    figure;
    plot(tsol,exp(-p*tsol).*cos(sqrt((-p^2+omega2))*tsol),'c-',...
        tsol,wsol(:,1),'b.',... % tsol,0,'g.',
        'MarkerSize',12,'LineWidth',2);
    legend('exact','RK2','Location','NorthEast');
    set(gca,'fontsize',24);
    xlim([t0 tf]);
    xlabel('t'); ylabel('w_1');

    numerr = sum(abs((exp(-p*tsol).*cos(sqrt((-p^2+omega2))*tsol))' - ...
        wsol(:,1)))/(n+1);
    fprintf('|err|_1 = %g\n',numerr);
elseif prob == 2
    figure;
    plot(tsol,wsol(:,1),'c-',tsol,wsol(:,1),'b.',... % tsol,0,'g.',
        'MarkerSize',12,'LineWidth',2);
    set(gca,'fontsize',24);
    xlim([t0 tf]);
    xlabel('t'); ylabel('w_1');
else
    figure;
    plot(tsol,wsol(:,1),'c-',tsol,wsol(:,1),'b.',...
        'MarkerSize',12,'LineWidth',2);
    set(gca,'fontsize',24);
    xlim([t0 tf]);
    xlabel('t'); ylabel('w_1');
end

if prob ~= 3
    figure;
    plot(wsol(:,1),wsol(:,2),'b-',wsol(1,1),wsol(1,2),'r.',...
        'MarkerSize',24,'LineWidth',2);
    if prob == 1
        axis equal;
    else
        axis square;
    end
    set(gca,'fontsize',24);
    xlabel('w_1'); ylabel('w_2');
    title('attractor with transients');

    if prob ~= 1
        N = size(tsol,2);
        M = round(0.5*N); % to discard transients in attractor for first M steps
        figure;
        plot(wsol(M:N,1),wsol(M:N,2),'b-','MarkerSize',24,'LineWidth',2);
        if prob == 1
            axis equal;
        else
            axis square;
        end
        set(gca,'fontsize',24);
        xlabel('w_1'); ylabel('w_2');
        title('attractor');
    end
end

    function wprime = f(w)
        if prob == 1 % damped harmonic oscillator
            wprime = [w(2); -omega2*w(1) - 2*p*w(2)];
        elseif prob == 2 % van der Pol oscillator
            wprime = [w(2); -w(1) + (4 - w(1)^2)*w(2)];
        elseif prob == 3 % y' = y^2, y(0) = 1: y(t) = 1/(1-t)
            wprime = w^2;
        elseif prob == 4 % Lorenz Eqs
            wprime = [10*(w(2) - w(1))
                28*w(1) - w(2) - w(1)*w(3)
                w(1)*w(2) - 8*w(3)/3];
        else
            error('incorrect type = %g',prob);
        end
    end

end









