function burgers1(N,steps,CFL,method,prob)
%burgers1(N,steps,CFL,method,prob) VERSION 8-24-2023
% solves the inviscid Burgers equation 
%     u_t + f_x = 0, f(u) = u^2/2
% with N+1 grid points for 0 <= x <= 1, h = 1/N, dt = CFL*h/umax, and
%     t_f = steps*dt.
% prob = 1 (shock), 2 (rarefaction), 3 (merging shocks), or 4 (N wave) 
% Uses conservative upwind method (if method = 1) and through-flow BCs.
% Methods: 1 conservative upwind
%          2 second upwind for u > 0
%          3 nonconservative upwind for u > 0
%          Lax-Friedrichs (four different versions):
%          4 original version of LF
%          5 third best version of LF
%          6 second best version of LF
%          7 best version of LF
% Try: burgers1(200,120,0.9,1,1) for Riemann problem shock
%      burgers1(200,100,0.9,1,2) for Riemann problem rarefaction
%      burgers1(200,180,0.9,1,3) for merging shocks
%      burgers1(200,90,0.9,1,4)  for N wave

tic; 

if method == 0
    fprintf('conservative upwind for u > 0\n');
elseif method == 1
    fprintf('conservative upwind\n');
elseif method == 2
    fprintf('second upwind for u > 0\n');
elseif method == 3
    fprintf('nonconservative upwind for u > 0\n');
elseif method == 4
    fprintf('original version of LF\n');
elseif method == 5
    fprintf('third best version of LF\n');
elseif method == 6
    fprintf('second best version of LF\n');
elseif method == 7
    fprintf('best version of LF\n');
else
    error('incorrect method type: %g',method);
end

fprintf('N = %g, steps = %g, CFL = %g\n',N,steps,CFL);

h = 1/N;
x = linspace(-0.5,0.5,N+1); t = 0;
u = zeros(N+1,1); % N+1 dimensional column vector
Fminus = zeros(N+1,1); % flux at j-1/2
Fplus = zeros(N+1,1);  % flux at j+1/2

% ICs: 
if prob == 1
    fprintf('Riemann problem shock wave\n');
    j = 1:N/2; u(j) = 2; 
    j = N/2+1:N+1; u(j) = 1;
elseif prob == 2
    fprintf('Riemann problem rarefaction wave\n');
    j = 1:N/2+1; u(j) = 1;  
	j = N/2+2:N+1; u(j) = 2;    
elseif prob == 3
    fprintf('merging shocks\n');
    j = 1:N/4; u(j) = 3;
    j = N/4+1:N/2; u(j) = 2;
    j = N/2+1:N+1; u(j) = 1;
elseif prob == 4
    fprintf('N wave\n');
    j = N/4+1:3*N/4+1; u(j) = 4*(4*(j-N/4-1)/N - 1);
else
    error('incorrect problem type: %g',prob);
end

figure; % for movie
for n = 1:steps % timestep loop
    umax = max(abs(u));
    dt = CFL*h/umax;
    t = t + dt;
    % with thru-flow BCs
    % grid points: 1 | 1 2 ... N N+1 | N+1
    uleft = [u(1);u(1:N)];      % j-1
    uright = [u(2:N+1);u(N+1)]; % j+1

    umid = [u(1);(u(1:N)+u(2:N+1))/2;u(N+1)];
    uminus = umid(1:N+1); % at j-1/2
    uplus = umid(2:N+2);  % at j+1/2

    if method == 0
        % conservative upwind for u > 0: correct shock speed
        u = u - dt*(u.^2 - uleft.^2)/(2*h);
    elseif method == 1
        % conservative upwind: correct shock speed
        for j = 1:N+1
            % at j+1/2
            if uplus(j) >= 0
                Fplus(j) = u(j)^2/2;
            else
                Fplus(j) = uright(j)^2/2;
            end
            % at j-1/2
            if uminus(j) >= 0
                Fminus(j) = uleft(j)^2/2;
            else
                Fminus(j) = u(j)^2/2;
            end
        end
        u = u - dt*(Fplus - Fminus)/h;
    elseif method == 2
        % second upwind (shock overshoot) for u > 0: not a good method for
        % Burgers eq -- need CFL <= 0.5 or too dispersive
        u = u - dt*(uplus.*u - uminus.*uleft)/(2*h);
    elseif method == 3
        % nonconservative upwind for u > 0: incorrect shock speed, slightly
        % incorrect rarefaction speed, and incorrect unphysical rarefaction
        % (backward shock) for CFL = 1
        u = u - dt*u.*(u - uleft)/h;
    elseif method == 4
        % Lax-Friedrichs (four different versions)
        % original version of LF
        u = (uleft + uright)/2 - dt*(uright.^2 - uleft.^2)/(4*h);
    elseif method == 5
        % used in weno3.m:
        % F_LF = 0.5*(f(1:N+2,:) + f(2:N+3,:) - alpha*(qghost(2:N+3,:) - ...
        %     qghost(1:N+2,:)));
        % third best version of LF
        r = umax*dt/h;
        u = (r*uleft + 2*(1 - r)*u + r*uright)/2 - ...
            dt*(uright.^2 - uleft.^2)/(4*h);
    elseif method == 6
        % second best version of LF
        rleft = max(abs(uleft),abs(u))*dt/h;
        rright = max(abs(uright),abs(u))*dt/h;
        u = (rleft.*uleft + (2 - rleft - rright).*u + rright.*uright)/2 - ...
            dt*(uright.^2 - uleft.^2)/(4*h);
    elseif method == 7
        % best version of LF
        rminus = abs(uminus)*dt/h;
        rplus = abs(uplus)*dt/h;
        u = (rminus.*uleft + (2 - rminus - rplus).*u + rplus.*uright)/2 - ...
            dt*(uright.^2 - uleft.^2)/(4*h);
    end

    if prob == 3 % merging shocks
        plot(x,uexact(t),'r-',x,u,'b-','LineWidth',2);
        axis([-0.5 0.5 0 3.2]);
    elseif prob == 4 % N wave
        plot(x,uexact(t),'r-',x,u,'b-','LineWidth',2);
        axis([-0.5 0.5 -4 4]);
    else % RP
        plot(x,uexact(t),'r-',x,u,'b-','LineWidth',2);
        axis([-0.5 0.5 0 2.2]); % RP
    end
    % legend('exact','computed','Location','South');
    set(gca,'fontsize',24);
    xlabel('x'); ylabel('u');
    getframe;
end
fprintf('tf = %g\n',t);

toc;

numerr = sum(abs(uexact(t) - u))/(N+1);
fprintf('|err|_1 = %g\n',numerr);

    function ue = uexact(t)
        ue = zeros(N+1,1);
        if prob == 1
            % shock wave
            s = 1.5;
            Ns = round(N/2 + s*t*N);
            Ns = min(N+1,Ns);
            i = 1:Ns; ue(i) = 2;
            i = Ns+1:N+1; ue(i) = 1;
        elseif prob == 2
            % rarefaction wave
            xleft = t;
            xright  = 2*t;
            for i = 1:N+1
                if x(i) <= xleft
                    ue(i) = 1;
                elseif x(i) >= xright
                    ue(i) = 2;
                elseif t > 0
                    ue(i) = x(i)/t;
                end
            end
        elseif prob == 3
            % ICs above were:
            %     j = 1:N/4; u(j) = 3;
            %     j = N/4+1:N/2; u(j) = 2;
            %     j = N/2+1:N+1; u(j) = 1;
            % merging shocks
            s1 = 2.5;
            Ns1 = round(N/4 + s1*t*N);
            Ns1 = min(N+1,Ns1);
            s2 = 1.5;
            Ns2 = round(N/2 + s2*t*N);
            Ns2 = min(N+1,Ns2);
            if Ns1 < Ns2
                i = 1:Ns1; ue(i) = 3;
                i = Ns1+1:Ns2; ue(i) = 2;
                i = Ns2+1:N+1; ue(i) = 1;
            else
                tm = 1/(4*(s1-s2));
                Nm = round(N/4 + s1*tm*N);
                s3 = 2;
                Ns3 = Nm + round(s3*(t-tm)*N);
                Ns3 = min(N+1,Ns3);
                i = 1:Ns3; ue(i) = 3;
                i = Ns3+1:N+1; ue(i) = 1;
            end
        elseif prob == 4
            % N wave
            t = t + 1/16;
            xleft = -sqrt(t);
            xright  = sqrt(t);
            for i = 1:N+1
                if x(i) < xleft
                    ue(i) = 0;
                elseif x(i) > xright
                    ue(i) = 0;
                else
                    ue(i) = x(i)/t;
                end
            end
        end
    end

end









