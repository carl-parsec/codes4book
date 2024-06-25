function diffusion1(N,steps)
%diffusion1(N,steps) VERSION 8-22-2023
% solves the diffusion (heat) equation u_t = u_xx 
% with N dx (N+1 grid points 1,...,N+1) to time t = steps*dt
% with Dirichlet BCs u(-1,t) = 0 = u(1,t) and ICs u(x,0)= 1 - x^2.
% Uses TRBDF2 method with fixed timestep.
% An "exact" Fourier solution is also provided.
% Try: diffusion1(100,100)

tic;

fprintf('N = %g, steps = %g\n',N,steps);
GAMMA = 2 - sqrt(2);
h = 2/N;
% dt = h^2/2; % forward Euler stability limit
dt = sqrt(N)*h^2/2;
fprintf('tf = %g\n',steps*dt);

% Note there are N-1 interior points
e = ones(N-1,1);
D2 = spdiags([e -2*e e],[-1 0 1],N-1,N-1)/h^2;

j = (1:N-1)';
x = -1 + h*j;
u = 1 - x.^2; % ICs and allocates memory for u

CONST = GAMMA/2;
CONST1 = (1 - GAMMA)/(2 - GAMMA);
CONST2 = 1/(GAMMA*(2 - GAMMA));
CONST3 = (1 - GAMMA)^2/(GAMMA*(2 - GAMMA));
t = 0;
figure; % for movie
for n = 1:steps % timestep loop
    t = t + dt;
    
    umid = (speye(N-1) - CONST*dt*D2)\((speye(N-1) + CONST*dt*D2)*u); % TR
    u = (speye(N-1) - CONST1*dt*D2)\(CONST2*umid - CONST3*u); % BDF2
    
    plot([-1; x; 1],[0; u; 0],'b-','LineWidth',2);
    set(gca,'fontsize',24);
    xlabel('x'); ylabel('u');
    axis([-1 1 0 1]);
    getframe;
end

toc;

uexact = fourier(x,t);
figure;
plot([-1; x; 1],[0; uexact; 0],'r.',[-1; x; 1],[0; u; 0],'b-',...
    'LineWidth',2,'MarkerSize',12);
set(gca,'fontsize',24);
xlabel('x'); ylabel('u');
axis([-1 1 0 1]);

numerr = sum(abs(uexact - u))/(N+1);
fprintf('|err|_1 = %g\n',numerr);

    function ue = fourier(x,t)
        Nfourier = 4096;
        m = (1:Nfourier)';
        f = 32./((2*m-1)*pi).^3;

        fsum = 0;
        for m = 1:Nfourier
            fsum = fsum + f(m)*sin(0.5*pi*(2*m-1)*(x+1))*...
                exp(-(2*m-1)^2*pi^2*t/4);
        end
        ue = fsum;
    end

end









