function lorenzeigs
%lorenzeigs VERSION 8-21-2023
% calculates the eigenvalues of the Lorenz Jacobian at the equilibria of 
% the Lorenz equations.  
% Note: Lorenz orbits are (mostly) stable for r <= 24.737

b = 8/3; sigma = 10;
r = 28;    % standard chaotic attractor
% r = 24.736 % Moler's value 
% r = 24.735 % Try values near Moler's value
% r = 24.74  % unstable
% r = 24.73  % stable
% r = 10;    % -> fixed point
% r = 99.65; % complicated periodic
% r = 100.5; % complicated periodic
% r = 160;   % complicated periodic
% r = 350;   % simplest periodic orbit

eta = sqrt(b*(r-1));

x = eta; y = eta; z = r - 1;
J = [-sigma sigma 0
     r-z    -1    -x
     y      x     -b];
fprintf('Eigenvalues of J for unstable equilibrium at center of butterfly');
fprintf(' wing:\n(x,y,z) = (%g,%g,%g)\n',x,y,z);
lambda_plus = eig(J)

x = -eta; y = -eta; z = r - 1;
J = [-sigma sigma 0
     r-z    -1    -x
     y      x     -b];
fprintf('Eigenvalues of J for unstable equilibrium at center of butterfly');
fprintf(' wing:\n(x,y,z) = (%g,%g,%g)\n',x,y,z);
lambda_minus = eig(J)

x = 0; y = 0; z = 0;
J = [-sigma sigma 0
     r-z    -1    -x
     y      x     -b];
fprintf('Eigenvalues of J for very unstable equilibrium at origin:\n');
fprintf('(x,y,z) = (%g,%g,%g)\n',x,y,z);
lambda0 = eig(J)

end









