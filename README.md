## Computer Codes for *Applied Numerical Methods for Partial Differential Equations* by Carl L. Gardner, Springer (November 2024)

Email comments to 
     carl.gardner@asu.edu

The computer programs are all self-contained for simplicity of use. For example, **chd.m** has a built-in WENO3 solver identical to the one in **weno3.m**.

The WENO3 solver was developed by Jeremiah Jones and Carl Gardner, SoMSS, Arizona State University.

**Errata and Clarifications** for the Book are in the file **errata.pdf**.

---

**burgers1.m** conservative upwind method and four versions of the Lax-Friedrichs method for the inviscid Burgers equation

**burgers3.m** WENO3 method for the inviscid Burgers equation

**bvp.m** solves a linear BVP using central differences and a tridiagonal direct solve

**chd.m** WENO3 method for the classical hydrodynamic model for semiconductor devices

**diffusion0.m** backward Euler method for the linear diffusion equation with fixed timestep and Dirichlet BCs

**diffusion1.m** TRBDF2 method for the linear diffusion equation with fixed timestep and Dirichlet BCs

**diffusion2.m** TRBDF2 method for the linear diffusion equation with dynamic timestep and Dirichlet BCs

**diffusionNBC.m** TRBDF2 method for the linear diffusion equation with dynamic timestep and homogeneous Neumann BCs

**diffusion2D.m** TRBDF2 method for the 2D linear diffusion equation with fixed timestep for movie and Dirichlet BCs

**ivp.m** solves the IVP 𝑑𝑦/𝑑𝑡 = −𝑦, 𝑦(0) = 𝑦0, using forward Euler, backward Euler, and TR

**jets.m** simulates 2D gas dynamical supersonic jets using the WENO3 method 

**laplace0.m** banded matrix direct solve for the 2D Laplace equation with Dirichlet BCs

**laplace1.m** Jacobi iteration for the 2D Laplace equation with Dirichlet BCs

**laplace2.m** Gauss-Seidel iteration for the 2D Laplace equation with Dirichlet BCs

**laplace3.m** SOR iteration for the 2D Laplace equation with Dirichlet BCs 

**laplace4.m** conjugate gradient method for the 2D Laplace equation with Dirichlet BCs

**laplace5.m** PCG method for the 2D Laplace equation with Dirichlet BCs 

**laplace6.m** MATLAB built-in CG, PCG, or GMRES method for the 2D Laplace equation with Dirichlet BCs

**laplaceNBC.m** SOR method for the 2D Laplace equation with a Neumann BC on one side

**layer.m** solves the nonlinear layer BVP using central differences by Newton iteration (with a tridiagonal direct solve)

**lorenz1.m** solves the Lorenz equations using fourth/fifth-order Runge-Kutta or TRBDF2

**lorenz2.m** solves the Lorenz equations with two different sets of initial conditions using fourth/fifth-order Runge-Kutta

**lorenz3.m** solves the Lorenz equations with three different methods or sets of initial conditions using fourth/fifth-order Runge-Kutta

**lorenzeigs.m** calculates the eigenvalues of the Lorenz Jacobian at the equilibria of the Lorenz equations

**nonlin_diffusion.m** solves the nonlinear diffusion equation with homogeneous Neumann BCs using TRBDF2 with dynamic timestep and Newton’s method 

**nonlin_laplace.m** solves the 2D nonlinear Laplace equation with Dirichlet BCs using Newton’s method with a banded matrix direct solve or GMRES 

**pendulum.m** plots the phase space diagram/direction field for the nonlinear pendulum

**rk2.m** solves the IVP 𝑑𝑤/𝑑𝑡 = 𝑓(𝑤) using second-order Runge-Kutta with fixed timestep

**rk2dyn.m** solves the IVP 𝑑𝑤/𝑑𝑡 = 𝑓(𝑤) using second-order Runge-Kutta with dynamic timestep

**shaw.m** solves the Shaw equations using fourth/fifth-order Runge-Kutta 

**wave1.m** solves the first-order wave equation using the original upwind method 

**wave2.m** solves the second-order wave equation as a first-order system using the Lax-Friedrichs or Lax-Wendroff method

**wave3.m** solves the first-order wave equation using the WENO3 method 

**weno3.m** solves 1D gas dynamical Riemann problems using the WENO3 method with Lax-Friedrichs flux splitting



