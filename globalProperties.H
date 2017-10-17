  // Dissipation rate normalized
  epsilon =  (0.5*nu*average( magSqr(fvc::grad(U) + fvc::grad(U)().T())))/(Uinit*Uinit*Uinit/L);
  // Kinetic enery normalized
  Ek = 0.5 * average(magSqr(U)) / (Uinit*Uinit);

