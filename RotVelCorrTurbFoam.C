/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    RotVelCorrFoam

Description
    Transient solver for incompressible flow using a rotational velocity-correction strategy. 
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    This version is intended to work for LES 
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "IFstream.H"
#include "OFstream.H"
#include "LESfilter.H"
#include "IOmanip.H" // for input/ouput format control
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"   
    #include "createFields.H"
    #include "initContinuityErrs.H"  
    #include "readTimeControls.H" 
    
    #include "readTransportProperties.H"
    #include "createOutFile.H" // Create output variables and outfile
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Calculate analytical field at time t=0;
    #include "initialize.H"
    
    // Calculate global properties
    #include "globalProperties.H"
    #include "writeOutput.H"  
    
    scalar gamma0; // temoral discretization coeffs
    dimensionedScalar dt = runTime.deltaT();
    
    U0.correctBoundaryConditions();     
       
    // effective viscosity    
    volScalarField nuEff(turbulence->nut() + nu);   
    volSymmTensorField SymmGradU     = twoSymm(fvc::grad(U));
    volSymmTensorField SymmGradUstar = twoSymm(fvc::grad(Ustar));
    while (runTime.loop())
    {   
        Info<< "Time = " << runTime.timeName()<< endl;
        #include "CourantNo.H"
        #include "setDeltaT.H" 
            
        // Start with one first order backward step
        if(runTime.value() == dt.value())
        {
            Info<< "\nU0 unavaiable: 1st order backward temporal discretization\n" << endl;
            gamma0 = 1.0;
            Uhat   = U;
            Ustar  = U;           
        } else {
            //Info<< "\nSecond-order backward temporal discretization\n" << endl;
            gamma0 = 1.5;
            Uhat   = 2.*U - 0.5*U0; 
            Ustar  = 2.*U-U0;      
        }
        
        // Update in every time step
        dt = runTime.deltaT();  // update deltaT
        volVectorField Utemp = U;  //hold a copy of current velocity           
        Info<< "Updating W^(n+1) velocity.... " << endl;
        W.correctBoundaryConditions(); // update W =U^(n+1) at dirichlet boundaries
        
        // Pressure Poisson for new pressure----------------------------------------//
        Info<< "Calcuting phi for pEqu:..." << endl;
        phi = (fvc::interpolate(U) & mesh.Sf());
        SymmGradU = twoSymm(fvc::grad(U));
        
        //Info<< "updates p^(n+1) at dirichlet boundary with W^(n+1)" << endl;
        //p.correctBoundaryConditions(); Will be performed automatically at pEqn.solve()
        
        fvScalarMatrix pEqn
        (
        fvm::laplacian(p) == fvc::div( Uhat/dt - (fvc::div(phi, U) - fvc::div(phi)*U) + (SymmGradU & fvc::grad(nuEff)))
  //    fvm::laplacian(p) == fvc::div( Uhat/dt - (fvc::div(phi, U) - fvc::div(phi)*U) + (twoSymm(fvc::grad(U)) & fvc::grad(nuEff)))
  //    fvm::laplacian(p) == fvc::div(Uhat/dt - (fvc::div(phi, U) - fvc::div(phi)*U))
        );
        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();
        
        // update p at all boundaries
        Info << "Final Pressure Update....." << endl;
        p.correctBoundaryConditions();
        
        // Helmholz equation for velocity ------------------------------------------//
        Info<< "Calcuting phi for UEqu:..." << endl;
        phi = (fvc::interpolate(Ustar) & mesh.Sf());
        SymmGradUstar = twoSymm(fvc::grad(Ustar));
        // update velocity at Dirichlet and open boundary
        //U.correctBoundaryConditions(); Will be performed automatically at UEqn.solve()
        fvVectorMatrix UEqn
        (
        fvm::Sp(gamma0/dt, U) - nuEff*fvm::laplacian(U)  == (Uhat/dt - (fvc::div(phi, Ustar) - fvc::div(phi)*Ustar) - fvc::grad(p) + (SymmGradUstar & fvc::grad(nuEff)))       
  //    fvm::Sp(gamma0/dt, U) - nuEff*fvm::laplacian(U)  == (Uhat/dt - (fvc::div(phi, Ustar) - fvc::div(phi)*Ustar) - fvc::grad(p) + (twoSymm(fvc::grad(Ustar)) & fvc::grad(nuEff)))       
  //    fvm::Sp(gamma0/dt, U) - nuEff * fvm::laplacian(U) == (Uhat/dt - (fvc::div(phi, Ustar) - fvc::div(phi)*Ustar) - fvc::grad(p))
        );
        UEqn.solve();
        // update velocity at every boundary
        //U.correctBoundaryConditions(); is automatically updated after UEqn.solve()
        
        // update U0 boundaries at previous time step 
        //(time step down from current time is addressed in the boundary condition)
        U0 = Utemp;
        Info<< "Updating U0^(n) velocity.... " << endl;
        U0.correctBoundaryConditions(); 
        
        #include "continuityErrs.H"
        
        // update turbulence veriables (such as nuSgs)
        turbulence->correct();
        nuEff = turbulence->nut() + nu;
        
        runTime.write();
            
        // Calculate global properties
        #include "globalProperties.H"
        #include "writeOutput.H"        
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
