// Intialize Velocity field
volScalarField x = mesh.C().component(vector::X);
volScalarField y = mesh.C().component(vector::Y);
volScalarField z = mesh.C().component(vector::Z); 

// Estimate the initial velocity field for error calculation
volVectorField* UstartPointer;

if (init && TaylorGreenDimension == "2D") {
    UstartPointer = new volVectorField (
                    Uinit* 
                    ( vector(1,0,0) * sin(x/L) * cos(y/L) 
                    - vector(0,1,0) * cos(x/L) * sin(y/L) 
                    + vector(0,0,1) * scalar(0.)) );
} else if (init && TaylorGreenDimension == "3D") {
    UstartPointer = new volVectorField (
                    Uinit* 
                    ( vector(1,0,0) * sin(x/L) * cos(y/L) * cos(z/L)
                    - vector(0,1,0) * cos(x/L) * sin(y/L) * cos(z/L)
                    + vector(0,0,1) * scalar(0.)) );
} //else { UstartPointer = new volVectorField (mesh, dimensionedVector("Ustart",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)); }

//volVectorField Ustart = *UstartPointer;

if(runTime.timeName() == "0") {
    // Reynolds number
    dimensionedScalar Re = Uinit*L/nu;
    Info<< "\nInitializing Taylor-Greem Vortex flow with Re = "<< Re.value() << "\n" << endl;
    
    Info<< "Writing initialized field U" << endl;
    U = *UstartPointer;
    U.correctBoundaryConditions();
    U.write();

    Uhat   = U;
    Ustar  = U;
    
    // initialize U0 as U
    U0 = U;
    U0.correctBoundaryConditions();
    Info<< "Writing field U0 = U" << endl;
    U0.write();
}

//Pressure-------------------------------------------------------------------- //
dimensionedScalar pRefValueAna("pRefValueAna", dimensionSet(0,2,-2,0,0), 0);
dimensionedScalar pRefValueAnaStart("pRefValueAna", dimensionSet(0,2,-2,0,0), 0);

// Estimate the initial pressure field for error calculation
volScalarField* pstartPointer;
if (init && TaylorGreenDimension == "2D") {
    pstartPointer = new volScalarField(sqr(Uinit)/4* (cos(2.*x/L) + cos(2.*y/L)));
} else if (init && TaylorGreenDimension == "3D"){
    pstartPointer = new volScalarField(sqr(Uinit)/4* (cos(2.*x/L) + cos(2.*y/L)));
}
//volScalarField pstart = *pstartPointer;


if(runTime.timeName() == "0") {
    p = *pstartPointer;
    
    if(pRefOn) {
        const scalar xPRefCell = p.mesh().C()[pRefCell].x();
        const scalar yPRefCell = p.mesh().C()[pRefCell].y();
        const scalar zPRefCell = p.mesh().C()[pRefCell].z();
    
        if (init && TaylorGreenDimension == "2D") {
            pRefValueAnaStart  = sqr(Uinit)* 0.25*(
                                Foam::cos(2.*xPRefCell/L.value())
                              + Foam::cos(2.*yPRefCell/L.value()) );
            
        } else if (init && TaylorGreenDimension == "3D") {
            pRefValueAnaStart  = sqr(Uinit)/16*
                                (   Foam::cos(2.*xPRefCell/L.value())
                                + Foam::cos(2.*yPRefCell/L.value()) 
                                )*( Foam::cos(2.*zPRefCell/L.value()) + 2.0);           
        }
        p = p - pRefValueAnaStart;
    }
    p.correctBoundaryConditions();
    Info<< "Writing field p" << endl;
    p.write();    
}



