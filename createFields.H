    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
  
    volVectorField U0
    (
        IOobject
        (
            "U0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );    
    //bool U0Exist = true;
    
    volVectorField W
    (
        IOobject
        (
            "W",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );    

//---------------------------------------------------------------------------------//    
//     IOobject U0Header
//     (
//         "U0",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ
//     );
//     
//     volVectorField* U0Pointer;
//     bool U0Exist = false;
// 
//     if (U0Header.headerOk())
//     {
//         U0Exist = true;
//         Info<< "Reading field U0\n" << endl;
//         
//         U0Pointer = new volVectorField
//         (
//             IOobject
//             (
//                 "U0",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh
//         );
//     } else {
//         Info<< "U0 unavaiable! Initializing U0 = U \n" << endl;
//         U0Pointer = new volVectorField
//         (
//             IOobject
//             (
//                 "U0",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             //U 
//             //mesh,
//             //dimensionedVector("U0",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero),
//             //U.boundaryField().types()
//         );  
//     }    
//     
//     volVectorField U0 = *U0Pointer;
//---------------------------------------------------------------------------------//
    
    Info<< "Creating field Ustar \n" << endl;
    volVectorField Ustar
    (
        IOobject
        (
            "Ustar",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("Ustar",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)
    );
    
    Info<< "Creating field Uhat \n" << endl;
    volVectorField Uhat
    (
        IOobject
        (
            "Uhat",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("Uhat",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)
    );
    
       
    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    
    #include "createPhi.H"
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

    const dictionary& pisoDict = mesh.solutionDict().subDict("PISO");
    const bool pRefOn = pisoDict.lookupOrDefault("pRefOn", false);
    Info << "pRefOn" << pRefOn << endl;
    
    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );
        
//     // Filter 
//     autoPtr<LESfilter> filterPtr(LESfilter::New(U.mesh(), transportProperties));
//     LESfilter& filter(filterPtr());
