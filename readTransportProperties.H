Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED, //MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

dimensionedScalar nu
(
    transportProperties.lookup("nu")
);
Info << "Viscosity: " << nu.value() << endl;

//  Read initial Velocity V_0
dimensionedScalar Uinit
(
    "Uinit",
    dimVelocity,    
    transportProperties.lookupOrDefault<scalar>("Uinit", 0)
);
Info << "Velocity Scale: " << Uinit.value() << endl;

//  Read length scale
dimensionedScalar L
(
    "L",
    dimLength,
    transportProperties.lookupOrDefault<scalar>("L", 1.0)
);
Info << "Length Scale: " << L.value() << endl;

const Switch init
(
    transportProperties.lookupOrDefault("init", false)
);

const word TaylorGreenDimension
(
    transportProperties.lookupOrDefault<word>("TaylorGreenDimension", "2D")   
);
// ************************************************************************* //
