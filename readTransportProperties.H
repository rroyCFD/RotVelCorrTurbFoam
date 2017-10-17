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
    transportProperties.lookup("Uinit")
);
Info << "Velocity Scale: " << Uinit.value() << endl;

//  Read length scale
dimensionedScalar L
(
    transportProperties.lookup("L")
);
Info << "Length Scale: " << L.value() << endl;

const Switch init
(
    transportProperties.lookup("init")
);

const word TaylorGreenDimension
(
    transportProperties.lookup("TaylorGreenDimension")   
);
// ************************************************************************* //
