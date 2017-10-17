if (Pstream::master() && Pstream::parRun())
{    
    //Info << "master processor " << endl;
    outFilePar << runTime.value()/tc.value() << "\t" << setprecision(10) << Ek.value() << "\t" << setprecision(10) << epsilon.value()<< endl;
} else {
    outFile    << runTime.value()/tc.value() << "\t" << setprecision(10) << Ek.value() << "\t" << setprecision(10) << epsilon.value()<< endl;
}

// outFile << runTime.value()/tc.value() << "\t" << Ek.value() << "\t" 
//         << epsilon.value()<< "\t" << UL2err.value() << "\t" << pL2err.value() 
//         << "\t" << ULinfErr.value() << "\t" << pLinfErr.value()<< endl;
