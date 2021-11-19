    // Momentum predictor
    Info<< "\ncorrectBoundaryVelocity\n" << endl;
    
    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        Info<< "Inside if(simple.momentumPredictor())" << endl;
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
