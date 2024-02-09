/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "zetaF.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
//#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(zetaF, 0);
addToRunTimeSelectionTable(RASModel, zetaF, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> zetaF::Tau() const
{
/*    
    volScalarField T_lb("T_lb", CTau_*sqrt(nu()/(epsilon_+epsilonMin_)));
    volScalarField T_ub("T_ub", a_/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_+mTSmall_));
    volScalarField T_nb("T_nb",  k_/(epsilon_+epsilonMin_));
    
    dimensionedScalar TSmall("TSmall", dimTime, 1e-15);
    
    volScalarField v_min = min(T_nb, T_ub);
    volScalarField I_min = 1.0*pos(T_nb - T_ub - TSmall); // = 1 wenn gebounded
    volScalarField I_max = 2.0*pos(T_lb - v_min - TSmall); // = 1 wenn gebounded
    volScalarField TInd("TInd", I_min + I_max);
    
    if (runTime_.outputTime())
    {
        TInd.write();
        T_lb.write();
        T_ub.write();
        T_nb.write();
    }
*/   
    return max(
                  min(
                           k_/(epsilon_+epsilonMin_),
                                 ( a_/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_+mTSmall_))
                       ),
                   CTau_*sqrt(nu()/(epsilon_+epsilonMin_))
               );
}

tmp<volScalarField> zetaF::L() const
{
/*
    volScalarField L_lb("L_lb", CEta_*pow( (pow(nu(),3)/(epsilon_+epsilonMin_)),0.25));
    volScalarField L_ub("L_ub",sqrt(k_)/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_));
    volScalarField L_nb("L_nb",pow(k_,1.5)/(epsilon_+epsilonMin_));
    
    dimensionedScalar LSmall("LSmall", dimLength, 1e-15);
    
    volScalarField v_min = min(L_nb, L_ub);
    volScalarField I_min = 1.0*pos(L_nb - L_ub - LSmall); // = 1 wenn gebounded
    volScalarField I_max = 2.0*pos(L_lb - v_min - LSmall); // = 1 wenn gebounded
    volScalarField LInd("LInd", I_min + I_max);
    
    if (runTime_.outputTime())
    {
        LInd.write();
        L_lb.write();
        L_ub.write();
        L_nb.write();
    }
*/
    return CL_*max(
                      min(
                              pow(k_,1.5)/(epsilon_+epsilonMin_),
                                sqrt(k_)/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_)
                          ),
                      CEta_*pow( (pow(nu(),3)/(epsilon_+epsilonMin_)),0.25)
                  );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

zetaF::zetaF
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.22
        )
    ),
    CEps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.9
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            0.4
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            0.65
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),
    sigmaZ_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaZ",
            coeffDict_,
            1.2
        )
    ),
    CTau_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTau",
            coeffDict_,
            6.0
        )
    ),
    CL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            coeffDict_,
            0.36
        )
    ),
    CEta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CEta",
            coeffDict_,
            85
        )
    ),
    a_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a",
            coeffDict_,
            0.6
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    zeta_
    (
        IOobject
        (
            "zeta",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    f_
    (
        IOobject
        (
            "f",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
   
    yr_(mesh_),

    mTSmall_
    (
	"mTSmall",	
	dimensionSet(0, 0, -1, 0, 0, 0, 0),
	1e-10
    ),
	
    zetaMin_
    (
        "zetaMin",
         dimless,
         SMALL
    ),
    fMin_
    (
        "fMin",
        dimless/dimTime,
        SMALL
    )
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);
    bound(zeta_, zetaMin_);
    printCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<volSymmTensorField> zetaF::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> zetaF::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> zetaF::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}

tmp<fvVectorMatrix> zetaF::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}

bool zetaF::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        CEps2_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaK_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        sigmaZ_.readIfPresent(coeffDict());
        CTau_.readIfPresent(coeffDict());
        CL_.readIfPresent(coeffDict());
        CEta_.readIfPresent(coeffDict());
        a_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void zetaF::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    volScalarField G("RASModel::G", nut_*2*magSqr(symm(fvc::grad(U_))));
    volScalarField CEps1_ = 1.4*(1.0+(0.012/(zeta_+zetaMin_)));

    volScalarField T_ = Tau();
    volScalarField L_ = L();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::Sp(fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        CEps1_*G/T_
      - fvm::Sp(CEps2_/T_, epsilon_)
/*
	CEps1_*G*epsilon_/k_
      - fvm::Sp(CEps2_*epsilon_/k_, epsilon_)
*/
    );

    epsEqn().relax();
    //epsEqn().boundaryManipulate(epsilon_.boundaryField());

    dimensionedScalar nu1 = nu()->internalField()[0];
    #include "BoundaryConditionEp.H"

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(epsilon_/k_, k_)
     );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);

    // f equation
    tmp<fvScalarMatrix> fEqn
    (
        - fvm::laplacian(f_)
     ==
        - fvm::Sp(1.0/sqr(L_),f_)
        - (C1_+(C2_*G/(epsilon_)))*((zeta_ - 2.0/3.0))/(sqr(L_)*T_)
    );
    fEqn().relax();
    solve(fEqn);
    bound(f_, fMin_);

    volScalarField fTilda = f_ - 2.0*nu()*zeta_/sqr(yr_);

    // zeta equation
    tmp<fvScalarMatrix> zetaEqn
    (
        fvm::ddt(zeta_)
      + fvm::div(phi_, zeta_)
      - fvm::Sp(fvc::div(phi_), zeta_)
      - fvm::laplacian(DzetaEff(), zeta_)
     ==
        fTilda
      - fvm::Sp(G/(k_), zeta_)
    );

    zetaEqn().relax();
    solve(zetaEqn);
    bound(zeta_, zetaMin_);

    // Re-calculate viscosity
    nut_ = Cmu_*zeta_*k_*T_;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
