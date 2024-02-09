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

#include "kZetaF.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "mappedWallFvPatch.H"

//#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kZetaF, 0);
addToRunTimeSelectionTable(RASModel, kZetaF, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> kZetaF::Tau() const
{
    return max
           (
               min
               (
                   k_/(epsilon_+epsilonMin_),
                   a_/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_+mTSmall_)    //BK: Durbin Original 1/(sqrt(6) zeta Cmu sqrt(Sij Sij)
               ),
               CTau_*sqrt((mu()/rho_)/(epsilon_+epsilonMin_))
           ); 
}

tmp<volScalarField> kZetaF::L() const
{
    return CL_*max
               (
                   min                                                              //BK: in Samules Version auskommentiert
                   (                                                                //BK: in Samules Version auskommentiert
                       pow(k_,1.5)/(epsilon_+epsilonMin_),
                       sqrt(k_)/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_)   //BK: in Samules Version auskommentiert
                   ),                                                               //BK: in Samules Version auskommentiert
                   CEta_*pow( (pow((mu()/rho_),3)/(epsilon_+epsilonMin_)),0.25)
               );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kZetaF::kZetaF
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const fluidThermo& thermophysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, rho, U, phi, thermophysicalModel, turbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.22
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.9
        )
    ),
    Ceps3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps3",
            coeffDict_,
            -0.33
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

    Prt_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt",
            coeffDict_,
            1.0
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
    
    mut_
    (
        IOobject
        (
            "mut",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    alphat_
    (
        IOobject
        (
            "alphat",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    yr_(mesh_),
    
    zetaMin_("zetaMin", dimless, SMALL),
    fMin_("fMin", dimless/dimTime, SMALL),
    mTSmall_("mTSmall", dimless/dimTime, SMALL)
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);
    bound(f_, fMin_);
    bound(zeta_, zetaMin_);
    
    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();
    
    mut_ = rho_*Cmu_*zeta_*k_*Tau();
    mut_.correctBoundaryConditions();
    
    printCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<volSymmTensorField> kZetaF::R() const
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
            ((2.0/3.0)*I)*k_ - (mut_/rho_) * twoSymm(fvc::grad(U_)), //BK: Scheint korrekt
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> kZetaF::devRhoReff() const
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
           -muEff()*dev(twoSymm(fvc::grad(U_))) //BK: Scheint korrekt
        )
    );
}


tmp<fvVectorMatrix> kZetaF::divDevRhoReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(muEff(), U)
      - fvc::div(muEff()*dev2(T(fvc::grad(U))))
    ); //BK: Scheint korrekt
}

bool kZetaF::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        Ceps2_.readIfPresent(coeffDict());
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


void kZetaF::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    volScalarField divU(fvc::div(phi_/fvc::interpolate(rho_)));
    
    if (mesh_.moving())
    {
        divU += fvc::div(mesh_.phi());
    }
    
    const volTensorField gradU(fvc::grad(U_));
    const volScalarField S2(2*magSqr(dev(symm(gradU)))); //BK: S2 = 2 Sij Sij - 2/3 Skk Skk; 2/3 k vernachlässigt
    const volScalarField G(GName(), mut_*S2);
    //volScalarField G("RASModel::G", mut_*2*magSqr(symm(fvc::grad(U_)))); //BK: symm(gradU) oder dev(symm(gradU))???
    
    volScalarField Ceps1_ = 1.4*(1.0+(0.012/zeta_));
    
    volScalarField T_ = Tau();
    volScalarField L_ = L();
    

    // Update epsilon (and possibly G) at the wall
    //epsilon_.boundaryField().updateCoeffs();
     
    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(rho_, epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        Ceps1_*G/T_
      - fvm::Sp(rho_*Ceps2_/T_, epsilon_) //BK: Dichte
/*
        Ceps1_*G*epsilon_/k_
      - fvm::Sp(Ceps2_*epsilon_/k_, epsilon_)
*/
      - fvm::SuSp(((2.0/3.0)*Ceps1_ + Ceps3_)*rho_*divU, epsilon_) //BK: Kompressibilität?!
    );

    epsEqn().relax();
    //epsEqn().boundaryManipulate(epsilon_.boundaryField());
   
    //dimensionedScalar rho1 = rho_[0];
    //dimensionedScalar mu1 = mu()[0];
    //dimensionedScalar nu1 = mu1/rho1; Info<<"WARNING: Only for rhoConst, ";
    //dimensionedScalar nu1("nu1", dimensionSet(0, 2, -2, 0, 0, 0, 0), 2e-5); 
    #include "BoundaryConditionEp.H"
    solve(epsEqn);
    bound(epsilon_, epsilonMin_);

    
    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rho_, k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(rho_*epsilon_/k_, k_) //BK: Dichte
      - fvm::SuSp((2.0/3.0)*rho_*divU, k_) //BK: Kompressibilität?!
    );

    kEqn().relax();
    solve(kEqn);
    //bound(k_, kMin_);
    dimensionedScalar kSmall_("kSmall", dimensionSet(0, 2, -2, 0, 0, 0, 0), 1e-10);
    bound(k_, kSmall_);


    // f equation
    tmp<fvScalarMatrix> fEqn
    (
        - fvm::laplacian(rho_, f_)
     ==
        - fvm::Sp(rho_/sqr(L_),f_) //BK: Dichte
        - (rho_*C1_+C2_*G/epsilon_)*(zeta_ - 2.0/3.0)/(sqr(L_)*T_) //BK: Dichte
    );
    
    fEqn().relax();
    solve(fEqn);
    bound(f_, fMin_);

    
    // Calculate f from the transformed fEqn
    volScalarField fTilda = f_ - 2.0*(mu()/rho_)*zeta_/sqr(yr_);

    
    // Zeta equation
    tmp<fvScalarMatrix> zetaEqn
    (
        fvm::ddt(rho_, zeta_)
      + fvm::div(phi_, zeta_)
      - fvm::laplacian(DzetaEff(), zeta_)
     ==
//      f_
        rho_*fTilda //BK: Dichte
      - fvm::Sp(G/(k_), zeta_)
//    - fvm::Sp(G/(k_+kSmall_), zeta_)
    );

    zetaEqn().relax();
    solve(zetaEqn);
    bound(zeta_, zetaMin_);
    zeta_ = min(zeta_, 2.0); //BK: Warum maximal 1.5?
    
    // Re-calculate viscosity
    mut_ = rho_*Cmu_*zeta_*k_*T_; //BK: Dichte
    //mut_.correctBoundaryConditions();

    alphat_ = mut_/Prt_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
