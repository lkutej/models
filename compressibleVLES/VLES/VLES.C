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

#include "VLES.H"
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

defineTypeNameAndDebug(VLES, 0);
addToRunTimeSelectionTable(RASModel, VLES, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> VLES::Tau() const
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
/*
    Info<<"T.Org, "; 
    return max
           (
               min
               (
                   k_/(epsilon_+epsilonMin_),
                   a_/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_+mTSmall_)    //BK: Durbin Original 1/(sqrt(6) zeta Cmu sqrt(Sij Sij)
               ),
               CTau_*sqrt(nu()/(epsilon_+epsilonMin_))
           ); 
*/
/*    
    Info<<"T.Org without realizability constraint, ";
    //Info<<"rho:"<<min(rho_)<<", epsilon:"<<min(epsilon_)<<endl;
    return max
           (
               k_/epsilon_,
               CTau_*sqrt((mu()/rho_)/epsilon_)
           );
*/

    Info<<"Tau = turbulent time scale, ";
    return k_/epsilon_;

/*    
    return max(k_/(epsilon_+epsilonMin_), sqrt(nu()/(epsilon_+epsilonMin_)));
*/
}

tmp<volScalarField> VLES::L() //const
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
/*
    Info<<"L.org"<<endl;
    return CL_*max
               (
                   min                                                              //BK: in Samules Version auskommentiert
                   (                                                                //BK: in Samules Version auskommentiert
                       pow(k_,1.5)/(epsilon_+epsilonMin_),
                       sqrt(k_)/((sqrt(6.0)*Cmu_*mag(symm(fvc::grad(U_))))*zeta_)   //BK: in Samules Version auskommentiert
                   ),                                                               //BK: in Samules Version auskommentiert
                   CEta_*pow( (pow(nu(),3)/(epsilon_+epsilonMin_)),0.25)
               );
*/

    //Info<<min(rho_)<<endl;
    //Lt_ = CL_ * pow(k_,1.5)/epsilon_;
    //Lk_ = CL_ * CEta_ * pow((pow((mu()/rho_),3)/epsilon_),0.25);
       
    Info<<"L.org without realizability constraint"<<endl;
    return CL_*max
               (
                   pow(k_,1.5)/epsilon_,
                   CEta_*pow((pow((mu()/rho_),3)/epsilon_),0.25)
               );   

/*	       
    return CL_*max(pow(k_,1.5)/(epsilon_+epsilonMin_), pow((pow(nu(),3)/(epsilon_+epsilonMin_)),0.25));
*/
}

void VLES::calculateDelta()
{
    /*
    // ~~~~~~~~~~~~~~~~~~~~~~~ max(dx,dy,dz) ~~~~~~~~~~~~~~~~~~~~~~~ //
  
    Info<<"delta=max(dx,dy,dz)"<<endl;

    tmp<volScalarField> hmax
    (
        new volScalarField
        (
            IOobject
            (
                "hmax",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zrero", dimLength, 0.0)
        )
    );

    const cellList& cells = mesh_.cells();
    const vectorField& cellC = mesh_.cellCentres();
    const vectorField& faceC = mesh_.faceCentres();
    const vectorField faceN(mesh_.faceAreas()/mag(mesh_.faceAreas()));

    forAll(cells, cellI)
    {
        scalar deltaMaxTmp = 0.0;
        const labelList& cFaces = cells[cellI];
        const point& cc = cellC[cellI];

        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];
            const point& fc = faceC[faceI];
            const vector& n = faceN[faceI];

            scalar tmp = magSqr(n*(n & (fc - cc)));
            if (tmp > deltaMaxTmp)
            {
                deltaMaxTmp = tmp;
            }
        }

        hmax()[cellI] = 2.0 * sqrt(deltaMaxTmp);
    }

    delta_.internalField() = hmax();
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    */
    
    
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~ mean  ~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    /* 
    Info<<"delta = 1/3*(dx+dy+dz)"<<endl;
    
    tmp<volScalarField> hsum
    (
        new volScalarField
        (
            IOobject
            (
                "hsum",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zrero", dimLength, 0.0)
        )
    );

    const cellList& cells = mesh_.cells();
    const vectorField& cellC = mesh_.cellCentres();
    const vectorField& faceC = mesh_.faceCentres();
    const vectorField faceN(mesh_.faceAreas()/mag(mesh_.faceAreas()));

    forAll(cells, cellI)
    {
        const labelList& cFaces = cells[cellI];
        const point& cc = cellC[cellI];

        int ii=0;
        forAll(cFaces, cFaceI)
        {
            label faceI = cFaces[cFaceI];
            const point& fc = faceC[faceI];
            const vector& n = faceN[faceI];

            scalar tmp = magSqr(n*(n & (fc - cc)));
            hsum()[cellI] += sqrt(tmp);
            ii++;
        }
	Info<<ii<<endl;
    }

    delta_.internalField() = 1.0/3.0 * hsum();
    */
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    
    
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~ cube root ~~~~~~~~~~~~~~~~~~~~~~~~~ //
    
    Info<<"delta = mesh.V^(1/3)"<<endl;
    delta_.internalField() = pow(mesh_.V(), 1.0/3.0);
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
    
}
	
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

VLES::VLES
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

    ksgs_
    (
        IOobject
        (
            "ksgs",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("ksgs", dimensionSet(0, 2, -2, 0, 0, 0, 0), 1.0e-10)
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
    
    fk_
    (
        IOobject
        (
            "fk",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_, 0.1
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
    
    delta_
    (
        IOobject
        (
            "delta",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_, 
        dimensionedScalar("delta", dimLength, 1.0e-10)
    ),

    Lmeso_
    (
        IOobject
        (
            "Lmeso",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pow(k_, 1.5)/epsilon_
    ),

    Lkom_
    (
        IOobject
        (
            "Lkom",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        CEta_*pow( (pow(mu()/rho_,3)/epsilon_),0.25)
    ),

    Lt_
    (  
        IOobject
        (  
            "Lt",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("LtInitial", dimLength, 0.0)
    ),

    Lk_
    (  
        IOobject
        (  
            "Lk",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("LkInitial", dimLength, 0.0)
    ),
   
    yr_(mesh_),
    
    mTSmall_("mTSmall", dimless/dimTime, SMALL),
    zetaMin_("zetaMin", dimless, SMALL),
    fMin_("fMin", dimless/dimTime, SMALL),
    TscMin_("TscMin", dimTime, SMALL),
    LscMin_("LscMin", dimLength, SMALL)
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);
    bound(f_, fMin_);
    bound(zeta_, zetaMin_);
    
    //alphat_ = mut_/Prt_;
    //alphat_.correctBoundaryConditions();
    
    //mut_ = rho_*Cmu_*zeta_*sqr(fk_)*k_*Tau();
    //mut_.correctBoundaryConditions();
    
    calculateDelta();    

    printCoeffs();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<volSymmTensorField> VLES::R() const
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
            ((2.0/3.0)*I)*ksgs_ - (mut_/rho_) * twoSymm(fvc::grad(U_)), //BK: Scheint korrekt
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> VLES::devRhoReff() const
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


tmp<fvVectorMatrix> VLES::divDevRhoReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(muEff(), U)
      - fvc::div(muEff()*dev2(T(fvc::grad(U))))
    ); //BK: Scheint korrekt
}

bool VLES::read()
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


void VLES::correct()
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
    
    //volScalarField Ceps1_ = 1.4*(1.0+(0.012/zeta_));
    volScalarField Ceps1_ = 1.4*(1.0+0.045/sqrt(zeta_)); //BK: Änderung zum Original
    
    volScalarField T_ = Tau() + TscMin_;
    volScalarField L_ = L() + LscMin_;
    

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

    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    if (mesh_.changing())
    {
	calculateDelta();
    };

    Lmeso_ = pow(k_, 1.5)/epsilon_ + LscMin_; Info<<"modified Lmeso"<<endl;
    Lkom_ = CEta_ * pow((pow((mu()/rho_),3)/epsilon_), 0.25) + LscMin_;
    
    fk_ = max(min(pow(delta_/Lmeso_, 2.0/3.0), 1.0), 1.0e-5);
    volScalarField Fr_ = sqr(fk_);
    
    ksgs_ = fk_*k_; 

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    
    // Re-calculate viscosity
    mut_ = rho_*Cmu_*zeta_*k_*T_*Fr_; //BK: Dichte
    //mut_.correctBoundaryConditions();

    alphat_ = mut_/Prt_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
