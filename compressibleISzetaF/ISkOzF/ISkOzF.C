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

#include "ISkOzF.H"
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

defineTypeNameAndDebug(ISkOzF, 0);
addToRunTimeSelectionTable(RASModel, ISkOzF, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> ISkOzF::Tau() const
{
    volScalarField Ti = k_/epsilon_;
    volScalarField Tk = TSwitch_*CTau_ * sqrt((mu()/rho_)/epsilon_);
    volScalarField TInd("TInd", pos(Ti-Tk));

    if (runTime_.outputTime())
    {
        TInd.write();
    }
    
    return max
           (
               k_/epsilon_,
               TSwitch_*CTau_*sqrt((mu()/rho_)/epsilon_)
           );
}

tmp<volScalarField> ISkOzF::L() const
{
    volScalarField Li = CL_ * pow(k_,1.5)/epsilon_;
    volScalarField Lk = CL_ * LSwitch_*CEta_ * pow(pow(mu()/rho_,3)/epsilon_,0.25);
    volScalarField LInd("LInd", pos(Li-Lk));

    if (runTime_.outputTime())
    {
        LInd.write();
    }

    if (LTypeSwitch_.value() == 1)
    {
        Info<<"USING 13*L_TAYLOR"<<endl;
        return 13.0*sqrt(10.0*mu()/rho_*k_/epsilon_);
    }
/*
    else if (LTypeSwitch_.value() == 2)
    {
        return CL_*max
                 (
                     pow(k_+kr_,1.5)/epsilon_,
                     LSwitch_*CEta_*pow(pow(mu()/rho_,3)/epsilon_,0.25)
                 );
    }
*/
    else
    {
        return CL_*max
                 (
                     pow(k_,1.5)/epsilon_,
                     LSwitch_*CEta_*pow(pow(mu()/rho_,3)/epsilon_,0.25)
                 );
    }
}

void ISkOzF::calculateDelta()
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
/*
void ISkOzF::writeAveragingProperties() const
{
    IOdictionary propsDict
    (
        IOobject
        (
            "VLESAveragingProperties",
            mesh_.time().timeName(),
            "uniform",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    propsDict.add("Dt", Dt_);
    propsDict.regIOobject::write();
}
*/	
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ISkOzF::ISkOzF
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
    sigmaCDv_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaCDv",
            coeffDict_,
            1.5
        )
    ),
    sigmaCDt_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaCDt",
            coeffDict_,
            4.0
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
    TSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "TSwitch",
            coeffDict_,
            1.0
        )
    ),
    LSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "LSwitch",
            coeffDict_,
            1.0
        )
    ),
    Csas_
    (   
        dimensioned<scalar>::lookupOrAddToDict
        (   
            "Csas",
            coeffDict_,
            3.0
        )
    ),
    CT2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CT2",
            coeffDict_,
            20.0
        )
    ),
    Clim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clim",
            coeffDict_,
            0.0
        )
    ),
    fwSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fwSwitch",
            coeffDict_,
            0.0
        )
    ), 
    fEqnPbykSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fEqnPbykSwitch",
            coeffDict_,
            0.0
        )
    ),
    DavidsonSwitch_
    (   
        dimensioned<scalar>::lookupOrAddToDict
        (   
            "DavidsonSwitch",
            coeffDict_,
            0.0
        )
    ),
    LTypeSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "LTypeSwitch",
            coeffDict_,
            0.0
        )
    ),
    SdiffSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "SdiffSwitch",
            coeffDict_,
            1.0
        )
    ),
    CDtboundSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDtboundSwitch",
            coeffDict_,
            0.0
        )
    ),
    Ceps1zetaSwitch_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1zetaSwitch",
            coeffDict_,
            1.0
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
/*
    kr_
    (
        IOobject
        (
            "kr",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("kr", dimensionSet(0, 2, -2, 0, 0, 0, 0), 1.0e-10)
    ),
*/
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
/* 
    epsilonr_
    (
        IOobject
        (
            "epsilonr",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("epsilonr", dimensionSet(0, 2, -3, 0, 0, 0, 0), 0.0)
    ),
*/
    omega_
    (
        IOobject
        (
            "omega",
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

    T1_
    (
        IOobject
        (
            "T1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("T1", dimless/sqr(dimTime), 0.0)
    ),

    T2_
    (
        IOobject
        (
            "T2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("T2", dimless/sqr(dimTime), 0.0)
    ),
/*
    VLESAveragingProperties_
    (
        IOobject
        (
            "VLESAveragingProperties",
            mesh_.time().timeName(),
            "uniform",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),

    Dt_("Dt", dimless, VLESAveragingProperties_.lookup("Dt")),

    UAvg_
    (
        IOobject
        (
            "UAvg",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("UAvg", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
    ), 
*/  
    zetaMin_("zetaMin", dimless, SMALL),
    fMin_("fMin", dimless/dimTime, SMALL)
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


tmp<volSymmTensorField> ISkOzF::R() const
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


tmp<volSymmTensorField> ISkOzF::devRhoReff() const
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


tmp<fvVectorMatrix> ISkOzF::divDevRhoReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(muEff(), U)
      - fvc::div(muEff()*dev2(T(fvc::grad(U))))
    ); //BK: Scheint korrekt
}

bool ISkOzF::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        Ceps2_.readIfPresent(coeffDict());
        Ceps3_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaK_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        sigmaCDv_.readIfPresent(coeffDict());
        sigmaCDt_.readIfPresent(coeffDict());
        sigmaZ_.readIfPresent(coeffDict());
        CTau_.readIfPresent(coeffDict());
        CL_.readIfPresent(coeffDict());
        CEta_.readIfPresent(coeffDict());
        a_.readIfPresent(coeffDict());
        TSwitch_.readIfPresent(coeffDict());
        LSwitch_.readIfPresent(coeffDict());
        Csas_.readIfPresent(coeffDict());
        CT2_.readIfPresent(coeffDict());
        Clim_.readIfPresent(coeffDict());
        fwSwitch_.readIfPresent(coeffDict());
        fEqnPbykSwitch_.readIfPresent(coeffDict());
        DavidsonSwitch_.readIfPresent(coeffDict());
        LTypeSwitch_.readIfPresent(coeffDict());
        SdiffSwitch_.readIfPresent(coeffDict());
        CDtboundSwitch_.readIfPresent(coeffDict());
        Ceps1zetaSwitch_.readIfPresent(coeffDict());
        Prt_.readIfPresent(coeffDict());
        
        return true;
    }
    else
    {
        return false;
    }
}


void ISkOzF::correct()
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
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ online averaging  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
/*
    dimensionedScalar dt = runTime_.deltaTValue();
    Dt_ += dt;

    Info<<"Dt_ is: "<<Dt_<<endl;

    dimensionedScalar alpha = (Dt_ - dt)/Dt_;
    dimensionedScalar beta  = dt/Dt_;

    //RAvg_ += sqr(UAvg_);
    
    UAvg_ = alpha*UAvg_ + beta*U_;

    //RAvg_ = alpha*RAvg_ + beta*sqr(U_) - sqr(UAvg_);

    kr_ = 0.5 * magSqr(U_-UAvg_);

    volTensorField graduPrime = fvc::grad(U_-UAvg_);
    volTensorField epsilonrTens_ = 2*mu()/rho_*(graduPrime.T() & graduPrime);
    epsilonr_ = 0.5*tr(epsilonrTens_);
*/
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    const volTensorField gradU(fvc::grad(U_));
    const volScalarField S2(2*magSqr(dev(symm(gradU)))); //BK: S2 = 2 Sij Sij - 2/3 Skk Skk; 2/3 k vernachlässigt
    const volScalarField G(GName(), mut_*S2);
    volScalarField Ceps1_ = 0.0*zeta_;
    if (Ceps1zetaSwitch_.value() == 1.0)
    {   
        Ceps1_ = 1.4*(1.0+(0.012/(zeta_+zetaMin_)));
        Info<<"Ceps1 with zeta"<<endl;
    }
    else if (Ceps1zetaSwitch_.value() == 0.0)
    {
        Ceps1_ = 1.4252+0.0*zeta_;
        Info<<"Ceps1 w/o zeta"<<endl;
    }
    //volScalarField Ceps1_ = 1.4252+0.0*zeta_;
    //Info<<"Ceps1_ = 1.4252"<<endl;

    volScalarField T_ = Tau();
    volScalarField L_ = L();
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SAS TERM MIT LIMITER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    //volScalarField S2(2.0*magSqr(symm(fvc::grad(U_))));
    dimensionedScalar Psaslim("Psaslim", dimensionSet(0, 0, -2, 0, 0), 0.0);
    T2_ = 3.0*k_*max(pow(omega_,-2.0)*(fvc::grad(omega_) & fvc::grad(omega_)), pow(k_,-2.0)*(fvc::grad(k_) & fvc::grad(k_)));
    //volScalarField L("L", 1/pow(0.09, 0.25)*max(10.0*pow(pow(nu(),3.0)/(k_*omega_),0.25), sqrt(k_)/omega_));
    //volScalarField Tstar("Tstar", 20.0*1.775*0.41*(1.0/0.09)*G*omega_/k_*pow(L/LvKdelta,0.5));

    volScalarField LvKdeltaTest("LvKdeltaTest", Clim_*0.11*sqrt(0.41*3.51/(0.8/0.09-0.44))*delta_);
    volScalarField Lderiv("Lderiv", mag(2.0*symm(fvc::grad(U_)))/mag(fvc::laplacian(U_)));
    Info<<"ATTENTION: Lderiv TEST WRITE!"<<endl;
    Lderiv.write();
    volScalarField T1lim("T1lim", 40.0*1.775*0.41*mag(2.0*symm(fvc::grad(U_)))*1.0/max(Lderiv, LvKdeltaTest)*sqrt(k_));
    volScalarField ISlim("ISlim", pos(Lderiv-LvKdeltaTest));
    volScalarField Psas("Psas", T2_*0.0);

    if (Clim_.value() == 0.0)
    {
        T1_ = 40.0*1.775*0.41*mag(fvc::laplacian(U_))*sqrt(k_);
        Psas = max(0.003*(T1_ - CT2_*T2_), Psaslim);
        Info<<"Psas limiter disabled"<<endl;
    }
    else
    {
        T1_ = 40.0*1.775*0.41*mag(2.0*symm(fvc::grad(U_)))*1.0/max(Lderiv, LvKdeltaTest)*sqrt(k_);
        Psas = max(0.003*(T1_ - CT2_*T2_), Psaslim);
        Info<<"Psas limiter enabled"<<endl;
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
   
    volScalarField CDv("CDv", (2.0/k_*mu()/sigmaCDv_*(fvc::grad(k_)&fvc::grad(omega_)))/omega_);
    volScalarField CDt("CDt", (2.0/k_*mut_/sigmaCDt_*(fvc::grad(k_)&fvc::grad(omega_)))/omega_);
    volScalarField CD("CD", CDv + CDt);
    if (CDtboundSwitch_.value() == 1.0)
    {
        CD=CDv+max(CDt,dimensionedScalar("0.0", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0));
        Info<<"Bounded turbulent cross diffusion"<<endl;
    }
 
    // omega equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(rho_, omega_)
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DepsilonEff(), omega_)
      ==
        (Ceps1_-1.0)*G/k_*omega_
      - fvm::SuSp(rho_*(Ceps2_-1.0)*omega_, omega_) //BK: Dichte
      + fvm::Sp(CD, omega_) // Cross diffusion
      + fvm::SuSp(SdiffSwitch_*(fvc::laplacian(DepsilonEff(), k_) - fvc::laplacian(DkEff(), k_))/k_, omega_) // Zero, if sigmaEps=sigmaK
      + Csas_*rho_*Psas //BK: Dichte
    );
    omegaEqn().relax();
    #include "BoundaryConditionOmega.H"
    solve(omegaEqn);
    bound(omega_, omegaMin_);
    //omega_.max(pol);

    epsilon_ = omega_ * k_;
    bound(epsilon_, epsilonMin_);

/*
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
      + Csas_*rho_*Psas*k_
      //- fvm::SuSp(((2.0/3.0)*Ceps1_ + Ceps3_)*rho_*divU, epsilon_) //BK: Kompressibilität?!
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
*/
    
    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rho_, k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(rho_*epsilon_/k_, k_) //BK: Dichte
      //- fvm::SuSp((2.0/3.0)*rho_*divU, k_) //BK: Kompressibilität?!
    );

    kEqn().relax();
    solve(kEqn);
    //bound(k_, kMin_);
    dimensionedScalar kSmall_("kSmall", dimensionSet(0, 2, -2, 0, 0, 0, 0), 1e-10);
    bound(k_, kSmall_);

    if (fwSwitch_.value() == 0.0)
    {
        Info<<"Using fTilda"<<endl;    
    
        T_ = Tau();
        L_ = L();    
        
        // f equation
        tmp<fvScalarMatrix> fEqn
        (
            - fvm::laplacian(rho_, f_)
         ==
            - fvm::Sp(rho_/sqr(L_),f_) //BK: Dichte
            - (rho_*C1_+C2_*G/epsilon_)*(zeta_ - 2.0/3.0)/(sqr(L_)*T_) //BK: Dichte
            + fEqnPbykSwitch_ * 1.0/sqr(L_)*2.0/3.0*G/k_
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
            rho_*fTilda //BK: Dichte
          - fvm::Sp(G/(k_), zeta_)
        );

        zetaEqn().relax();
        solve(zetaEqn);
        bound(zeta_, zetaMin_);
        zeta_ = min(zeta_, 2.0); //BK: Warum maximal 1.5?
    }
    else if (fwSwitch_.value() == 1.0)
    {
        Info<<"Using f"<<endl;
      
        T_ = Tau();
        L_ = L();    

        // f equation
        tmp<fvScalarMatrix> fEqn
        (
            - fvm::laplacian(rho_, f_)
         ==
            - fvm::Sp(rho_/sqr(L_),f_) //BK: Dichte
            - (rho_*C1_+C2_*G/epsilon_)*(zeta_ - 2.0/3.0)/(sqr(L_)*T_) //BK: Dichte
            + fEqnPbykSwitch_ * 1.0/sqr(L_)*2.0/3.0*G/k_
        );
    
        fEqn().relax();
        #include "BoundaryConditionf.H"
        solve(fEqn);
    
        // Zeta equation
        tmp<fvScalarMatrix> zetaEqn
        (
            fvm::ddt(rho_, zeta_)
          + fvm::div(phi_, zeta_)
          - fvm::laplacian(DzetaEff(), zeta_)
         ==
//            min(rho_*f_, (rho_*C1_+(C2_*G/(epsilon_)))*((zeta_ - 2.0/3.0))/T_)
            rho_*f_
          - fvm::Sp(G/(k_), zeta_)
        );

        zetaEqn().relax();
        solve(zetaEqn);
        bound(zeta_, zetaMin_);
        zeta_ = min(zeta_, 2.0); //BK: Warum maximal 1.5?
    }
    
    // Re-calculate viscosity
    if (DavidsonSwitch_.value() == 0.0)
    {
        mut_ = rho_*Cmu_*zeta_*k_*T_; //BK: Dichte
    }
    else if (DavidsonSwitch_.value() == 1.0)
    {
        mut_ = min(rho_*Cmu_*zeta_*k_*T_, rho_*0.09*sqr(k_)/epsilon_);
    }  
    else if (DavidsonSwitch_.value() == 2.0)
    {
        mut_ = rho_*0.09*sqr(k_)/epsilon_;
        Info<<"mut kE-mode"<<endl;
    } 
    //mut_.correctBoundaryConditions();

    alphat_ = mut_/Prt_;
    
    Info<<"k, min: "<<min(k_).value()<<" max: "<<max(k_).value()<<" average: "<<k_.weightedAverage(mesh_.V()).value()<<endl;
    Info<<"epsilon, min: "<<min(epsilon_).value()<<" max: "<<max(epsilon_).value()<<" average: "<<epsilon_.weightedAverage(mesh_.V()).value()<<endl;

    if(runTime_.outputTime())
    {  
        Psas.write();
//      writeAveragingProperties();
        ISlim.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
