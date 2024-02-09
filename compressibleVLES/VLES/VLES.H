/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

Class
    Foam::incompressible::RASModels::VLES

Description
    Zeta-effe turbulence model for incompressible flows.
    
    Reference article: 

    "A robust near-wall elliptic-relaxation eddy-viscosity turbulence model for CFD"
    M. Popovac et al. (2004), International Journal of Heat and Fluid Flow 25 (2004)

    The default model coefficients correspond to the following:
    \verbatim
        VLESCoeffs
        {
             Cmu_         0.22,
             Ceps2_       1.90,
             C1_          0.40,
             C2_          0.65,
             sigmaK_      1.00,
             sigmaEps_    1.30,
             sigmaZ_      1.20,
             CTau_        6.00,
             CL_          0.36,
             CEta_        85.0
        }
    \endverbatim

SourceFiles
    VLES.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleVLES_H
#define compressibleVLES_H

#include "RASModel.H"
#include "wallDistReflection.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class VLES Declaration
\*---------------------------------------------------------------------------*/

class VLES
:
    public RASModel
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar Ceps2_;
            dimensionedScalar Ceps3_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaZ_;
            dimensionedScalar CTau_;
            dimensionedScalar CL_;
            dimensionedScalar CEta_;
            dimensionedScalar a_;
            dimensionedScalar Prt_;

        // Fields

            volScalarField k_;
            volScalarField ksgs_;
            volScalarField epsilon_;
            volScalarField zeta_;
            volScalarField f_;
            volScalarField fk_;
            volScalarField mut_;
            volScalarField alphat_;

        // testFields

            volScalarField delta_; 
            volScalarField Lmeso_;
            volScalarField Lkom_;
            volScalarField Lt_;
            volScalarField Lk_;	

        // wall distance reflection

            wallDistReflection yr_;
       
        // local limiters

            dimensionedScalar mTSmall_;
            dimensionedScalar zetaMin_;
            dimensionedScalar fMin_;
            dimensionedScalar TscMin_;
            dimensionedScalar LscMin_;

    // protected member functions

            tmp<volScalarField> Tau() const;
            tmp<volScalarField> L(); //const;
            void calculateDelta();
    
public:

    //- Runtime type information
    TypeName("VLES");

    // Constructors

        //- Construct from components
        VLES
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const fluidThermo& thermophysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~VLES()
    {}


    // Member Functions
     
        //- Return the lower allowable limit for zeta (default: SMALL)
        const dimensionedScalar& zetaMin() const
        {
            return zetaMin_;
        }

        //- Allow zetaMin to be changed
        dimensionedScalar& zetaMin()
        {
            return zetaMin_;
        }

        //- Return the lower allowable limit for f (default: SMALL)
        const dimensionedScalar& fMin() const
        {
            return fMin_;
        }

        //- Allow fMin to be changed
        dimensionedScalar& fMin()
        {
            return fMin_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", mut_/sigmaK_ + mu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", mut_/sigmaEps_ + mu())
            );
        }

        //- Return the effective diffusivity for zeta
        tmp<volScalarField> DzetaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DzetaEff", mut_/sigmaZ_ + mu())
            );
        }
        
        //- Return the turbulence viscosity
        virtual tmp<volScalarField> mut() const
        {
            return mut_;
        }

        //- Return the turbulence thermal diffusivity
        virtual tmp<volScalarField> alphat() const
        {
            return alphat_;
        }
        
        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the v2-k ratio (zeta)
        virtual tmp<volScalarField> zeta() const
        {
            return zeta_;
        }

        virtual tmp<volScalarField> f() const
        {
            return f_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
