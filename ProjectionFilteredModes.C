/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 AUTHOR,AFFILIATION
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

Application
    ProjectionFilteringModes.C

Description

\*---------------------------------------------------------------------------*/

// #include "fvCFD.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "primitiveFields.H"
#include "FieldFields.H"
#include "scalarMatrices.H"
#include "SortableList.H"
#include "volFieldsFwd.H"
#include "forces.H"
#include "forceCoeffs.H"
#include "volFields.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
#include "ITHACAstream.H"
#include "ITHACAPOD.H"
#include "ITHACAparameters.H"
#include "ITHACAutilities.H"
//#include "fvPatch.H"
//#include "Modes.H"
#include "FilteredProblem.H"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/SparseLU>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

 //  FilteredProblem FP;   
  // PtrList<volVectorField> Vmodes;
  // PtrList<volScalarField> Smodes;
  //
  ITHACAparameters* para = ITHACAparameters::getInstance(mesh,runTime);


  #include "CheckPOD.H"

    label nSnapshots;
    label startTime;
    label endTime;
    
    IOdictionary ITHACAPODdict
    (
        IOobject
        (
            "ITHACAPODdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    List<instant> Times = runTime.times();


 // PtrList<volVectorField>& Vmodes = FP.Umodes;
 //  PtrList<volScalarField>& Smodes = FP.Hemodes;

  #include "Manager.H"
  #include "POD_part.H"
  
  FilteredProblem
   FP(
     Vfield,
     Sfield,
     Vmodes,
     Smodes
   );

  scalar vis_nu = 5e4;

  FP.closureViscosity = vis_nu;
  
  FP.FilteringStepU();
  FP.FilteringStepHe();
  
  FP.WriteFilteredSolutions();
  
  
}


// ************************************************************************* //
