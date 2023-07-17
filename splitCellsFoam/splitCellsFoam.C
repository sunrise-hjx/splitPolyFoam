/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "Time.H"
#include "Pair.H"
#include "fvCFD.H"
#include "scalar.H"
#include "stdFoam.H"
#include "argList.H"
#include "topoSet.H"
#include "polyMesh.H"
#include "cellShape.H"
#include "meshTools.H"
#include "pyrMatcher.H"
#include "tetMatcher.H"
#include "hexMatcher.H"
#include "polyAddCell.H"
#include "polyAddFace.H"
#include "globalIndex.H"
#include "mapPolyMesh.H"
#include "prismMatcher.H"
#include "IOdictionary.H"
#include "polyAddPoint.H"
#include "edgeCollapser.H"
#include "polyTopoChange.H"
#include "boundaryCutter.H"
#include "MYcellSplitter.H"
#include "polyModifyFace.H"
#include "processorMeshes.H"
#include "polyModifyPoint.H"
#include "pointBoundaryMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Manipulate mesh elements.\n"
        "For example, moving points, splitting/collapsing edges etc."
    );
    #include "addOverwriteOption.H"
    argList::addOption("dict", "file", "Use alternative modifyMeshDict");

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"
    #include "createTopoMesh.H"

    const word oldInstance = mesh.pointsInstance();
    const bool overwrite = args.found("overwrite");
    const cellList& cells = mesh.cells();
    const labelListList& faceEdges = mesh.faceEdges();
    const vectorField &faceCentres = mesh.faceCentres();

    labelList polyLabel;

    pyrMatcher pyr;
    tetMatcher tet;
    hexMatcher hex;
    prismMatcher prism;
    cellShape shape;
    forAll(mesh.cells(), cellI)
    {
        bool isHex = hex.matches(mesh, cellI, shape);
        bool isTet = tet.matches(mesh, cellI, shape);
        bool isPyr = pyr.matches(mesh, cellI, shape);
        bool isPrism = prism.matches(mesh, cellI, shape);
        if (!isHex && !isTet && !isPyr && !isPrism)
        {
            polyLabel.append(cellI);
        }
    }

    label polyNum = polyLabel.size();

    Map<point> cellToPyrCentre(polyNum);
    forAll(polyLabel,i)
    {
        label cellI = polyLabel[i];
        vector cellC = topoMesh.C()[cellI];
        cellToPyrCentre.insert(cellI, cellC);
    }

    Info<< nl << "Start decomposing the polyhedron !" << endl;

    // Mesh change engine
    MYcellSplitter cutter(mesh);

    // Topo change container
    polyTopoChange meshMod(mesh);

    // Insert commands into meshMod
    cutter.setRefinement(cellToPyrCentre, meshMod);

    // Do changes
    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }

    cutter.updateMesh(morphMap());

    if (!overwrite)
    {
        ++runTime;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    // Write resulting mesh
    Info << "Writing modified mesh to time " << runTime.timeName() << endl;
    mesh.write();
    topoSet::removeFiles(mesh);
    processorMeshes::removeFiles(mesh);
    
    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
