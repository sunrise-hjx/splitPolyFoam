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
#include "polyModifyFace.H"
#include "processorMeshes.H"
#include "polyModifyPoint.H"
#include "pointBoundaryMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void MYgetFaceInfo
(
    const label facei,
    label& patchID,
    label& zoneID,
    label& zoneFlip,
    const polyMesh& mesh
)
{
    patchID = -1;

    if (!mesh.isInternalFace(facei))
    {
        patchID = mesh.boundaryMesh().whichPatch(facei);
    }

    zoneID = mesh.faceZones().whichZone(facei);

    zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
    }
}

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

    const word oldInstance = mesh.pointsInstance();
    const bool overwrite = args.found("overwrite");
    const cellList& cells = mesh.cells();
    const labelListList& faceEdges = mesh.faceEdges();
    const vectorField &faceCentres = mesh.faceCentres();

    labelList polyFaceLabel;
    forAll(faceEdges,faceI)
    {
        if(faceEdges[faceI].size() > 4)
        {
            polyFaceLabel.append(faceI);
        }
    }

    label polyFaceNum = polyFaceLabel.size();

    Map<point> faceToDecompose(polyFaceNum);
    Map<label> faceAddedPoint(polyFaceNum);
    forAll(polyFaceLabel, i)
    {
        label faceI = polyFaceLabel[i];
        vector faceC = faceCentres[faceI];
        faceToDecompose.insert(faceI, faceC);
    }
    
    Info << nl << "Start triangulating polygonal faces !" << nl;

    // Topo change container
    polyTopoChange meshMod(mesh);
    forAllConstIters(faceToDecompose, iter)
    {
        const label facei = iter.key();
        const face& f = mesh.faces()[facei];

        label addedPointi =
            meshMod.setAction
            (
                polyAddPoint
                (
                    iter.val(), // point
                    f[0],   // master point
                    -1,     // zone for point
                    true    // supports a cell
                )
            );
        faceAddedPoint.insert(facei, addedPointi);
    }

    forAllConstIters(faceAddedPoint, iter)
    {
        const label facei = iter.key();
        const label addedPointi = iter.val();

        // Get face with new points on cut edges.
        face newFace(mesh.faces()[facei]);

        // Information about old face
        label patchID, zoneID, zoneFlip;
        MYgetFaceInfo(facei, patchID, zoneID, zoneFlip, mesh);
        label own = mesh.faceOwner()[facei];
        label nei = -1;
        if(mesh.isInternalFace(facei)) { nei = mesh.faceNeighbour()[facei]; }
        label masterPoint = mesh.faces()[facei][0];

        // Triangulate face around mid point
        face tri(3);
        forAll(newFace, fp)
        {
            label nextV = newFace.nextLabel(fp);

            tri[0] = newFace[fp];
            tri[1] = nextV;
            tri[2] = addedPointi;

            if (fp == 0)
            {
                // Modify the existing face.
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        tri,                        // face
                        facei,
                        own,                        // owner
                        nei,                        // neighbour
                        false,                      // flux flip
                        patchID,                    // patch for face
                        false,                      // remove from zone
                        zoneID,                     // zone for face
                        zoneFlip                    // face zone flip
                    )
                );
            }
            else
            {
                // Add additional faces
                meshMod.setAction
                (
                    polyAddFace
                    (
                        tri,                        // face
                        own,                        // owner
                        nei,                         // neighbour
                        masterPoint,                // master point
                        -1,                         // master edge
                        -1,                         // master face for addition
                        false,                      // flux flip
                        patchID,                    // patch for face
                        zoneID,                     // zone for face
                        zoneFlip                    // face zone flip
                    )
                );
            }
        }
    }

    // Do changes
    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);
    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }
    const mapPolyMesh &morphMapNew = morphMap;

    Map<label> newAddedPoints(faceAddedPoint.size());
    forAllConstIters(faceAddedPoint, iter)
    {
        const label oldFacei = iter.key();
        const label oldPointi = iter.val();

        const label newFacei = morphMapNew.reverseFaceMap()[oldFacei];
        const label newPointi = morphMapNew.reversePointMap()[oldPointi];

        if (newFacei >= 0 && newPointi >= 0)
        {
            newAddedPoints.insert(newFacei, newPointi);
        }
    }
    faceAddedPoint.transfer(newAddedPoints);

    if (!overwrite)
    {
        ++runTime;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    // Write resulting mesh
    Info<< "Writing modified mesh to time " << runTime.timeName() << endl;
    mesh.write();
    topoSet::removeFiles(mesh);
    processorMeshes::removeFiles(mesh);

    Info<< "\nEnd\n" << endl;
    return 0;
}

// ************************************************************************* //