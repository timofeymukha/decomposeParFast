/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::domainDecomposition

Description
    Automatic domain decomposition class for finite-volume meshes

SourceFiles
    domainDecomposition.C
    domainDecompositionDistribute.C
    domainDecompositionMesh.C

\*---------------------------------------------------------------------------*/

#ifndef domainDecomposition_H
#define domainDecomposition_H

#include "fvMesh.H"
#include "labelList.H"
#include "point.H"
#include "Time.H"
#include "volFields.H"

namespace Foam
{

// Forward Declarations
class decompositionModel;

/*---------------------------------------------------------------------------*\
                     Class domainDecomposition Declaration
\*---------------------------------------------------------------------------*/

class domainDecomposition
:
    public fvMesh
{
    // Private Data

        //- Optional: points at the facesInstance
        autoPtr<pointIOField> facesInstancePointsPtr_;

        //- Optional non-standard file for decomposeParDict
        const fileName decompDictFile_;

        //- Number of processors in decomposition
        label nProcs_;

        //- Is the decomposition data to be distributed for each processor
        bool distributed_;

        //- Processor label for each cell
        labelList cellToProc_;

        //- Labels of points for each processor
        labelListList procPointAddressing_;

        //- Labels of faces for each processor
        // Note: Face turning index is stored as the sign on addressing
        // Only the processor boundary faces are affected: if the sign of the
        // index is negative, the processor face is the reverse of the
        // original face. In order to do this properly, all face
        // indices will be incremented by 1 and the decremented as
        // necessary to avoid the problem of face number zero having no
        // sign.
        List<DynamicList<label>> procFaceAddressing_;

        //- Labels of cells for each processor
        labelListList procCellAddressing_;

        //- Sizes for processor mesh patches
        // Excludes inter-processor boundaries
        labelListList procPatchSize_;

        //- Start indices for processor patches
        // Excludes inter-processor boundaries
        labelListList procPatchStartIndex_;


        // Per inter-processor patch information

        //- Neighbour processor ID for inter-processor boundaries
        labelListList procNeighbourProcessors_;

        //- Sizes for inter-processor patches
        labelListList procProcessorPatchSize_;

        //- Start indices (in procFaceAddressing_) for inter-processor patches
        labelListList procProcessorPatchStartIndex_;

        //- Sub patch IDs for inter-processor patches
        List<labelListList> procProcessorPatchSubPatchIDs_;

        //- Sub patch sizes for inter-processor patches
        List<labelListList> procProcessorPatchSubPatchStarts_;


        // Work data for constructing processor meshes

        //- Lookups/zone markers prepared
        mutable bool procMeshBuildReady_ = false;

        //- Full-mesh point/face to processor-local index (-1 if unused).
        //  Allocated once, only touched entries are reset per processor.
        mutable labelList pointLookup_;
        mutable labelList faceLookup_;

        //- Single zone per point/face/cell (-1: none, -2: multiple)
        mutable labelList pointToZone_;
        mutable labelList faceToZone_;
        mutable labelList cellToZone_;


    // Private Member Functions

        void distributeCells();

        //- Mark all elements with value or -2 if occur twice
        static void mark
        (
            const labelList& zoneElems,
            const label zoneI,
            labelList& elementToZone
        );

        //- Add face to inter-processor patch
        void addInterProcFace
        (
            const label facei,
            const label ownerProc,
            const label nbrProc,

            List<Map<label>>&,
            List<DynamicList<DynamicList<label>>>&
        ) const;

        //- Generate sub patch info for processor cyclics
        template<class BinaryOp>
        void processInterCyclics
        (
            const polyBoundaryMesh& patches,
            List<DynamicList<DynamicList<label>>>& interPatchFaces,
            List<Map<label>>& procNbrToInterPatch,
            List<labelListList>& subPatchIDs,
            List<labelListList>& subPatchStarts,
            bool owner,
            BinaryOp bop
        ) const;

        //- Prepare lookups and zone markers (once)
        void initProcMeshBuild() const;

        //- Reset the lookup entries touched by processor proci
        void resetProcMeshBuild(const label proci) const;

        //- Construct processor mesh (points, faces, cells, patches, zones).
        //  Leaves the lookups set for proci; call resetProcMeshBuild after.
        template<class MeshType>
        autoPtr<MeshType> buildProcMesh
        (
            const label proci,
            const Time& processorDb
        ) const;


public:

    // Constructors

        //- Construct from components.
        //  \param io the IOobject for mesh
        //  \param decompDictFile optional non-standard location for the
        //      decomposeParDict file
        explicit domainDecomposition
        (
            const IOobject& io,
            const fileName& decompDictFile = ""
        );


    //- Destructor
    ~domainDecomposition() = default;


    // Member Functions

        //- The mesh
        const fvMesh& mesh() const noexcept
        {
            return *this;
        }


    // Settings

        //- Number of processor in decomposition
        label nProcs() const noexcept
        {
            return nProcs_;
        }

        //- Is decomposition data to be distributed for each processor
        bool distributed() const noexcept
        {
            return distributed_;
        }

        //- Change distributed flag
        bool distributed(const bool on) noexcept
        {
            bool old(distributed_);
            distributed_ = on;
            return old;
        }

        //- Return decomposition model used
        const decompositionModel& model() const;

        //- Update flags based on the decomposition model settings
        //  Sets "distributed"
        void updateParameters(const dictionary& params);


    // Mappings

        //- Cell-processor decomposition labels
        const labelList& cellToProc() const noexcept
        {
            return cellToProc_;
        }


    // Decompose

        //- Decompose mesh
        void decomposeMesh();

        //- Write decomposition
        bool writeDecomposition(const bool decomposeSets);

        //- True if the processor meshes can be rebuilt in memory
        //- (decomposition done, mesh points not from a different instance)
        bool canBuildProcMeshes() const noexcept
        {
            return
            (
                !facesInstancePointsPtr_
             && procCellAddressing_.size() == nProcs_
             && procPatchSize_.size() == nProcs_
             && (nProcs_ == 0 || procPatchSize_[0].size())
            );
        }

        //- In-memory processor fvMesh, identical to the one written by
        //- writeDecomposition (avoids reading it back from disk)
        autoPtr<fvMesh> processorFvMesh
        (
            const label proci,
            const Time& processorDb
        ) const;

        //- Point addressing (processor to complete mesh) for proci
        const labelList& procPointAddressing(const label proci) const
        {
            return procPointAddressing_[proci];
        }

        //- Face addressing (with turning index) for proci
        const labelUList& procFaceAddressing(const label proci) const
        {
            return procFaceAddressing_[proci];
        }

        //- Cell addressing for proci
        const labelList& procCellAddressing(const label proci) const
        {
            return procCellAddressing_[proci];
        }

        //- Boundary addressing for proci (identity, -1 for processor patches)
        labelList procBoundaryAddressing(const label proci) const;


    // Write

        //- Write decomposition as volScalarField for visualization
        void writeVolField(const word& timeName) const;

        //- Write cell distribution as VTU file (binary)
        void writeVTK(const fileName& file) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "domainDecompositionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
