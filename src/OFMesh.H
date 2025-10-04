/*---------------------------------------------------------------------------*\

    pyOFM  : Python interface for OpenFOAM mesh
    Version : v1.2

    Description:
        Functions to read OpenFOAM mesh information

\*---------------------------------------------------------------------------*/

#ifndef OFMesh_H
#define OFMesh_H

#include "fvCFD.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class OFMesh Declaration
\*---------------------------------------------------------------------------*/

class OFMesh
{

private:
    /// Disallow default bitwise copy construct
    OFMesh(const OFMesh&);

    /// Disallow default bitwise assignment
    void operator=(const OFMesh&);

    /// all the arguments
    char* argsAll_;

    autoPtr<fvMesh> meshPtr_;

    autoPtr<Time> runTimePtr_;

    autoPtr<argList> argsPtr_;

    label nLocalPoints_;

    label nLocalCells_;

    label nLocalFaces_;

    label nLocalInternalFaces_;

    label nLocalBoundaryPatches_;

    pointField pointField_;

public:
    /// Constructors
    OFMesh(char* argsAll);

    /// Destructor
    virtual ~OFMesh();

    void readMesh();

    double getMeshPointCoord(
        const label pointI,
        const label compI) const
    {
#if defined(CODI_ADR) || defined(CODI_ADF)
        return meshPtr_().points()[pointI][compI].getValue();
#else
        return meshPtr_().points()[pointI][compI];
#endif
    }

    void setMeshPointCoord(
        const label pointI,
        const label compI,
        const double value)
    {
        pointField_[pointI][compI] = value;
    }

    label getMeshFacePointIndex(
        const label faceI,
        const label pointI) const
    {
        const label& facePointSize = meshPtr_().faces()[faceI].size();
        if (pointI > facePointSize - 1)
        {
            FatalErrorIn("getMeshFacePointIndex") << "Face " << faceI << " has " << facePointSize << " points."
                                                  << " while the require pointI is " << pointI << abort(FatalError);
        }
        return meshPtr_().faces()[faceI][pointI];
    }

    label getNLocalPoints() const
    {
        return nLocalPoints_;
    }

    label getNLocalCells() const
    {
        return nLocalCells_;
    }

    label getNLocalFaces() const
    {
        return nLocalFaces_;
    }

    label getNLocalInternalFaces() const
    {
        return nLocalInternalFaces_;
    }

    label getNFacePoints(const label faceI) const
    {
        return meshPtr_().faces()[faceI].size();
    }

    void writeMesh() const
    {
        meshPtr_().write();
    }

    void updateMesh()
    {
        meshPtr_().movePoints(pointField_);
    }

    label getNLocalBoundaryPatches() const
    {
        return nLocalBoundaryPatches_;
    }

    word getLocalBoundaryName(const label patchI) const
    {
        return meshPtr_().boundaryMesh()[patchI].name();
    }

    word getLocalBoundaryType(const label patchI) const
    {
        return meshPtr_().boundaryMesh()[patchI].type();
    }

    label getLocalBoundaryStartFace(const label patchI) const
    {
        return meshPtr_().boundaryMesh()[patchI].start();
    }

    label getLocalBoundaryNFaces(const label patchI) const
    {
        return meshPtr_().boundaryMesh()[patchI].size();
    }

    label getLocalFaceOwner(const label faceI) const
    {
        return meshPtr_().owner()[faceI];
    }

    label getLocalFaceNeighbour(const label faceI) const
    {
        return meshPtr_().neighbour()[faceI];
    }

    label getLocalBoundaryFaceOwner(
        const label patchI,
        const label faceI) const
    {
        const UList<label>& pFaceCells = meshPtr_().boundaryMesh()[patchI].faceCells();
        return pFaceCells[faceI];
    }

    void readField(
        const word fieldName,
        const word fieldType,
        const word timeName,
        double* field)
    {
        if (fieldType == "volScalarField")
        {
            volScalarField state(
                IOobject(
                    fieldName,
                    timeName,
                    meshPtr_(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false),
                meshPtr_());

            forAll(state, cellI)
            {
#if defined(CODI_ADR) || defined(CODI_ADF)
                field[cellI] = state[cellI].getValue();
#else
                field[cellI] = state[cellI];
#endif
            }
        }
        else if (fieldType == "volVectorField")
        {
            volVectorField state(
                IOobject(
                    fieldName,
                    timeName,
                    meshPtr_(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false),
                meshPtr_());

            label counterI = 0;
            forAll(state, cellI)
            {
                for (label i = 0; i < 3; i++)
                {
#if defined(CODI_ADR) || defined(CODI_ADF)
                    field[counterI] = state[cellI][i].getValue();
#else
                    field[counterI] = state[cellI][i];
#endif
                    counterI++;
                }
            }
        }
        else
        {
            FatalErrorIn("readField") << "fieldType " << fieldType
                                      << " not supported! Options are volScalariField or volVectorField" << abort(FatalError);
        }
    }

    void writeField(
        const word fieldName,
        const word fieldType,
        const double* field)
    {
        if (fieldType == "volScalarField")
        {
            volScalarField state(
                IOobject(
                    fieldName,
                    "0",
                    meshPtr_(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false),
                meshPtr_(),
                dimensionedScalar(fieldName, dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
                "zeroGradient");

            forAll(state, cellI)
            {
                state[cellI] = field[cellI];
            }
            state.correctBoundaryConditions();
            state.write();
        }
        else if (fieldType == "volVectorField")
        {
            volVectorField state(
                IOobject(
                    fieldName,
                    "0",
                    meshPtr_(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false),
                meshPtr_(),
                dimensionedVector(fieldName, dimensionSet(0, 0, 0, 0, 0, 0, 0), vector::zero),
                "zeroGradient");

            label counterI = 0;
            forAll(state, cellI)
            {
                for (label i = 0; i < 3; i++)
                {
                    state[cellI][i] = field[counterI];
                    counterI++;
                }
            }
            state.correctBoundaryConditions();
            state.write();
        }
        else
        {
            FatalErrorIn("readField") << "fieldType " << fieldType
                                      << " not supported! Options are volScalariField or volVectorField" << abort(FatalError);
        }
    }

    // ------------------------------------------------------------------
    // Added geometry accessor methods to expose OpenFOAM cached geometric
    // primitives for golden reference generation in py-ofmesh project.
    // These simply forward to fvMesh on-demand computed fields. Any
    // deviation from OpenFOAM algorithms would happen there; here we only
    // expose read-only scalars.
    // ------------------------------------------------------------------

    double getFaceCentre(const label faceI, const label compI) const
    {
        // mesh.Cf(): face centres (vectorField)
        return meshPtr_().Cf()[faceI][compI];
    }

    double getFaceAreaVector(const label faceI, const label compI) const
    {
        // mesh.Sf(): face area vectors (vectorField)
        return meshPtr_().Sf()[faceI][compI];
    }

    double getFaceAreaMag(const label faceI) const
    {
        // NOTE (py-ofmesh parity fix): We deliberately avoid using mesh.magSf()[faceI]
        // here because in our golden-reference extraction workflow we observed
        // inconsistent/stale magnitudes for boundary faces (internal faces were fine).
        // Instead we derive the magnitude directly from the area vector Sf to
        // guarantee consistency: mag(Sf) == ||Sf||.
        //
        // OpenFOAM reference implementation for face centres and area vectors
        // (and implicitly the magnitudes) lives in:
        //   ./submodules/openfoam/src/OpenFOAM/meshes/primitiveMesh/
        //       primitiveMeshGeometry.C
        //   ./submodules/openfoam/src/OpenFOAM/meshes/primitiveMesh/
        //       primitiveMeshTools.C (makeFaceCentresAndAreas)
        // which computes both Sf and magSf from the same underlying data.
        // If magSf is ever guaranteed fresh this can revert to magSf(). Until
        // then we enforce correctness deterministically.
        const vector& sf = meshPtr_().Sf()[faceI];
        return Foam::mag(sf);
    }

    double getCellCentre(const label cellI, const label compI) const
    {
        // mesh.C(): cell centres (vectorField)
        return meshPtr_().C()[cellI][compI];
    }

    double getCellVolume(const label cellI) const
    {
        // mesh.V(): cell volumes (scalarField)
        return meshPtr_().V()[cellI];
    }

    // ------------------------------------------------------------------
    // Diagnostic: expose face pyramid volumes (primitiveMeshTools::facePyramidVolume)
    // Returns owner pyramid volumes (size nFaces) and neighbour pyramid volumes
    // (size nInternalFaces). Matches OpenFOAM sign convention in that both arrays
    // store positive magnitudes: primitiveMeshTools negates owner internally when
    // forming the scalarField for checks, but the underlying mag() call produces
    // positive volume. We reproduce the direct output of facePyramidVolume here.
    // ------------------------------------------------------------------
    void getFacePyramidVolumes(scalarField& ownPyr, scalarField& neiPyr) const;
    // Convenience helper returning std::vector for Cython (avoids exposing Foam::scalarField in Cython layer)
    void getFacePyramidVolumesStd(std::vector<double>& ownPyr, std::vector<double>& neiPyr) const;

    // ------------------------------------------------------------------
    // Additional checkMesh primitive accessors (augmentation for py-ofmesh parity work)
    // ------------------------------------------------------------------
    // Expose face centres array pointer (component-wise) not needed; keep scalar access pattern.

    // Cell determinant (primitiveMeshTools::cellDeterminant) - compute on demand.
    // We return a std::vector<double> for Python consumption.
    void getCellDeterminantStd(std::vector<double>& det) const;

    // Min pyramid volume written by checkMesh (already available on disk as volScalarField
    // when -writeAllFields used). Provided here for direct in-memory capture to remove file IO.
    bool getMinPyrVolumeStd(std::vector<double>& minPyr) const;

    // Min tet volume (if written). Returns false if field absent.
    bool getMinTetVolumeStd(std::vector<double>& minTet) const;

    // Wrong oriented faces: replicate logic in primitiveMeshGeometry::checkFacePyramids
    // We collect indices of faces with wrong orientation into vector. Returns count.
    int getWrongOrientedFacesStd(std::vector<int>& faceIds) const;

    // Cell negative / inverted flag: detect any cell with negative volume or negative tet volume
    // (simple heuristic) and return mask (0/1).
    void getCellInvertedMaskStd(std::vector<int>& inverted) const;

    // Bulk dump convenience: gather a suite of diagnostics in one call. Any optional field
    // not found is returned empty.
    void dumpDiagnostics(
        std::vector<double>& ownPyr3,
        std::vector<double>& neiPyr3,
        std::vector<double>& minPyr,
        std::vector<double>& minTet,
        std::vector<double>& cellDet,
        std::vector<int>& wrongFaces,
        std::vector<int>& invertedMask) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
