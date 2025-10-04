/*---------------------------------------------------------------------------*\

    pyOFM  : Python interface for OpenFOAM mesh
    Version : v1.2

\*---------------------------------------------------------------------------*/
#include "OFMesh.H"
#include "primitiveMeshTools.H" // For primitiveMeshTools::facePyramidVolume
#include "cellQuality.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Constructors
OFMesh::OFMesh(
    char* argsAll)
    : meshPtr_(nullptr),
      runTimePtr_(nullptr),
      argsPtr_(nullptr)
{
    argsAll_ = argsAll;
}

OFMesh::~OFMesh()
{
}

void OFMesh::readMesh()
{
// read the mesh and compute some sizes
#include "setArgs.H"
#include "setRootCasePython.H"

    Info << "Reading the OpenFOAM mesh.." << endl;

    runTimePtr_.reset(
        new Time(
            Time::controlDictName,
            args));

    Time& runTime = runTimePtr_();

    word regionName = fvMesh::defaultRegion;
    meshPtr_.reset(
        new fvMesh(
            IOobject(
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ)));

    fvMesh& mesh = meshPtr_();

    nLocalPoints_ = mesh.nPoints();

    nLocalCells_ = mesh.nCells();

    nLocalFaces_ = mesh.nFaces();

    nLocalInternalFaces_ = mesh.nInternalFaces();

    // initialize pointField and assign its values based on the initial mesh
    pointField_.setSize(nLocalPoints_);

    forAll(pointField_, pointI)
    {
        for (label compI = 0; compI < 3; compI++)
        {
            pointField_[pointI][compI] = mesh.points()[pointI][compI];
        }
    }

    nLocalBoundaryPatches_ = 0;
    forAll(mesh.boundaryMesh(), patchI)
    {
        nLocalBoundaryPatches_++;
    }

    return;
}

void OFMesh::getFacePyramidVolumes(scalarField& ownPyr, scalarField& neiPyr) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getFacePyramidVolumes") << "Mesh not loaded" << abort(FatalError);
    }
    const fvMesh& mesh = meshPtr_();
    Foam::primitiveMeshTools::facePyramidVolume(
        mesh,
        mesh.points(),
        mesh.cellCentres(),
        ownPyr,
        neiPyr);
}

void OFMesh::getFacePyramidVolumesStd(std::vector<double>& ownPyr, std::vector<double>& neiPyr) const
{
    scalarField own;
    scalarField nei;
    getFacePyramidVolumes(own, nei);
    ownPyr.resize(own.size());
    for (size_t i = 0; i < own.size(); ++i)
    {
        ownPyr[i] = own[i];
    }
    neiPyr.resize(nei.size());
    for (size_t i = 0; i < nei.size(); ++i)
    {
        neiPyr[i] = nei[i];
    }
}

void OFMesh::getCellDeterminantStd(std::vector<double>& det) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getCellDeterminantStd") << "Mesh not loaded" << abort(FatalError);
    }
    const fvMesh& mesh = meshPtr_();
    const vectorField& Sf = mesh.Sf();
    const label nCells = mesh.nCells();
    det.assign(nCells, 0.0);
    // Parity-inspired reimplementation of primitiveMeshTools::cellDeterminant (3D case)
    // Reference: primitiveMeshTools.C (consulted externally). We compute area tensor
    // A = sum_{internal faces} (Sf/avgA) outer (Sf/avgA) where avgA is average |Sf| of internal faces.
    for (label c = 0; c < nCells; ++c)
    {
        const cell& cellFaces = mesh.cells()[c];
        scalar avgA = 0.0;
        label nInt = 0;
        forAll(cellFaces, i)
        {
            label f = cellFaces[i];
            if (f < mesh.nInternalFaces())
            {
                const vector& sf = Sf[f];
                avgA += mag(sf);
                ++nInt;
            }
        }
        if (nInt == 0 || avgA <= VSMALL)
        {
            det[c] = 0.0;
            continue;
        }
        avgA /= static_cast<scalar>(nInt);
        scalar a00 = 0, a01 = 0, a02 = 0, a11 = 0, a12 = 0, a22 = 0; // symmetric
        forAll(cellFaces, i)
        {
            label f = cellFaces[i];
            if (f < mesh.nInternalFaces())
            {
                const vector& sf = Sf[f];
                scalar sfx = sf.x() / avgA;
                scalar sfy = sf.y() / avgA;
                scalar sfz = sf.z() / avgA;
                a00 += sfx * sfx;
                a01 += sfx * sfy;
                a02 += sfx * sfz;
                a11 += sfy * sfy;
                a12 += sfy * sfz;
                a22 += sfz * sfz;
            }
        }
        // Determinant of symmetric matrix
        scalar detA = a00 * (a11 * a22 - a12 * a12) - a01 * (a01 * a22 - a12 * a02) + a02 * (a01 * a12 - a11 * a02);
        det[c] = std::fabs(detA) / 8.0; // normalization per OpenFOAM code path
    }
}

static bool readVolFieldIfPresent(const fvMesh& mesh, const word& name, scalarField& out)
{
    IOobject io(
        name,
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE);
    if (!io.typeHeaderOk<volScalarField>(true))
        return false; // also tests existence
    volScalarField fld(io, mesh);
    out.setSize(fld.internalField().size());
    forAll(out, i) out[i] = fld[i];
    return true;
}

bool OFMesh::getMinPyrVolumeStd(std::vector<double>& minPyr) const
{
    if (!meshPtr_.valid())
        return false;
    const fvMesh& mesh = meshPtr_();
    scalarField fld;
    if (!readVolFieldIfPresent(mesh, "minPyrVolume", fld))
        return false;
    minPyr.resize(fld.size());
    forAll(fld, i) minPyr[i] = fld[i];
    return true;
}

bool OFMesh::getMinTetVolumeStd(std::vector<double>& minTet) const
{
    if (!meshPtr_.valid())
        return false;
    const fvMesh& mesh = meshPtr_();
    scalarField fld;
    if (!readVolFieldIfPresent(mesh, "minTetVolume", fld))
        return false;
    minTet.resize(fld.size());
    forAll(fld, i) minTet[i] = fld[i];
    return true;
}

int OFMesh::getWrongOrientedFacesStd(std::vector<int>& faceIds) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getWrongOrientedFacesStd") << "Mesh not loaded" << abort(FatalError);
    }
    const fvMesh& mesh = meshPtr_();
    // Following primitiveMeshGeometry::checkFacePyramids logic (simplified): a face is wrong
    // oriented if its owner pyramid volume < -SMALL tolerance when using Sf·(fc - Cc)
    const vectorField& fc = mesh.Cf();
    const vectorField& sf = mesh.Sf();
    const vectorField& cc = mesh.C();
    faceIds.clear();
    forAll(mesh.faces(), fI)
    {
        label own = mesh.owner()[fI];
        const vector d = fc[fI] - cc[own];
        scalar pyr3 = sf[fI] & d; // 3 * volume (signed)
        if (pyr3 < -SMALL)
        {
            faceIds.push_back(fI);
        }
    }
    return static_cast<int>(faceIds.size());
}

void OFMesh::getCellInvertedMaskStd(std::vector<int>& inverted) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getCellInvertedMaskStd") << "Mesh not loaded" << abort(FatalError);
    }
    const fvMesh& mesh = meshPtr_();
    const scalarField& vols = mesh.V();
    inverted.assign(vols.size(), 0);
    forAll(vols, cI)
    {
        if (vols[cI] < 0)
            inverted[cI] = 1;
    }
}

void OFMesh::dumpDiagnostics(
    std::vector<double>& ownPyr3,
    std::vector<double>& neiPyr3,
    std::vector<double>& minPyr,
    std::vector<double>& minTet,
    std::vector<double>& cellDet,
    std::vector<int>& wrongFaces,
    std::vector<int>& invertedMask) const
{
    // Face pyramid volumes
    getFacePyramidVolumesStd(ownPyr3, neiPyr3);
    // Cell determinant
    getCellDeterminantStd(cellDet);
    // Optional fields (only filled if present)
    getMinPyrVolumeStd(minPyr);
    getMinTetVolumeStd(minTet);
    // Wrong orientation & inverted mask
    getWrongOrientedFacesStd(wrongFaces);
    getCellInvertedMaskStd(invertedMask);
}

void OFMesh::getOwnersStd(std::vector<int>& owners) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getOwnersStd") << "Mesh not loaded" << abort(FatalError);
    }
    const fvMesh& mesh = meshPtr_();
    const labelUList& own = mesh.owner();
    owners.resize(own.size());
    forAll(own, i) owners[i] = own[i];
}

void OFMesh::getNeighboursStd(std::vector<int>& neighbours) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getNeighboursStd") << "Mesh not loaded" << abort(FatalError);
    }
    const fvMesh& mesh = meshPtr_();
    const labelUList& nei = mesh.neighbour();
    neighbours.resize(nei.size());
    forAll(nei, i) neighbours[i] = nei[i];
}

void OFMesh::getFacePointsCsrStd(std::vector<int>& offsets, std::vector<int>& points) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getFacePointsCsrStd") << "Mesh not loaded" << abort(FatalError);
    }
    const fvMesh& mesh = meshPtr_();
    const faceList& fcs = mesh.faces();
    const label nF = fcs.size();
    offsets.resize(nF + 1);
    // First pass: prefix sums
    std::size_t running = 0;
    for (label f = 0; f < nF; ++f)
    {
        offsets[f] = static_cast<int>(running);
        running += fcs[f].size();
    }
    offsets[nF] = static_cast<int>(running);
    points.resize(running);
    // Second pass: copy point indices preserving original order
    for (label f = 0; f < nF; ++f)
    {
        const face& fc = fcs[f];
        const int base = offsets[f];
        forAll(fc, i)
        {
            points[base + i] = fc[i];
        }
    }
}

void OFMesh::getCellFacesCsrStd(std::vector<int>& offsets, std::vector<int>& faces) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getCellFacesCsrStd") << "Mesh not loaded" << abort(FatalError);
    }
    const fvMesh& mesh = meshPtr_();
    const cellList& cls = mesh.cells();
    const label nC = cls.size();
    offsets.resize(nC + 1);
    std::size_t running = 0;
    for (label c = 0; c < nC; ++c)
    {
        offsets[c] = static_cast<int>(running);
        running += cls[c].size();
    }
    offsets[nC] = static_cast<int>(running);
    faces.resize(running);
    for (label c = 0; c < nC; ++c)
    {
        const cell& cf = cls[c];
        const int base = offsets[c];
        forAll(cf, i)
        {
            faces[base + i] = cf[i];
        }
    }
}

// ------------------------------------------------------------------
// Bulk geometry exports (py-ofmesh parity diagnostics)
// ------------------------------------------------------------------
void OFMesh::getFaceAreaVectorsStd(std::vector<double>& data) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getFaceAreaVectorsStd") << "Mesh not loaded" << abort(FatalError);
    }
    const fvMesh& mesh = meshPtr_();
    const vectorField& Sf = mesh.Sf(); // primitiveMeshGeometry.C
    const label nF = Sf.size();
    data.resize(static_cast<std::size_t>(3 * nF));
    for (label f = 0; f < nF; ++f)
    {
        const vector& v = Sf[f];
        data[3 * f + 0] = v.x();
        data[3 * f + 1] = v.y();
        data[3 * f + 2] = v.z();
    }
}

void OFMesh::getFaceCentresStd(std::vector<double>& data) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getFaceCentresStd") << "Mesh not loaded" << abort(FatalError);
    }
    const fvMesh& mesh = meshPtr_();
    const vectorField& Cf = mesh.Cf(); // primitiveMeshGeometry.C
    const label nF = Cf.size();
    data.resize(static_cast<std::size_t>(3 * nF));
    for (label f = 0; f < nF; ++f)
    {
        const vector& v = Cf[f];
        data[3 * f + 0] = v.x();
        data[3 * f + 1] = v.y();
        data[3 * f + 2] = v.z();
    }
}

void OFMesh::getCellCentresStd(std::vector<double>& data) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getCellCentresStd") << "Mesh not loaded" << abort(FatalError);
    }
    const fvMesh& mesh = meshPtr_();
    const vectorField& Cc = mesh.C(); // primitiveMeshGeometry.C
    const label nC = Cc.size();
    data.resize(static_cast<std::size_t>(3 * nC));
    for (label c = 0; c < nC; ++c)
    {
        const vector& v = Cc[c];
        data[3 * c + 0] = v.x();
        data[3 * c + 1] = v.y();
        data[3 * c + 2] = v.z();
    }
}

void OFMesh::getCellVolumesStd(std::vector<double>& data) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getCellVolumesStd") << "Mesh not loaded" << abort(FatalError);
    }
    const fvMesh& mesh = meshPtr_();
    const scalarField& vols = mesh.V(); // primitiveMeshGeometry.C
    const label nC = vols.size();
    data.resize(static_cast<std::size_t>(nC));
    for (label c = 0; c < nC; ++c)
    {
        data[c] = vols[c];
    }
}

void OFMesh::getOwnerPyr3Std(std::vector<double>& data) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getOwnerPyr3Std") << "Mesh not loaded" << abort(FatalError);
    }
    // pyr3 = Sf · (Cf - C_owner) (signed 3 * pyramid volume from owner side)
    const fvMesh& mesh = meshPtr_();
    const vectorField& Sf = mesh.Sf();
    const vectorField& Cf = mesh.Cf();
    const vectorField& Cc = mesh.C();
    const labelUList& own = mesh.owner();
    const label nF = Sf.size();
    data.resize(static_cast<std::size_t>(nF));
    for (label f = 0; f < nF; ++f)
    {
        const vector d = Cf[f] - Cc[own[f]];
        data[f] = Sf[f] & d; // primitiveMeshTools::facePyramidVolume internal form
    }
}

void OFMesh::getNeighbourPyr3Std(std::vector<double>& data) const
{
    if (!meshPtr_.valid())
    {
        FatalErrorIn("OFMesh::getNeighbourPyr3Std") << "Mesh not loaded" << abort(FatalError);
    }
    // pyr3 (neighbour) = -Sf · (Cf - C_nei)
    const fvMesh& mesh = meshPtr_();
    const vectorField& Sf = mesh.Sf();
    const vectorField& Cf = mesh.Cf();
    const vectorField& Cc = mesh.C();
    const labelUList& nei = mesh.neighbour();
    const label nIF = mesh.nInternalFaces();
    data.resize(static_cast<std::size_t>(nIF));
    for (label f = 0; f < nIF; ++f)
    {
        const vector d = Cf[f] - Cc[nei[f]];
        data[f] = -(Sf[f] & d); // sign per neighbour definition
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
