
# distutils: language = c++
# distutils: sources = OFMesh.C

'''

    DAFoam  : Python interface for OpenFOAM mesh
    Version : v1.2

    Description:
        Cython wrapper functions that call OpenFOAM libraries defined
        in the *.C and *.H files. The python naming convention is to 
        add "py" before the C++ class name

'''

from libcpp.vector cimport vector
from libcpp.string cimport string
cimport numpy as np

# declear cpp functions
cdef extern from "OFMesh.H" namespace "Foam":
    cppclass OFMesh:
        OFMesh(char*) except +
        void readMesh()
        double getMeshPointCoord(int,int)
        void setMeshPointCoord(int,int,double)
        int getNLocalPoints()
        int getNLocalCells()
        int getNLocalFaces()
        int getNLocalInternalFaces()
        int getNFacePoints(int)
        int getMeshFacePointIndex(int,int)
        void writeMesh()
        void updateMesh()
        int getNLocalBoundaryPatches()
        string getLocalBoundaryName(int)
        string getLocalBoundaryType(int)
        int getLocalBoundaryStartFace(int)
        int getLocalBoundaryNFaces(int)
        int getLocalFaceOwner(int)
        int getLocalBoundaryFaceOwner(int,int)
        int getLocalFaceNeighbour(int)
        void readField(char* , char *, char *, double *)
        void writeField(char* , char *, double *)
        double getFaceCentre(int,int)
        double getFaceAreaVector(int,int)
        double getFaceAreaMag(int)
        double getCellCentre(int,int)
        double getCellVolume(int)
        void getFacePyramidVolumesStd(vector[double]& , vector[double]&)
        void getCellDeterminantStd(vector[double]&)
        bint getMinPyrVolumeStd(vector[double]&)
        bint getMinTetVolumeStd(vector[double]&)
        int getWrongOrientedFacesStd(vector[int]&)
        void getCellInvertedMaskStd(vector[int]&)
        void dumpDiagnostics(vector[double]&, vector[double]&, vector[double]&, vector[double]&, vector[double]&, vector[int]&, vector[int]&)

# create python wrappers that call cpp functions
cdef class pyOFMesh:

    # define a class pointer for cpp functions
    cdef:
        OFMesh * _thisptr

    # initialize this class pointer with NULL
    def __cinit__(self):
        self._thisptr = NULL

    # deallocate the class pointer, and
    # make sure we don't have memory leak
    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    # point the class pointer to the cpp class constructor
    def __init__(self, argsAll):
        '''
        argsAll: string that contains all the arguments
        for running OpenFOAM solvers, including
        the name of the solver.

        For example, in OpenFOAM, if we run the following:

        mpirun -np 2 simpleFoam -parallel

        Then, the corresponding call in pySimpleFoam is:

        pySimpleFoam(b"simpleFoam -parallel")
        '''
        self._thisptr = new OFMesh(argsAll)
    
    # wrap all the other memeber functions in the cpp class
    def readMesh(self):
        self._thisptr.readMesh()
    
    def getMeshPointCoord(self, pointI, compI):
        return self._thisptr.getMeshPointCoord(pointI,compI)
    
    def setMeshPointCoord(self, pointI, compI, value):
        self._thisptr.setMeshPointCoord(pointI,compI,value)

    def getNLocalPoints(self):
        return self._thisptr.getNLocalPoints()
    
    def getNLocalCells(self):
        return self._thisptr.getNLocalCells()
    
    def getNLocalFaces(self):
        return self._thisptr.getNLocalFaces()
    
    def getNLocalInternalFaces(self):
        return self._thisptr.getNLocalInternalFaces()

    def getNFacePoints(self, faceI):
        return self._thisptr.getNFacePoints(faceI)
    
    def getMeshFacePointIndex(self, faceI, pointI):
        return self._thisptr.getMeshFacePointIndex(faceI,pointI)
    
    def writeMesh(self):
        self._thisptr.writeMesh()
    
    def updateMesh(self):
        self._thisptr.updateMesh()
    
    def getNLocalBoundaryPatches(self):
        return self._thisptr.getNLocalBoundaryPatches()
    
    def getLocalBoundaryName(self, patchI):
        return self._thisptr.getLocalBoundaryName(patchI)
    
    def getLocalBoundaryType(self, patchI):
        return self._thisptr.getLocalBoundaryType(patchI)
    
    def getLocalBoundaryStartFace(self, patchI):
        return self._thisptr.getLocalBoundaryStartFace(patchI)
    
    def getLocalBoundaryNFaces(self, patchI):
        return self._thisptr.getLocalBoundaryNFaces(patchI)

    def getLocalFaceOwner(self, faceI):
        return self._thisptr.getLocalFaceOwner(faceI)
    
    def getLocalBoundaryFaceOwner(self, patchI, faceI):
        return self._thisptr.getLocalBoundaryFaceOwner(patchI,faceI)
    
    def getLocalFaceNeighbour(self, faceI):
        return self._thisptr.getLocalFaceNeighbour(faceI)
    
    def readField(self, fieldName, fieldType, timeName, np.ndarray[double, ndim=1, mode="c"] field):
        if fieldType == "volScalarField":
            assert len(field) == self.getNLocalCells(), "invalid array size!"
        elif fieldType == "volVectorField":
            assert len(field) == self.getNLocalCells() * 3, "invalid array size!"
        else:
            print("fieldType invalid!")
            exit(1)
        
        cdef double *field_data = <double*>field.data
        self._thisptr.readField(fieldName.encode(), fieldType.encode(), timeName.encode(), field_data)
    
    def writeField(self, fieldName, fieldType, np.ndarray[double, ndim=1, mode="c"] field):
        if fieldType == "volScalarField":
            assert len(field) == self.getNLocalCells(), "invalid array size!"
        elif fieldType == "volVectorField":
            assert len(field) == self.getNLocalCells() * 3, "invalid array size!"
        else:
            print("fieldType invalid!")
            exit(1)
        
        cdef double *field_data = <double*>field.data
        self._thisptr.writeField(fieldName.encode(), fieldType.encode(), field_data)

    # ------------------------------------------------------------------
    # Geometry accessors (added for golden reference extraction)
    # ------------------------------------------------------------------
    def getFaceCentre(self, faceI, compI):
        return self._thisptr.getFaceCentre(faceI, compI)

    def getFaceAreaVector(self, faceI, compI):
        return self._thisptr.getFaceAreaVector(faceI, compI)

    def getFaceAreaMag(self, faceI):
        return self._thisptr.getFaceAreaMag(faceI)

    def getCellCentre(self, cellI, compI):
        return self._thisptr.getCellCentre(cellI, compI)

    def getCellVolume(self, cellI):
        return self._thisptr.getCellVolume(cellI)

    def getFacePyramidVolumes(self):
        cdef vector[double] own
        cdef vector[double] nei
        self._thisptr.getFacePyramidVolumesStd(own, nei)
        # Convert C++ vectors to numpy arrays (copies)
        import numpy as np
        a = np.empty(len(own), dtype=float)
        for i in range(len(own)):
            a[i] = own[i]
        b = np.empty(len(nei), dtype=float)
        for i in range(len(nei)):
            b[i] = nei[i]
        return a, b

    def getCellDeterminant(self):
        cdef vector[double] det
        self._thisptr.getCellDeterminantStd(det)
        import numpy as np
        arr = np.empty(len(det), dtype=float)
        for i in range(len(det)):
            arr[i] = det[i]
        return arr

    def getMinPyrVolume(self):
        cdef vector[double] v
        ok = self._thisptr.getMinPyrVolumeStd(v)
        if not ok:
            return None
        import numpy as np
        arr = np.empty(len(v), dtype=float)
        for i in range(len(v)):
            arr[i] = v[i]
        return arr

    def getMinTetVolume(self):
        cdef vector[double] v
        ok = self._thisptr.getMinTetVolumeStd(v)
        if not ok:
            return None
        import numpy as np
        arr = np.empty(len(v), dtype=float)
        for i in range(len(v)):
            arr[i] = v[i]
        return arr

    def getWrongOrientedFaces(self):
        cdef vector[int] f
        count = self._thisptr.getWrongOrientedFacesStd(f)
        import numpy as np
        arr = np.empty(len(f), dtype=np.int32)
        for i in range(len(f)):
            arr[i] = f[i]
        return arr

    def getCellInvertedMask(self):
        cdef vector[int] m
        self._thisptr.getCellInvertedMaskStd(m)
        import numpy as np
        arr = np.empty(len(m), dtype=np.int32)
        for i in range(len(m)):
            arr[i] = m[i]
        return arr

    def dumpDiagnostics(self):
        cdef vector[double] own
        cdef vector[double] nei
        cdef vector[double] minP
        cdef vector[double] minT
        cdef vector[double] det
        cdef vector[int] wrong
        cdef vector[int] inv
        self._thisptr.dumpDiagnostics(own, nei, minP, minT, det, wrong, inv)
        import numpy as np
        def to_np_d(v):
            a = np.empty(len(v), dtype=float)
            for i in range(len(v)):
                a[i] = v[i]
            return a
        def to_np_i(v):
            a = np.empty(len(v), dtype=np.int32)
            for i in range(len(v)):
                a[i] = v[i]
            return a
        return {
            'face_owner_pyr3': to_np_d(own),
            'face_nei_pyr3': to_np_d(nei),
            'cell_minPyrVolume': None if len(minP)==0 else to_np_d(minP),
            'cell_minTetVolume': None if len(minT)==0 else to_np_d(minT),
            'cell_determinant': to_np_d(det),
            'wrong_oriented_faces': to_np_i(wrong),
            'cell_inverted_mask': to_np_i(inv),
        }
