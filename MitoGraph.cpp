// ==============================================================
// MitoGraph: Quantifying Mitochondrial Content in Living Cells
// Written by Matheus P. Viana - vianamp@gmail.com - 2014.05.28
//
// Susanne Rafelski Lab, University of California Irvine
//
// The official documentation will be soon available at
//
//      The official software website
//      - https://github.com/vianamp/MitoGraph
//
//      or
//
//      - http://www.rafelski.com/susanne/Home.html
//
// A protocol paper describing how to use MitoGraph is also
// coming soon.
// ==============================================================

#include <list>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <vtkMath.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkLongArray.h>
#include <vtkImageData.h>
#include <vtkDataArray.h>
#include <vtkTransform.h>
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkVectorText.h>
#include <vtkInformation.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPoints.h>
#include <vtkInformationVector.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageExtractComponents.h>
#include <vtkStructuredPointsReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkKdTreePointLocator.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkTIFFReader.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkImageRFFT.h>
#include <vtkImageCast.h>
#include <vtkPoints.h>

#include "MitoThinning.h"

    double _rad = 0.150;
    double _dxy, _dz = -1.0;
    bool _export_graph_files = true;
    bool _export_image_resampled = false;
    bool _adaptive_threshold = false;
    bool _scale_polydata_before_save = true;
    bool _export_nodes_label = true;
    double _div_threshold = 0.1666667;
    bool _checkonly = false;
    bool _width = false;
    bool _improve_skeleton_quality = true; // when this is true nodes with degree zero
                                           // expanded and detected. Additional checking
                                           // is also done to garantee that all non-zero
                                           // voxels were analysized.


    //               |------06------|
    //               |------------------------18------------------------|
    //               |---------------------------------------26----------------------------------|
    int ssdx_sort[26] = { 0,-1, 0, 1, 0, 0,-1, 0, 1, 0,-1, 1, 1,-1,-1, 0, 1, 0, -1, 1, 1,-1,-1, 1, 1,-1};
    int ssdy_sort[26] = { 0, 0,-1, 0, 1, 0, 0,-1, 0, 1,-1,-1, 1, 1, 0,-1, 0, 1, -1,-1, 1, 1,-1,-1, 1, 1};
    int ssdz_sort[26] = {-1, 0, 0, 0, 0, 1,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1, -1,-1,-1,-1, 1, 1, 1, 1};

// In order to acess the voxel (x,y,z) from ImageJ, I should use
// GetId(x,(Dim[1]-1)-y,z,Dim);
// Or change the volume orientation...

/**========================================================
 Auxiliar functions
 =========================================================*/


// This routine returns the x of the id-th point of a 3D volume
// of size Dim[0]xDim[1]xDim[2]
int  GetX(vtkIdType id, int *Dim);

// This routine returns the y of the id-th point of a 3D volume
// of size Dim[0]xDim[1]xDim[2]
int  GetY(vtkIdType id, int *Dim);

// This routine returns the z of the id-th point of a 3D volume
// of size Dim[0]xDim[1]xDim[2]
int  GetZ(vtkIdType id, int *Dim);

// This routine returns the id of a point located at coordinate
// (x,y,z) of a 3D volume of size Dim[0]xDim[1]xDim[2]
vtkIdType GetId(int x, int y, int z, int *Dim);

// Swap values
void Swap(double *x, double *y);

// Simple sorting algorithm and the output is such that
// l3 >= l2 >= l3
void Sort(double *l1, double *l2, double *l3);

// Calculate the Frobenius norm of a given 3x3 matrix
// http://mathworld.wolfram.com/FrobeniusNorm.html
double FrobeniusNorm(double M[3][3]);

// This routine scales the polydata points to the correct dimension
// given by parameters _dxy and _dz.
void ScalePolyData(vtkSmartPointer<vtkPolyData> PolyData);

/* ================================================================
   IMAGE TRANSFORM
=================================================================*/

// This routine converts 16-bit volumes into 8-bit volumes by
// linearly scaling the original range of intensities [min,max]
// in [0,255] (http://rsbweb.nih.gov/ij/docs/guide/146-28.html)
vtkSmartPointer<vtkImageData> Convert16To8bit(vtkSmartPointer<vtkImageData> Image);

// Apply a threshold to a ImageData and converts the result in
// a 8-bit ImageData.
vtkSmartPointer<vtkImageData> BinarizeAndConvertDoubleToChar(vtkSmartPointer<vtkImageData> Image, double threshold);

// Fill holes in the 3D image
void FillHoles(vtkSmartPointer<vtkImageData> ImageData);

/* ================================================================
   I/O ROUTINES
=================================================================*/

// Export maximum projection of a given ImageData as a
// PNG file.
void ExportMaxProjection(vtkSmartPointer<vtkImageData> Image, const char FileName[], bool binary);

// Export maximum projection of bottom and top parts of
// a given Tiff image as a PNG file as well as the polydata
// surface points
void ExportDetailedMaxProjection(const char FileName[]);

/* ================================================================
   ROUTINES FOR VESSELNESS CALCUATION VIA DISCRETE APPROCH
=================================================================*/

// This routine uses a discrete differential operator to
// calculate the derivatives of a given 3D volume
void GetImageDerivativeDiscrete(vtkSmartPointer<vtkDataArray> Image, int *dim, char direction, vtkSmartPointer<vtkFloatArray> Derivative);

// This routine calculate the Hessian matrix for each point
// of a 3D volume and its eigenvalues (Discrete Approach)
void GetHessianEigenvaluesDiscrete(double sigma, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkDoubleArray> L1, vtkSmartPointer<vtkDoubleArray> L2, vtkSmartPointer<vtkDoubleArray> L3);
void GetHessianEigenvaluesDiscreteZDependentThreshold(double sigma, vtkImageData *Image, vtkSmartPointer<vtkDoubleArray> L1, vtkSmartPointer<vtkDoubleArray> L2, vtkSmartPointer<vtkDoubleArray> L3);

// Calculate the vesselness at each point of a 3D volume based
// based on the Hessian eigenvalues
void GetVesselness(double sigma, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkDoubleArray> L1, vtkSmartPointer<vtkDoubleArray> L2, vtkSmartPointer<vtkDoubleArray> L3);

// Calculate the vesselness over a range of different scales
int MultiscaleVesselness(const char FileName[], double _sigmai, double _dsigma, double _sigmaf, double *attributes);

/* ================================================================
   DIVERGENCE FILTER
=================================================================*/

// This routine calculates the divergence filter of a 3D volume
// based on the orientation of the gradient vector field
void GetDivergenceFilter(int *Dim, vtkSmartPointer<vtkDoubleArray> Scalars);

/* ================================================================
   WIDTH ANALYSIS
=================================================================*/

// Approximation of tubule width by the distance of the skeleton
// to the closest point over the surface.
void EstimateTubuleWidth(vtkSmartPointer<vtkPolyData> Skeleton, vtkSmartPointer<vtkPolyData> Surface, double *attributes);

/* ================================================================
   INTENSITY MAPPING
=================================================================*/

// Intensities of the original TIFF image is mapped into a scalar
// component of the skeleton.
void MapImageIntensity(vtkSmartPointer<vtkPolyData> Skeleton, vtkSmartPointer<vtkImageData> ImageData, int nneigh);

/**========================================================
 Auxiliar functions
 =========================================================*/

int GetX(vtkIdType id, int *Dim) {
    return (int) ( id % (vtkIdType)Dim[0] );
}

int GetY(vtkIdType id, int *Dim) {
    return (int) (  ( id % (vtkIdType)(Dim[0]*Dim[1]) ) / (vtkIdType)Dim[0]  );
}

int GetZ(vtkIdType id, int *Dim) {
    return (int) ( id / (vtkIdType)(Dim[0]*Dim[1]) );
}

vtkIdType GetId(int x, int y, int z, int *Dim) {
    return (vtkIdType)(x+y*Dim[0]+z*Dim[0]*Dim[1]);
}

vtkIdType GetReflectedId(int x, int y, int z, int *Dim) {
    int rx = ceil(0.5*Dim[0]);
    int ry = ceil(0.5*Dim[1]);
    int rz = ceil(0.5*Dim[2]);
    int sx = (x-(rx-0.5)<0) ? -rx : Dim[0]-rx;
    int sy = (y-(ry-0.5)<0) ? -ry : Dim[1]-ry;
    int sz = (z-(rz-0.5)<0) ? -rz : Dim[2]-rz;
    return GetId(x-sx,y-sy,z-sz,Dim);
}

void Swap(double *x, double *y) {
    double t = *y; *y = *x; *x = t;
}

void Sort(double *l1, double *l2, double *l3) {
    if (fabs(*l1) > fabs(*l2)) Swap(l1,l2);
    if (fabs(*l2) > fabs(*l3)) Swap(l2,l3);
    if (fabs(*l1) > fabs(*l2)) Swap(l1,l2);
}

double FrobeniusNorm(double M[3][3]) {
    double f = 0.0;
    for (int i = 3;i--;)
        for (int j = 3;j--;)
            f += M[i][j]*M[i][j];
    return sqrt(f);
}

void ScalePolyData(vtkSmartPointer<vtkPolyData> PolyData) {
    double r[3];
    vtkPoints *Points = PolyData -> GetPoints();
    for (vtkIdType id = 0; id < Points -> GetNumberOfPoints(); id++) {
        Points -> GetPoint(id,r);
        Points -> SetPoint(id,_dxy*r[0],_dxy*r[1],_dz*r[2]);
    }
    Points -> Modified();
}

/* ================================================================
   I/O ROUTINES
=================================================================*/

void ExportMaxProjection(vtkSmartPointer<vtkImageData> Image, const char FileName[]) {

    #ifdef DEBUG
        printf("Saving Max projection...\n");
    #endif

    int *Dim = Image -> GetDimensions();
    vtkSmartPointer<vtkImageData> MaxP = vtkSmartPointer<vtkImageData>::New();
    MaxP -> SetDimensions(Dim[0],Dim[1],1);
    vtkIdType N = Dim[0] * Dim[1];

    vtkSmartPointer<vtkUnsignedCharArray> MaxPArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
    MaxPArray -> SetNumberOfComponents(1);
    MaxPArray -> SetNumberOfTuples(N);

    int x, y, z;
    double v, vproj;
    for (x = Dim[0]; x--;) {
        for (y = Dim[1]; y--;) {
            vproj = 0;
            for (z = Dim[2]; z--;) {
                v = Image -> GetScalarComponentAsFloat(x,y,z,0);
                vproj = (v > vproj) ? v : vproj;
            }
            MaxPArray -> SetTuple1(MaxP->FindPoint(x,y,0),(unsigned char)vproj);
        }
    }
    MaxPArray -> Modified();

    MaxP -> GetPointData() -> SetScalars(MaxPArray);

    vtkSmartPointer<vtkPNGWriter> PNGWriter = vtkSmartPointer<vtkPNGWriter>::New();
    PNGWriter -> SetFileName(FileName);
    PNGWriter -> SetFileDimensionality(2);
    PNGWriter -> SetCompressionLevel(0);
    PNGWriter -> SetInputData(MaxP);
    PNGWriter -> Write();

    #ifdef DEBUG
        printf("File Saved!\n");
    #endif

}

void ExportDetailedMaxProjection(const char FileName[]) {

    // =======================================================================================
    //                 |                |                 |                |                 |
    //  Total MaxProj  |      First     |       Top       |       Top      |       Top       |
    // (original tiff) |      Slice     |     8 slices    |     surface    |     skeleton    |
    //                 |                |                 |                |                 |
    // =======================================================================================
    //                 |                |                 |                |                 |
    //  Total MaxProj  |      Last      |      Bottom     |     Bottom     |     Bottom      |
    //    (surface)    |      Slice     |     8 slices    |     surface    |     skeleton    |
    //                 |                |                 |                |                 |
    // =======================================================================================

    #ifdef DEBUG
        printf("Saving Detailed Max projection...\n");
    #endif

    char _fullpath[256];
    sprintf(_fullpath,"%s.tif",FileName);

    vtkSmartPointer<vtkTIFFReader> TIFFReader = vtkSmartPointer<vtkTIFFReader>::New();
    int errlog = TIFFReader -> CanReadFile(_fullpath);
    // File cannot be opened
    if (!errlog) {

        printf("File %s connot be opened.\n",_fullpath);

    } else {

        // Loading TIFF File
        TIFFReader -> SetFileName(_fullpath);
        TIFFReader -> Update();
        vtkSmartPointer<vtkImageData> Image = TIFFReader -> GetOutput();

        //16bit to 8bit Conversion
        Image = Convert16To8bit(Image);

        // Loading PolyData Surface
        sprintf(_fullpath,"%s_surface.vtk",FileName);
        vtkSmartPointer<vtkPolyDataReader> PolyDaTaReader = vtkSmartPointer<vtkPolyDataReader>::New();
        PolyDaTaReader -> SetFileName(_fullpath);
        PolyDaTaReader -> Update();
        vtkSmartPointer<vtkPolyData> Surface = PolyDaTaReader -> GetOutput();

        // Loading Skeleton
        sprintf(_fullpath,"%s_skeleton.vtk",FileName);
        vtkSmartPointer<vtkPolyDataReader> PolyDaTaReaderSkell = vtkSmartPointer<vtkPolyDataReader>::New();
        PolyDaTaReaderSkell -> SetFileName(_fullpath);
        PolyDaTaReaderSkell -> Update();
        vtkSmartPointer<vtkPolyData> Skeleton = PolyDaTaReaderSkell -> GetOutput();

        #ifdef DEBUG
            printf("\t#Points [%s] = %d\n",_fullpath,(int)Surface->GetNumberOfPoints());
            printf("\t#Points [%s] = %d\n",_fullpath,(int)Skeleton->GetNumberOfPoints());
        #endif

        // Stack Dimensions
        int *Dim = Image -> GetDimensions();
        
        #ifdef DEBUG
            printf("Dim = [%d,%d,%d]\n",Dim[0],Dim[1],Dim[2]);
        #endif

        // Surface Bounds
        double *Bounds = Surface -> GetBounds();

        int zi = round(Bounds[4]/_dz); zi -= (zi>1) ? 1 : 0;
        int zf = round(Bounds[5]/_dz); zf += (zi<Dim[2]-1) ? 1 : 0;

        #ifdef DEBUG
            printf("Z = [%d,%d]\n",zi,zf);
        #endif

        // Plane
        vtkSmartPointer<vtkImageData> Plane = vtkSmartPointer<vtkImageData>::New();
        Plane -> SetDimensions(5*Dim[0],2*Dim[1],1);
        vtkIdType N = 10 * Dim[0] * Dim[1];

        // Scalar VEctor
        vtkSmartPointer<vtkUnsignedCharArray> MaxPArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
        MaxPArray -> SetNumberOfComponents(1);
        MaxPArray -> SetNumberOfTuples(N);
        double range[2];
        Image -> GetScalarRange(range);
        MaxPArray -> FillComponent(0,range[0]);

        #ifdef DEBUG
            printf("Range = [%f,%f]\n",range[0],range[1]);
        #endif

        //Max Projection bottom
        int x, y, z;
        double v, vproj;
        for (x = Dim[0]; x--;) {
            for (y = Dim[1]; y--;) {
                vproj = 0;
                for (z = zi; z <= zi+8; z++) {
                    v = Image -> GetScalarComponentAsFloat(x,y,z,0);
                    vproj = (v > vproj) ? v : vproj;
                }
                MaxPArray -> SetTuple1(Plane->FindPoint(x+2*Dim[0],y,0),(unsigned char)vproj);
            }
        }
        //Max Projection top
        for (x = Dim[0]; x--;) {
            for (y = Dim[1]; y--;) {
                vproj = 0;
                for (z = zf-8; z <= zf; z++) {
                    v = Image -> GetScalarComponentAsFloat(x,y,z,0);
                    vproj = (v > vproj) ? v : vproj;
                }
                MaxPArray -> SetTuple1(Plane->FindPoint(x+2*Dim[0],y+Dim[1],0),(unsigned char)vproj);
            }
        }

        // Partial surface Projection
        double r[3];
        for (vtkIdType id=0; id < Surface -> GetPoints() -> GetNumberOfPoints(); id++) {
            Surface -> GetPoints() -> GetPoint(id,r);
            x = round(r[0]/_dxy);
            y = round(r[1]/_dxy);
            z = round(r[2]/_dz);
            if ( z >= zi && z <= zi+8 ) {
                MaxPArray -> SetTuple1(Plane->FindPoint(x+3*Dim[0],y,0),255);
            }
            if ( z >= zf-8 && z <= zf ) {
                MaxPArray -> SetTuple1(Plane->FindPoint(x+3*Dim[0],y+Dim[1],0),255);
            }
        }

        // Complete surface Projection
        for (vtkIdType id=0; id < Surface -> GetPoints() -> GetNumberOfPoints(); id++) {
            Surface -> GetPoints() -> GetPoint(id,r);
            x = round(r[0]/_dxy);
            y = round(r[1]/_dxy);
            z = round(r[2]/_dz);
            MaxPArray -> SetTuple1(Plane->FindPoint(x,y,0),255);
        }

        // Partial skeleton Projection
        for (vtkIdType id=0; id < Skeleton -> GetPoints() -> GetNumberOfPoints(); id++) {
            Skeleton -> GetPoints() -> GetPoint(id,r);
            x = round(r[0]/_dxy);
            y = round(r[1]/_dxy);
            z = round(r[2]/_dz);
            if ( z >= zi && z <= zi+8 ) {
                MaxPArray -> SetTuple1(Plane->FindPoint(x+4*Dim[0],y,0),255);
            }
            if ( z >= zf-8 && z <= zf ) {
                MaxPArray -> SetTuple1(Plane->FindPoint(x+4*Dim[0],y+Dim[1],0),255);
            }
        }

        //Complete max projection
        for (x = Dim[0]; x--;) {
            for (y = Dim[1]; y--;) {
                vproj = 0;
                for (z = 0; z < Dim[2]; z++) {
                    v = Image -> GetScalarComponentAsFloat(x,y,z,0);
                    vproj = (v > vproj) ? v : vproj;
                }
                MaxPArray -> SetTuple1(Plane->FindPoint(x,y+Dim[1],0),(unsigned char)vproj);
            }
        }

        //First and last slice
        for (x = Dim[0]; x--;) {
            for (y = Dim[1]; y--;) {
                v = Image -> GetScalarComponentAsFloat(x,y,0,0);
                MaxPArray -> SetTuple1(Plane->FindPoint(x+Dim[1],y,0),(unsigned char)v);
                v = Image -> GetScalarComponentAsFloat(x,y,Dim[2]-1,0);
                MaxPArray -> SetTuple1(Plane->FindPoint(x+Dim[1],y+Dim[1],0),(unsigned char)v);
            }
        }


        MaxPArray -> Modified();

        Plane -> GetPointData() -> SetScalars(MaxPArray);

        sprintf(_fullpath,"%s_detailed.png",FileName);

        // Saving PNG File
        vtkSmartPointer<vtkPNGWriter> PNGWriter = vtkSmartPointer<vtkPNGWriter>::New();
        PNGWriter -> SetFileName(_fullpath);
        PNGWriter -> SetFileDimensionality(2);
        PNGWriter -> SetCompressionLevel(0);
        PNGWriter -> SetInputData(Plane);
        PNGWriter -> Write();

    }

    #ifdef DEBUG
        printf("File Saved!\n");
    #endif

}

/* ================================================================
   IMAGE TRANSFORM
=================================================================*/

vtkSmartPointer<vtkImageData> Convert16To8bit(vtkSmartPointer<vtkImageData> Image) {

    // 8-Bit images
    if (Image -> GetScalarType() == VTK_UNSIGNED_CHAR) {

        return Image;

    // 16-Bit images
    } else if (Image -> GetScalarType() == VTK_UNSIGNED_SHORT) {

        #ifdef DEBUG
            printf("Converting from 16-bit to 8-bit...\n");
        #endif

        vtkSmartPointer<vtkImageData> Image8 = vtkImageData::New();
        Image8 -> ShallowCopy(Image);

        vtkDataArray *ScalarsShort = Image -> GetPointData() -> GetScalars();
        unsigned long int N = ScalarsShort -> GetNumberOfTuples();
        double range[2];
        ScalarsShort -> GetRange(range);

        #ifdef DEBUG
            printf("\tOriginal intensities range: [%d-%d]\n",(int)range[0],(int)range[1]);
        #endif

        vtkSmartPointer<vtkUnsignedCharArray> ScalarsChar = vtkSmartPointer<vtkUnsignedCharArray>::New();
        ScalarsChar -> SetNumberOfComponents(1);
        ScalarsChar -> SetNumberOfTuples(N);
        
        double x, y;
        vtkIdType register id;
        for ( id = N; id--; ) {
            x = ScalarsShort -> GetTuple1(id);
            y = 255.0 * (x-range[0]) / (range[1]-range[0]);
            ScalarsChar -> SetTuple1(id,(unsigned char)y);
        }
        ScalarsChar -> Modified();

        Image8 -> GetPointData() -> SetScalars(ScalarsChar);
        return Image8;

    // Other depth
    } else {
        return NULL;
    }
}

vtkSmartPointer<vtkImageData> BinarizeAndConvertDoubleToChar(vtkSmartPointer<vtkImageData> Image, double threshold) {

    vtkSmartPointer<vtkImageData> Image8 = vtkImageData::New();
    Image8 -> ShallowCopy(Image);

    vtkDataArray *ScalarsDouble = Image -> GetPointData() -> GetScalars();
    unsigned long int N = ScalarsDouble -> GetNumberOfTuples();
    double range[2];
    ScalarsDouble -> GetRange(range);

    vtkSmartPointer<vtkUnsignedCharArray> ScalarsChar = vtkSmartPointer<vtkUnsignedCharArray>::New();
    ScalarsChar -> SetNumberOfComponents(1);
    ScalarsChar -> SetNumberOfTuples(N);
        
    double x;

    if (threshold > 0) {
        for ( vtkIdType id = N; id--; ) {
            x = ScalarsDouble -> GetTuple1(id);
            if (x<=threshold) {
                ScalarsChar -> SetTuple1(id,0);
            } else {
                ScalarsChar -> SetTuple1(id,255);
            }
        }
    } else {
        for ( vtkIdType id = N; id--; ) {
            x = ScalarsDouble -> GetTuple1(id);
            ScalarsChar -> SetTuple1(id,(int)(255*(x-range[0])/(range[1]-range[0])));
        }        
    }
    ScalarsChar -> Modified();

    Image8 -> GetPointData() -> SetScalars(ScalarsChar);
    return Image8;

}

void FillHoles(vtkSmartPointer<vtkImageData> ImageData) {

    #ifdef DEBUG
        printf("\tSearching for holes in the image...\n");
    #endif

    int *Dim = ImageData -> GetDimensions();
    vtkIdType N = ImageData -> GetNumberOfPoints();

    vtkIdType i, s, ido, id;

    int x, y, z;
    double v, r[3];
    bool find = true;
    long long int ro[3];
    long int scluster, label;
    ro[0] = Dim[0] * Dim[1] * Dim[2];
    ro[1] = Dim[0] * Dim[1];

    vtkSmartPointer<vtkIdList> CurrA = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> NextA = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkLongArray> CSz = vtkSmartPointer<vtkLongArray>::New();
    vtkSmartPointer<vtkLongArray> Volume = vtkSmartPointer<vtkLongArray>::New();
    Volume -> SetNumberOfComponents(1);
    Volume -> SetNumberOfTuples(N);
    Volume -> FillComponent(0,0);

    for (x = 1; x < Dim[0]-1; x++) {
        for (y = 1; y < Dim[1]-1; y++) {
            for (z = 1; z < Dim[2]-1; z++) {
                id = ImageData -> FindPoint(x,y,z);
                if ((unsigned short int)ImageData->GetScalarComponentAsDouble(x,y,z,0)) {
                    Volume -> SetTuple1(id,0);
                } else {
                    Volume -> SetTuple1(id,1);
                }
            }
        }
    }

    Volume -> Modified();

    label = 0;
    while (find) {
        for (s = 0; s < CurrA->GetNumberOfIds(); s++) {
            ido = CurrA -> GetId(s);
            ImageData -> GetPoint(ido,r);
            x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
            for (i = 0; i < 6; i++) {
                id = ImageData -> FindPoint(x+ssdx_sort[i],y+ssdy_sort[i],z+ssdz_sort[i]);
                v = Volume -> GetTuple1(id);
                if ((long int)v > 0) {
                    NextA -> InsertNextId(id);
                    Volume -> SetTuple1(id,-label);
                    scluster++;
                }
            }
        }
        if (!NextA->GetNumberOfIds()) {
            find = false;
            for (id=ro[0]; id--;) {
                v = Volume -> GetTuple1(id);
                if ((long int)v > 0) {
                    find = true;
                    ro[0] = id;
                    break;
                }
            }
            if (label) {
                CSz -> InsertNextTuple1(scluster);
            }
            if (find) {
                label++;
                scluster = 1;
                Volume -> SetTuple1(id,-label);
                CurrA -> InsertNextId(id);
            }
        } else {
            CurrA -> Reset();
            CurrA -> DeepCopy(NextA);
            NextA -> Reset();
        }
    }

    for (id = N; id--;) {
        if ((long int)Volume->GetTuple1(id)<-1) {
            ImageData -> GetPointData() -> GetScalars() -> SetTuple1(id,255);
        }
    }
    ImageData -> GetPointData() -> GetScalars() -> Modified();

    #ifdef DEBUG
        printf("\tNumber of filled holes: %ld\n",(long int)CSz->GetNumberOfTuples()-1);
    #endif

}


/* ================================================================
   ROUTINES FOR VESSELNESS CALCUATION VIA DISCRETE APPROCH
=================================================================*/

void GetImageDerivativeDiscrete(vtkSmartPointer<vtkDataArray> Image, int *dim, char direction, vtkSmartPointer<vtkFloatArray> Derivative) {
    #ifdef DEBUG
        printf("Calculating Image Derivatives (Discrete)...\n");
    #endif

    double d, f1, f2;
    vtkIdType i, j, k;
    if (direction=='x') {
        for (i = (vtkIdType)dim[0]; i--;) {
            for (j = (vtkIdType)dim[1]; j--;) {
                for (k = (vtkIdType)dim[2]; k--;) {
                    if (i==0) {
                        Image -> GetTuple(GetId(1,j,k,dim),&f1);
                        Image -> GetTuple(GetId(0,j,k,dim),&f2);
                    }
                    if (i==dim[0]-1) {
                        Image -> GetTuple(GetId(dim[0]-1,j,k,dim),&f1);
                        Image -> GetTuple(GetId(dim[0]-2,j,k,dim),&f2);
                    }
                    if (i>0&i<dim[0]-1) {
                        Image -> GetTuple(GetId(i+1,j,k,dim),&f1);
                        Image -> GetTuple(GetId(i-1,j,k,dim),&f2);
                        f1 /= 2.0; f2 /= 2.0;
                    }
                    d = f1 - f2;
                    Derivative -> SetTuple(GetId(i,j,k,dim),&d);
                }
            }
        }
    }

    if (direction=='y') {
        for (i = (vtkIdType)dim[0]; i--;) {
            for (j = (vtkIdType)dim[1]; j--;) {
                for (k = (vtkIdType)dim[2]; k--;) {
                    if (j==0) {
                        Image -> GetTuple(GetId(i,1,k,dim),&f1);
                        Image -> GetTuple(GetId(i,0,k,dim),&f2);
                    }
                    if (j==dim[1]-1) {
                        Image -> GetTuple(GetId(i,dim[1]-1,k,dim),&f1);
                        Image -> GetTuple(GetId(i,dim[1]-2,k,dim),&f2);
                    }
                    if (j>0&j<dim[1]-1) {
                        Image -> GetTuple(GetId(i,j+1,k,dim),&f1);
                        Image -> GetTuple(GetId(i,j-1,k,dim),&f2);
                        f1 /= 2.0; f2 /= 2.0;
                    }
                    d = f1 - f2;
                    Derivative -> SetTuple(GetId(i,j,k,dim),&d);
                }
            }
        }
    }

    if (direction=='z') {
        for (i = (vtkIdType)dim[0]; i--;) {
            for (j = (vtkIdType)dim[1]; j--;) {
                for (k = (vtkIdType)dim[2]; k--;) {
                    if (k==0) {
                        Image -> GetTuple(GetId(i,j,1,dim),&f1);
                        Image -> GetTuple(GetId(i,j,0,dim),&f2);
                    }
                    if (k==dim[2]-1) {
                        Image -> GetTuple(GetId(i,j,dim[2]-1,dim),&f1);
                        Image -> GetTuple(GetId(i,j,dim[2]-2,dim),&f2);
                    }
                    if (k>0&k<dim[2]-1) {
                        Image -> GetTuple(GetId(i,j,k+1,dim),&f1);
                        Image -> GetTuple(GetId(i,j,k-1,dim),&f2);
                        f1 /= 2.0; f2 /= 2.0;
                    }
                    d = f1 - f2;
                    Derivative -> SetTuple(GetId(i,j,k,dim),&d);
                }
            }
        }
    }

}

void GetHessianEigenvaluesDiscrete(double sigma, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkDoubleArray> L1, vtkSmartPointer<vtkDoubleArray> L2, vtkSmartPointer<vtkDoubleArray> L3) {
    #ifdef DEBUG
        printf("Calculating Hessian Eigeinvalues (Discrete)...\n");
    #endif

    int *Dim = Image -> GetDimensions();
    vtkIdType id, N = Image -> GetNumberOfPoints();
    double H[3][3], Eva[3], Eve[3][3], dxx, dyy, dzz, dxy, dxz, dyz, l1, l2, l3, frobnorm;

    #ifdef DEBUG
        printf("Calculating Gaussian Convolution...\n");
    #endif

    vtkSmartPointer<vtkImageGaussianSmooth> Gauss = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    Gauss -> SetInputData(Image);
    Gauss -> SetDimensionality(3);
    Gauss -> SetRadiusFactors(10,10,10);
    Gauss -> SetStandardDeviations(sigma,sigma,sigma);
    Gauss -> Update();

    vtkFloatArray *ImageG = (vtkFloatArray*) Gauss -> GetOutput() -> GetPointData() -> GetScalars();

    vtkSmartPointer<vtkFloatArray> Dx = vtkSmartPointer<vtkFloatArray>::New(); Dx -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dy = vtkSmartPointer<vtkFloatArray>::New(); Dy -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dz = vtkSmartPointer<vtkFloatArray>::New(); Dz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dxx = vtkSmartPointer<vtkFloatArray>::New(); Dxx -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dyy = vtkSmartPointer<vtkFloatArray>::New(); Dyy -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dzz = vtkSmartPointer<vtkFloatArray>::New(); Dzz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dxy = vtkSmartPointer<vtkFloatArray>::New(); Dxy -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dxz = vtkSmartPointer<vtkFloatArray>::New(); Dxz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dyz = vtkSmartPointer<vtkFloatArray>::New(); Dyz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Fro = vtkSmartPointer<vtkFloatArray>::New();
    Fro -> SetNumberOfComponents(1);
    Fro -> SetNumberOfTuples(N);

    GetImageDerivativeDiscrete(ImageG,Dim,'x',Dx);
    GetImageDerivativeDiscrete(ImageG,Dim,'y',Dy);
    GetImageDerivativeDiscrete(ImageG,Dim,'z',Dz);
    GetImageDerivativeDiscrete(Dx,Dim,'x',Dxx);
    GetImageDerivativeDiscrete(Dy,Dim,'y',Dyy);
    GetImageDerivativeDiscrete(Dz,Dim,'z',Dzz);
    GetImageDerivativeDiscrete(Dy,Dim,'x',Dxy);
    GetImageDerivativeDiscrete(Dz,Dim,'x',Dxz);
    GetImageDerivativeDiscrete(Dz,Dim,'y',Dyz);

    for ( id = N; id--; ) {
        l1 = l2 = l3 = 0.0;
        H[0][0]=Dxx->GetTuple1(id); H[0][1]=Dxy->GetTuple1(id); H[0][2]=Dxz->GetTuple1(id);
        H[1][0]=Dxy->GetTuple1(id); H[1][1]=Dyy->GetTuple1(id); H[1][2]=Dyz->GetTuple1(id);
        H[2][0]=Dxz->GetTuple1(id); H[2][1]=Dyz->GetTuple1(id); H[2][2]=Dzz->GetTuple1(id);
        frobnorm = FrobeniusNorm(H);
        if (H[0][0]+H[1][1]+H[2][2]<0.0) {
            vtkMath::Diagonalize3x3(H,Eva,Eve);
            l1 = Eva[0]; l2 = Eva[1]; l3 = Eva[2];
            Sort(&l1,&l2,&l3);
        }
        L1 -> SetTuple1(id,l1);
        L2 -> SetTuple1(id,l2);
        L3 -> SetTuple1(id,l3);
        Fro -> SetTuple1(id,frobnorm);
    }
    double ftresh,frobenius_norm_range[2];
    Fro -> GetRange(frobenius_norm_range);
    ftresh = sqrt(frobenius_norm_range[1]);

    for ( id = N; id--; ) {
        if ( Fro->GetTuple1(id) < ftresh) {
            L1 -> SetTuple1(id,0.0);
            L2 -> SetTuple1(id,0.0);
            L3 -> SetTuple1(id,0.0);
        }
    }
    L1 -> Modified();
    L2 -> Modified();
    L3 -> Modified();

}

void GetHessianEigenvaluesDiscreteZDependentThreshold(double sigma, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkDoubleArray> L1, vtkSmartPointer<vtkDoubleArray> L2, vtkSmartPointer<vtkDoubleArray> L3) {
    #ifdef DEBUG
        printf("Calculating Hessian Eigeinvalues (Discrete)...\n");
    #endif

    int *Dim = Image -> GetDimensions();
    vtkIdType id, N = Image -> GetNumberOfPoints();
    double H[3][3], Eva[3], Eve[3][3], dxx, dyy, dzz, dxy, dxz, dyz, l1, l2, l3, frobnorm;

    #ifdef DEBUG
        printf("Calculating Gaussian Convolution...\n");
    #endif

    vtkSmartPointer<vtkImageGaussianSmooth> Gauss = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    Gauss -> SetInputData(Image);
    Gauss -> SetDimensionality(3);
    Gauss -> SetRadiusFactors(10,10,10);
    Gauss -> SetStandardDeviations(sigma,sigma,sigma);
    Gauss -> Update();

    vtkFloatArray *ImageG = (vtkFloatArray*) Gauss -> GetOutput() -> GetPointData() -> GetScalars();

    vtkSmartPointer<vtkFloatArray> Dx = vtkSmartPointer<vtkFloatArray>::New(); Dx -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dy = vtkSmartPointer<vtkFloatArray>::New(); Dy -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dz = vtkSmartPointer<vtkFloatArray>::New(); Dz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dxx = vtkSmartPointer<vtkFloatArray>::New(); Dxx -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dyy = vtkSmartPointer<vtkFloatArray>::New(); Dyy -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dzz = vtkSmartPointer<vtkFloatArray>::New(); Dzz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dxy = vtkSmartPointer<vtkFloatArray>::New(); Dxy -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dxz = vtkSmartPointer<vtkFloatArray>::New(); Dxz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Dyz = vtkSmartPointer<vtkFloatArray>::New(); Dyz -> SetNumberOfTuples(N);
    vtkSmartPointer<vtkFloatArray> Fro = vtkSmartPointer<vtkFloatArray>::New();
    Fro -> SetNumberOfComponents(1);
    Fro -> SetNumberOfTuples(N);

    GetImageDerivativeDiscrete(ImageG,Dim,'x',Dx);
    GetImageDerivativeDiscrete(ImageG,Dim,'y',Dy);
    GetImageDerivativeDiscrete(ImageG,Dim,'z',Dz);
    GetImageDerivativeDiscrete(Dx,Dim,'x',Dxx);
    GetImageDerivativeDiscrete(Dy,Dim,'y',Dyy);
    GetImageDerivativeDiscrete(Dz,Dim,'z',Dzz);
    GetImageDerivativeDiscrete(Dy,Dim,'x',Dxy);
    GetImageDerivativeDiscrete(Dz,Dim,'x',Dxz);
    GetImageDerivativeDiscrete(Dz,Dim,'y',Dyz);

    int x, y, z;
    double *FThresh = new double[Dim[2]];
    for ( id = Dim[2]; id--; ) FThresh[id] = 0.0;

    for ( id = N; id--; ) {
        l1 = l2 = l3 = 0.0;
        H[0][0]=Dxx->GetTuple1(id); H[0][1]=Dxy->GetTuple1(id); H[0][2]=Dxz->GetTuple1(id);
        H[1][0]=Dxy->GetTuple1(id); H[1][1]=Dyy->GetTuple1(id); H[1][2]=Dyz->GetTuple1(id);
        H[2][0]=Dxz->GetTuple1(id); H[2][1]=Dyz->GetTuple1(id); H[2][2]=Dzz->GetTuple1(id);
        frobnorm = FrobeniusNorm(H);
        if (H[0][0]+H[1][1]+H[2][2]<0.0) {
            vtkMath::Diagonalize3x3(H,Eva,Eve);
            l1 = Eva[0]; l2 = Eva[1]; l3 = Eva[2];
            Sort(&l1,&l2,&l3);
        }
        L1 -> SetTuple1(id,l1);
        L2 -> SetTuple1(id,l2);
        L3 -> SetTuple1(id,l3);
        Fro -> SetTuple1(id,frobnorm);
        z = GetZ(id,Dim);
        FThresh[z] = (frobnorm > FThresh[z]) ? frobnorm : FThresh[z];
    }

    for ( z = Dim[2]; z--; ) FThresh[z] = sqrt(FThresh[z]);

    int j;
    double frobneigh;
    for ( id = N; id--; ) {
        x = GetX(id,Dim);
        y = GetY(id,Dim);
        z = GetZ(id,Dim);
        frobneigh = 0.0;
        if (x>0&&x<Dim[0]-1&&y>0&&y<Dim[1]-1&&z>0&&z<Dim[2]-1) {
            for (j = 0; j < 6; j++) {
                frobneigh += Fro -> GetTuple1(GetId(x+ssdx_sort[j],y+ssdy_sort[j],z+ssdz_sort[j],Dim));
            }
            frobneigh /= 6.0;
        }
        if ( frobneigh < FThresh[z] ) {
            L1 -> SetTuple1(id,0.0);
            L2 -> SetTuple1(id,0.0);
            L3 -> SetTuple1(id,0.0);
        }
    }
    L1 -> Modified();
    L2 -> Modified();
    L3 -> Modified();

    delete[] FThresh;
}

/* ================================================================
   VESSELNESS ROUTINE
=================================================================*/

void GetVesselness(double sigma, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkDoubleArray> L1, vtkSmartPointer<vtkDoubleArray> L2, vtkSmartPointer<vtkDoubleArray> L3) {

    #ifdef DEBUG
        printf("Calculating Vesselness...\n");
    #endif

    double c = 500.0;
    double beta = 0.5;
    double alpha = 0.5;
    double std = 2 * c * c;
    double rbd = 2 * beta * beta;
    double rad = 2 * alpha * alpha;
    double l1, l2, l3, ra, ran, rb, rbn, st, stn, ft_old, ft_new;
    vtkIdType id, N = Image -> GetNumberOfPoints();

    if (_adaptive_threshold) {
        GetHessianEigenvaluesDiscreteZDependentThreshold(sigma,Image,L1,L2,L3);
    } else {
        GetHessianEigenvaluesDiscrete(sigma,Image,L1,L2,L3);
    }
    
    for ( id = N; id--; ) {
        l1 = L1 -> GetTuple1(id);
        l2 = L2 -> GetTuple1(id);
        l3 = L3 -> GetTuple1(id);
        if (l2<0&&l3<0) {

            ra = fabs(l2) / fabs(l3);
            ran = -ra * ra;

            rb = fabs(l1) / sqrt(l2*l3);
            rbn = -rb * rb;

            st = sqrt(l1*l1+l2*l2+l3*l3);
            stn = -st * st;

            ft_new = (1-exp(ran/rad)) * exp(rbn/rbd) * (1-exp(stn/std));

            //L1 is used to return vesselness values
            L1 -> SetTuple1(id,ft_new);

        } else L1 -> SetTuple1(id,0.0);
    }
    L1 -> Modified();
}

/* ================================================================
   DIVERGENCE FILTER
=================================================================*/

void GetDivergenceFilter(int *Dim, vtkSmartPointer<vtkDoubleArray> Scalars) {

    #ifdef DEBUG
        printf("Calculating Divergent Filter...\n");
    #endif

    vtkIdType id;
    int register j, i;
    int x, y, z, s = 2;
    double v, norm, V[6][3];
    int Dx[6] = {1,-1,0,0,0,0};
    int Dy[6] = {0,0,1,-1,0,0};
    int Dz[6] = {0,0,0,0,1,-1};
    int MI[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

    vtkSmartPointer<vtkDoubleArray> Div = vtkSmartPointer<vtkDoubleArray>::New();
    Div -> SetNumberOfComponents(1);
    Div -> SetNumberOfTuples(Scalars->GetNumberOfTuples());
    Div -> FillComponent(0,0.0);

    for (z = s+1; z < Dim[2]-s-1; z++) {
        for (y = s+1; y < Dim[1]-s-1; y++) {
            for (x = s+1; x < Dim[0]-s-1; x++) {
                v = 0.0;
                id = GetId(x,y,z,Dim);
                if (Scalars->GetTuple1(id)) {
                    for (i = 0; i < 6; i++) {
                        for (j = 0; j < 3; j++) {
                            V[i][j]  = Scalars -> GetTuple1(GetId(x+s*Dx[i]+MI[j][0],y+s*Dy[i]+MI[j][1],z+s*Dz[i]+MI[j][2],Dim));
                            V[i][j] -= Scalars -> GetTuple1(GetId(x+s*Dx[i]-MI[j][0],y+s*Dy[i]-MI[j][1],z+s*Dz[i]-MI[j][2],Dim));
                        }
                        norm = sqrt(pow(V[i][0],2)+pow(V[i][1],2)+pow(V[i][2],2));
                        if (norm) {V[i][0]/=norm; V[i][1]/=norm; V[i][2]/=norm; }
                    }
                    v = (V[0][0]-V[1][0])+(V[2][1]-V[3][1])+(V[4][2]-V[5][2]);
                    v = (v<0) ? -v / 6.0 : 0.0;
                }
                Div -> InsertTuple1(id,v);
            }
        }
    }
    Div -> Modified();
    Scalars -> DeepCopy(Div);
    Scalars -> Modified();

}

/* ================================================================
   WIDTH ANALYSIS
=================================================================*/

void EstimateTubuleWidth(vtkSmartPointer<vtkPolyData> Skeleton, vtkSmartPointer<vtkPolyData> Surface, double *attributes) {

    #ifdef DEBUG
        printf("Calculating tubules width...\n");
    #endif

    vtkIdType N = Skeleton -> GetNumberOfPoints();
    vtkSmartPointer<vtkDoubleArray> Width = vtkSmartPointer<vtkDoubleArray>::New();
    Width -> SetName("Width");
    Width -> SetNumberOfComponents(1);
    Width -> SetNumberOfTuples(N);

    #ifdef DEBUG
        printf("\tGenerating point locator...\n");
    #endif

    vtkSmartPointer<vtkKdTreePointLocator> Tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
    Tree -> SetDataSet(Surface);
    Tree -> BuildLocator();

    int k, n = 3;
    vtkIdType id, idk;
    double r[3], rk[3], w;
    vtkSmartPointer<vtkIdList> List = vtkSmartPointer<vtkIdList>::New();
    
    double av_w = 0.0, sd_w = 0.0;

    for (id = 0; id < N; id++) {
        w = 0;
        Skeleton -> GetPoint(id,r);
        Tree -> FindClosestNPoints(n,r,List);
        for (k = 0; k < n; k++) {
            idk = List -> GetId(k);
            Surface -> GetPoint(idk,rk);
            w += 2.0*sqrt( pow(r[0]-rk[0],2) + pow(r[1]-rk[1],2) + pow(r[2]-rk[2],2) );
        }
        w /= n;
        av_w += w;
        sd_w += w*w;
        Width -> SetTuple1(id,w);

    }
    Width -> Modified();
    Skeleton -> GetPointData() -> SetScalars(Width);

    attributes[3] = av_w / N;
    attributes[4] = sqrt(sd_w/N - (av_w/N)*(av_w/N));

}

/* ================================================================
   INTENSITY MAPPING
=================================================================*/

void MapImageIntensity(vtkSmartPointer<vtkPolyData> Skeleton, vtkSmartPointer<vtkImageData> ImageData, int nneigh) {

    vtkIdType N = Skeleton -> GetNumberOfPoints();
    vtkSmartPointer<vtkDoubleArray> Intensity = vtkSmartPointer<vtkDoubleArray>::New();
    Intensity -> SetName("Intensity");
    Intensity -> SetNumberOfComponents(1);
    Intensity -> SetNumberOfTuples(N);

    vtkIdType id;
    int x, y, z, k;
    double v, r[3];
    for (id = 0; id < N; id++) {
        v = 0;
        Skeleton -> GetPoint(id,r);
        x = round(r[0] / _dxy);
        y = round(r[1] / _dxy);
        z = round(r[2] / _dz);
        for (k = 0; k < nneigh; k++) {
            v += ImageData -> GetScalarComponentAsDouble(x+ssdx_sort[k],y+ssdy_sort[k],z+ssdz_sort[k],0);
        }
        Intensity -> SetTuple1(id,v/((double)nneigh));
    }
    Intensity -> Modified();
    Skeleton -> GetPointData() -> AddArray(Intensity);

}

/* ================================================================
   WIDTH-CORRECTED VOLUME
=================================================================*/

void GetVolumeFromSkeletonLengthAndWidth(vtkSmartPointer<vtkPolyData> PolyData, double *attributes) {
    double h, r1[3], r2[3], R1, R2, volume = 0.0, length = 0.0;
    vtkPoints *Points = PolyData -> GetPoints();
    for (vtkIdType edge=PolyData->GetNumberOfCells();edge--;) {
        for (vtkIdType n = 1; n < PolyData->GetCell(edge)->GetNumberOfPoints(); n++) {
            PolyData -> GetPoint(PolyData->GetCell(edge)->GetPointId(n-1),r1);
            PolyData -> GetPoint(PolyData->GetCell(edge)->GetPointId(n  ),r2);
            R1 = 0.5*PolyData -> GetPointData() -> GetArray("Width") -> GetTuple1(PolyData->GetCell(edge)->GetPointId(n-1));
            R2 = 0.5*PolyData -> GetPointData() -> GetArray("Width") -> GetTuple1(PolyData->GetCell(edge)->GetPointId(n  ));
            h = sqrt(pow(r2[0]-r1[0],2)+pow(r2[1]-r1[1],2)+pow(r2[2]-r1[2],2));
            length += h;
            volume += 3.141592 / 3.0 * h * (R1*R1+R2*R2+R1*R2);
        }
    }
    attributes[1] = length;
    attributes[2] = length * (acos(-1.0)*pow(_rad,2));
    //attributes[3] = volume; // Validation needed
}


/* ================================================================
   TOPOLOGICAL ATTRIBUTES FROM SKELETON
=================================================================*/

void GetTopologicalAttributes(vtkSmartPointer<vtkPolyData> PolyData, double *attributes) {
    long int n, k;
    std::list<vtkIdType> Endpoints;
    for (vtkIdType edge=PolyData->GetNumberOfCells();edge--;) {
        n = PolyData->GetCell(edge)->GetNumberOfPoints();
        Endpoints.push_back(PolyData->GetCell(edge)->GetPointId(0));
        Endpoints.push_back(PolyData->GetCell(edge)->GetPointId(n-1));
    }
    std::list<vtkIdType> Nodes = Endpoints;
    Nodes.sort();
    Nodes.unique();
    long int ne = 0, nb = 0;
    for (std::list<vtkIdType>::iterator it = Nodes.begin(); it!=Nodes.end(); it++) {
        k = std::count(Endpoints.begin(),Endpoints.end(),*it);
        ne += (k==1) ? 1 : 0;
        nb += (k>=3) ? 1 : 0;
    }
    attributes[5] = ne;
    attributes[6] = nb;
}

/* ================================================================
   MULTISCALE VESSELNESS
=================================================================*/

int MultiscaleVesselness(const char FileName[], double _sigmai, double _dsigma, double _sigmaf, double *attributes) {     

    // Loading multi-paged TIFF file (Supported by VTK 6.2 and higher)
    char _fullpath[256];
    sprintf(_fullpath,"%s.tif",FileName);
    vtkSmartPointer<vtkTIFFReader> TIFFReader = vtkSmartPointer<vtkTIFFReader>::New();
    int errlog = TIFFReader -> CanReadFile(_fullpath);
    // File cannot be opened
    if (!errlog) {
        printf("File %s cannnot be opened.\n",_fullpath);
        return -1;
    }
    TIFFReader -> SetFileName(_fullpath);
    TIFFReader -> Update();

    int *Dim = TIFFReader -> GetOutput() -> GetDimensions();

    #ifdef DEBUG
        printf("MitoGraph V2.0 [DEBUG mode]\n");
        printf("File name: %s\n",_fullpath);
        printf("Volume dimensions: %dx%dx%d\n",Dim[0],Dim[1],Dim[2]);
        printf("Scales to run: [%1.3f:%1.3f:%1.3f]\n",_sigmai,_dsigma,_sigmaf);
        printf("Threshold: %1.5f\n",_div_threshold);
    #endif

    // Exporting resampled images

    if (_export_image_resampled) {
        sprintf(_fullpath,"%s_resampled.vtk",FileName);
        SaveImageData(TIFFReader->GetOutput(),_fullpath,true);
    }

    // Conversion 16-bit to 8-bit
    vtkSmartPointer<vtkImageData> Image = Convert16To8bit(TIFFReader->GetOutput());

    if (!Image) printf("Format not supported.\n");

    //VESSELNESS
    //----------

    vtkIdType id, N = Image -> GetNumberOfPoints();

    vtkSmartPointer<vtkDoubleArray> AUX1 = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> AUX2 = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> AUX3 = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> VSSS = vtkSmartPointer<vtkDoubleArray>::New();

    AUX1 -> SetNumberOfTuples(N);
    AUX2 -> SetNumberOfTuples(N);
    AUX3 -> SetNumberOfTuples(N);
    VSSS -> SetNumberOfTuples(N);
    AUX1 -> FillComponent(0,0.0);
    AUX2 -> FillComponent(0,0.0);
    AUX3 -> FillComponent(0,0.0);
    VSSS -> FillComponent(0,0.0);

    double sigma, vn, vo;

    for ( sigma = _sigmai; sigma <= _sigmaf+0.5*_dsigma; sigma += _dsigma ) {
        
        #ifdef DEBUG
            printf("Running sigma = %1.3f\n",sigma);
        #endif
        
        GetVesselness(sigma,Image,AUX1,AUX2,AUX3);
        
        for ( id = N; id--; ) {
            vn = AUX1 -> GetTuple1(id);
            vo = VSSS -> GetTuple1(id);
            if ( vn > vo ) {
                VSSS -> SetTuple1(id,vn);
            }
        }

    }
    VSSS -> Modified();

    //DIVERGENCE FILTER
    //-----------------

    GetDivergenceFilter(Dim,VSSS);

    vtkSmartPointer<vtkImageData> ImageEnhanced = vtkSmartPointer<vtkImageData>::New();
    ImageEnhanced -> GetPointData() -> SetScalars(VSSS);
    ImageEnhanced -> SetDimensions(Dim);

    #ifdef DEBUG
        sprintf(_fullpath,"%s_div.vtk",FileName);
        SaveImageData(BinarizeAndConvertDoubleToChar(ImageEnhanced,-1),_fullpath);
    #endif

    #ifdef DEBUG
        printf("Clear boundaries and removing tiny components...\n");
    #endif

    CleanImageBoundaries(ImageEnhanced);

    long int cluster;
    std::vector<long int> CSz;
    vtkSmartPointer<vtkDoubleArray> Volume = vtkSmartPointer<vtkDoubleArray>::New();
    Volume -> SetNumberOfComponents(0);
    Volume -> SetNumberOfTuples(N);
    Volume -> FillComponent(0,0);
    long int ncc = LabelConnectedComponents(ImageEnhanced,Volume,CSz,6,_div_threshold);

    if (ncc > 1) {
        for (id = N; id--;) {
            cluster = (long int)Volume -> GetTuple1(id);
            if (cluster < 0) {
                if (CSz[-cluster-1] <= 5) {
                    ImageEnhanced -> GetPointData() -> GetScalars() -> SetTuple1(id,0);
                }
            }
        }
    }

    //CREATING SURFACE POLYDATA
    //-------------------------

    vtkSmartPointer<vtkContourFilter> Filter = vtkSmartPointer<vtkContourFilter>::New();
    Filter -> SetInputData(ImageEnhanced);
    Filter -> SetValue(1,_div_threshold);
    Filter -> Update();

    vtkSmartPointer<vtkPolyData> Surface = Filter -> GetOutput();
    ScalePolyData(Surface);

    //SAVING SURFACE
    //--------------

    sprintf(_fullpath,"%s_surface.vtk",FileName);
    SavePolyData(Filter->GetOutput(),_fullpath);

    //BINARIZATION
    //------------

    vtkSmartPointer<vtkImageData> Binary = BinarizeAndConvertDoubleToChar(ImageEnhanced,_div_threshold);

    //FILLING HOLES
    //-------------
    if (_improve_skeleton_quality) FillHoles(Binary);

    //MAX PROJECTION
    //--------------

    sprintf(_fullpath,"%s.png",FileName);
    ExportMaxProjection(Binary,_fullpath);

    //SKELETONIZATION
    //---------------

    vtkSmartPointer<vtkPolyData> Skeleton = Thinning3D(Binary,FileName,attributes);
    ScalePolyData(Skeleton);

    //TUBULES WIDTH
    //-------------

    EstimateTubuleWidth(Skeleton,Filter->GetOutput(),attributes);

    //INTENSITY PROFILE ALONG THE SKELETON
    //------------------------------------

    sprintf(_fullpath,"%s.tif",FileName);
    TIFFReader = vtkSmartPointer<vtkTIFFReader>::New();
    TIFFReader -> SetFileName(_fullpath);
    TIFFReader -> Update();
    vtkSmartPointer<vtkImageData> ImageData = TIFFReader -> GetOutput();

    MapImageIntensity(Skeleton,ImageData,6);

    vtkDataArray *W = Skeleton -> GetPointData() -> GetArray("Width");
    vtkDataArray *I = Skeleton -> GetPointData() -> GetArray("Intensity");

    vtkIdType p;
    double length;
    sprintf(_fullpath,"%s.width",FileName);
    FILE *fw = fopen(_fullpath,"w");
    for (vtkIdType edge = 0; edge < Skeleton -> GetNumberOfCells(); edge++) {
        length = GetEdgeLength(edge,Skeleton);
        for (vtkIdType id = 0; id < Skeleton -> GetCell(edge) -> GetNumberOfPoints(); id++) {
            p = Skeleton -> GetCell(edge) -> GetPointId(id);
            fprintf(fw,"%d\t%1.5f\t%1.5f\t%1.5f\n",(int)edge,length,W->GetTuple1(p),I->GetTuple1(p));
        }
    }
    fclose(fw);

    GetVolumeFromSkeletonLengthAndWidth(Skeleton,attributes); //Validation needed

    GetTopologicalAttributes(Skeleton,attributes);

    //SAVING SKELETON
    //---------------

    sprintf(_fullpath,"%s_skeleton.vtk",FileName);
    SavePolyData(Skeleton,_fullpath);

    return 0;
}

/* ================================================================
   MAIN ROUTINE
=================================================================*/

int main(int argc, char *argv[]) {     

    int i;
    char _prefix[64];
    char _impath[128];
    sprintf(_impath,"");
    double _sigmai = 1.00;
    double _sigmaf = 1.50;
       int _nsigma = 6;

    // Collecting input parameters
    for (i = 0; i < argc; i++) {
        if (!strcmp(argv[i],"-path")) {
            sprintf(_impath,"%s//",argv[i+1]);
        }
        if (!strcmp(argv[i],"-xy")) {
            _dxy = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-z")) {
            _dz = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-rad")) {
            _rad = atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-scales")) {
            _sigmai = atof(argv[i+1]);
            _sigmaf = atof(argv[i+2]);
            _nsigma = atoi(argv[i+3]);
        }
        if (!strcmp(argv[i],"-adaptive")) {
            _adaptive_threshold = true;
        }
        if (!strcmp(argv[i],"-threshold")) {
            _div_threshold = (double)atof(argv[i+1]);
        }
        if (!strcmp(argv[i],"-scale_off")) {
            _scale_polydata_before_save = false;
        }        
        if (!strcmp(argv[i],"-graph_off")) {
            _export_graph_files = false;
        }
        if (!strcmp(argv[i],"-labels_off")) {
            _export_nodes_label = false;
        }
        if (!strcmp(argv[i],"-checkonly")) {
            _checkonly = true;
        }
        if (!strcmp(argv[i],"-width")) {
            _width = true;
        }
        if (!strcmp(argv[i],"-precision_off")) {
            _improve_skeleton_quality = false;
        }
        if (!strcmp(argv[i],"-export_image_resampled")) {
            _export_image_resampled = true;
        }
    }

    if (_dz<0) {
        printf("Please, use -dxy and -dz to provide the pixel size.\n");
        return -1;
    }

    // Generating list of files to run
    char _cmd[256];
    sprintf(_cmd,"ls %s*.tif | sed -e 's/.tif//' > %smitograph.files",_impath,_impath);
    system(_cmd);

    char _tifffilename[256];
    char _tifflistpath[128];
    sprintf(_tifflistpath,"%smitograph.files",_impath);
    FILE *f = fopen(_tifflistpath,"r");

    if (_checkonly) {

        while (fgets(_tifffilename,256, f) != NULL) {
            _tifffilename[strcspn(_tifffilename, "\n" )] = '\0';

            ExportDetailedMaxProjection(_tifffilename);

            printf("%s\n",_tifffilename);
        }

    } else {

        double _dsigma = (_sigmaf-_sigmai) / (_nsigma-1);

        // Generating summary file and writing the header
        char _summaryfilename[256];
        char _individfilename[256];
        sprintf(_summaryfilename,"%ssummary.txt",_impath);
        FILE *fsummary = fopen(_summaryfilename,"w");
        if (_adaptive_threshold) {
            fprintf(fsummary,"MitoGraph V2.1Beta [Adaptive Algorithm]\n");
        } else {
            fprintf(fsummary,"MitoGraph V2.0\n");
        }
        fprintf(fsummary,"Folder: %s\n",_impath);
        fprintf(fsummary,"Pixel size: -xy %1.4fum, -z %1.4fum\n",_dxy,_dz);
        fprintf(fsummary,"Average tubule radius: -r %1.4fum\n",_rad);
        fprintf(fsummary,"Scales: -scales %1.2f",_sigmai);
        for ( double sigma = _sigmai+_dsigma; sigma < _sigmaf+0.5*_dsigma; sigma += _dsigma )
            fprintf(fsummary," %1.2f",sigma);
        fprintf(fsummary,"\nPost-divergence threshold: -threshold %1.5f\n",_div_threshold);
        time_t now = time(0);
        fprintf(fsummary,"%s\n",ctime(&now));
        fprintf(fsummary,"Image\tsurface-volume_(um3)\ttotal_length_(um)\tskeleton-volume_(um3)\tavg_width_(um)\tstd_width_(um)\t#Endpoints\t#Branches\n");
        fclose(fsummary);

        // Multiscale vesselness
        double *attributes = new double[7];
        sprintf(_tifflistpath,"%smitograph.files",_impath);
        FILE *f = fopen(_tifflistpath,"r");
        while (fgets(_tifffilename,256, f) != NULL) {
            _tifffilename[strcspn(_tifffilename, "\n" )] = '\0';

            MultiscaleVesselness(_tifffilename,_sigmai,_dsigma,_sigmaf,attributes);
            
            // Saving network attributes in the group file
            fsummary = fopen(_summaryfilename,"a");
            fprintf(fsummary,"%s\t%1.5f\t%1.5f\t%1.5f\t%1.5f\t%1.5f\t%d\t%d\n",_tifffilename,attributes[0],attributes[1],attributes[2],attributes[3],attributes[4],(int)attributes[5],(int)attributes[6]);
            fclose(fsummary);

            // Saving network attributes in the individual file
            sprintf(_individfilename,"%s.mitograph",_tifffilename);
            FILE *findv = fopen(_individfilename,"w");
            fprintf(findv,"surface-volume_(um3)\ttotal_length_(um)\tskeleton-volume_(um3)\tavg_width_(um)\tstd_width (um)\n");
            fprintf(findv,"%1.5f\t%1.5f\t%1.5f\t%1.5f\t%1.5f\n",attributes[0],attributes[1],attributes[2],attributes[3],attributes[4]);
            fclose(findv);

            // Also printing on the screen
            printf("%s\t[done]\n",_tifffilename);

        }

    }

    fclose(f);
    return 0;
}
