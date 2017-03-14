// =====================================================================================================
// MitoThinning: This routine receives as input an ImageData binary from MitoGraph main routine a and
// apply a thinning process over it, resulting in a new but topologically equivalent image.
//
// Matheus P. Viana - vianamp@gmail.com - 2014.06.10
// Susanne Rafelski Lab, University of California Irvine
// =====================================================================================================

#include "MitoThinning.h"

    int ssdx[26] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    int ssdy[26] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int ssdz[26] = { 1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0, 0,-1,-1,-1};

// Routine used to save .gnet and .coo files representing
// the skeleton of the mitochondrial network.
void ExportGraphFiles(vtkSmartPointer<vtkPolyData> PolyData, long int nnodes, const char Prefix[]);

// Routine to create the file _nodes.vtk. This file contains
// little speres located at the junctions (nodes) coordinates
// and might have the nodes label depending on if the variable
// export_node_labels is true or false.
void ExportNodes(vtkSmartPointer<vtkPolyData> PolyData, long int nnodes, long int *ValidId, _mitoObject *mitoObject);

// Routine used to check whether a voxel is locates at the border
// of the binary image or not. A border voxel is defined as those
// that have at least one of their 6-neighbors equal to zero.
bool IsBorder(vtkIdType id, vtkSmartPointer<vtkImageData> Image);

// Routine used to smooth the parametric curves that describe the
// skeleton edges. A simple average filter is used to do the job.
// Sigma represents the number of times the filter is applied.
void SmoothEdgesCoordinates(vtkSmartPointer<vtkPolyData> PolyData, double sigma);

// Estimate the mitochondrial volume by counting the number of
// pixels in the binary image used as input for thinning.
void GetVolumeFromVoxels(vtkSmartPointer<vtkImageData> Image, std::vector<attribute> *Atts);

// Estimate the mitochondrial volume by using the skeleton
// total length and assuming constant radius.
//void GetVolumeFromSkeletonLength(vtkSmartPointer<vtkPolyData> PolyData, double *attributes);

// Calculate the length of a given edge.
//double GetEdgeLength(vtkIdType edge, vtkSmartPointer<vtkPolyData> PolyData);

// Returns the number of voxels around the voxel (x,y,z) with
// value given different of "value".
char GetNumberOfNeighborsWithoutValue(vtkSmartPointer<vtkImageData> Image, int x, int y, int z, long int value);

// Returns the number of voxels around the voxel (x,y,z) with
// value different of "value" in the vector "Volume".
char GetNumberOfNeighborsWithoutValue(vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, int x, int y, int z, long int value);

// Returns the number of voxels around the voxel (x,y,z) with
// value given by "value".
char GetNumberOfNeighborsWithValue(vtkSmartPointer<vtkImageData> Image, int x, int y, int z, long int value);

// Returns the number of voxels around the voxel (x,y,z) with
// value "value" in the vector "Volume".
char GetNumberOfNeighborsWithValue(vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, int x, int y, int z, long int value);

// Returns one neighbor of (x,y,z) with value "value" in the
// vector "Volume".
vtkIdType GetOneNeighborWithValue(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, long int value);

// Returns one neighbor of (x,y,z) with value different of
// "value" in the vector "Volume".
vtkIdType GetOneNeighborWithoutValue(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, long int value);

// Returns one neighbor of (x,y,z) with value different of
// "value".
vtkIdType GetOneNeighborWithoutValue(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, long int value);

// Merge together junctions that touch each other.
bool JunctionsMerge(std::list<vtkIdType> Junctions, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume);

// Track an edge starting at voxel (x,y,z) in the volume "Volume".
std::list<vtkIdType> GetEdgeStartingAt(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume);

// Returns an adjacency edge "edge_label".
long int GetOneAdjacentEdge(vtkSmartPointer<vtkPolyData> PolyData, long int edge_label, long int junction_label, bool *found_on_left);

// Track all the nodes and edges of a 3D structured thinned by
// the routine Thinning3D.
vtkSmartPointer<vtkPolyData> Skeletonization(vtkSmartPointer<vtkImageData> Image, _mitoObject *mitoObject);

// Replace two edges A and B attached to a node of degree 2
// with one edge C that corresponds to A + B. The degree of the
// node is set to -1 and A and B are moreved from PolyData.
bool MergeEdgesOfDegree2Nodes(vtkSmartPointer<vtkPolyData> PolyData, long int nedges_before_filtering, int *K);

// After merging edges attached to nodes of degree 2, the original
// edges must be deleted throught this routine.
void RemoveCellsWithInvalidNodes(vtkSmartPointer<vtkPolyData> PolyData, int *K);

/* ================================================================
   I/O ROUTINES
=================================================================*/

void SaveImageData(vtkSmartPointer<vtkImageData> Image, const char FileName[], bool _resample) {
    #ifdef DEBUG
        printf("Saving ImageData File...\n");
    #endif

    vtkSmartPointer<vtkImageFlip> Flip = vtkSmartPointer<vtkImageFlip>::New();
    Flip -> SetInputData(Image);
    Flip -> SetFilteredAxis(1);
    Flip -> PreserveImageExtentOn();
    Flip -> Update();

    vtkSmartPointer<vtkTIFFWriter> writer = vtkSmartPointer<vtkTIFFWriter>::New();
    writer -> SetFileName(FileName);

    // TIF writer does not support data of type double
    if (Image -> GetScalarType() == VTK_DOUBLE) {

        double range[2];
        Image -> GetScalarRange( range );
        vtkSmartPointer<vtkImageShiftScale> ShiftFilter = vtkSmartPointer<vtkImageShiftScale>::New();
        ShiftFilter -> SetInputData(Flip->GetOutput());
        ShiftFilter -> SetScale( 65535./(range[1]-range[0]));
        ShiftFilter -> SetShift( -range[0] );
        ShiftFilter -> SetOutputScalarTypeToUnsignedShort();
        ShiftFilter -> Update();
        
        writer -> SetInputData(ShiftFilter->GetOutput());

    } else {

        if (_resample) {
            #ifdef DEBUG
                printf("\tResampling data...%f\t%f\n",_dxy,_dz);
            #endif
            vtkSmartPointer<vtkImageResample> Resample = vtkSmartPointer<vtkImageResample>::New();
            Resample -> SetInterpolationModeToLinear();
            Resample -> SetDimensionality(3);
            Resample -> SetInputData(Flip->GetOutput());
            Resample -> SetAxisMagnificationFactor(0,1.0);
            Resample -> SetAxisMagnificationFactor(1,1.0);
            Resample -> SetAxisMagnificationFactor(2,_dz/_dxy);
            Resample -> Update();

            vtkSmartPointer<vtkImageData> ImageResampled = Resample -> GetOutput();
            ImageResampled -> SetSpacing(1,1,1);

            writer -> SetInputData(ImageResampled);
        } else {

            writer -> SetInputData(Flip->GetOutput());

        }
        
    }
    
    writer -> Write();

    #ifdef DEBUG
        printf("\tFile Saved!\n");
    #endif
}

void SavePolyData(vtkSmartPointer<vtkPolyData> PolyData, const char FileName[]) {

    #ifdef DEBUG
        printf("Saving PolyData...\n");
    #endif

    #ifdef DEBUG
        printf("\t#Points in PolyData file: %llu.\n",(vtkIdType)PolyData->GetNumberOfPoints());
    #endif

    vtkSmartPointer<vtkPolyDataWriter> Writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    #ifndef DEBUG
        Writer -> SetFileType(VTK_BINARY);
    #endif
    Writer -> SetFileName(FileName);
    Writer -> SetInputData(PolyData);
    Writer -> Write();

    #ifdef DEBUG
        printf("\tFile Saved!\n");
    #endif
}

void ExportGraphFiles(vtkSmartPointer<vtkPolyData> PolyData, long int nnodes, long int *ValidId, const char Prefix[]) {
    
    #ifdef DEBUG
        printf("Saving .coo file...\n");
    #endif

    char _fullpath[256];

    vtkPoints *Points = PolyData -> GetPoints();

    double r[3];
    long int node, exact_nnodes = 0;
    for (node = 0; node < nnodes; node++) {
        if (ValidId[node]>=0) exact_nnodes++;
    }

    sprintf(_fullpath,"%s.coo",Prefix);
    FILE *fcoo = fopen(_fullpath,"w");
    for (node = 0; node < nnodes; node++) {
        if (ValidId[node]>=0) {
            Points -> GetPoint(node,r);
            fprintf(fcoo,"%1.4f\t%1.4f\t%1.4f\n",_dxy*r[0],_dxy*r[1],_dz*r[2]);
        }
    }
    fclose(fcoo);

    double length;
    vtkIdType edge, npoints, i, j;
    sprintf(_fullpath,"%s.gnet",Prefix);
    FILE *fgnet = fopen(_fullpath,"w");
    fprintf(fgnet,"%ld\n",exact_nnodes);
    for (edge = 0; edge < PolyData -> GetNumberOfCells(); edge++) {
        npoints = PolyData -> GetCell(edge) -> GetNumberOfPoints();
        i = PolyData -> GetCell(edge) -> GetPointId(0);
        j = PolyData -> GetCell(edge) -> GetPointId(npoints-1);
        if ( ValidId[i] >= 0 && ValidId[j] >= 0 ) {
            length = GetEdgeLength(edge,PolyData);
            fprintf(fgnet,"%ld\t%ld\t%1.5f\n",ValidId[i],ValidId[j],length);
        }
    }
    fclose(fgnet);

}

void ExportNodes(vtkSmartPointer<vtkPolyData> PolyData, long int nnodes, long int *ValidId, _mitoObject *mitoObject) {
    
    #ifdef DEBUG
        if (_export_nodes_label) {
            printf("Exporting nodes and their labels...\n");
        } else {
            printf("Exporting nodes...\n");
        }
    #endif

    vtkPoints *Points = PolyData -> GetPoints();
    vtkSmartPointer<vtkAppendPolyData> Append = vtkSmartPointer<vtkAppendPolyData>::New();
    Append -> SetOutputPointsPrecision(vtkAlgorithm::DEFAULT_PRECISION);

    double r[3];
    long int node;
    char node_txt[16];
    for (node = 0; node < nnodes; node++) {
        if ( ValidId[node] >= 0 ) {
            Points -> GetPoint(node,r);
        
            vtkSmartPointer<vtkSphereSource> Node = vtkSmartPointer<vtkSphereSource>::New();
            if (_scale_polydata_before_save) {
                Node -> SetRadius(2*_dxy);
                Node -> SetCenter(_dxy*(r[0]+mitoObject->Ox),_dxy*(r[1]+mitoObject->Oy),_dz*(r[2]+mitoObject->Oz));
            } else {
                Node -> SetRadius(2.0);
                Node -> SetCenter(r[0]+mitoObject->Ox,r[1]+mitoObject->Oy,r[2]+mitoObject->Oz);
            }

            Node -> SetThetaResolution(36);
            Node -> SetPhiResolution(36);
            Node -> Update();
            Append -> AddInputData(Node->GetOutput());
            Append -> Update();

            if (_export_nodes_label) {
                sprintf(node_txt,"%ld",ValidId[node]);
                vtkSmartPointer<vtkVectorText> PolyText = vtkSmartPointer<vtkVectorText>::New();
                PolyText -> SetText(node_txt);
                PolyText -> Update();

                vtkSmartPointer<vtkTransform> T = vtkSmartPointer<vtkTransform>::New();
                if (_scale_polydata_before_save) {
                    T -> Translate(_dxy*(r[0]+1+mitoObject->Ox),_dxy*(r[1]+1+mitoObject->Oy),_dz*(r[2]+mitoObject->Oz));
                    T -> Scale(2*_dxy,2*_dxy,1);
                } else {
                    T -> Translate(r[0]+2+mitoObject->Ox,r[1]+mitoObject->Oy,r[2]+mitoObject->Oz);
                    T -> Scale(2,2,1);
                }

                vtkSmartPointer<vtkTransformPolyDataFilter> Trans = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
                Trans -> SetInputData(PolyText->GetOutput());
                Trans -> SetTransform(T);
                Trans -> Update();

                Append -> AddInputData(Trans -> GetOutput());
                Append -> Update();
            }
        }

    }
    
    SavePolyData(Append->GetOutput(),(mitoObject->FileName+"_nodes.vtk").c_str());
}


/* ================================================================
   IMAGE TRANSFORMATION
=================================================================*/

void CleanImageBoundaries(vtkSmartPointer<vtkImageData> ImageData) {
    #ifdef DEBUG
        printf("Cleaning the image boundaries...\n");
    #endif
    int p, q;
    int *Dim = ImageData -> GetDimensions();
    for ( p = Dim[0]; p--; ) {
        for ( q = Dim[1]; q--; ) {
            ImageData -> SetScalarComponentFromDouble(p,q,0,0,0);
            ImageData -> SetScalarComponentFromDouble(p,q,Dim[2]-1,0,0);
        }
    }
    for ( p = Dim[0]; p--; ) {
        for ( q = Dim[2]; q--; ) {
            ImageData -> SetScalarComponentFromDouble(p,0,q,0,0);
            ImageData -> SetScalarComponentFromDouble(p,Dim[1]-1,q,0,0);
        }
    }
    for ( p = Dim[1]; p--; ) {
        for ( q = Dim[2]; q--; ) {
            ImageData -> SetScalarComponentFromDouble(0,p,q,0,0);
            ImageData -> SetScalarComponentFromDouble(Dim[0]-1,p,q,0,0);
        }
    }
}

/* ================================================================
   AUXILIAR ROUTINES
=================================================================*/


bool IsBorder(vtkIdType id, vtkSmartPointer<vtkImageData> Image) {
    double r[3];
    int x, y, z;
    Image -> GetPoint(id,r);
    x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
    if(!Image->GetScalarComponentAsDouble(x,y,z,0))   return false;
    if(!Image->GetScalarComponentAsDouble(x+1,y,z,0)) return true;
    if(!Image->GetScalarComponentAsDouble(x-1,y,z,0)) return true;
    if(!Image->GetScalarComponentAsDouble(x,y+1,z,0)) return true;
    if(!Image->GetScalarComponentAsDouble(x,y-1,z,0)) return true;
    if(!Image->GetScalarComponentAsDouble(x,y,z+1,0)) return true;
    if(!Image->GetScalarComponentAsDouble(x,y,z-1,0)) return true;
    return false;
}

void SmoothEdgesCoordinates(vtkSmartPointer<vtkPolyData> PolyData, double sigma){

    double r1[3], r2[3], r3[3];
    long int id, line, nlines, n, s;
    vtkPoints *Points = PolyData -> GetPoints();

    vtkCell *Line;
    nlines = PolyData -> GetNumberOfCells();

    for (line = 0; line < nlines; line++){
        Line = PolyData -> GetCell(line);
        n = Line -> GetNumberOfPoints();
        double *X = new double[n];
        double *Y = new double[n];
        double *Z = new double[n];
        for (s = 0; s < int(sigma); s++) {
            for (id = 1; id < n-1; id++ ) {
                Points -> GetPoint(Line->GetPointId(id-1),r1);
                Points -> GetPoint(Line->GetPointId(id+0),r2);
                Points -> GetPoint(Line->GetPointId(id+1),r3);
                X[id] = 0.125*(r1[0]+6*r2[0]+r3[0]);
                Y[id] = 0.125*(r1[1]+6*r2[1]+r3[1]);
                Z[id] = (1.0/3.0)*(r1[2]+r2[2]+r3[2]);
            }
            for (id = 1; id < n-1; id++ ) {
                Points -> SetPoint(Line->GetPointId(id),X[id],Y[id],Z[id]);
            }
        }
        delete[] X; delete[] Y; delete[] Z;
    }
    PolyData -> Modified();
}

void GetVolumeFromVoxels(vtkSmartPointer<vtkImageData> Image, std::vector<attribute> *Atts) {
    double v;
    unsigned long int nv = 0;
    for (vtkIdType id=Image->GetNumberOfPoints();id--;) {
        v = Image -> GetPointData() -> GetScalars() -> GetTuple1(id);
        if (v) nv++;
    }
    attribute newAtt = {"Volume from voxels",nv * (_dxy * _dxy * _dz)};
    Atts -> push_back(newAtt);
}
/*
void GetVolumeFromSkeletonLength(vtkSmartPointer<vtkPolyData> PolyData, double *attributes) {
    double r1[3], r2[3], length = 0.0;
    vtkPoints *Points = PolyData -> GetPoints();
    for (vtkIdType edge=PolyData->GetNumberOfCells();edge--;) {
        length += GetEdgeLength(edge,PolyData);
    }
    attributes[1] = length;
    attributes[2] = length * (acos(-1.0)*pow(_rad,2));
}
*/
double GetEdgeLength(vtkIdType edge, vtkSmartPointer<vtkPolyData> PolyData) {
    double r1[3], r2[3];
    double length = 0.0;
    for (vtkIdType n = 1; n < PolyData->GetCell(edge)->GetNumberOfPoints(); n++) {
        PolyData -> GetPoint(PolyData->GetCell(edge)->GetPointId(n-1),r1);
        PolyData -> GetPoint(PolyData->GetCell(edge)->GetPointId(n  ),r2);
        length += sqrt(pow(_dxy*(r2[0]-r1[0]),2)+pow(_dxy*(r2[1]-r1[1]),2)+pow(_dz*(r2[2]-r1[2]),2));
    }
    return length;
}

/* ================================================================
   THINNING 3D
=================================================================*/

vtkSmartPointer<vtkPolyData> Thinning3D(vtkSmartPointer<vtkImageData> ImageData, _mitoObject *mitoObject) {

    #ifdef DEBUG
        SaveImageData(ImageData,(mitoObject->FileName+"_binary.tif").c_str());
    #endif

    vtkIdType N = ImageData -> GetNumberOfPoints();
    GetVolumeFromVoxels(ImageData,&mitoObject->attributes);

    int x, y, z;
    int count_thin = 0;
    double r[3], v, vl;
    vtkIdType ndels, id;
    ssThinVox *STV = new ssThinVox();

    CleanImageBoundaries(ImageData);

    int i, j, ***Vol = new int**[3];
    for (i = 0; i < 3; i++) {
        Vol[i] = new int*[3];
            for (j = 0; j < 3; j++) {
                Vol[i][j] = new int[3];
            }
    }

    #ifdef DEBUG
        printf("Starting thinning process...\n");
    #endif

    std::list<vtkIdType> ToBeDeleted;
    std::list<vtkIdType> OnTheSurface;
    std::list<vtkIdType>::iterator itId;

    do {
        count_thin = count_thin + 1;
        ndels = 0;
        for (id = N; id--;) {
            if (IsBorder(id, ImageData)) {
                OnTheSurface.insert(OnTheSurface.begin(),id);
            }
        }

        for (int direction = 6; direction--;) {
            for ( itId=OnTheSurface.begin(); itId!=OnTheSurface.end(); itId++) {
                id = *itId;
                ImageData -> GetPoint(id,r);
                x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
                v = ImageData -> GetScalarComponentAsDouble(x,y,z,0);
                if (v) {
                    for (i = 26; i--;) {
                        vl = ImageData -> GetScalarComponentAsDouble(x+ssdx[i],y+ssdy[i],z+ssdz[i],0);
                        Vol[1+ssdx[i]][1+ssdy[i]][1+ssdz[i]] = (vl) ? 1 : 0;
                    }
                    if ( STV -> match(direction,Vol) ) ToBeDeleted.insert(ToBeDeleted.begin(),id);
                }
            }
            itId = ToBeDeleted.begin();
            while (itId != ToBeDeleted.end()) {
                ndels++;
                ImageData -> GetPoint(*itId,r);
                x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
                ImageData -> SetScalarComponentFromDouble(x,y,z,0,0);
                ToBeDeleted.erase(itId++);
            }
        }

        #ifdef DEBUG
            printf("\t#Surface = %llu / #Deletions = %llu\n",(long long int)OnTheSurface.size(),ndels);
        #endif

        OnTheSurface.clear();

    } while(ndels>0 & count_thin<1);
    
    delete STV;

    #ifdef DEBUG
        printf("Thinning done!\n");
    #endif

    OnTheSurface.clear();
    ToBeDeleted.clear();

    return Skeletonization(ImageData,mitoObject);

}

/* ================================================================
   SKELETONIZATION
=================================================================*/

char GetNumberOfNeighborsWithoutValue(vtkSmartPointer<vtkImageData> Image, int x, int y, int z, long int value) {
    double r[3];
    char nn = 0;
    for (char k = 26; k--;) {
        if (Image->GetScalarComponentAsDouble(x+ssdx[k],y+ssdy[k],z+ssdz[k],0) != value) nn++;
    }
    return nn;
}

char GetNumberOfNeighborsWithoutValue(vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, int x, int y, int z, long int value) {
    double r[3];
    char nn = 0;
    vtkIdType idk;
    for (char k = 26; k--;) {
        idk = Image->FindPoint(x+ssdx[k],y+ssdy[k],z+ssdz[k]);
        if (Volume->GetTuple1(idk) != value) nn++;
    }
    return nn;
}

char GetNumberOfNeighborsWithValue(vtkSmartPointer<vtkImageData> Image, int x, int y, int z, long int value) {
    double r[3];
    char nn = 0;
    for (char k = 26; k--;) {
        if (Image->GetScalarComponentAsDouble(x+ssdx[k],y+ssdy[k],z+ssdz[k],0) == value) nn++;
    }
    return nn;
}

char GetNumberOfNeighborsWithValue(vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, int x, int y, int z, long int value) {
    double r[3];
    char nn = 0;
    vtkIdType idk;
    for (char k = 26; k--;) {
        idk = Image -> FindPoint(x+ssdx[k],y+ssdy[k],z+ssdz[k]);
        if (Volume->GetTuple1(idk) == value) nn++;
    }
    return nn;
}

vtkIdType GetOneNeighborWithValue(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, long int value) {
    vtkIdType idk;
    for (char k = 26; k--;) {
        idk = Image -> FindPoint(x+ssdx[k],y+ssdy[k],z+ssdz[k]);
        if (Volume -> GetTuple1(idk) == value) {
            return idk;
        }
    }
    return 0;  // We can do it because, by construction the voxel at id 0 should always be empty
}

vtkIdType GetOneNeighborWithoutValue(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume, long int value) {
    vtkIdType idk;
    for (char k = 26; k--;) {
        idk = Image -> FindPoint(x+ssdx[k],y+ssdy[k],z+ssdz[k]);
        if (Volume -> GetTuple1(idk) && Volume -> GetTuple1(idk) != value) {
            return idk;
        }
    }
    return 0;  // We can do it because, by construction the voxel at id 0 should always be empty
}

vtkIdType GetOneNeighborWithoutValue(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, long int value) {
    vtkIdType idk;
    for (char k = 26; k--;) {
        idk = Image -> FindPoint(x+ssdx[k],y+ssdy[k],z+ssdz[k]);
        if (Image->GetScalarComponentAsDouble(x+ssdx[k],y+ssdy[k],z+ssdz[k],0)!=value) return idk;
    }
    return 0;  // We can do it because, by construction the voxel at id 0 should always be empty
}

bool JunctionsMerge(std::list<vtkIdType> Junctions, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume) {
    char nk;
    double r[3];
    vtkIdType id;
    int x, y, z, k;
    bool _has_changed = false;
    std::list<vtkIdType>::iterator itId;
    long int junction_label, neigh_junction_label;
    for (itId=Junctions.begin(); itId!=Junctions.end(); itId++) {
        Image -> GetPoint(*itId,r);
        junction_label = Volume -> GetTuple1(*itId);
        x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
        for (k = 0; k < 26; k++) {
            id = Image->FindPoint(x+ssdx[k],y+ssdy[k],z+ssdz[k]);
            neigh_junction_label = Volume -> GetTuple1(id);
            if (junction_label < neigh_junction_label) {
                _has_changed = true;
                Volume -> SetTuple1(id,junction_label);
            }
        }
    }
    Volume -> Modified();
    return _has_changed;
}

std::list<vtkIdType> GetEdgeStartingAt(int x, int y, int z, vtkSmartPointer<vtkImageData> Image, vtkSmartPointer<vtkLongArray> Volume) {
    double r[3];
    vtkIdType idk;
    std::list<vtkIdType> Edge;
    do {
        idk = GetOneNeighborWithValue(x,y,z,Image,Volume,-1);
        if (idk) {
            Volume -> SetTuple1(idk,0);
            Edge.insert(Edge.end(),idk);
            Image -> GetPoint(idk,r);
            x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
        }
    } while(idk);
    return Edge;
}

long int GetOneAdjacentEdge(vtkSmartPointer<vtkPolyData> PolyData, long int edge_label, long int original_source, bool *_common_source) {
    vtkCell *Edge;
    long int source, target;
    for (long int edge = PolyData->GetNumberOfCells();edge--;) {
        if (edge != edge_label) {
            Edge = PolyData -> GetCell((vtkIdType)edge);
            source = Edge -> GetPointId(0);
            target = Edge -> GetPointId(Edge->GetNumberOfPoints()-1);
            if (original_source==source) {
                *_common_source = true;
                return edge;
            }
            if (original_source==target) {
                *_common_source = false;
                return edge;
            }
        }
    }
    return -1;
}

vtkSmartPointer<vtkPolyData> Skeletonization(vtkSmartPointer<vtkImageData> Image, _mitoObject *mitoObject) {

    #ifdef DEBUG
        SaveImageData(Image,(mitoObject->FileName+"_thinned.tif").c_str());
    #endif

    // Creating raw polyData
    vtkSmartPointer<vtkPolyData> PolyData = vtkPolyData::New();
    return PolyData;
}

bool MergeEdgesOfDegree2Nodes(vtkSmartPointer<vtkPolyData> PolyData, long int nedges_before_filtering, int *K) {

    // Given this edge: (source) o---->----o (target), it's necessary
    // to check whether either source or target are nodes of degree 2.

    std::list<vtkIdType> Merg;
    std::list<vtkIdType>::iterator itId;

    vtkCell *Edge, *NeighEdge;
    long int i, source, target, edge, edge_length, neigh_edge, neigh_edge_length;
    bool _common_source, _merged = false;

    for (edge = 0; edge < PolyData -> GetNumberOfCells(); edge++) {

        Edge = PolyData -> GetCell(edge);
        edge_length = Edge -> GetNumberOfPoints();
        source = (long int)Edge -> GetPointId(0);
        target = (long int)Edge -> GetPointId(Edge->GetNumberOfPoints()-1);

        if (K[source]==2 && source!=target) {

            //  [  neigh  | original ]
            //  o----?----o---->-----o
            //            s          t

            _merged = true;
            neigh_edge = GetOneAdjacentEdge(PolyData,edge,source,&_common_source);
            NeighEdge = PolyData -> GetCell((vtkIdType)neigh_edge);
            neigh_edge_length = NeighEdge -> GetNumberOfPoints();

            if (_common_source) {//  o----<----o---->----o
                for (i=0;i<edge_length;i++)
                    Merg.insert(Merg.begin(),PolyData->GetCell(edge)->GetPointId(i));
                for (i=1;i<neigh_edge_length;i++)
                    Merg.insert(Merg.end(),PolyData->GetCell(neigh_edge)->GetPointId(i));
            } else { //  o---->----o---->----o
                for (i=0;i<neigh_edge_length;i++)
                    Merg.insert(Merg.end(),PolyData->GetCell(neigh_edge)->GetPointId(i));                
                for (i=1;i<edge_length;i++)
                    Merg.insert(Merg.end(),PolyData->GetCell(edge)->GetPointId(i));
            }
            K[source] = -1; // Tagged as not valid node

        } else if (K[target]==2 && source!=target) {

            //  [ original |  neigh  ]
            //  o---->-----o----?----o
            //  s          t

            _merged = true;
            neigh_edge = GetOneAdjacentEdge(PolyData,edge,target,&_common_source);
            NeighEdge = PolyData -> GetCell((vtkIdType)neigh_edge);
            neigh_edge_length = NeighEdge -> GetNumberOfPoints();

            if (_common_source) {//  o---->----o---->----o
                for (i=0;i<edge_length;i++)
                    Merg.insert(Merg.end(),PolyData->GetCell(edge)->GetPointId(i));
                for (i=1;i<neigh_edge_length;i++)
                    Merg.insert(Merg.end(),PolyData->GetCell(neigh_edge)->GetPointId(i));
            } else {//  o---->----o----<----o
                for (i=0;i<edge_length;i++)
                    Merg.insert(Merg.end(),PolyData->GetCell(edge)->GetPointId(i));
                for (i=neigh_edge_length-1;i--;)
                    Merg.insert(Merg.end(),PolyData->GetCell(neigh_edge)->GetPointId(i));                
            }
            K[target] = -1; // Tagged as not valid node
        }

        if (_merged) {
            i = 0;
            // Removing the two pieces that were merged together
            PolyData -> BuildLinks();
            PolyData -> DeleteCell(edge);
            PolyData -> DeleteCell(neigh_edge);
            PolyData -> RemoveDeletedCells();
            // Adding the new "merged" edge
            vtkIdType IdList[Merg.size()];
            for (itId=Merg.begin(); itId!=Merg.end(); itId++) {
                IdList[i] = *itId;
                i++;
            }
            PolyData -> InsertNextCell(VTK_POLY_LINE,Merg.size(),IdList);
            Merg.clear();
            return true;
        }

    }
    return false;
}

long int LabelConnectedComponents(vtkSmartPointer<vtkImageData> ImageData, vtkSmartPointer<vtkDataArray> Volume, std::vector<long int> &CSz, int ngbh, double threshold) {

    #ifdef DEBUG
        printf("\tCalculating connected components...\n");
    #endif

    int *Dim = ImageData -> GetDimensions();

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

    Volume -> CopyComponent(0,ImageData->GetPointData()->GetScalars(),0);
    Volume -> Modified();

    label = 0;
    while (find) {
        for (s = 0; s < CurrA->GetNumberOfIds(); s++) {
            ido = CurrA -> GetId(s);
            ImageData -> GetPoint(ido,r);
            x = (int)r[0]; y = (int)r[1]; z = (int)r[2];
            for (i = 0; i < ngbh; i++) {
                id = ImageData -> FindPoint(x+ssdx_sort[i],y+ssdy_sort[i],z+ssdz_sort[i]);
                v = Volume -> GetTuple1(id);
                if (v > threshold) {
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
                if (v > threshold) {
                    find = true;
                    ro[0] = id;
                    break;
                }
            }
            if (label) {
                CSz.push_back(scluster);
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

    #ifdef DEBUG
        printf("\tNumber of detected components: %ld\n",(long int)CSz.size());
    #endif

    //return (long int)CSz->GetNumberOfTuples();
    return (long int)CSz.size();
}
