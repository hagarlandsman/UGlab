#ifndef UGBURDEN_H
#define UGBURDEN_H

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <iostream>
#include <cmath>
#include <vector>
#include <iostream>
#include "vcpkg/installed/x64-linux/include/nanoflann.hpp"
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <TGraph.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraph2D.h>
#include <TRandom3.h>
#include <TApplication.h>
#include "TPolyLine3D.h"
#include <sstream>
#include "TH2F.h"
#include <vector>
#include <string>



using namespace std;
using namespace nanoflann;



class UGburden : public TObject {

struct PointCloud {
    std::vector<std::array<double, 2>> xy_points;  // Store (x, y) points
    std::vector<double> z_values;                  // Corresponding z values

    // Required by nanoflann: returns the number of data points
    inline size_t kdtree_get_point_count() const { return xy_points.size(); }

    // Returns a specific dimension value for a point
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return xy_points[idx][dim];  // dim: 0 for x, 1 for y
    }

    // Bounding box optimization not required
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, PointCloud>,  // Distance metric
    PointCloud,                                        // Point cloud data
    2                                                   // Dimensionality (x, y)
> KDTree;


private:
    double refEast = 0.0;
    double refNorth = 0.0;
    PointCloud cloud;
    TString filename;
    KDTree* kd_tree = nullptr;  // KDTree as a pointer
    int n_cols;
    double refLat, refLong;

public:
    UGburden();
    UGburden(TString filename) ;
                             // Default constructor
           // Constructor with parameters
    virtual ~UGburden();                   // Destructor
    int load(TString filename_) ;

    void buildKDTree();  // Method to build KD-Tree
    double get_one(double x_in, double y_in) ;
    double get_one() ;
    void draw_surface();


    void setrefLong(double val) {refLong = val;};
void setrefLat(double val) {refLat = val;};
double getrefLong() { return refLong;};
double getrefLat() { return refLat;};


    TGraph2D* draw_azimuth(TVector3 start_point, double zenith, double azimuth) ;

    TVector3 propagateUntilZExceeds(
        TVector3 start_point,
        TVector3 direction);
    double getRburden(
        TVector3 start_point,
        TVector3 direction);
   double getRburden(
        TVector3 start_point,
        double zenith, double azimuth);


    int getMap();
    void latLonToNorthingEasting(
        double refLat, double refLon, double lat, double lon,
        double& northing, double& easting
    );



    int check_n_columns(); // Add logic for this helper

     // Print method

    ClassDef(UGburden, 1);                 // ROOT I/O macro
};

double NorthingEastingToangleNorth(double northing, double easting);

double NorthingEastingToangleEast(double northing, double easting);


#endif // UGBURDEN_H
