// g++ -o knn_search read_geometry.c -Ivcpkg/installed/x64-linux/include -pthread `root-config --cflags --libs`

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



int get_map(TString file_name,  PointCloud& cloud) ;
int check_n_columns(TString filename);
void  latLonToNorthingEasting(double refLat, double refLon, double lat, double lon, double &northing, double &easting) ;
// Constants
const double EARTH_RADIUS_METERS = 6378137.0; // WGS-84 radius
const double DEG_TO_RAD = M_PI / 180.0;

// Converts latitude and longitude (in degrees) to northing and easting

void  latLonToNorthingEasting(double refLat, double refLon, double lat, double lon, double &northing, double &easting) {
    // Convert degrees to radians
    double refLatRad = refLat * DEG_TO_RAD;
    double refLonRad = refLon * DEG_TO_RAD;
    double latRad = lat * DEG_TO_RAD;
    double lonRad = lon * DEG_TO_RAD;
    // Differences in coordinates
    double dLat = latRad - refLatRad;
    double dLon = lonRad - refLonRad;
    // Compute northing and easting using a simple equirectangular projection
    northing = dLat * EARTH_RADIUS_METERS;
    easting = dLon * EARTH_RADIUS_METERS * cos(refLatRad);
}


void read_geometry(){

   // latLonToNorthingEasting(refLat, refLon, lat, lon, northing, easting);




}
/*

void get_ray() {
    TRandom3 rand;
    rand.SetSeed(0);
    Double_t x,y,z;
    rand.Sphere(x,y,z,1);
    TVector3 v3(x,y,z);
    Double_t theta = v3.Theta();
    Double_t phi = v3.Phi();

    double u = 0.0;
    TVector3 current_point = start_point;
    TVector3 nearest_point = findNearestPoint(start_point, points);

    double step_size = 0.1;
     while (current_point.Z() <= nearest_point.Z()) {
        u += step_size;
        current_point = start_point + u * direction;  // Parametric equation: r(u) = start_point + u * direction
    }


*/



TVector3 generateRandomVectorOnUpperHemisphere() {
    // Generate random zenith angle (theta) between 0 and pi/2
    double theta = acos(2.0 * rand() / RAND_MAX - 1.0);  // Uniform distribution in [0, π/2]

    // Generate random azimuthal angle (phi) between 0 and 2π
    double phi = 2.0 * M_PI * rand() / RAND_MAX;  // Uniform distribution in [0, 2π]

    // Convert spherical coordinates to Cartesian coordinates (x, y, z)
    double x = sin(theta) * cos(phi);
    double y = sin(theta) * sin(phi);
    double z = cos(theta);

    return TVector3(x, y,fabs(z));
}

//TVector3 propagateUntilZExceeds(const TVector3& start_point,const TVector3& direction, const nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2>& kdtree, const PointCloud& cloud, std::vector<TVector3>& path) {

TVector3 propagateUntilZExceeds( TVector3 start_point, TVector3 direction, const nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2>& kdtree, const PointCloud& cloud, std::vector<TVector3>& path) {
    // Initialize the parametric variable u
    double u = 0.0;
    TVector3 current_point = start_point;

    // Propagate in the direction of the vector until z exceeds nearest point's z value
    double step_size = 0.1;  // Adjust step size for finer control



     while (true) {
        // Query KD-Tree for nearest neighbor to current point (x, y only)
        double query_pt[2] = {current_point.X(), current_point.Y()};  // 2D query (x, y)
        unsigned int nearest_idx;
        double out_dist_sqr;
        kdtree.knnSearch(query_pt, 1, &nearest_idx, &out_dist_sqr);

        // Get the nearest point's z value
        const double nearest_z = cloud.z_values[nearest_idx];
        //path.push_back(current_point);
        //printf("Current point coordinates: (%f, %f, %f)\n", current_point.X(), current_point.Y(), current_point.Z());

        // Break loop if z value exceeds nearest point's z value
        if (current_point.Z() > nearest_z) {
         //   printf ("Positive \n");
            break;
        }
        if (std::sqrt(out_dist_sqr) > 10000) {
        //    printf ("too far \n");
            break;
        }
        // Increment the parametric variable and propagate the point
        u += step_size;
        current_point = start_point + u * direction;  // Parametric equation: r(u) = start_point + u * direction
    }

    return current_point;


}
int get_map(TString file_name,  PointCloud& cloud) {
//    double refLat = 31.8944;    // rehovot
//    double refLon = 34.8115;    // rehovot

    double refEast =  0; //249943.68;
    double refNorth = 0; //722580.77;

    double refLat = 32.597179;  // Kochav
    double refLon = 35.529270; // Kochav
    int n_col = check_n_columns(file_name);
    TFile *output_file = new TFile("output_tree.root", "RECREATE");

    double northing, easting;
    float lat2,lon2,elav2, dx,dy,land;
    TTree *T10 = new TTree("T10", "T10");
//    T10->ReadFile(file_name, "lat2:lon2:elav2:dx:dy:land");
    cout<<"ncol="<<n_col<<endl;
    if (n_col == 4) {
        T10->ReadFile(file_name, "lat2:lon2:elav2:dx:dy");
    }
    else if (n_col == 5) {
     T10->ReadFile(file_name, "lat2:lon2:elav2:dx:dy:land");

    }
    else {
    printf ("ileggal file structure %s ",file_name.Data());
    return -1;

    }
    T10->SetBranchAddress("lat2", &lat2);
    T10->SetBranchAddress("lon2", &lon2);
    T10->SetBranchAddress("elav2", &elav2);
    T10->SetBranchAddress("land", &land);
    TBranch *branch_easting = T10->Branch("northing", &northing, "northing/D");
    TBranch *branch_northing = T10->Branch("easting", &easting, "easting/D");

    Long64_t nEntries = T10->GetEntries();
    double x_min=99999;
    double y_min=99999;

    for (Long64_t i = 0; i < nEntries; i++) {
        T10->GetEntry(i);
        latLonToNorthingEasting(refLat, refLon, lat2, lon2, northing, easting);
        branch_easting->Fill();
        branch_northing->Fill();
       // printf (" %f %f, %f northing,  %feasting \n",elav2, lat2,(northing-refNorth)/1000,(easting-refEast)/1000);
        cloud.xy_points.push_back({northing-refNorth, easting-refEast});
            cloud.z_values.push_back(elav2);
        if (fabs(northing-refNorth)<y_min and fabs(easting-refEast)<x_min) {
            x_min = easting - refEast;
            y_min = northing - refNorth;

        }
        }
        printf("Min= %f %f \n",x_min,y_min);
        T10->Write();
    output_file->Close();
    printf ("Got %lld entries \n",nEntries);

    return(1);
}

int check_n_columns(TString filename){
    std::ifstream file(filename);
    std::string line;
    int comma_count = 0;
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return 0;
    }


        std::getline(file, line);  // Skip the first line
        std::getline(file, line);  // Read the second line

        cout<<"line = "<<line<<endl;
        for (char ch : line) {
            if (ch == ',') {
                comma_count++;
            }
        }
        // Check number of columns


        return (comma_count);

}


void get_one(const nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2>& kdtree, const PointCloud& cloud) {
    double x_in, y_in;
    std::cout << "Enter x, y coordinates to find nearest z: ";
    std::cin >> x_in >> y_in;
    printf("Read: %f, %f\n", x_in, y_in);

    // Check if point cloud is empty
    if (cloud.xy_points.empty()) {
        std::cout << "No data points available.\n";
        return;
    }

    // Perform nearest neighbor search
    const size_t num_results = 1;
    unsigned int ret_index;
    double out_dist_sqr;

    double query_pt[2] = {x_in, y_in};

    kdtree.knnSearch(&query_pt[0], num_results, &ret_index, &out_dist_sqr);

    // Output nearest point and corresponding z value
    std::cout << "Nearest point: (" << cloud.xy_points[ret_index][0] << ", "
              << cloud.xy_points[ret_index][1] << ")\n";
    std::cout << "Distance squared: " << out_dist_sqr << std::endl;
    std::cout << "Distance: " << std::sqrt(out_dist_sqr) << std::endl;
    std::cout << "Corresponding z value: " << cloud.z_values[ret_index] << std::endl;
}



//TVector3 propagateUntilZExceeds(const TVector3& start_point, const TVector3& direction, const nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2>& kdtree) {
  //  TVector3 propagated_point = propagateUntilZExceeds(start_point, direction, kdtree);

void generateAndDrawPaths(int nPaths, const nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2>& kd_tree, const PointCloud& cloud) {
    // Create a new canvas for 3D plotting
    std::ofstream outFile("output.txt",std::ios::app);


    // Plot the 3D scatter plot of the points in the cloud
/*
    TGraph2D* gr = new TGraph2D();
    // Add data points from the KD-Tree (x, y, z) to the graph
    for (size_t i = 0; i < cloud.xy_points.size(); ++i) {
        gr->SetPoint(i, cloud.xy_points[i][0], cloud.xy_points[i][1], cloud.z_values[i]);
    }

      TCanvas* canvas = new TCanvas("canvas", "Generated Paths", 800, 600);
      gr->Draw("surf3");
     Create a TGraph2D to hold the 3D paths (x, y, z)
  */
    TGraph2D* pathGraph = new TGraph2D();

    TH2F *h2 = new TH2F("h2", "Azimuth vs Zenith", 100, -180, 180, 100, 0, 90);

    // Define starting point and KD-Tree (for simplicity, I'm assuming kd_tree and cloud are already populated)
    TVector3 start_point(0, 0, -274.505);  // Starting at some point (x, y, z)

    // Generate 100 random paths
    std::vector<TVector3> path;
    TVector3 direction(0,0,-100);

    for (int i = 0; i < nPaths; ++i) {
        //printf ("N= %d \n",i);
        // Generate random direction for each path
        direction = generateRandomVectorOnUpperHemisphere();

        // Propagate the path
        TVector3 final_point = propagateUntilZExceeds(start_point, direction, kd_tree, cloud,path);


        double distance = (start_point-final_point).Mag();
        h2->Fill(direction.Phi()*180/TMath::Pi(),direction.Theta()*180/TMath::Pi(),distance);  // Phi is azimuth
        //outFile << direction.Phi()*180/TMath::Pi()<< "\t" << direction.Theta()*180/TMath::Pi() << "\t" << distance << std::endl;  // Tab-separated values

        // Add final point to graph
        pathGraph->SetPoint(i, final_point.X(), final_point.Y(), final_point.Z());
        if (i%100 == 0) printf ("%d). Final point %f,%f,%f,  R=%f \n",i,final_point.X(), final_point.Y(), final_point.Z(),distance);

        /*
        TPolyLine3D* line = new TPolyLine3D(2);  // 2 points for each line


        line->SetPoint(0,start_point.X(),start_point.Y(),start_point.Z());
        line->SetPoint(1,final_point.X(), final_point.Y(), final_point.Z());
        line->SetLineColor(kRed);
        line->SetLineWidth(3);

        line->Draw();
        */
    }

    // Draw the 3D paths using the TGraph2D object
    pathGraph->SetMarkerSize(1);
    pathGraph->SetMarkerStyle(1);

    pathGraph->SetMarkerColor(kRed);

    pathGraph->Draw("P same");

    // Display the canvas
//    canvas->Update();
    TCanvas* canvas2 = new TCanvas("canvas2", "Generated Paths", 800, 600);

    h2->Draw("colz");
    outFile.close();

}

void get_one_and_draw(const nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2>& kd_tree, const PointCloud& cloud) {
// Container to hold the path of propagated points
    std::vector<TVector3> path;

    TVector3 start_point(0, 0, -274.505);  // Starting at (1.0, 1.0, 0.5)
    TVector3 direction(0.0, 10.0, 1.0);    // Propagate along z-axis (z direction)
    direction = generateRandomVectorOnUpperHemisphere();

    TVector3 final_point = propagateUntilZExceeds(start_point, direction, kd_tree, cloud,path);

    // Input point for which we need the closest neighbor
    std::cout << "Final point after propagation: ("
              << final_point.X() << ", " << final_point.Y() << ", " << final_point.Z() << ")\n";

    TCanvas* canvas = new TCanvas("canvas", "3D Path Plot", 800, 600);
    TGraph2D* gr = new TGraph2D();

    // Add data points from the KD-Tree (x, y, z) to the graph
    for (size_t i = 0; i < cloud.xy_points.size(); ++i) {
        gr->SetPoint(i, cloud.xy_points[i][0], cloud.xy_points[i][1], cloud.z_values[i]);
    }

    // Plot the 3D scatter plot of the points in the cloud
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kRed);
    gr->Draw("surf3");

    // Add the propagated path points to the graph
    TGraph2D* pathGraph = new TGraph2D();
    for (size_t i = 0; i < path.size(); ++i) {
        pathGraph->SetPoint(i, path[i].X(), path[i].Y(),path[i].Z());
        //printf ("%f %f %f\n",path[i].X(), path[i].Y(),path[i].Z());
    }

    pathGraph->SetMarkerColor(kRed);
    pathGraph->SetMarkerSize(20);
    pathGraph->Draw("P same");

    canvas->Update();



}
int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);

    PointCloud cloud;
//   get_map("R10km_0.005.csv",cloud);
   get_map("R60km_0.1.csv",cloud);

   //get_map("R10km_0.1.csv",cloud);
    // Build KD-Tree
    typedef KDTreeSingleIndexAdaptor<
        L2_Simple_Adaptor<double, PointCloud>,  // Distance metric
        PointCloud,                             // Point cloud data
        2                                       // Dimensionality (x, y)
    > KDTree;
    KDTree kd_tree(2, cloud, KDTreeSingleIndexAdaptorParams(10));
    kd_tree.buildIndex();

          //  get_one(kd_tree, cloud);
//    get_one_and_draw(kd_tree,cloud);
    generateAndDrawPaths(10000,kd_tree,cloud);
    app.Run();  // Start the ROOT application (GUI)


    return 0;
}


