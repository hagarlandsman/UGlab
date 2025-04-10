#include "UGburden.h"


#include <TClass.h>
#include <TObject.h>





/*


In LinkDef.cpp put this line:
#pragma link C++ class UGburden;


Then, do:
rootcling -f LinkDef_rdict.cxx -rml -c LinkDef.cpp UGburden.h

and compile with:
g++ -o main main.cpp LinkDef_rdict.cxx UGburden.cpp -Ivcpkg/installed/x64-linux/include -pthread  `root-config --cflags --libs`
*/


ClassImp(UGburden)


/*
 .L UGburden.cpp+

TString map_file = "R60km_0.1.csv";
UGburden burden(map_file);
TVector3 start_point(0, 0, -274.505);
burden.draw_surface();
burden.draw_azimuth(start_point, TMath::Pi()/2, 0);

burden.get_one(10,10)


*/

UGburden::UGburden()
{
    //ctor
}


UGburden::UGburden(TString filename_) {

    refLat = 32.597179;  // kochav - based on east/north 249943.68; 722580.77;
    refLong = 35.529270;   // kochav  -   based on east/north 249943.68; 722580.77;


    kd_tree = nullptr;
    filename = filename_;
    if (getMap() != 1) {
        std::cerr << "Error loading map from file: " << filename << std::endl;
    } else {
        std::cout << "Successfully loaded map from file: " << filename << std::endl;
    }
    buildKDTree();

    typedef KDTreeSingleIndexAdaptor<
        L2_Simple_Adaptor<double, PointCloud>,  // Distance metric
        PointCloud,                             // Point cloud data
        2                                       // Dimensionality (x, y)
    > KDTree;
    KDTree kd_tree(2, cloud, KDTreeSingleIndexAdaptorParams(10));
    kd_tree.buildIndex();

    // Optional: Perform operations on cloud here if needed

 printf ("===> Units of coordinates are  {northing , easting } w.r.t refernce  Lat= %f, Long=%f, \n",refLat,refLong);

}



UGburden::~UGburden()
{
    if (kd_tree) {
        delete kd_tree;  // Clean up the KD-Tree
    }
}




void UGburden::buildKDTree() {
    if (kd_tree) {
        delete kd_tree;  // Clean up previous KD-Tree if necessary
    }

    kd_tree = new KDTree(2, cloud, KDTreeSingleIndexAdaptorParams(10));
    kd_tree->buildIndex();  // Build the index
}



TGraph2D* UGburden::draw_azimuth(TVector3 start_point, double zenith, double azimuth) {
        TGraph2D* pathGraph = new TGraph2D();
        pathGraph->SetPoint(0, start_point.X(), start_point.Y(), start_point.Z());

        TVector3 direction(1,1,1);
        direction.SetMag(10000);                // Set vector magnitude
        direction.SetTheta(zenith); // Set the zenith (polar) angle in radians
        direction.SetPhi(azimuth);
        pathGraph->SetPoint(1, direction.X(), direction.Y(), direction.Z());
        pathGraph->SetLineColor(kRed);
        pathGraph->SetLineWidth(3);
        pathGraph->SetMarkerStyle(20);
        pathGraph->Draw("pl same");
        //canvas->Update();
        return pathGraph;

}
void UGburden::draw_surface(){
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
    canvas->Update();
    canvas->SaveAs("canvas1.png");


}
double UGburden::get_one() {
    double x_in, y_in;
    std::cout << "Enter x, y coordinates to find nearest z: ";
    std::cin >> x_in >> y_in;
    printf("Read: %f, %f\n", x_in, y_in);
    return get_one(x_in,y_in);

}


double UGburden::get_one(double x_in, double y_in) {

    // Check if point cloud is empty
    if (cloud.xy_points.empty()) {
        std::cout << "No data points available.\n";
        return -1;
    }

    // Perform nearest neighbor search
    const size_t num_results = 1;
    unsigned int ret_index;
    double out_dist_sqr;

    double query_pt[2] = {x_in, y_in};

    kd_tree->knnSearch(&query_pt[0], num_results, &ret_index, &out_dist_sqr);

    // Output nearest point and corresponding z value
    std::cout << "Nearest point: (" << cloud.xy_points[ret_index][0] << ", "
              << cloud.xy_points[ret_index][1] << ")\n";
    std::cout << "Distance squared: " << out_dist_sqr << std::endl;
    std::cout << "Distance: " << std::sqrt(out_dist_sqr) << std::endl;
    std::cout << "Corresponding z value: " << cloud.z_values[ret_index] << std::endl;
    return (std::sqrt(out_dist_sqr) );
}

    double UGburden::getRburden(
        TVector3 start_point,
        TVector3 direction)
        {
            printf("Start Vector components: X = %.2f, Y = %.2f, Z = %.2f\n", start_point.X(), start_point.Y(), start_point.Z());
            printf("direction Vector components: X = %.2f, Y = %.2f, Z = %.2f\n", direction.X(), direction.Y(), direction.Z());

            TVector3 final_point = UGburden::propagateUntilZExceeds(start_point, direction);
            printf("End Vector components: X = %.2f, Y = %.2f, Z = %.2f\n", final_point.X(), final_point.Y(), final_point.Z());

            double distance = (start_point-final_point).Mag();
            return distance;


        }

    double UGburden::getRburden(
        TVector3 start_point,
        double zenith, double azimuth)
        {
            TVector3 direction(1,1,1);
            direction.SetMag(1);                // Set vector magnitude
            direction.SetTheta(zenith); // Set the zenith (polar) angle in radians
            direction.SetPhi(azimuth);  // Set the azimuthal angle in radians


            return getRburden(start_point, direction);


        }


TVector3 UGburden::propagateUntilZExceeds(
    TVector3 start_point, TVector3 direction
) {
    double u = 0.0;
    TVector3 current_point = start_point;
    double step_size = 0.1;

    while (true) {
        double query_pt[2] = {current_point.X(), current_point.Y()};
      //  printf ("x= %f y= %f, angle = %f \n",current_point.X(), current_point.Y(),NorthingEastingToangleNorth(current_point.X(),current_point.Y())* 180.0 / M_PI);
        unsigned int nearest_idx;
        double out_dist_sqr;
        kd_tree->knnSearch(query_pt, 1, &nearest_idx, &out_dist_sqr);

        const double nearest_z = cloud.z_values[nearest_idx];

        if (current_point.Z() > nearest_z) break;
        if (std::sqrt(out_dist_sqr) > 10000) break;

        u += step_size;
        current_point = start_point + u * direction;
    }

    return current_point;
}

int UGburden::getMap() {


    int n_col = check_n_columns();

    TFile* output_file = new TFile("output_tree.root", "RECREATE");
    double northing, easting;
    double angle_from_north_deg, angle_from_east_deg ;
    float lat2, lon2, elav2, dx, dy, land;
    TTree* T10 = new TTree("T10", "T10");

    if (n_col == 4) {
        T10->ReadFile(filename, "lat2:lon2:elav2:dx:dy");
    } else if (n_col == 5) {
        T10->ReadFile(filename, "lat2:lon2:elav2:dx:dy:land");
    } else {
        std::cerr << "Illegal file structure: " << filename.Data() << std::endl;
        return -1;
    }

    T10->SetBranchAddress("lat2", &lat2);
    T10->SetBranchAddress("lon2", &lon2);
    T10->SetBranchAddress("elav2", &elav2);
    TBranch* branch_northing = T10->Branch("northing", &northing, "northing/D");
    TBranch* branch_easting = T10->Branch("easting", &easting, "easting/D");
    TBranch* branch_fromnorth = T10->Branch("from_north_deg", &angle_from_north_deg, "from_north_deg/D");
    TBranch* branch_fromeast = T10->Branch("from_east_deg", &angle_from_east_deg, "from_east_deg/D");

    Long64_t nEntries = T10->GetEntries();
    double x_min = 99999, y_min = 99999;

    for (Long64_t i = 0; i < nEntries; i++) {
        T10->GetEntry(i);
        latLonToNorthingEasting(refLat, refLong, lat2, lon2, northing, easting);
        double angle_from_north_rad = std::atan2(easting, northing);
        double angle_from_north_deg = angle_from_north_rad * 180.0 / M_PI;
        double angle_from_east_rad = std::atan2(northing, easting);
        double angle_from_east_deg = angle_from_east_rad * 180.0 / M_PI;

        branch_easting->Fill();
        branch_northing->Fill();
        cloud.xy_points.push_back({northing , easting });
        cloud.z_values.push_back(elav2);
        if (fabs(northing ) < y_min && fabs(easting) < x_min) {
            x_min = easting ;
            y_min = northing ;
        }
    }

    std::cout << "Min = " << x_min << ", " << y_min << std::endl;
    T10->Write();
    output_file->Close();

    std::cout << "Got " << nEntries << " entries." << std::endl;
    return 1;
}

const double EARTH_RADIUS_METERS = 6378137.0; // WGS-84 radius
const double DEG_TO_RAD = M_PI / 180.0;

// Converts latitude and longitude (in degrees) to northing and easting

void  UGburden::latLonToNorthingEasting(double refLat, double refLon, double lat, double lon, double &northing, double &easting) {
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

double  NorthingEastingToangleNorth(double northing, double easting) {
    // Angle east of north (measured clockwise from north)
    double angle_from_north_rad = std::atan2(easting, northing);
    double angle_from_north_deg = angle_from_north_rad * 180.0 / M_PI;
    return angle_from_north_rad;
}

double  NorthingEastingToangleEast(double northing, double easting) {

    // Angle north of east (measured counter-clockwise from east)
    double angle_from_east_rad = std::atan2(northing, easting);
    double angle_from_east_deg = angle_from_east_rad * 180.0 / M_PI;
    return angle_from_east_rad;
}
int UGburden::check_n_columns(){
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

    return comma_count;
}



