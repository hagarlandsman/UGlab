#include <iostream>
#include <cmath>


void  latLonToNorthingEasting(double refLat, double refLon, double lat, double lon, double &northing, double &easting);



int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <refLat> <refLon> <lat> <lon>\n";
        return 1;
    }

    double refLat = std::atof(argv[1]);
    double refLon = std::atof(argv[2]);
    double lat    = std::atof(argv[3]);
    double lon    = std::atof(argv[4]);

    double northing, easting;
    latLonToNorthingEasting(refLat, refLon, lat, lon, northing, easting);

    std::cout << "Northing: " << northing << " meters\n";
    std::cout << "Easting: " << easting << " meters\n";




        // Angle east of north (measured clockwise from north)
        double angle_from_north_rad = std::atan2(easting, northing);
        double angle_from_north_deg = angle_from_north_rad * 180.0 / M_PI;

        // Angle north of east (measured counter-clockwise from east)
        double angle_from_east_rad = std::atan2(northing, easting);
        double angle_from_east_deg = angle_from_east_rad * 180.0 / M_PI;

        std::cout << "Angle east of north: " << angle_from_north_deg << "°" << std::endl;
        std::cout << "Angle north of east: " << angle_from_east_deg << "°" << std::endl;



    return 0;
}
