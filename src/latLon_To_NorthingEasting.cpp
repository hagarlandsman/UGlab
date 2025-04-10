#include <iostream>
#include <cmath>

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

int latLon_To_NorthingEasting() {
    //double refLat = 40.7128;   // Reference latitude (e.g., New York City)
    //double refLon = -74.0060;  // Reference longitude
    //double lat = 34.0522;      // Target latitude (e.g., Los Angeles)
    //double lon = -118.2437;    // Target longitude

    double refLat = 31.8944;    // rehovot
    double refLon = 34.8115;    // rehovot

    double lat = 32.7940; // haifa
    double lon = 34.9896; //haifa

    double northing, easting;
    latLonToNorthingEasting(refLat, refLon, lat, lon, northing, easting);

    std::cout << "Northing: " << northing << " meters" << std::endl;
    std::cout << "Easting: " << easting << " meters" << std::endl;

    return 0;
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