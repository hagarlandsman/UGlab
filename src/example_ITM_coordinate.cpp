#include "ITM_coordinate.h"
#include <iostream>
#include <proj.h>
#include <cmath>

int main() {
    std::cout<<" Code not working well, giving strange results \n";

    double lat = 31.73436;
    double lon = 35.20452;
    double refLat = 31.73200;
    double refLon = 35.20000;

    double deltaE, deltaN;
    RelativeOffsetFromReference(lat, lon, refLat, refLon, deltaE, deltaN);

    std::cout << "Offset from reference: Easting = " << deltaE << " m, Northing = " << deltaN << " m\n";


//    print_axis_units("EPSG:2039"); // ITM
//    print_axis_units("EPSG:4326"); // WGS84


    double northing = 722580.77 /1000.;  // Kochav site = Lat=26.1673 Lon=33.0432
    double easting = 249943.68/1000.; // kochav site

    ITMtoAbsWGS84( easting,  northing,  lat,  lon);
    std::cout<<" ("<<northing<<","<<easting<<") translates to Lat="<<lat<<" Lon="<<lon<<std::endl;

    ITMtoAbsWGS84( 0,  0,  lat,  lon);
    std::cout<<" ("<<0<<","<<0<<") translates to lat="<<lat<<" lon="<<lon<<std::endl;


    WGS84toITM( lat,  lon, easting, northing);
    std::cout<<" Bacward:"<<lat<<","<<lon<<"==> "<<easting<<","<<northing<<std::endl;



PJ_CONTEXT *ctx = proj_context_create();
// Create transformation from EPSG:2039 (ITM) to EPSG:4326 (WGS84)
PJ *proj = proj_create_crs_to_crs(ctx, "EPSG:2039", "EPSG:4326", NULL);
PJ_COORD itm_coord = proj_coord(249943.68/100, 722580.77/100, 0, 0);     // E,N in ITM (meters)
PJ_COORD geo_coord = proj_trans(proj, PJ_FWD, itm_coord);        // transform forward to WGS84
double lat_deg = geo_coord.lp.phi * 180.0/M_PI;                  // radians to degrees
double lon_deg = geo_coord.lp.lam * 180.0/M_PI;
std::cout << "Lat: " << lat_deg << "°, Lon: " << lon_deg << "°\n";
proj_destroy(proj);
proj_context_destroy(ctx);

    return 0;
}
