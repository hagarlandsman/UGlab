#include "ITM_coordinate.h"
#include <proj.h>
#include <iostream>


/*
ITM (0,0) =  Latitude = 26.061987 , Longitude= 33.011727


*/


void OffsetFromReferenceToWGS84(double refLat, double refLon,double deltaE, double deltaN, double& newLat, double& newLon) {
// Step 1: Convert reference WGS84 to ITM
double e0, n0;
WGS84toITM(refLat, refLon, e0, n0);

// Step 2: Apply offset in meters
double e1 = e0 + deltaE;
double n1 = n0 + deltaN;

// Step 3: Convert back to WGS84
ITMtoWGS84(e1, n1, newLat, newLon);
}

void ITMtoAbsWGS84( double deltaE, double deltaN,    double& newLat, double& newLon) {
    double e0, n0;
    double refLat = 26.1673; // ITM
    double refLon = 33.0432;
    WGS84toITM(refLat, refLon, e0, n0);
    std::cout<<"ITMtoAbsWGS84:: When converting Lat="<<refLat<<" and Long="<<refLon<<", we get in IDT e0="<<e0<<" and n0="<<n0<<std::endl;
    // Step 2: Apply offset in meters
    double e1 = e0 + deltaE;
    double n1 = n0 + deltaN;
    // Step 3: Convert back to WGS84
    ITMtoWGS84(e1, n1, newLat, newLon);

}
void ITMtoWGS84(double easting, double northing, double& lat, double& lon) {
    PJ_CONTEXT* ctx = proj_context_create();
    PJ* transform = proj_create_crs_to_crs(ctx, "EPSG:2039", "EPSG:4326", NULL); //https://epsg.io/?q=4326
    if (!transform) {
        std::cerr << "Failed to create transformation from EPSG:2039 to EPSG:4326\n";
        return;
    }

    PJ_COORD input = proj_coord(easting, northing, 0, 0);
    PJ_COORD output = proj_trans(transform, PJ_FWD, input);

    lon = proj_todeg(output.lp.lam);
    lat = proj_todeg(output.lp.phi);
    std::cout<<" ITMtoWGS84: easting="<<easting<<"\tNorthing="<<northing<<"===> lon="<<lon<<"  lat="<<lat<<std::endl;
    proj_destroy(transform);
    proj_context_destroy(ctx);
}

void WGS84toITM(double lat, double lon, double& easting, double& northing) {
    PJ_CONTEXT* ctx = proj_context_create();
    PJ* transform = proj_create_crs_to_crs(ctx, "EPSG:4326", "EPSG:2039", NULL);
    if (!transform) {
        std::cerr << "Failed to create transformation from EPSG:4326 to EPSG:2039\n";
        return;
    }

    PJ_COORD input = proj_coord(proj_torad(lon), proj_torad(lat), 0, 0);
    PJ_COORD output = proj_trans(transform, PJ_FWD, input);

    easting = output.xy.x;
    northing = output.xy.y;
    std::cout<<" WGS84toITM:  lon="<<lon<<"  lat="<<lat<<"===> easting="<<easting<<"\tNorthing="<<northing<<std::endl;

    proj_destroy(transform);
    proj_context_destroy(ctx);
}

void RelativeOffsetFromReference(double lat, double lon,
                                 double refLat, double refLon,
                                 double& deltaEasting, double& deltaNorthing) {
    double e1, n1, e0, n0;
    WGS84toITM(lat, lon, e1, n1);
    WGS84toITM(refLat, refLon, e0, n0);
    deltaEasting = e1 - e0;
    deltaNorthing = n1 - n0;
}
void print_axis_units(const char* crs_string) {
    PJ_CONTEXT* ctx = proj_context_create();
    PJ* crs = proj_create(ctx, crs_string);
    if (!crs) {
        std::cerr << "Failed to load CRS: " << crs_string << "\n";
        return;
    }

    const PJ* cs = proj_crs_get_coordinate_system(ctx, crs);
    if (!cs) {
        std::cerr << "Failed to get coordinate system.\n";
        proj_destroy(crs);
        proj_context_destroy(ctx);
        return;
    }

    for (int i = 0; i < 2; ++i) {
        const char* unit_name = nullptr;
        double conversion_factor = 0.0;

        const char* name = nullptr;
const char* direction = nullptr;
const char* unit_auth_name = nullptr;
const char* unit_code = nullptr;
const char* axis_orientation = nullptr;

if (proj_cs_get_axis_info(
        ctx, cs, i,
        &name, &direction,
        &unit_name, &conversion_factor,
        &unit_auth_name, &unit_code, &axis_orientation)) {
    std::cout << "Axis " << i << " unit: " << unit_name
              << " (conversion to meters: " << conversion_factor << ")\n";
}
    }

    proj_destroy(crs);
    proj_context_destroy(ctx);
}
