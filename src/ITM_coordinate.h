#ifndef ITM_COORDINATE_H
#define ITM_COORDINATE_H



// Convert ITM to absolute long, lat
void ITMtoAbsWGS84( double deltaE, double deltaN,    double& newLat, double& newLon);

// Helper to caluclate offset from some reference point in northing, easting
void OffsetFromReferenceToWGS84(double refLat, double refLon,
    double deltaE, double deltaN,
    double& newLat, double& newLon) ;


// Convert ITM to WGS84
void ITMtoWGS84(double easting, double northing, double& lat, double& lon);

// Convert WGS84 to ITM
void WGS84toITM(double lat, double lon, double& easting, double& northing);

// Compute offset from reference WGS84 point (in ITM meters)
void RelativeOffsetFromReference(double lat, double lon,
                                 double refLat, double refLon,
                                 double& deltaEasting, double& deltaNorthing);

void print_axis_units(const char* crs_string);
#endif
