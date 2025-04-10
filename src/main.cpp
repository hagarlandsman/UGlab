#include "UGburden.h"
#include "TVector3.h"

int main() {
   // TString map_file = "../../../works/UGlab/kochav_Hayarden/sim/R60km_0.1.csv";
    TString map_file = "../../../works/UGlab/kochav_Hayarden/sim/R10km_0.005.csv";
    UGburden burden(map_file);
    burden.draw_surface();


    burden.get_one(10,10);  // The units are : northing - refNorth, easting - refEast
    TVector3 start_point(0, 0, -274.505);  // Starting at (1.0, 1.0, 0.5)
    double dx,dy,dz;
    dx = 0;
    dy = 1;
    dz = 10;
    TVector3 direction(dx, dy, dz);    // Propagate along z-axis (z direction)
// (1,0,1) = north
// (0,1,0) = east


    TVector3 up(0, 0, 10);    // Propagate along z-axis (z direction)
    printf ("up: theta=%f, phi=%f \n ", up.Theta()* 180.0 / M_PI,up.Phi() *180.0/ M_PI);
    TVector3 side(1, 0, 0);    // Propagate along z-axis (z direction)
    printf ("side: theta=%f, phi=%f \n", side.Theta()* 180.0 / M_PI,side.Phi() *180.0/ M_PI);

    TVector3 direction0(1, 0, 20);    // Propagate along z-axis (z direction)

    TVector3 direction1(0, 1, 2);    // Propagate along z-axis (z direction)
    TVector3 direction2(-1, 0, 2);    // Propagate along z-axis (z direction)
    TVector3 direction3(0, -1, 2);    // Propagate along z-axis (z direction)

    printf ("===\n");
    double a0 = burden.getRburden(start_point, direction0);
    double a1 = burden.getRburden(start_point, direction1);
    double a2 = burden.getRburden(start_point, direction2);
    double a3 = burden.getRburden(start_point, direction3);

   // double b = burden.getRburden(start_point, 1.2,0.3);
   printf ("\n L=%f, Angle north = %f , zenith = %f  \n",a0,direction0.Theta()* 180.0 / M_PI,direction0.Phi() *180.0/ M_PI);

    printf ("\n L=%f, Angle north = %f , zenith = %f  \n",a1,direction1.Theta()* 180.0 / M_PI,direction1.Phi() *180.0/ M_PI);
    printf ("\n L=%f, Angle north = %f , zenith = %f  \n",a2,direction2.Theta()* 180.0 / M_PI,direction2.Phi() *180.0/ M_PI);
    printf ("\n L=%f, Angle north = %f , zenith = %f  \n",a3,direction3.Theta()* 180.0 / M_PI,direction3.Phi() *180.0/ M_PI);

    return 0;
}
