#include "UGburden.h"
#include "TVector3.h"

int main() {
    TString map_file = "R60km_0.1.csv";
    UGburden burden(map_file);
    burden.draw_surface();
    burden.get_one(10,10);
        TVector3 start_point(0, 0, -274.505);  // Starting at (1.0, 1.0, 0.5)
    TVector3 direction(0.0, 10.0, 1.0);    // Propagate along z-axis (z direction)

    double a = burden.getRburden(start_point, direction);
    double b = burden.getRburden(start_point, 1.2,0.3);

    cout<<b<<endl;
    return 0;
}
