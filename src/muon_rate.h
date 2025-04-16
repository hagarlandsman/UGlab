#ifndef MUONRATE_H
#define MUONRATE_H

double depth_to_mwe_rock(double r);
double depth_to_mwe(double r,double rho);
double depth_to_gram_cm2(double r, double rho) ;
double depth_to_gram_cm2_rock(double r);
double depth_to_gram_cm2_water(double r);
double get_flux_ug(double x);
double get_flux_ug(double x, TString map_file, TVector3 start_point);
double energy_min(double X) ;
double surface_muon_flux(double E, double theta_rad) ;

#endif // MUONRATE_H