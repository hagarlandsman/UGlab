#include <iostream>
#include <cmath>
#include <vector>
#include "UGburden.h"
#include "muon_rate.h"

// Physical constants
// const double a = 2.0;  // MeV / (g/cm^2)
const double a = 2.0e-3;  // GeV / (g/cm^2)
const double b = 4.0e-6;  // 1 / (g/cm^2)
const double rho_rock = 2.65;  // g/cm³ (rock)
const double rho_water = 1.00;  // g/cm³ (water)
const double pi = M_PI;

double depth_to_mwe_rock(double r){
    /*
        r - depth in some units
        output: depth in similiar units
    */
    return depth_to_mwe(r,rho_rock);
}

double depth_to_mwe(double r,double rho){
    /*
        r  - depth in some units
        rho - density in g/cm^3
        return - m.w.e depth in similiar units to input
    */
    return (r * rho / rho_water);
 }

double depth_to_gram_cm2(double r, double rho) {
    /*
     r  - in [m]
     100  - in [cm/m]
     rho  - in [g/cm^3]
     out = r * 100 * rho = [g/cm^2]

    */
    return r * 100.0 * rho;

}

double depth_to_gram_cm2_rock(double r){
    /*
     r  - in [m]
     out = r in [g/cm^2] for rock.
    */
    return (depth_to_gram_cm2(r,2.65));
 }

double depth_to_gram_cm2_water(double r){
       /*
     r  - in [m]
     out = r in [g/cm^2] for water.
    */
    return (depth_to_gram_cm2(r,1.00));
 }

// Gaisser's surface muon flux parameterization
double surface_muon_flux(double E, double theta_rad) {
    /*
    See equation 30.4, in PDG Cosmic rays review:
     theta_rad - in [radians]
     E - in [GeV]
     Output - dN/dE/dOmega = 1/m^2/s/Sr

    */
    double cos_theta = std::cos(theta_rad);
    double f = 0.14 * std::pow(E, -2.7) *
           (1.0 / (1.0 + 1.1 * E * cos_theta / 115.0) +
            0.054 / (1.0 + 1.1 * E * cos_theta / 850.0));
    return f;

}

// Minimum energy to survive depth X [g/cm^2]
double energy_min(double X) {
    /*
         See equation 30.5,30.6 in PDG Cosmic rays review:
         a/b ~= 500 GeV in standard rock
         X = Range  = g/cm^2
         a = GeV * cm^2 / g
         b = cm^2 / g
         output - in GeV

    */
   // printf ("%f vs. %f \n",(a / b) * (std::exp(b * X) - 1.0) / 1000.0, a * X);
    return a*X;  // This approximation is true when X<<1/b~2.5 km.w.e
    return (a / b) * (std::exp(b * X) - 1.0);  // in GeV
}
double get_flux_ug(double x) {
    TVector3 start_point(0, 0, -274.505);
    TString map_file = "../mytools/scan_1km_step10m.csv";
    return (get_flux_ug(x,map_file,start_point));
}



double get_flux_ug(double x, TString map_file, TVector3 start_point){ // depth in m real (not mwe). If x is negative, do the ray tracing depth from a google api file.
    //TVector3 start_point(0, 0, -274.505);
//    TString map_file = "../mytools/scan_1km_step10m.csv";

    UGburden burden;
    const int N_theta = 90;
    const int N_phi = 120;
    const int N_E = 100;
    if (x<0){
        burden.load(map_file);
    }

    // Energy bins [1 GeV to 10 TeV]
    std::vector<double> E_vals(N_E);
    for (int i = 0; i < N_E; ++i) {
        E_vals[i] = std::pow(10.0, 0.0 + 6.0 * i / double(N_E - 1));
       // cout<<E_vals[i]<<"\t";
    }
    cout<<endl;
    // Integration
    double total_rate = 0.0;
    double total_rate_surface = 0.0;
    bool print_vertical = false;
    for (int i = 0; i < N_theta; ++i) {   // From 0 to pi/2
        if (x<0) printf ("==>\t %d // %d \n",i,N_theta);
        if (i>0 && print_vertical == false ){
            print_vertical = true;
        //    printf ("vertical = %e (%e) \n", total_rate, total_rate_surface);

        }
        double theta = (i + 0.5) * (pi / 2.0) / N_theta;
        double dtheta = (pi / 2.0) / N_theta;
        double sin_theta = std::sin(theta);
        for (int j = 0; j < N_phi; ++j) {   // From 0 to 2pi
            double dphi = 2.0 * pi / N_phi;
            double phi = dphi * j;
            double rock_m = 0;
            double depth = 0;
            if (x > 0) {  // Use constant depth model
                depth = x / std::cos(theta);  // in m, for rock
            }
            else {
                TVector3 dir;
                dir.SetMagThetaPhi(1.0, theta, phi);
                depth = burden.getRburden(start_point, dir); // in m
                          }
            double depth_mwe = depth_to_mwe_rock(depth); // in m for water
  //          printf ("depth = %f theta=%f, l=%f %f \n",x,theta,depth,depth_mwe);
            double X = depth_to_gram_cm2_rock(depth); // in g/cm^2
            double E_threshold = energy_min(X); // in GeV
          //  printf ("X=%f, depth=%f, depth_mwe=%f, thereshold=%e \n",X,depth,depth_to_mwe_rock(depth)/1000.,E_threshold );

            for (int k = 0; k < N_E; ++k) {
                double E = E_vals[k];
                double dE = (k < N_E - 1) ? (E_vals[k + 1] - E_vals[k]) : (E_vals[k] - E_vals[k - 1]);
                double flux = surface_muon_flux(E, theta);
                double dr = flux * dE * sin_theta * dtheta * dphi;
               // printf ("flux= %f\t dE=%f\t dt=%f\t dphi=%f \tdr=%f  ",flux,dE, dtheta, dphi,dr );

                if (E > E_threshold) {
                    total_rate += flux * dE * sin_theta * dtheta * dphi;
                }
                total_rate_surface += flux * dE * sin_theta * dtheta * dphi;
              //  printf ("E=%f, E_thr=%f ,dr=%e,  tot=%e, surface=%e\n",E,E_threshold,flux * dE * sin_theta * dtheta * dphi,total_rate,total_rate_surface);

            }
        }
    }
   /* std::cout << "Estimated underground muon rate: "
              << total_rate << "ug  cm^-2 s^-1" << std::endl
              << total_rate_surface << "surface  cm^-2 s^-1" << std::endl
             <<"ratio = "<<total_rate / total_rate_surface<<std::endl;
            printf ("ug flux = %e \t surface = %e \n",total_rate*100*100,total_rate_surface*100*100);
    */
 //  printf ("x= %f\t xmwe= %f\t  g/cm^2=%e\t f=%e \t surface=%e\n",x,depth_to_mwe_rock(x),depth_to_gram_cm2_rock(x),total_rate/4/M_PI,total_rate_surface/4/M_PI);
   return total_rate;
}



int main_like() {

    TVector3 start_point(0, 0, -274.505);
    TString map_file = "../mytools/scan_1km_step10m.csv";

    double h = get_flux_ug(-4,map_file,start_point);  // Some negative value as depth, so the function will use the 3d surface map for ray tracing.
    printf ("Kochav = %f\n",h);

    TGraph *g1=new TGraph();
    TGraph *g1i=new TGraph();

    for (float x=40;x<1000; x=x+10){
    double r= get_flux_ug(x);
        g1->SetPoint(g1->GetN(),depth_to_mwe_rock(x),r);
        g1i->SetPoint(g1i->GetN(),r,depth_to_mwe_rock(x));

        printf ("x= %f m \t xmwe= %f m\t  f=%e \n",x,depth_to_mwe_rock(x),r);
    }
    printf ("Kochav depth like %f \n",g1i->Eval(h));
    TGraph *g2 = new TGraph();
    g2->SetPoint(0,g1i->Eval(h),h);
    g2->SetMarkerColor(2);
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(3);

    TCanvas *c1=new TCanvas("c1","c1",800,800);
   // 7.640681e-08
    g1->SetTitle("Flat depth; Depth [m.w.e]; m^{-2}s^{-1}");
    g2->SetTitle("Kochav; Depth [m.w.e]; m^{-2}s^{-1}");
    g1->Draw("A*l");
    g2->Draw("*");
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetGrid();
    c1->SaveAs("d2.png");
    TFile* fout = new TFile("coszenith_phi_hist.root", "RECREATE");
    g1->Write();
    g2->Write();
    fout->Close();

//double kk = get_flux_kochav();
  //  printf ("kk= %e \n",kk);
    return -1;
    // depth in mwe



}
