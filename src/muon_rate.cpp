#include <iostream>
#include <cmath>
#include <vector>
#include "UGburden.h"

// Physical constants
const double a = 2.0e-3;  // MeV / (g/cm^2)
const double b = 4.0e-6;  // 1 / (g/cm^2)
const double rho = 2.65;  // g/cm³ (rock)
const double pi = M_PI;

// Convert depth (m.w.e.) to g/cm²
double depth_to_gram_cm2(double mwe) {
    return mwe * 100.0 * rho;
}

// Gaisser's surface muon flux parameterization
double surface_muon_flux(double E, double theta_rad) {
    double cos_theta = std::cos(theta_rad);
    double f = 0.14 * std::pow(E, -2.7) *
           (1.0 / (1.0 + 1.1 * E * cos_theta / 115.0) +
            0.054 / (1.0 + 1.1 * E * cos_theta / 850.0));
    return f;

}

// Minimum energy to survive depth X [g/cm^2]
double energy_min(double X) {
   // printf ("%f vs. %f \n",(a / b) * (std::exp(b * X) - 1.0) / 1000.0, a * X);
    return a*X;
    return (a / b) * (std::exp(b * X) - 1.0) / 1000.0;  // in GeV
}

double get_flux_ug_deep(double x){ // depth in mwe
    const int N_theta = 90;
    const int N_phi = 180;
    const int N_E = 100;


    // Energy bins [1 GeV to 10 TeV]
    std::vector<double> E_vals(N_E);
    for (int i = 0; i < N_E; ++i) {
        E_vals[i] = std::pow(10.0, 0.0 + 4.0 * i / double(N_E - 1));
       // cout<<E_vals[i]<<"\t";
    }
    cout<<endl;
    // Integration
    double total_rate = 0.0;
    double total_rate_surface = 0.0;

    for (int i = 0; i < N_theta; ++i) {
        double theta = (i + 0.5) * (pi / 2.0) / N_theta;
        double dtheta = (pi / 2.0) / N_theta;
        double sin_theta = std::sin(theta);
  //      printf ("theta= %f,  %d/%d \n ",theta,i,N_theta);
        for (int j = 0; j < N_phi; ++j) {
            double dphi = 2.0 * pi / N_phi;
            double phi = dphi * j;
          //  printf ("theta= %f, phi=%f \n",theta,phi);
            TVector3 dir;
            dir.SetMagThetaPhi(1.0, theta, phi);
/*
            double rock_m = burden.getRburden(start_point, dir); // in m

            double rock_mwe = rock_m *2.6;
            double X = depth_to_gram_cm2(rock_mwe);
            double E_threshold = energy_min(X);
            printf ("rock = %f, E_thr=%f \n",rock_mwe,E_threshold);
*/          //  printf ("E_threshold = %f \n",E_threshold);
            // Example overburden model (simplified): deeper at vertical

             double depth_mwe = x / std::cos(theta);
             double X = depth_to_gram_cm2(depth_mwe);
             double E_threshold = energy_min(X);

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
   return total_rate/4/M_PI;
    }

    double get_flux_kochav(){ // depth in mwe
        const int N_theta = 90;
        const int N_phi = 180;
        const int N_E = 100;
        TString map_file = "../mytools/scan_1km_step10m.csv";
        TVector3 start_point(0, 0, -274.505);
        UGburden burden(map_file);


        // Energy bins [1 GeV to 10 TeV]
        std::vector<double> E_vals(N_E);
        for (int i = 0; i < N_E; ++i) {
            E_vals[i] = std::pow(10.0, 0.0 + 4.0 * i / double(N_E - 1));
           // cout<<E_vals[i]<<"\t";
        }
        cout<<endl;
        // Integration
        double total_rate = 0.0;
        double total_rate_surface = 0.0;
        TVector3 dir;

        for (int i = 0; i < N_theta; ++i) {
            double theta = (i + 0.5) * (pi / 2.0) / N_theta;
            double dtheta = (pi / 2.0) / N_theta;
            double sin_theta = std::sin(theta);
           printf ("theta= %f,  %d/%d \n ",theta,i,N_theta);
            for (int j = 0; j < N_phi; ++j) {
                double dphi = 2.0 * pi / N_phi;
                double phi = dphi * j;
              //  printf ("theta= %f, phi=%f \n",theta,phi);
                dir.SetMagThetaPhi(1.0, theta, phi);

                double rock_m = burden.getRburden(start_point, dir); // in m

                double rock_mwe = rock_m *2.6;
                double X = depth_to_gram_cm2(rock_mwe);
                double E_threshold = energy_min(X);
     //           printf ("rock = %f, E_thr=%f \n",rock_mwe,E_threshold);
             //  printf ("E_threshold = %f \n",E_threshold);
                // Example overburden model (simplified): deeper at vertical

          /*       double depth_mwe = x / std::cos(theta);
                 double X = depth_to_gram_cm2(depth_mwe);
                 double E_threshold = energy_min(X);
    */
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
        std::cout << "Estimated underground muon rate: "
                  << total_rate << "ug  cm^-2 s^-1" << std::endl
                  << total_rate_surface << "surface  cm^-2 s^-1" << std::endl
                 <<"ratio = "<<total_rate / total_rate_surface<<std::endl;
                printf ("ug flux = %e \t surface = %e \n",total_rate*100*100,total_rate_surface*100*100);

       return total_rate/4/M_PI;
        }



int main() {

    TGraph *g1=new TGraph();
    for (float x=700;x<900; x=x+2){
    double r= get_flux_ug_deep(x);
        g1->SetPoint(g1->GetN(),x,r);
        printf ("%f %e \n",x,r);
    }
    TCanvas *c1=new TCanvas("c1","c1",800,800);
   // 7.640681e-08
    g1->Draw("A*l");
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetGrid();
    c1->SaveAs("d2.png");
    TFile* fout = new TFile("coszenith_phi_hist.root", "RECREATE");
    g1->Write();
    fout->Close();

    double kk = get_flux_kochav();
    printf ("kk= %e \n",kk);
    return -1;
    // depth in mwe



    const int N_theta = 90;
    const int N_phi = 180;
    const int N_E = 100;
    TString map_file = "../mytools/scan_1km_step10m.csv";
    TVector3 start_point(0, 0, -274.505);
    UGburden burden(map_file);

    // Energy bins [1 GeV to 10 TeV]
    std::vector<double> E_vals(N_E);
    for (int i = 0; i < N_E; ++i) {
        E_vals[i] = std::pow(10.0, 0.0 + 4.0 * i / double(N_E - 1));
        cout<<E_vals[i]<<"\t";
    }
    cout<<endl;
    // Integration
    double total_rate = 0.0;
    double total_rate_surface = 0.0;

    for (int i = 0; i < N_theta; ++i) {
        double theta = (i + 0.5) * (pi / 2.0) / N_theta;
        double dtheta = (pi / 2.0) / N_theta;
        double sin_theta = std::sin(theta);
      //  printf ("theta= %f,  %d/%d \n ",theta,i,N_theta);
        for (int j = 0; j < N_phi; ++j) {
            double dphi = 2.0 * pi / N_phi;
            double phi = dphi * j;
          //  printf ("theta= %f, phi=%f \n",theta,phi);
            TVector3 dir;
            dir.SetMagThetaPhi(1.0, theta, phi);
/*
            double rock_m = burden.getRburden(start_point, dir); // in m

            double rock_mwe = rock_m *2.6;
            double X = depth_to_gram_cm2(rock_mwe);
            double E_threshold = energy_min(X);
            printf ("rock = %f, E_thr=%f \n",rock_mwe,E_threshold);
*/          //  printf ("E_threshold = %f \n",E_threshold);
            // Example overburden model (simplified): deeper at vertical

             double depth_mwe = 2000 / std::cos(theta);
             double X = depth_to_gram_cm2(depth_mwe);
             double E_threshold = energy_min(X);

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

    std::cout << "Estimated underground muon rate: "
              << total_rate << "ug  cm^-2 s^-1" << std::endl
              << total_rate_surface << "surface  cm^-2 s^-1" << std::endl
             <<"ratio = "<<total_rate / total_rate_surface<<std::endl;
            printf ("ug flux = %e \t surface = %e \n",total_rate*100*100,total_rate_surface*100*100);
    return 0;
}
