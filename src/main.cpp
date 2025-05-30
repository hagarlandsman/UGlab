#include "UGburden.h"
#include "TVector3.h"
#include "TH2F.h"
#include "TStyle.h"
#include "muon_rate.h"

// Plots for the 04/2025 cosmic rays paper
void make_1d_slunt_plot(UGburden *burden) ;
void make_2d_slunt_plot(UGburden *burden) ;
int get_flux();
int main() {
   // TString map_file = "../../../works/UGlab/kochav_Hayarden/sim/R60km_0.1.csv";
    //TString map_file = "../../../works/UGlab/kochav_Hayarden/sim/R10km_0.005.csv";



    printf ("\n\n=====> \t\t Loading map \n\n");
   // TString map_file = "../mytools/scan_10km_step10m.csv";

    TString map_file = "../mytools/scan_10km_step10m.csv";
//    UGburden burden(map_file);
    UGburden *burden=new UGburden(map_file);

    burden->draw_surface();

    //some_examples(burden);
    printf ("\n\n=====> \t\t 1D slunt plot\n\n");
    make_1d_slunt_plot(burden);
    printf ("\n\n=====> \t\t 2D slunt plot\n\n");
    make_2d_slunt_plot(burden);

    printf ("\n\n=====> \t\t flux plot\n\n");
    get_flux();
    return (0);
}


int get_flux(){

        TVector3 start_point(0, 0, -274.505);
        TString map_file = "../mytools/scan_10km_step10m.csv";

        double h = get_flux_ug(-4,map_file,start_point);  // Some negative value as depth, so the function will use the 3d surface map for ray tracing.
        printf ("Kochav = %f\n",h);

        TGraph *g1=new TGraph();
        TGraph *g1i=new TGraph();

        for (float x=40;x<2000; x=x+20){
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
        c1->SaveAs("get_flux.png");
        TFile* fout = new TFile("get_flux.root", "RECREATE");
        g1->Write();
        g2->Write();
        fout->Close();

    //double kk = get_flux_kochav();
      //  printf ("kk= %e \n",kk);
        return -1;
        // depth in mwe






}

void some_examples(UGburden *burden) {
    burden->get_one(-929.994,939.698);//10.2376, 5.02477);
    burden->get_one(500,500);  // The units are : northing - refNorth, easting - refEast in m
    burden->get_one(1,1);  // The units are : northing - refNorth, easting - refEast
    burden->get_one(0,0);  // The units are : northing - refNorth, easting - refEast
    TVector3 start_point(0, 0, -274.505);  // Starting at (1.0, 1.0, 0.5)


    TVector3 up(0, 0, 10);    // Propagate along z-axis (z direction)
    printf ("up: theta=%f, phi=%f \n ", up.Theta()* 180.0 / M_PI,up.Phi() *180.0/ M_PI);
    TVector3 side(1, 0, 0);    // Propagate along z-axis (z direction)
    printf ("side: theta=%f, phi=%f \n", side.Theta()* 180.0 / M_PI,side.Phi() *180.0/ M_PI);

    TVector3 vertical(0,0,10); //
    double b_vertical = burden->getRburden(start_point, vertical);
    printf ( " Vertical up: burden=%f ,  Angle north = %f , zenith = %f  \n",b_vertical,vertical.Theta()* 180.0 / M_PI,vertical.Phi() *180.0/ M_PI);

    TVector3 north(1, 0, 1);    // Propagate along z-axis (z direction)
    north.SetMagThetaPhi(1.0, 10/180.*M_PI, 0.0);
    double b_north = burden->getRburden(start_point, north);
    printf ( " north: burden=%f ,  Angle north = %f , zenith = %f  \n",b_vertical,vertical.Theta()* 180.0 / M_PI,vertical.Phi() *180.0/ M_PI);

    TVector3 east(0, 1, 1);    // Propagate along z-axis (z direction)
    east.SetMagThetaPhi(1.0, 10/180.*M_PI,  0.5 * M_PI);
    double b_east = burden->getRburden(start_point, east);
    printf ( " east: burden=%f ,  Angle north = %f , zenith = %f  \n",b_east,east.Theta()* 180.0 / M_PI,east.Phi() *180.0/ M_PI);


}
void make_1d_slunt_plot(UGburden *burden) {
    TVector3 start_point(0, 0, -274.505);  // Starting at (1.0, 1.0, 0.5)

    TVector3 east(0,1,20);
    TGraph *g1 = new TGraph();
    TVector3 north(1, 0, 20);    // Propagate along z-axis (z direction)

    for (float theta=-60; theta<80; theta=theta+2){

        north.SetMagThetaPhi(1.0, theta*M_PI/180., 0.0);  // φ = 0 = North direction in XZ plane

        double a0 = burden->getRburden(start_point, north);
        if (theta==30 || theta==-30)  printf (" %f north (%f,%f,%f), theta= %f, phi=%f \n",a0,north.X(),north.Y(),north.Z(),north.Theta(),north.Phi());

        g1->SetPoint(g1->GetN(),theta,a0);
    }
    printf (" East \n");
    TGraph *g2 = new TGraph();

    for (float theta=-60; theta<80; theta=theta+2){
        double phi = M_PI / 2.0;  // 90°, horizontal plane
        east.SetMagThetaPhi(1.0, theta*M_PI/180.,phi);  // φ = 0 = North direction in XZ plane




        double a0 = burden->getRburden(start_point, east);
        if (theta==30 || theta==-30) printf ("%f east (%f,%f,%f), theta= %f, phi=%f \n",a0,east.X(),east.Y(),east.Z(),east.Theta(),east.Phi());
        g2->SetPoint(g2->GetN(),theta,a0);

    }
        TCanvas *c1=new TCanvas("c1","c1",800,800);
        g1->Draw("A*l");
        g2->SetLineColor(2);
        g2->SetMarkerColor(2);

        g2->Draw("*l");
        c1->SaveAs("slunt.png");




}
void make_2d_slunt_plot(UGburden *burden) {
    TVector3 start_point(0, 0, -274.505);  // Starting at (1.0, 1.0, 0.5)

    const int NcosTheta = 100;  // bins in cos(θ), from 0 to 1
    const int Nphi = 180;       // bins in φ (azimuth), 0 to 360°
    TCanvas *c1=new TCanvas("c1","c1",800,800);

    TH2F* h = new TH2F("h_coszen_phi", "cos(zenith) vs Azimuth",
                       Nphi, 0, 360,            // φ in degrees
                       NcosTheta, 0.2, 1.0);      // cos(θ) from 0 (horizon) to 1 (zenith)
    for (int i = 0; i <= NcosTheta; ++i) {
       // double cos_theta = static_cast<double>(i) / NcosTheta;
        double cos_theta = h->GetYaxis()->GetBinCenter(i);
        printf (" %d / %d costheta= %f \n ",i,NcosTheta, cos_theta);
        double theta_rad = std::acos(cos_theta);  // convert to θ
        double theta_deg = theta_rad * 180. / M_PI;
        if (theta_deg>70 or theta_deg<-70 ) continue;
        // double theta_deg = theta_rad * 180.0 / M_PI;

        for (int j = 0; j < Nphi; ++j) {
//            double phi_deg = 360.0 * j / Nphi;
            double phi_deg = h->GetXaxis()->GetBinCenter(j);

            double phi_rad = phi_deg * M_PI / 180.0;

            TVector3 dir;
            dir.SetMagThetaPhi(1.0, theta_rad, phi_rad);
            double a1 = burden->getRburden(start_point, dir);
           // printf(" %f, %f %f\n",phi_deg,theta_deg,a1);
            h->Fill(phi_deg, cos_theta,a1);  // y-axis is cos(zenith)
            }
    }
    gStyle->SetOptStat(0);
    h->GetXaxis()->SetTitle("Azimuth #phi [deg]");
    h->GetYaxis()->SetTitle("cos(#theta)");
    h->GetZaxis()->SetRangeUser(300,1400);
    h->GetZaxis()->SetTitle("Overburden [m]");
    h->SetTitle("Slant Depth");
    gStyle->SetPalette(kViridis);  // or any other from above
    gStyle->SetNumberContours(100);  // smooth gradient

    h->Draw("colz");
    c1->SaveAs("slunt_2d.png");
    c1->SaveAs("slunt_2d.pdf");
    TFile* fout = new TFile("slunt_2d.root", "RECREATE");
    h->Write();
    fout->Close();

    return ;
}

//
