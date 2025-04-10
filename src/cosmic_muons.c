#include <TMath.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <iostream>
#include "UGburden.h"
/*

g++ -o muon cosmic_muons.c LinkDef_rdict.cxx UGburden.cpp -Ivcpkg/installed/x64-linux/include -pthread  `root-config --cflags --libs`


*/

TH1D* strawmanMuonSimulationWithAngles(double azimuth0=0, double azimuth_open=0) {
    const int nMuons = 1000;             // Number of simulated muons
    const double minEnergy = 1.0;          // Minimum energy in GeV
    const double maxEnergy = 1000.0;       // Maximum energy in GeV
    const double characteristicDepth = 500.0;  // Decay depth in meters
    TRandom3 randGen(0);                   // Random number generator
    TVector3 start_point(0, 0, -274.505);  // Starting at (1.0, 1.0, 0.5)
    TVector3 direction(0, 0,1);  // Starting at (1.0, 1.0, 0.5)

    // Histograms
    TH1D *hEnergy = new TH1D("hEnergy", "Energy Distribution (GeV)", 100, minEnergy, maxEnergy);
    TH2D *hAngleDist = new TH2D("hAngleDist", "Zenith vs Azimuth Distribution;Azimuth (deg);Zenith (deg)",
                                100, 0, 360, 100, 0, 90);
    TH1D *hDepth = new TH1D("hDepth", "Depth", 500, 0, 5000);
    TH1D *hDepthSurvival = new TH1D("hDepthSurvival", "Depth Survival Probability", 1000, 0, 1);

    TString map_file = "R60km_0.1.csv";
    UGburden burden(map_file);

    for (int i = 0; i < nMuons; i++) {
        if (i%1000 == 0)
            printf ("%d .... \n",i);
        // Generate random energy with E^-2.7 distribution
        double E = minEnergy * pow((maxEnergy / minEnergy), randGen.Rndm());
        double weight = pow(E, -2.7);
        hEnergy->Fill(E, weight);

        // Random zenith and azimuth generation
        double zenith = acos(sqrt(randGen.Rndm())) * 180.0 / TMath::Pi(); // Follows cos^2(theta) distribution
        //double azimuth = randGen.Uniform(0, 360);                         // Uniform azimuth [0, 360] degrees
        //double azimuth0 = 90;
        double azimuth;
        if (azimuth_open > 0)
                azimuth = randGen.Uniform(-azimuth_open, azimuth_open) + azimuth0;                         // Uniform azimuth [0, 360] degrees
        else
                azimuth  = randGen.Uniform(0, 360);
        double sign = randGen.Uniform(-1,1);

        if (sign<0) { azimuth = azimuth + 180; }

        hAngleDist->Fill(azimuth, zenith);
        double depth = burden.getRburden(start_point,zenith* TMath::DegToRad(),azimuth* TMath::DegToRad());
       // double depth = burden.getRburden(start_point,direction);
       // cout<<zenith<<"\t"<<azimuth<<"\t"<<depth<<endl;

        // Propagate muon to random depth and simulate survival
//       double depth = randGen.Uniform(0, 1000); // Depth in meters
        double survivalProb = exp(-depth / (characteristicDepth * cos(zenith * TMath::Pi() / 180.0)));
       // double E_in_MeV;
       // double maxRange = 2298.2 * TMath::Log( 0.001920 * E_in_MeV + 0.998)
       // if (randGen.Rndm() < survivalProb) {
            hDepth->Fill(depth);
            hDepthSurvival->Fill(survivalProb);

       // }
    }

    // Draw results
    TCanvas *c1 = new TCanvas("c1", "Cosmic Muons Simulation with Angles", 1200, 800);
    c1->Divide(2, 2);

    c1->cd(1);
    hEnergy->Draw("hist");

    c1->cd(2);
    hAngleDist->Draw("colz");

    c1->cd(3);
    hDepthSurvival->Draw("hist");

    c1->cd(4);
    hDepth->Draw("hist");

    c1->SaveAs(Form("cosmic_muon_simulation_%.2f_%.2f.pdf",azimuth0,azimuth_open));
    return hDepth;
}


int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);

    TH1D *h_north = strawmanMuonSimulationWithAngles(0,0.05);
    TH1D *h_east = strawmanMuonSimulationWithAngles(90,0.05);
    TH1D *h_all = strawmanMuonSimulationWithAngles();

        app.Run();  // Start the ROOT application (GUI)

    return 1;
    }
