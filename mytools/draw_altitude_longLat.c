
struct coordinate {
    Double_t  longitude;
    Double_t  latitude;
    Double_t  elevation;
};

class coord {       // The class
  public:             // Access specifier
    Double_t longitude;
    Double_t latitude;
    Double_t elevation;

    Double_t Rlati, Rlong;

    coord(Double_t long_, Double_t lat_, Double_t z_) { // Constructor with parameters
     longitude = long_;
     latitude = lat_;
     elevation = z_;

     Rlati = 110.574;
     Rlong = 111.320 * cos(longitude/180.*  TMath::Pi());
    }

    void add_to_graph2d(TGraph2D *g){
        g->AddPoint(  longitude,latitude,elevation);
    }

    void add_to_graph(TGraph *g){
        g->SetPoint(g->GetN(),longitude, latitude);
    }

     TGraph *g_dxdy(TGraph *g){
        int N=g->GetN();
        Double_t *xv=g->GetX();
        Double_t *yv=g->GetY();
        Double_t zFactor = 1;
        TGraph *gdxdy = new TGraph();
        for (int i=0; i<N; i++){
            Double_t xKm = (xv[i] - longitude) * Rlong ;
            Double_t yKm = (yv[i] - latitude) * Rlati ;
            gdxdy->SetPoint(gdxdy->GetN(),xKm,yKm);
        }
        gdxdy->SetTitle(";km from Longtitude  [east][Km]; Km from Latitude  [north][Km];Altitude[m]") ;

        return (gdxdy);

        }
    TGraph2D *g_dxdy(TGraph2D *g){
        int N=g->GetN();
        Double_t *xv=g->GetX();
        Double_t *yv=g->GetY();
        Double_t *zv=g->GetZ();
        Double_t zFactor = 1;
        TGraph2D *gdxdy = new TGraph2D();
        for (int i=0; i<N; i++){
            Double_t xKm = (xv[i] - longitude) * Rlong ;
            Double_t yKm = (yv[i] - latitude) * Rlati ;
            Double_t zKm = zv[i] * zFactor;
            gdxdy->AddPoint(xKm,yKm,zKm);
        }
        gdxdy->SetTitle(";km from Longtitude  [east][Km]; Km from Latitude  [north][Km];Altitude[m]") ;

        return (gdxdy);

        }


};

void append2D(TGraph2D *g1, TGraph2D *g2){
        int N=g2->GetN();
        Double_t *xv=g2->GetX();
        Double_t *yv=g2->GetY();
        Double_t *zv=g2->GetZ();
        for (int i=0; i<N; i++){
            g1->AddPoint(xv[i],yv[i],zv[i]);
        }

}


TGraph2D *get_map(TString file_name) {
    TTree *T10 = new TTree("T10", "T10");
    T10->ReadFile(file_name, "lat2:lon2:elav2:dx:dy:land");
    int N = T10->GetEntries();
    printf ("got %d Entries \n",N);
    T10->Draw("lat2:lon2:elav2", "", "surf1");
    Double_t *xv = T10->GetV2();
    Double_t *yv = T10->GetV1();
    Double_t *zv = T10->GetV3();
    TGraph2D *g2_long_lat = new TGraph2D(N, xv,yv,zv );
    g2_long_lat->SetTitle(";Longtitude  [east]; Latitude  [north];Altitude[m]") ;
    return (g2_long_lat);
}

TGraph2D *get_path3D(){
    coord *tunnel_entrance=new coord(35.539444, 32.59722, -171.43);
    coord *pipe_top=new coord(35.51976, 32.599451, 223);
    coord *pipe_bot=new coord(35.51976, 32.599451, -280);
    coord *site = new coord(35.52881, 32.598255, -280 );
    TGraph2D *g_path2d = new TGraph2D();
    tunnel_entrance->add_to_graph2d(g_path2d);
    site->add_to_graph2d(g_path2d);
    pipe_bot->add_to_graph2d(g_path2d);
    pipe_top->add_to_graph2d(g_path2d);
    g_path2d->SetLineWidth(3);

    return (g_path2d);
}
TGraph *get_path2D(){
    coord *tunnel_entrance=new coord(35.539444, 32.59722, -171.43);
    coord *pipe_top=new coord(35.51976, 32.599451, 223);
    coord *pipe_bot=new coord(35.51976, 32.599451, -280);
    coord *site = new coord(35.52881, 32.598255, -280 );
    TGraph *g_path = new TGraph();
    tunnel_entrance->add_to_graph(g_path);
    site->add_to_graph(g_path);
    pipe_bot->add_to_graph(g_path);
    pipe_top->add_to_graph(g_path);
    g_path->SetLineWidth(3);

    return (g_path);
}


void draw_altitude_longLat(){
    //  .L draw_altitude_longLat.c

  //TGraph2D *g10 = get_map("R10km_0.1.csv");

    TGraph2D *g = get_map("R2km_0.005.csv");
    g->SetNpx(500);
    g->SetNpy(500);
    TH2D *h2 = g->GetHistogram();

    TGraph2D *g5 = get_map("R5km_0.01.csv");
    append2D(g5,g);
    g5->SetNpx(500);
    g5->SetNpy(500);
    TH2D *h5 = g5->GetHistogram();

    TGraph2D *g30 = get_map("R30km_0.1.csv");
    append2D(g30,g5);
    g30->SetNpx(500);
    g30->SetNpy(500);
    TH2D *h30 = g30->GetHistogram();

    TGraph2D *g100 = get_map("R100km_1.csv");
    append2D(g100,g30);
    g100->SetNpx(500);
    g100->SetNpy(500);
    TH2D *h100 = g100->GetHistogram();

    coord *site = new coord(35.52881, 32.598255, -280 );
    TGraph2D *gxy = site->g_dxdy(g);
    TGraph2D *g5xy = site->g_dxdy(g5);
    TGraph2D *g30xy = site->g_dxdy(g30);
    TGraph2D *g100xy = site->g_dxdy(g100);
    gxy->SetNpx(500);
    gxy->SetNpy(500);
    g5xy->SetNpx(500);
    g5xy->SetNpy(500);
    g30xy->SetNpx(500);
    g30xy->SetNpy(500);
    g100xy->SetNpx(500);
    g100xy->SetNpy(500);
    TH2D *hxy = gxy->GetHistogram();
    TH2D *h5xy = g5xy->GetHistogram();
    TH2D *h30xy = g30xy->GetHistogram();
    TH2D *h100xy = g100xy->GetHistogram();

    TGraph2D *path3D = get_path3D();
    TGraph *path2D = get_path2D();
    TGraph2D *path3d_xy = site->g_dxdy(path3D);
    TGraph *path2d_xy = site->g_dxdy(path2D);
    path3D->SetLineWidth(6);
    path2D->SetLineWidth(6);
    path3d_xy->SetLineWidth(6);
    path2d_xy->SetLineWidth(6);

gStyle->SetPalette(kCividis);
   // gStyle->SetPalette(1); //kInvertedDarkBodyRadiation);

    TCanvas *c4 = new TCanvas("c4","c4",2200,500);
    c4->Divide(3,1);
//    c4->cd(1);
 //   gxy->Draw("colz");
  //  path2d_xy->Draw("sameline");
    c4->cd(1);
    h5xy->Draw("colz");
   h5xy->GetZaxis()->SetLabelSize(0.035);
   h5xy->GetZaxis()->SetTitleSize(0.035);
   h5xy->Draw("cont4z");
   gPad->Update();
   gPad->SetRightMargin(0.15);
   TPaletteAxis *palette = (TPaletteAxis*)h5xy->GetListOfFunctions()->FindObject("palette");
   palette->SetX1NDC(0.85);
   palette->SetX2NDC(0.9);
   palette->SetY1NDC(0.2);
   palette->SetY2NDC(0.8);
   gPad->Modified();
  gPad->Update();
    path2d_xy->Draw("sameline");
    c4->cd(2);
    g30xy->Draw("colz");
    path2d_xy->Draw("sameline");
    gPad->Update();
   gPad->SetRightMargin(0.15);
   palette = (TPaletteAxis*)h5xy->GetListOfFunctions()->FindObject("palette");
   palette->SetX1NDC(0.85);
   palette->SetX2NDC(0.9);
   palette->SetY1NDC(0.1);
   palette->SetY2NDC(0.9);
   gPad->Modified();
  gPad->Update();

    c4->cd(3);
    g100xy->Draw("colz");
    path2d_xy->Draw("sameline");


    TCanvas *c5 = new TCanvas("c5","c5",1000,1000);
    c5->Divide(2,2);
    c5->cd(1);
    g->Draw("surf2");
    path3d_xy->Draw("sameline");
    c5->cd(2);
    g5xy->Draw("surf3");
    path3d_xy->Draw("sameline");
    c5->cd(3);
    g30xy->Draw("surf3");
    path3d_xy->Draw("sameline");
    c5->cd(4);
    g100xy->Draw("surf3");
    path3d_xy->Draw("sameline");





  //  g->Draw("TRI2");




}

/*

from numpy import pi, cos, sin, arccos, arange
import mpl_toolkits.mplot3d
import matplotlib.pyplot as pp

num_pts = 1000
indices = arange(0, num_pts, dtype=float) + 0.5

phi = arccos(1 - 2*indices/num_pts)
theta = pi * (1 + 5**0.5) * indices

x, y, z = cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi);

pp.figure().add_subplot(111, projection='3d').scatter(x, y, z);
pp.show()


*/