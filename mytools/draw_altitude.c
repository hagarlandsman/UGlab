{

    Double_t pi = TMath::Pi();
    Float_t LatKochav = 32.598255; // 33.194;
    Float_t LonKochav = 35.52881; // 35.547;      722700 / 249900     // https://zvikabenhaim.appspot.com/software/ITM/
    Float_t Depth1 = -280;
    //  32.599451, 35.51976
    Float_t Lat_top = 32.599451;
    Float_t Lon_top = 35.51976 ;  //    722830 / 249050  // https://zvikabenhaim.appspot.com/software/ITM/
    Float_t H_top = 223; //222.57;   // from google api

    // 32.59733, 35.538392
    // 35 32' 22'' E 32 35'50''N  32.59722222  35.53944444
    //Float_t Lat_entrance = 32.59733;
    //Float_t Lon_entrance = 35.538392;   //    722600 / 250800

    Float_t Lat_entrance = 32.59722;
    Float_t Lon_entrance = 35.539444;   //    722600 / 250800
    Float_t H_entrance = -171.43; // -150.54; // from google api

    Float_t Rlati = 110.574;
    Float_t Rlong = 111.320 * cos(LonKochav/180.*pi); // .837;
    Float_t zFactor = 1.; //e-3;

    TGraph2D *gtunnel_ll=new TGraph2D();
    gtunnel_ll->SetPoint(0,Lat_top,Lon_top,H_top);
    gtunnel_ll->SetPoint(1,Lat_top,Lon_top,Depth1);
    gtunnel_ll->SetPoint(2,Lat_entrance,Lon_entrance,H_entrance);


    TGraph2D *gStation=new TGraph2D();
    gStation->SetPoint(0,(LatKochav-LatKochav)*Rlati,(LonKochav-LonKochav)*Rlong , Depth1);
    gStation->SetMarkerStyle(34);

    TGraph2D *gPeer=new TGraph2D();
    gPeer->SetPoint(0,(Lat_top-LatKochav)*Rlati,(Lon_top-LonKochav)*Rlong , H_top);
    gPeer->SetPoint(1,(Lat_top-LatKochav)*Rlati,(Lon_top-LonKochav)*Rlong , Depth1);
    gPeer->SetPoint(2,(LatKochav-LatKochav)*Rlati,(LonKochav-LonKochav)*Rlong , Depth1);
    gPeer->SetPoint(3,(Lat_entrance-LatKochav)*Rlati,(Lon_entrance-LonKochav)*Rlong , H_entrance);
    gPeer->SetLineWidth(3);
    gPeer->SetLineColor(2);
    TGraph2D *gPeer1=new TGraph2D();
    gPeer1->SetPoint(0,(Lat_top-LatKochav)*Rlati,(Lon_top-LonKochav)*Rlong , H_top);
    gPeer1->SetPoint(1,(Lat_top-LatKochav)*Rlati,(Lon_top-LonKochav)*Rlong, Depth1);
    TGraph2D *gPeer2=new TGraph2D();
     gPeer2->SetPoint(0,(Lat_top-LatKochav)*Rlati,(Lon_top-LonKochav)*Rlong , Depth1);
     gPeer2->SetPoint(1,(LatKochav-LatKochav)*Rlati,(LonKochav-LonKochav)*Rlong , Depth1);
    gPeer2->SetPoint(2,(Lat_entrance-LatKochav)*Rlati,(Lon_entrance-LonKochav)*Rlong , H_entrance);
    gPeer1->SetLineWidth(2);
    gPeer1->SetLineColor(2);
    gPeer2->SetLineWidth(2);
    gPeer2->SetLineColor(1);

    TGraph *gHere=new TGraph();

    gHere->SetPoint(0,(Lat_top-LatKochav)*Rlati,(Lon_top-LonKochav)*Rlong);
    gHere->SetPoint(1,(Lat_entrance-LatKochav)*Rlati,(Lon_entrance-LonKochav)*Rlong );
    gHere->SetLineWidth(3);
    gHere->SetLineColor(2);

    TGraph *gEntrance=new TGraph();
    gEntrance->SetPoint(0,(Lat_entrance-LatKochav)*Rlati,(Lon_entrance-LonKochav)*Rlong );
    gEntrance->SetMarkerStyle(21);
    gEntrance->SetMarkerSize(1.5);
    TGraph *gUp=new TGraph();
    gUp->SetPoint(0,(Lat_top-LatKochav)*Rlati,(Lon_top-LonKochav)*Rlong);
    gUp->SetMarkerStyle(23);
    gUp->SetMarkerSize(1.5);
    Double_t *xv50;
    Double_t *yv50;
    Double_t *zv50;
    // 1 degree of latitude is 111km, and one degree of longitude is 111km * cos(latitude):
    TTree *T10 = new TTree("T10", "T10");
//    T10->ReadFile("kochav_0006_a.txt", "lat:longi:elav");

//    T10->ReadFile("R10km_0.1.csv", "lat2:lon2:elav2:dx:dy:land"); //:land");
 //  T10->ReadFile("R5km_0.01.csv", "lat2:lon2:elav2:dx:dy:land"); //:land");
    //  T10->ReadFile("R20km_0.01.csv", "lat2:lon2:elav2:dx:dy:land"); //:land");
//    T10->ReadFile("kochav.txt", "type/C:lat/F:longi/F:elav/F:dist/F");
    T10->ReadFile("R5km_0.01.csv", "lat2:lon2:elav2:dx:dy:land");
   // T10->ReadFile("R2km_0.005.csv", "lat2:lon2:elav2:dx:dy:land");

  //  Double_t lat1 = 33.162712565  32.60199361926758   # lab_lat power room
  //  Double_t lon1 = 35.525451159   35.53016839929137   # lab_lon power room

    int N = T10->GetEntries();
    printf ("got %d Entries \n",N);
    T10->Draw("lat2:lon2:elav2", "", "surf1");

    //int N = T10->Draw("(32.60199361926758-lat2)*111:(35.530168-lon2)*111*cos(32.6/180*pi):elav2", "", "surf1");
   // xv50 = T10->GetV2();
   // yv50 = T10->GetV1();
    xv50 = T10->GetV2();
    yv50 = T10->GetV1();
    zv50 = T10->GetV3();
    TGraph2D *g2_50 = new TGraph2D();
    TGraph2D *g2_50_long_lati = new TGraph2D();

    g2_50->SetTitle(";Longtitude to km [east][km]; Latitude to km [north][km];Altitude[km]") ;
    Float_t Clati = LatKochav * Rlati;
    Float_t Clong = LonKochav * Rlong;
    TCanvas *c0=new TCanvas("c0","c0",800,800);
    g2_50_long_lati->Draw("colz");
//    gtunnel_ll->Draw("lsame");
    for (int i = 0; i < N-10; i++)
    {
//        Float_t xKm = xv50[i] * Rlong - Clong;
 //       Float_t yKm = yv50[i] * Rlati - Clati;
        Float_t xKm = xv50[i] * Rlong - Clong;
        Float_t yKm = yv50[i] * Rlati - Clati;

        Float_t zKm = zv50[i] * zFactor;
	printf ("%d/%d  %f %f %f \n",i,N,xKm,yKm,zKm);
        g2_50->SetPoint(g2_50->GetN(), xKm, yKm, zKm);
        g2_50_long_lati->SetPoint(g2_50_long_lati->GetN(),xv50[i],yv50[i],zv50[i]);

    }
    Float_t xKm = Lon_entrance * Rlong - Clong;
    Float_t yKm = Lat_entrance * Rlati - Clati;
    Float_t zKm = H_entrance * zFactor;
        g2_50->SetPoint(g2_50->GetN(), xKm, yKm, zKm);
        g2_50_long_lati->SetPoint(g2_50_long_lati->GetN(),Lon_entrance,Lat_entrance,H_entrance);

    printf ("g2_50=%f",g2_50->GetN());

        TCanvas *c1=new TCanvas("c1","c1",800,800);

    g2_50->Draw("surf1");
    gPeer1->Draw("sameline");
    gPeer2->Draw("sameline");

    gStation->SetMarkerColor(2);
    gStation->SetMarkerSize(2);
    gStation->Draw("samep");

    TCanvas *c2=new TCanvas("c2","c2",800,800);
        g2_50->Draw("colz");
	gEntrance->Draw("p");
	gUp->Draw("p");
	gHere->Draw("l");

}
