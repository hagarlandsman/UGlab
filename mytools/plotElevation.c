

void makeData();
// Reads the txt files holding the longtiude/latitude/elevation data, puts them in 2Dhistograms and in TGraph2Ds and saves it to  output file called ManaraData.root .
// Offsets are being done here on X,Y.

void draw() ;
// Reads the ManaraData.root and generates plots.

TH1F *rays(Int_t N);
// Wrapper to call rays with the default location   0,0,7./1000

TH1F *rays(Int_t N,Double_t x0, Double_t y0, Double_t z0, Int_t color) ;
// Generates N rays distriubuted on the upper sphere. and checks the paths lengths to (x0,y0,z0).
// Saves the output to a .root file and also draws on screen

void go();
// Wrapper that calls "rays" with various geometry and draws a summary.

void drawFile(); // Reads the root files created by "rays", The list and number of files to be read are hardcoded in the function.
void drawSummary(); // Reads the root files created by "rays", and generates a summary plot of all "geometries" on a single figure.

void makeData()
{

  /*
    Latitude: 1 deg = 110.574 km
    Longitude: 1 deg = 111.320*cos(latitude) =  km
  */

  // /pm 10 km with 1000  units

  int N;
  Double_t *xv1;
  Double_t *yv1;
  Double_t *zv1;
  Double_t *xv2;
  Double_t *yv2;
  Double_t *zv2;
  Double_t *xv10;
  Double_t *yv10;
  Double_t *zv10;
  Double_t *xv50;
  Double_t *yv50;
  Double_t *zv50;

  Float_t Rlati=110.574;
  Float_t Rlong=111.320*0.837;

  Float_t  LatManara=33.1808 ; //33.194;
  Float_t LonManara=35.5458; //35.547;
  Float_t Clati=LatManara*Rlati;
  Float_t Clong=LonManara*Rlong;
  Float_t zFactor=1.e-3;

  TH2F *h2_50=new TH2F("h2_50","Manara;km;km;km",1000-1,-50,50,1000-1,-50,50);
  TH2F *h2_10=new TH2F("h2_10","Manara;km;km;km",1000-1,-10,10,1000-1,-10,10);
  TH2F *h2_2=new TH2F("h2_2","Manara;km;km;km",500-1,-2,2,500-1,-2,2);
  TH2F *h2_1=new TH2F("h2_1","Manara;km;km;km",100-1,-2,2,100-1,-2,2);
  TGraph2D *g2=new TGraph2D();
  TGraph2D *g2_50=new TGraph2D();
  TGraph2D *g2_10=new TGraph2D();
  TGraph2D *g2_2=new TGraph2D();
  TGraph2D *g2_1=new TGraph2D();
  g2->SetTitle(";Longtitude to km [east][km]; Latitude to km [north][km];Altitude[km]") ;
  g2_1->SetTitle(";Longtitude to km [east][km]; Latitude to km [north][km];Altitude[km]") ;
  g2_2->SetTitle(";Longtitude to km [east][km]; Latitude to km [north][km];Altitude[km]") ;
  g2_10->SetTitle(";Longtitude to km [east][km]; Latitude to km [north][km];Altitude[km]") ;
  g2_50->SetTitle(";Longtitude to km [east][km]; Latitude to km [north][km];Altitude[km]") ;
  printf ("starting T1 \n");
  TTree *T1=new TTree("T1","T1");
  T1->ReadFile("R30km_0.1.csv","type/C:lat/F:longi/F:elav/F");
  N=T1->Draw("lat:longi:elav","");
  xv1=T1->GetV2();
  yv1=T1->GetV1();
  zv1=T1->GetV3();
  for (int i=0; i<N; i++) {
    Float_t xKm=xv1[i]*Rlong-Clong;
    Float_t yKm=yv1[i]*Rlati-Clati;
    Float_t zKm=zv1[i]*zFactor;
    g2_1->SetPoint(g2_1->GetN(),xKm,yKm,zKm);
    h2_1->SetBinContent(h2_1->FindBin(xKm, yKm),zKm);
  }



  printf ("starting T50 \n");
  TTree *T50=new TTree("T50","T50");
  T50->ReadFile("grid50km_1000_withAl.txt","type/C:lat/F:longi/F:elav/F");
  N=T50->Draw("lat:longi:elav","");
  xv50=T50->GetV2();
  yv50=T50->GetV1();
  zv50=T50->GetV3();
  for (int i=0; i<N; i++) {
    Float_t xKm=xv50[i]*Rlong-Clong;
    Float_t yKm=yv50[i]*Rlati-Clati;
    Float_t zKm=zv50[i]*zFactor;
    g2_50->SetPoint(g2_50->GetN(),xKm,yKm,zKm);
    if (fabs(xKm)>10 || fabs(yKm)>10)
      {
	g2->SetPoint(g2->GetN(),xKm,yKm,zKm);
	h2_50->SetBinContent(h2_50->FindBin(xKm, yKm),zKm);
      }
  }
  printf ("starting T10 \n");

  TTree *T10=new TTree("T10","T10");
  T10->ReadFile("grid10km_1000_withAl.txt","type/C:lat/F:longi/F:elav/F");
  N=T10->Draw("lat:longi:elav","");

  xv10=T10->GetV2();
  yv10=T10->GetV1();
  zv10=T10->GetV3();
  if (N>1e6) N=1e6;
  for (int i=0; i<N; i++) {
    Float_t xKm=xv10[i]*Rlong-Clong;
    Float_t yKm=yv10[i]*Rlati-Clati;
    Float_t zKm=zv10[i]*zFactor;

    //    printf ("T10:: %d %f %f %f \n",g2_10->GetN(),xKm,yKm,zKm);

    g2_10->SetPoint(g2_10->GetN(),xKm,yKm,zKm);
    if (fabs(xKm)>2 || fabs(yKm)>2)
      {
	g2->SetPoint(g2->GetN(),xKm,yKm,zKm);
	h2_10->SetBinContent(h2_10->FindBin(xKm, yKm),zKm);
      }
 }

  printf ("starting T2 \n");

  TTree *T2=new TTree("T2","T2");
  T2->ReadFile("grid2km_500_withAl.txt","type/C:lat/F:longi/F:elav/F");
  N=T2->Draw("lat:longi:elav","");
  xv2=T2->GetV2();
  yv2=T2->GetV1();
  zv2=T2->GetV3();
  for (int i=0; i<N; i++) {
    Float_t xKm=xv2[i]*Rlong-Clong;
    Float_t yKm=yv2[i]*Rlati-Clati;
    Float_t zKm=zv2[i]*zFactor;
    h2_2->SetBinContent(h2_2->FindBin(xKm, yKm),zKm);
    g2->SetPoint(g2->GetN(),xKm,yKm,zKm);
    g2_2->SetPoint(g2_2->GetN(),xKm,yKm,zKm);
  }
  printf ("Done Ts\n");
  TH2F *h2road=(TH2F*) h2_2->Clone("h2road");
 TGraph2D *g2Water=new TGraph2D();
 g2Water->SetPoint(0,35+32./60+45./3600,33+10./60+51./3600,737.1*zFactor );
 g2Water->SetPoint(1,35+32./60+45./3600,33+10./60+51./3600, 66*zFactor);
 g2Water->SetPoint(2,35+32./60+58./3600,33+10./60+51./3600, 6*zFactor);
 g2Water->SetPoint(3,35+32./60+60./3600,33+10./60+51./3600, 6*zFactor);
 for (int i=0; i<4; i++) {
   g2Water->GetX()[i]=g2Water->GetX()[i]*Rlong-Clong;
   g2Water->GetY()[i]=g2Water->GetY()[i]*Rlati-Clati;

 }
 g2Water->SetLineWidth(3);
 g2Water->SetLineColor(4);


 TTree *Tline=new TTree("Tline","T");
 N=Tline->ReadFile("powerLine.txt","type/C:lat/F:longi/F:elav/F");

 TGraph2D *g2line=new TGraph2D(); // not filled yet
 TTree *Troad=new TTree("Troad","T");
 Troad->ReadFile("road2.txt","type/C:lat/F:longi/F:elav/F");
 N=Troad->Draw("lat:longi:elav","");
 TGraph2D *g2Road=new TGraph2D();
 Double_t *xv=Troad->GetV2();
 Double_t *yv=Troad->GetV1();
 Double_t *zv=Troad->GetV3();
 for (int i=0; i<N; i++) {
   g2Road->SetPoint(i,xv[i]*Rlong-Clong,yv[i]*Rlati-Clati,zv[i]*zFactor);
   h2road->Fill(xv[i]*Rlong-Clong,yv[i]*Rlati-Clati,zv[i]*zFactor);
 }
 TTree *Tyetzoor=new TTree("Tyetzoor","T");
 Tyetzoor->ReadFile("yetzoor.txt","type/C:lat/F:longi/F:elav/F");
 N=Tyetzoor->Draw("lat:longi:elav","");
 TGraph2D *g2Yetzoor=new TGraph2D();
 xv=Tyetzoor->GetV2();
 yv=Tyetzoor->GetV1();
 zv=Tyetzoor->GetV3();
 for (int i=0; i<N; i++) g2Yetzoor->SetPoint(i,xv[i]*Rlong-Clong,yv[i]*Rlati-Clati,zv[i]*zFactor);
 g2->SetMaximum(1);
 // g2->SetMinimum(0);
 /* g2->GetXaxis()->SetRangeUser(-4,4);
 g2->GetYaxis()->SetRangeUser(-4,4);
 g2->GetXaxis()->SetRangeUser(-1,2);
 g2->GetYaxis()->SetRangeUser(-1,2);
 g2->GetXaxis()->SetRangeUser(-11,11);
 g2->GetYaxis()->SetRangeUser(-11,11);
 */

 //g2->Draw("surf1z");
 //g2->Draw("pcolz");
 // g2->Draw("TRI1z");
 g2->Draw("TRI2z");
 g2Road->SetLineColor(2);
 g2Yetzoor->SetLineColor(3);
 g2Road->Draw("same:line");

 TFile *f=new TFile("ManaraData.root","recreate");
 g2->Write("g2");
 g2_50->Write("g2_50");
 g2_10->Write("g2_10");
 g2_2->Write("g2_2");
 g2_1->Write("g2_1");
 g2Water->Write("g2Water");
 g2Yetzoor->Write("g2Yetzoor");
 g2Road->Write("g2Road");
 h2_50->Write("h2_50");
 h2_10->Write("h2_10");
 h2_2->Write("h2_2");
 h2_1->Write("h2_1");
 f->Close();
}

void draw() {
  TFile *f = new TFile("ManaraData.root");
  gStyle->SetPalette(1);
  //gStyle->SetPalette(54);

  TGraph2D *g2=(TGraph2D*) f->Get("g2");
  g2->SetMaximum(2.9);
  g2->SetMinimum(-0.2);
  TGraph2D *g2_10=(TGraph2D*) f->Get("g2_10");
  TGraph2D *g2_2=(TGraph2D*) f->Get("g2_2");
  TGraph2D *g2_1=(TGraph2D*) f->Get("g2_1");
  TGraph2D *g2Yetzoor=(TGraph2D*) f->Get("g2Yetzoor");
  TGraph2D *g2Road=(TGraph2D*) f->Get("g2Road");
  TGraph2D *g2Water=(TGraph2D*) f->Get("g2Water");
  //   g2_2->Draw("TRI2z");
  g2Road->SetLineColor(2);
  g2Yetzoor->SetLineColor(3);

  TGraph2D *g22=new TGraph2D();
  for (int i=0; i<g2->GetN(); i++) {
    Double_t x=g2->GetX()[i];
    Double_t y=g2->GetY()[i];
    Double_t z=g2->GetZ()[i];
    if (x<1.9 && y<1.9 && x>-1.9 && y>-1.9) g22->SetPoint(g22->GetN(),x,y,z);

  }

  TCanvas *c1=new TCanvas("c1","c1",1500,500);
  c1->Divide(3,1);

  c1->cd(1);
  gPad->SetTheta(90);
  gPad->SetPhi(0);
  gPad->SetMargin(0.1,0.2,0.13,0.1);
  g2->SetTitle("R=50Km;Longtitude to km [east][km]; Latitude to km [north][km];Altitude[km]") ;
  // g2->Draw("surf2z");
  g2->GetYaxis()->SetRangeUser(-50,30);

  g2->GetYaxis()->SetTitleOffset(-1);
  g2->GetXaxis()->SetTitleOffset(1);
  g2->GetZaxis()->SetLabelOffset(-1);
  g2->GetXaxis()->SetTitleOffset(1.4);
  g2->GetXaxis()->SetLabelOffset(-0.1);
  g2->GetZaxis()->SetLabelOffset(0.01);
  g2->GetZaxis()->SetTitleOffset(1.3);

  g2->Draw("TRI2z");
  g2Road->Draw("same:line");
  g2Water->Draw("same:line");
  g2Yetzoor->Draw("same:line");


  gPad->Update();


  c1->cd(2);
  gPad->SetTheta(90);
  gPad->SetPhi(0);
  gPad->SetMargin(0.1,0.2,0.13,0.1);
  g2_10->SetTitle("R=10Km;Longtitude to km [east][km]; Latitude to km [north][km];Altitude[km]") ;
  g2_10->GetHistogram()->GetYaxis()->SetTitleOffset(-1);
  g2_10->GetHistogram()->GetXaxis()->SetTitleOffset(1);
  g2_10->GetHistogram()->GetZaxis()->SetLabelOffset(-1);
  g2_10->GetHistogram()->GetXaxis()->SetTitleOffset(1.4);
  g2_10->GetHistogram()->GetXaxis()->SetLabelOffset(-0.1);
  g2_10->GetHistogram()->GetZaxis()->SetLabelOffset(0.01);
  g2_10->GetHistogram()->GetZaxis()->SetTitleOffset(1.3);
  //g2_10->Draw("SURF2Z");
  g2_10->Draw("TRI2Z");
  g2Road->Draw("same:line");
  g2Water->Draw("same:line");
  g2Yetzoor->Draw("same:line");
  gPad->Update();


  c1->cd(3);
  gPad->SetTheta(90);
  gPad->SetPhi(0);
  gPad->SetMargin(0.1,0.2,0.13,0.1);
  g2_2->SetTitle("R=2Km;Longtitude to km [east][km]; Latitude to km [north][km];Altitude[km]") ;
  //  g22->SetFillColor(0);
  //  g22->Draw("SURF2Z");
  gPad->SetGrid();
  g2_2->Draw("TRI2Z");

  g2Road->Draw("same:line");
  g2Water->Draw("same:line");
  g2Yetzoor->Draw("same:line");
  gPad->SetGrid();

  g2_2->GetHistogram()->GetYaxis()->SetTitleOffset(-1);
  g2_2->GetHistogram()->GetXaxis()->SetTitleOffset(1);
  g2_2->GetHistogram()->GetZaxis()->SetLabelOffset(-1);
  g2_2->GetHistogram()->GetXaxis()->SetTitleOffset(1.4);
  g2_2->GetHistogram()->GetXaxis()->SetLabelOffset(-0.1);
  g2_2->GetHistogram()->GetZaxis()->SetLabelOffset(0.01);
  g2_2->GetHistogram()->GetZaxis()->SetTitleOffset(1.3);
    gPad->Update();

  c1->SaveAs("Manara_xy_geography.png");


  TCanvas *c2=new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);

  c2->cd(1);
  gPad->SetMargin(0.1,0.1,0.1,0.1);

  gPad->SetTheta(10);
  gPad->SetPhi(2);
 // g22->Draw("SURF2");
  g2_1->SetTitle("R=1-2Km;east[km]; Latitude to km north[km];Altitude[km]") ;

  g2_1->Draw("surf1");
  g2Road->Draw("same:line");
  g2Water->Draw("same:line");
  g2Yetzoor->Draw("same:line");
  g2_1->GetHistogram()->GetXaxis()->SetRangeUser(-1.,2);

  g2_1->GetHistogram()->GetXaxis()->SetTitleOffset(1.4);
  g2_1->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
  g2_1->GetHistogram()->GetZaxis()->SetTitleOffset(1.4);

  gPad->Modified(); gPad->Update();

   c2->cd(2);
  gPad->SetMargin(0.1,0.1,0.1,0.1);

  gPad->SetTheta(10);
  gPad->SetPhi(2);
 // g22->Draw("SURF2");
  g2_10->Draw("surf1");
  g2Road->Draw("same:line");
  g2Water->Draw("same:line");
  g2Yetzoor->Draw("same:line");
  g2_10->SetTitle("R=10Km;east[km]; Latitude to km north[km];Altitude[km]") ;

  g2_10->GetHistogram()->GetXaxis()->SetTitleOffset(1.4);
  g2_10->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
  g2_10->GetHistogram()->GetZaxis()->SetTitleOffset(1.4);

  gPad->Modified(); gPad->Update();
  c2->SaveAs("Manara_xyz_geography.png");

  // g22->GetYaxis()->SetTitleOffset(-1);
 // g22->GetXaxis()->SetTitleOffset(1);
//  g22->GetZaxis()->SetTitleOffset(-1);





}

TH1F *rays(Int_t N) {
  return rays(N,0,0,7./1000, 1);
  }
TH1F *rays(Int_t N,Double_t x0, Double_t y0, Double_t z0, Int_t color) {
  TFile *f = new TFile("ManaraData.root");
  TGraph2D *g2=(TGraph2D*) f->Get("g2");
  TGraph2D *g2_10=(TGraph2D*) f->Get("g2_10");
  TGraph2D *g2_2=(TGraph2D*) f->Get("g2_2");
  TGraph2D *g2_1=(TGraph2D*) f->Get("g2_1");

  TH2D *h2_50=(TH2D*) f->Get("h2_50");
  TH2D *h2_10=(TH2D*) f->Get("h2_10");
  TH2D *h2_2=(TH2D*) f->Get("h2_2");
  TH1F *hR=new TH1F("hR","hR;R [km]",1000,0,10);
  hR->SetLineColor(color);
  TH1F *hRcos=new TH1F("hRcos","hRcos;R [km]",1000,0,10);
  hRcos->SetLineColor(2);
  hRcos->SetLineColor(color);

  TRandom3 rand;
  rand.SetSeed(0);
  /* int N=100;
  Double_t x0=0;
  Double_t y0=0;
  Double_t z0=7./1000.;
  */
  TGraph2D *gs[1000];
  TGraph2D *gThetaPhi=new TGraph2D();
  TGraph2D *gCThetaPhi=new TGraph2D();
  TGraph2D *gShadow=new TGraph2D();
  TGraph *gPattern=new TGraph();
  int ok=0;
  int nok=0;
  TH2F *hout=new TH2F("hout","h;cos(theta);phi",100,-1,1,100,-TMath::Pi(),TMath::Pi());
  TH2F *houtc=new TH2F("houtc","h;cos(theta);phi",100,-1,1,100,-TMath::Pi(),TMath::Pi());
  for (int i=0; i<N; i++){
    if (i<1000) {
      gs[i]=new TGraph2D();
      gs[i]->SetPoint(0,x0,y0,z0);
      }
//   Double_t r1=rand.Uniform(0,1);
//   Double_t r2=rand.Uniform(0,1);
//   Double_t phi = acos(2*r2-1);
//   Double_t theta = 2*TMath::Pi()*r1;
 //  Double_t x=1*cos(theta)*sin(phi);
 //  Double_t y=1*sin(theta)*sin(phi);
//   Double_t z=1*cos(phi);
   Double_t x,y,z;
   rand.Sphere(x,y,z,1);
   TVector3 v3(x,y,z);
   Double_t theta = v3.Theta();
   Double_t phi = v3.Phi();
   if (z<0) {
//     printf ("z<0 phi=%f r2=%f \n",phi,r2);
     continue;
     }
   Double_t x1=x+x0;
   Double_t y1=y+y0;
   Double_t z1=z+z0;

	  /*
	    Two points on a line:  x1,x0:

	    x=(1-u)*x0 + u*x1
	    y=(1-u)*y0 + u*y1
	    z=(1-u)*z0 + u*z1

	    when u=0, r=(x0,y0,z0), when u=1 r=(x1,y1,z1)

	    To scan up to altitude of 10 km, u should go up to u=(z-z0)/(z1-z0)=(10-z0)/(z1-z0)

	  */

   Float_t u;
   //   printf ("%f %f %f  R= %f \n ",x1Test,y1Test,z1Test,sqrt(x1Test*x1Test+y1Test*y1Test+z1Test*z1Test));
   Double_t pzG=0;
   Double_t zG=0;
   Float_t umax=(10-z0)/(z1-z0);
   //      printf ("umax= %f \n",umax);
   Int_t gotit=0;

   for (u=0; u<=umax; u=u+(umax-u)/1000) {
      Double_t xt=(1-u)*x0+x1*u;
      Double_t yt=(1-u)*y0+y1*u;
      Double_t zt=(1-u)*z0+z1*u;
      pzG=zG;

      if (fabs(xt)>50 || fabs(yt)>50) break;
      Float_t hillH=0;
      hillH=h2_2->GetBinContent(h2_2->FindBin(xt,yt));
//      printf ("r=%f by 2 = %f \n",sqrt(x*x+y*y),hillH);
      if (hillH==0) 	{
	hillH=h2_10->GetBinContent(h2_10->FindBin(xt,yt));
//	printf ("\t by 10 = %f \n",hillH);
      }
      if (hillH==0) 	{
	hillH=h2_50->GetBinContent(h2_50->FindBin(xt,yt));
	//	printf ("\t\t by 50 = %f \n",hillH);
	}
      if (hillH==0)     {
	hillH=g2->Interpolate(xt,yt);
	//printf ("\t\t\t by g = %f \n",hillH);
	}
      if (hillH==0) break;
      zG=zt-hillH;
      //      printf ("%f %f - ",z,hillH);

//      if (pzG<0 && zG>0) {
      if (zG>0) {
	gotit=1;
	ok++;
	Double_t Rhill=sqrt((hillH-z0)*(hillH-z0)+(xt-x0)*(xt-x0)+(yt-y0)*(yt-y0));
	gShadow->SetPoint(gShadow->GetN(),x,y,z);
	//	gShadow->SetPoint(gShadow->GetN(),x,y,hillH);
	//	if (Rhill>10 ) printf ("%d u=%f theta=%f pho=%f  \t %f %f hillH=%f R=%f (%f,%f,%f)\n",i,u/umax,theta,phi,zG,pzG,hillH,Rhill,xt,yt,zt);
	gThetaPhi->SetPoint(gThetaPhi->GetN(),phi,theta,Rhill);
	gThetaPhi->SetTitle("gThetaPhi;phi;theta;Rhill");
	gCThetaPhi->SetPoint(gThetaPhi->GetN(),phi,cos(theta),Rhill);
	gCThetaPhi->SetTitle("gThetaPhi;phi;cos(theta);Rhill");
	hR->Fill(Rhill);
	hRcos->Fill(Rhill,cos(theta)*cos(theta));
	gPattern->SetPoint(gPattern->GetN(),theta,phi);
	hout->Fill(cos(theta),phi,Rhill);
	houtc->Fill(cos(theta),phi,1);

    //	gs[i]->SetPoint(2,x1,y1,z1);
      if (int(i/1000)==i/1000.) printf ("aa=%d \n",i);

	if (i<100) {
	  gs[i]->SetPoint(1,xt,yt,zt);
	  //	  gs[i]->Draw("line:same:*");
	}
//	printf ("done\n");
	break;}
  // printf ("\tu=%f (%f %f %f) \t zG=%f \t surface=%f \n",u,x,y,z,zG,z-zG);
  }

if (gotit==0) {
    nok++;
    if (i<100) {
      gs[i]->SetPoint(1,x1,y1,x1);
      gs[i]->SetLineColor(2);

    }
}
// printf ("\n");



  }
  printf ("done loop\n");
  TGraph2D *g2Yetzoor=(TGraph2D*) f->Get("g2Yetzoor");
  TGraph2D *g2Road=(TGraph2D*) f->Get("g2Road");
  TGraph2D *g2Water=(TGraph2D*) f->Get("g2Water");
  //   g2_2->Draw("TRI2z");
  g2Road->SetLineColor(2);
  g2Yetzoor->SetLineColor(3);

  TCanvas *cS=new TCanvas("cS","cS",1200,400);
  cS->Divide(3,1);
  cS->cd(1);
  g2_2->GetHistogram()->GetXaxis()->SetRangeUser(-1.,2);
  gPad->SetTheta(10);
  gPad->SetPhi(2);
  g2_2->Draw("TRI2z");
  g2Road->Draw("same:line");
  g2Water->Draw("same:line");
  g2Yetzoor->Draw("same:line");
  Int_t imax=0;
  if (nok<imax) imax=nok;
  for (int i=0; i<imax; i++) {
    gs[i]->Draw("same:line");
   }
  TGraph2D *gDetector=new TGraph2D();
  gDetector->SetPoint(0,-10,-10,-10);
  gDetector->SetPoint(1,-10,-10,-20);
  gDetector->SetPoint(2,x0,y0,z0);
  gDetector->SetMarkerColor(color);


  gDetector->Draw("p0:same");
  g2_2->GetHistogram()->GetYaxis()->SetTitleOffset(-1);
  g2_2->GetHistogram()->GetXaxis()->SetTitleOffset(1);
  g2_2->GetHistogram()->GetZaxis()->SetLabelOffset(-1);
  g2_2->GetHistogram()->GetXaxis()->SetTitleOffset(1.4);
  g2_2->GetHistogram()->GetXaxis()->SetLabelOffset(-0.1);
  g2_2->GetHistogram()->GetZaxis()->SetLabelOffset(0.01);
  g2_2->GetHistogram()->GetZaxis()->SetTitleOffset(1.3);

  cS->cd(2);
    gPad->SetTheta(90);
  gPad->SetPhi(0);
  /*g2_1->GetHistogram()->GetYaxis()->SetTitleOffset(-1);
  g2_1->GetHistogram()->GetXaxis()->SetTitleOffset(1);
  g2_1->GetHistogram()->GetZaxis()->SetLabelOffset(-1);
  g2_1->GetHistogram()->GetXaxis()->SetTitleOffset(1.4);
  g2_1->GetHistogram()->GetXaxis()->SetLabelOffset(-0.1);
  g2_1->GetHistogram()->GetZaxis()->SetLabelOffset(0.01);
  g2_1->GetHistogram()->GetZaxis()->SetTitleOffset(1.3);
  */

  h2_2->Draw("colz");
  g2Road->Draw("same:line");
  g2Water->Draw("same:line");
  g2Yetzoor->Draw("same:line");
  for (int i=0; i<imax; i++) {
    gs[i]->Draw("same:line");
   }
  gDetector->Draw("p0:same");

  cS->cd(3);
  gPad->SetLogy();

	hRcos->Draw("hist");
	hR->Draw("same");
	hRcos->Draw("same");

  TCanvas *c4=new TCanvas("c4","c4",1000,600);

  c4->cd(1);
  //  gShadow->Draw("PCOL");
  gShadow->Draw("P");
  c4->cd(3);
   gPad->SetLogz();
   gThetaPhi->Draw("SURF");
   hout->Draw("colz");
   c4->cd(4);
   houtc->Draw("colz");
   printf ("ok=%d, nok=%d \n",ok,nok);


   TCanvas *c5=new TCanvas("c5","c5",1200,1200);
   c5->Divide(2,2);
   c5->cd(1);
   g2_1->Draw("surf1");
   g2Road->Draw("same:line");
   g2Water->Draw("same:line");
   g2Yetzoor->Draw("same:line");
   g2_1->GetHistogram()->GetXaxis()->SetRangeUser(-1.,2);
   if (nok<imax) imax=nok-1;

   for (int i=0; i<imax; i++) {
     gs[i]->Draw("same:line");
   }
   c5->cd(2);
   hout->Draw("colz");
   c5->cd(3);
   gThetaPhi->Draw("SURF");
   c5->cd(4);
   houtc->Draw("colz");
   printf ("ok=%d, nok=%d \n",ok,nok);
   TFile *ff=new TFile(Form("a_%f_%f_%f_%d.root",x0,y0,z0,N),"recreate");
   hout->Write("hout");
   houtc->Write("houtc");
   gThetaPhi->Write("gThetaPhi");
   gCThetaPhi->Write("gCThetaPhi");
   hRcos->Write("hRcos");
   hR->Write("hR");
   gShadow->Write("gShadow");
   gPattern->Write("gPattern");
   gDetector->Write("gDetector");
   ff->Close();

   return hRcos;

 }

void go() {
  int N=1000000;
  TFile *f = new TFile("ManaraData.root");
  TGraph2D *g2_1=(TGraph2D*) f->Get("g2_1");
  TH2F *h2_1=(TH2F*) f->Get("h2_1");

  TGraph2D *g2_2=(TGraph2D*) f->Get("g2_2");
  TGraph2D *g2Yetzoor=(TGraph2D*) f->Get("g2Yetzoor");
  TGraph2D *g2Road=(TGraph2D*) f->Get("g2Road");
  TGraph2D *g2Water=(TGraph2D*) f->Get("g2Water");

  g2_1->GetHistogram()->GetXaxis()->SetRangeUser(-1.,2);



  TGraph2D *gDetector0=new TGraph2D();
  gDetector0->SetPoint(0,-10,-10,-10);
  gDetector0->SetPoint(1,-10,-10,-20);
  gDetector0->SetPoint(2,0,0,0.06);
  gDetector0->SetMarkerColor(1);
  gDetector0->SetMarkerSize(1.1);
  TH1F *h0=rays(N,0,0,0.06,1);
  h0->SetLineColor(1);


  TGraph2D *gDetector1=new TGraph2D();
  gDetector1->SetPoint(0,-10,-10,-10);
  gDetector1->SetPoint(1,-10,-10,-20);
  gDetector1->SetPoint(2,0,0,0);
  gDetector1->SetMarkerColor(2);
  gDetector1->SetMarkerSize(1.1);

  TH1F *h1=rays(N,0,0,0,2);
  h1->SetLineColor(2);


  TGraph2D *gDetector2=new TGraph2D();
  gDetector2->SetPoint(0,-10,-10,-10);
  gDetector2->SetPoint(1,-10,-10,-20);
  gDetector2->SetPoint(2,0.4,-0.04,0);
  gDetector2->SetMarkerColor(3);
  TH1F *h2=rays(N,0.4,-0.04,0,3);
  h2->SetLineColor(3);

  TGraph2D *gDetector3=new TGraph2D();
  gDetector3->SetPoint(0,-10,-10,-10);
  gDetector3->SetPoint(1,-10,-10,-20);
  gDetector3->SetPoint(2,0,0,-0.07);
  gDetector3->SetMarkerColor(4);
//  gDetector3->SetMarkerSize(0.8);
  TH1F *h3=rays(N,0.,0,-0.07,4);
  h3->SetLineColor(4);

  TGraph2D *gDetector4=new TGraph2D();
  gDetector4->SetPoint(0,-10,-10,-10);
  gDetector4->SetPoint(1,-10,-10,-20);
  gDetector4->SetPoint(2, 0.4,-0.04,-0.07);
  gDetector4->SetMarkerColor(5);
  TH1F *h4=rays(N ,0.4,-0.04,-0.07,5);
  h4->SetLineColor(5);



  TGraph2D *gDetector5=new TGraph2D();
  gDetector5->SetPoint(0,-10,-10,-10);
  gDetector5->SetPoint(1,-10,-10,-20);
  gDetector5->SetPoint(2,1.5,-0.18,0.00);
  gDetector5->SetMarkerColor(6);
  TH1F *h5=rays(N,1.5,-0.18,0.00,6);
  h5->SetLineColor(6);
  TGraph2D *gDetector6=new TGraph2D();
  gDetector6->SetPoint(0,-10,-10,-10);
  gDetector6->SetPoint(1,-10,-10,-20);
  gDetector6->SetPoint(2,1.5,-0.18,0.08);
  gDetector6->SetMarkerColor(7);
  TH1F *h6=rays(N,1.5,-0.18,0.08,7);
  h6->SetLineColor(7);

  TCanvas *cSg=new TCanvas("cSg","cSg",1200,400);
  cSg->Divide(3,1);
  cSg->cd(1);
  g2_2->GetHistogram()->GetXaxis()->SetRangeUser(-1.,2);
  gPad->SetTheta(10);
  gPad->SetPhi(2);
  h2_1->SetMinimum(-0.1);
  h2_1->Draw("SURF2");
  g2Road->Draw("same:line");
  g2Water->Draw("same:line");
  g2Yetzoor->Draw("same:line");
  gDetector0->Draw("p0:same");
  gDetector1->Draw("p0:same");
  gDetector2->Draw("p0:same");
  gDetector3->Draw("p0:same");
  gDetector4->Draw("p0:same");
  gDetector5->Draw("p0:same");
  gDetector6->Draw("p0:same");

  g2_2->GetHistogram()->GetYaxis()->SetTitleOffset(-1);
  g2_2->GetHistogram()->GetXaxis()->SetTitleOffset(1);
  g2_2->GetHistogram()->GetZaxis()->SetLabelOffset(-1);
  g2_2->GetHistogram()->GetXaxis()->SetTitleOffset(1.4);
  g2_2->GetHistogram()->GetXaxis()->SetLabelOffset(-0.1);
  g2_2->GetHistogram()->GetZaxis()->SetLabelOffset(0.01);
  g2_2->GetHistogram()->GetZaxis()->SetTitleOffset(1.3);

  cSg->cd(2);
  gPad->SetTheta(90);
  gPad->SetPhi(0);

  h2_1->Draw("surf2");
  g2Road->Draw("same:line");
  g2Water->Draw("same:line");
  g2Yetzoor->Draw("same:line");
  g2_1->GetHistogram()->GetYaxis()->SetTitleOffset(-1);
  g2_1->GetHistogram()->GetXaxis()->SetTitleOffset(1);
  g2_1->GetHistogram()->GetZaxis()->SetLabelOffset(-1);
  g2_1->GetHistogram()->GetXaxis()->SetTitleOffset(1.4);
  g2_1->GetHistogram()->GetXaxis()->SetLabelOffset(-0.1);
  g2_1->GetHistogram()->GetZaxis()->SetLabelOffset(0.01);
  g2_1->GetHistogram()->GetZaxis()->SetTitleOffset(1.3);
  gDetector0->Draw("p0:same");
  gDetector1->Draw("p0:same");
  gDetector2->Draw("p0:same");
  gDetector3->Draw("p0:same");
  gDetector4->Draw("p0:same");
  gDetector5->Draw("p0:same");
  gDetector6->Draw("p0:same");


  cSg->cd(3);
 // gPad->SetLogy();
  h0->Draw("hist");
  h1->Draw("same:hist");
  h2->Draw("same:hist");
  h3->Draw("same:hist");
  h4->Draw("same:hist");
  h5->Draw("same:hist");
  h6->Draw("same:hist");
  TFile *ff=new TFile("a.root","recreate");
  h0->Write("h0");
  h1->Write("h1");
  h2->Write("h2");
  h3->Write("h3");
  h4->Write("h4");
  h5->Write("h5");
  h6->Write("h6");
  gDetector0->Write("gDetector0");
  gDetector1->Write("gDetector1");
  gDetector2->Write("gDetector2");
  gDetector3->Write("gDetector3");
  gDetector4->Write("gDetector4");
  gDetector5->Write("gDetector5");
  gDetector6->Write("gDetector6");
  ff->Close();
  cSg->SaveAs("cSg.png");


}


void drawSummary() {

  TFile *fM = new TFile("ManaraData.root");
 TGraph2D *g2_1=(TGraph2D*) fM->Get("g2_1");
  TH2F *h2_1=(TH2F*) fM->Get("h2_1");
  TGraph2D *g2Yetzoor=(TGraph2D*) fM->Get("g2Yetzoor");
  TGraph2D *g2Road=(TGraph2D*) fM->Get("g2Road");
  TGraph2D *g2Water=(TGraph2D*) fM->Get("g2Water");
  g2_1->GetHistogram()->GetXaxis()->SetRangeUser(-1.,2);
  g2_1->SetNpx(30);
  g2_1->SetNpy(30);
  g2_1->SetDirectory(0);
  g2_1->GetHistogram()->GetYaxis()->SetTitleOffset(-1);
  g2_1->GetHistogram()->GetXaxis()->SetTitleOffset(1);
  g2_1->GetHistogram()->GetZaxis()->SetLabelOffset(-1);
  g2_1->GetHistogram()->GetXaxis()->SetTitleOffset(1.4);
  g2_1->GetHistogram()->GetXaxis()->SetLabelOffset(-0.1);
  g2_1->GetHistogram()->GetZaxis()->SetLabelOffset(0.01);
  g2_1->GetHistogram()->GetZaxis()->SetTitleOffset(1.3);
  g2_1->SetTitle("; Longitude [Km east]; Latitude [Km North]; altitude [Km]");
  g2Yetzoor->SetLineColor(1);
  g2Road->SetDirectory(0);
  g2Yetzoor->SetDirectory(0);
  g2Water->SetDirectory(0);
  const int n=5;
  TH1F *hRs[n];
  TGraph2D *gDet[n];
  TFile *fdata = new TFile("a.root");
  //  TExec *ex1 = new TExec("ex1","gStyle->SetPalette(kLightTemperature);");
//  TExec *ex2 = new TExec("ex2","gStyle->SetPalette(56);");
  gStyle->SetPalette(kLightTemperature);
  for (int i=0; i<n; i++) {
     hRs[i]=(TH1F*) fdata->Get(Form("h%d",i));
     gDet[i]=(TGraph2D*) fdata->Get(Form("gDetector%d",i));
     hRs[i]->SetDirectory(0);
     }
  hRs[2]->SetLineColor(6);
  gDet[2]->SetMarkerColor(6);

  hRs[4]->SetLineColor(3);
  gDet[4]->SetMarkerColor(3);
 // hRs[6]->SetLineColor(5);
//  gDet[6]->SetMarkerColor(5);

  TCanvas *cc=new TCanvas("cc","cc",1200,500);
  cc->Divide(3,1);
  //cc->Divide(2,1);
  cc->cd(1);
  gPad->SetTheta(0); //10
  gPad->SetPhi(0); //2
  g2_1->SetMinimum(-0.1);
  g2_1->Draw("SURF1:X+");
  g2Road->Draw("same:line");
  g2Water->Draw("same:line");
  g2Yetzoor->Draw("same:line");
  for (int i=0; i<n;i++) {
    gDet[i]->Draw("p0:same");
  }
  cc->cd(2);
  gPad->SetGrid();
  gPad->SetTheta(90);
  gPad->SetPhi(0);
  g2_1->SetMinimum(-0.1);
  g2_1->Draw("SURF2");
  g2Road->Draw("same:line");
  g2Water->Draw("same:line");
  g2Yetzoor->Draw("same:line");
  for (int i=0; i<n;i++) {
    gDet[i]->Draw("p0:same");
  }

  cc->cd(3);
  gStyle->SetOptStat(0);
  gPad->SetLogy();
  gPad->SetGrid();
  hRs[0]->SetTitle(";R_{rock} [Km];frac/bin");
  hRs[0]->GetXaxis()->SetRangeUser(0.4,2.3);
  hRs[0]->Draw("hist");
  Double_t Norm;
  for (int i=0; i<n;i++) {
    Norm=hRs[i]->Integral();
    hRs[i]->Scale(1./Norm);
    hRs[i]->Draw("hist:same");
    hRs[i]->SetLineWidth(2);
  }

}

void drawFile(){
  TFile *fM = new TFile("ManaraData.root");
  TGraph2D *g2_1=(TGraph2D*) fM->Get("g2_1");
  TGraph2D *g2Yetzoor=(TGraph2D*) fM->Get("g2Yetzoor");
  TGraph2D *g2Road=(TGraph2D*) fM->Get("g2Road");
  TGraph2D *g2Water=(TGraph2D*) fM->Get("g2Water");
  g2_1->GetHistogram()->GetXaxis()->SetRangeUser(-1.,2);
  g2_1->SetNpx(20);
  g2_1->SetNpy(20);
  g2_1->SetMinimum(-0.1);

  //Float_t label=0.1;
  Float_t label=0.03;



  const int N=7;//7;
  TString files[N]={"a_0.000000_0.000000_0.000000_1000000.root","a_0.400000_-0.040000_-0.070000_1000000.root","a_0.000000_0.000000_0.060000_1000000.root",  "a_0.000000_0.000000_-0.070000_1000000.root",  "a_1.500000_-0.180000_0.000000_1000000.root","a_0.400000_-0.040000_0.000000_1000000.root" , "a_1.500000_-0.180000_0.080000_1000000.root"};


  TGraph2D *gDetectors[N];
  TGraph2D *gThetaPhis[N];
  TGraph2D *gCThetaPhis[N];
  TH1F *hRs[N];
  TH1F *hRcoss[N];
  g2_1->SetMinimum(-0.1);
  //  TCanvas *c=new TCanvas("c","c",1000,1000);
  // c->Divide(3,N);
  TCanvas *c[N];
  TExec *ex1 = new TExec("ex1","gStyle->SetPalette(kLightTemperature);");
  TExec *ex2 = new TExec("ex2","gStyle->SetPalette(56);");


  for (int i=0; i<N;i++) {
    c[i]=new TCanvas(Form("c%d",i),Form("c%d",i),1200,400);
    c[i]->Divide(3,1);

    TString file=files[i];
    printf ("Looking at %s \n",file.Data());
    TFile *f = new TFile(file.Data());
    gDetectors[i]=(TGraph2D*) f->Get("gDetector")->Clone(Form("gDetector%d",i));
    gDetectors[i]->SetDirectory(0);

    Double_t x0= gDetectors[i]->GetX()[2];
    Double_t y0= gDetectors[i]->GetY()[2];
    Double_t z0= gDetectors[i]->GetZ()[2];

    //    gThetaPhis[i]=(TGraph2D*) f->Get("gThetaPhi");
    gThetaPhis[i]=(TGraph2D*) f->Get("gThetaPhi")->Clone(Form("gThetaPhi%d",i));
    gThetaPhis[i]->SetDirectory(0);
    gThetaPhis[i]->GetHistogram()->GetXaxis()->SetLabelSize(label);
    gThetaPhis[i]->GetHistogram()->GetYaxis()->SetLabelSize(label);
    gThetaPhis[i]->GetHistogram()->GetZaxis()->SetLabelSize(label);
    gThetaPhis[i]->GetHistogram()->SetTitleSize(label*2);
    gThetaPhis[i]->GetHistogram()->GetXaxis()->SetTitleSize(label);
    gThetaPhis[i]->GetHistogram()->GetYaxis()->SetTitleSize(label);
    gThetaPhis[i]->GetHistogram()->GetZaxis()->SetTitleSize(label);
    gThetaPhis[i]->GetHistogram()->SetTitle(Form("Lab at %.2f,%.2f,%.2f; Phi; Theta; R_rock [km]",x0,y0,z0));
    gThetaPhis[i]->SetTitle(Form("Lab at %.2f,%.2f,%.2f; Phi; Theta; R_rock [km]",x0,y0,z0));

    //    gCThetaPhis[i]=(TGraph2D*) f->Get("gCThetaPhi");
    gCThetaPhis[i]=(TGraph2D*) f->Get("gCThetaPhi")->Clone(Form("gCThetaPhi%d",i));
    gCThetaPhis[i]->SetDirectory(0);
    gCThetaPhis[i]->GetXaxis()->SetLabelSize(label);
    gCThetaPhis[i]->GetYaxis()->SetLabelSize(label);
    gCThetaPhis[i]->GetZaxis()->SetLabelSize(label);

    gCThetaPhis[i]->GetHistogram()->SetTitleSize(label*2);
    gCThetaPhis[i]->GetHistogram()->GetXaxis()->SetTitleSize(label);
    gCThetaPhis[i]->GetHistogram()->GetYaxis()->SetTitleSize(label);
    gCThetaPhis[i]->GetHistogram()->GetZaxis()->SetTitleSize(label);
    gCThetaPhis[i]->GetHistogram()->SetTitle(Form("Lab at %.2f,%.2f,%.2f; Phi; cos(Theta); R_rock [km]",x0,y0,z0));
    gCThetaPhis[i]->SetTitle(Form("Lab at %.2f,%.2f,%.2f; Phi; cos(Theta); R_rock [km]",x0,y0,z0));
    gCThetaPhis[i]->SetNpx(100);
    gCThetaPhis[i]->SetNpy(100);

    hRs[i]=(TH1F*) (f->Get("hR")->Clone(Form("hRs%d",i)));
    hRs[i]->SetDirectory(0);
    hRs[i]->GetXaxis()->SetLabelSize(label);
    hRs[i]->GetXaxis()->SetTitleSize(label);
    hRs[i]->GetYaxis()->SetLabelSize(label);
    hRs[i]->GetXaxis()->SetRangeUser(0,5);
    hRs[i]->SetTitle(Form("%.2f,%.2f,%.2f (dash=uniform,solid=cos2theta; R_rock [Km]",x0,y0,z0));
    hRcoss[i]=(TH1F*) f->Get("hRcos");
    hRcoss[i]->SetDirectory(0);
    ex1->Draw();


    g2_1->SetTitle(Form("Lab at %.2f,%.2f,%.2f; east [Km]; North[km]; altitude [km]",x0,y0,z0));
//    c->cd(3*i+1);
    c[i]->cd(1);
    gPad->SetTheta(10);
    gPad->SetPhi(2);
    //gPad->SetMargin(0.1,0.1,0.2,0.1);

    g2_1->Draw("surf1");
    g2Yetzoor->Draw("line:same");
    g2Road->Draw("line:same");
    g2Water->Draw("line:same");
    gDetectors[i]->Draw("P0:same");

    /*
    c->cd(3*i+2);
    gPad->SetTheta(90);
    gPad->SetPhi(0);

    g2_1->Draw("surf2z");
    g2Yetzoor->Draw("line:same");
    g2Road->Draw("line:same");
    g2Water->Draw("line:same");
    gDetectors[i]->Draw("P0:same");

    */

    ex2->Draw();

//    gStyle->SetPalette(56);
//    c->cd(3*i+2);
    c[i]->cd(2);
    //gPad->SetMargin(0.1,0.12,0.2,0.1);

    gPad->SetLogz();
    gCThetaPhis[i]->GetHistogram()->GetXaxis()->SetLabelSize(label);
    gCThetaPhis[i]->GetHistogram()->GetYaxis()->SetLabelSize(label);
    gCThetaPhis[i]->GetHistogram()->GetZaxis()->SetLabelSize(label);
    gCThetaPhis[i]->GetHistogram()->GetXaxis()->SetTitleSize(label);
    gCThetaPhis[i]->GetHistogram()->GetYaxis()->SetTitleSize(label);
    gCThetaPhis[i]->GetHistogram()->GetZaxis()->SetTitleSize(label);
    //    gCThetaPhis[i]->GetHistogram()->GetYaxis()->SetTitleOffset(label);
    // gCThetaPhis[i]->GetHistogram()->GetZaxis()->SetTitleOffset(label);

//    gThetaPhis[i]->Draw("colz");
    gCThetaPhis[i]->GetHistogram()->Draw();

    gCThetaPhis[i]->Draw("colz:same");
    gPad->Update();

    //    c->cd(3*i+3);
    c[i]->cd(3);

    //    gPad->SetMargin(0.1,0.,0.25,0.2);

    hRs[i]->SetLineStyle(2);
    hRs[i]->Draw("h");
    if (hRcoss[i]->GetMaximum()>hRs[i]->GetMaximum()) hRs[i]->SetMaximum(hRcoss[i]->GetMaximum());
    gPad->SetLogy();
    gPad->SetGrid();
    hRcoss[i]->Draw("h:same");
    gPad->Update();
    c[i]->SaveAs(Form("Station%d.png",i));
    f->Close();
    delete f;
  }
}
/*
void getRays() {
  TFile *f = new TFile("ManaraData.root");
  TGraph2D *g2_2=(TGraph2D*) f->Get("g2_2");

  Double_t x0Test=0;
  Double_t y0Test=0;
  Double_t z0Test=7/1000.;
  Double_t t=0;
  TH1F *hR=new TH1F("hR","hR;R [km]",1000,0,50);
  TRandom *ra=new TRandom();
  Double_t x1Test,y1Test,z1Test;
  TGraph2D *gs[10000];
  int ok=0;
  int nok=0;
  TCanvas *c4=new TCanvas("c4","c4",1000,600);
  c4->Divide(2,1);
  c4->cd(1);
  g2->Draw("tri2z");
 for (int i=0; i<1000; i++){
   gs[i]=new TGraph2D();
   gs[i]->SetPoint(0,x0Test,y0Test,z0Test);
  ra->Sphere(x1Test,y1Test,z1Test,1);
  x1Test=x1Test+x0Test;
  y1Test=y1Test+y0Test;
  z1Test=z1Test+z0Test;
  if (z1>0) {
    Float_t u;
    printf ("%f %f %f  R= %f \n ",x1Test,y1Test,z1Test,sqrt(x1Test*x1Test+y1Test*y1Test+z1Test*z1Test));
    Double_t pzG=0;
    Double_t zG=0;
    Float_t umax=(10-z0)/(z1-z0);
    for (u=0; u<=umax; u=u+0.1) {
      Double_t x=(1-u)*x0Test+x1Test*u;
      Double_t y=(1-u)*y0Test+y1Test*u;
      Double_t z=(1-u)*z0Test+z1Test*u;
      // ga->SetPoint(ga->GetN(),x,y,z);
      pzG=zG;
      Float_t hillH=h2->GetBinContent(h2->FindBin(x,y));
      zG=z-hillH;
      if (pzG<0 && zG>0) { ok++; printf ("%f %f z=%f \n",zG,pzG,z); hR->Fill(sqrt(hillH*hillH+x*x+y*y));
	gs[i]->SetPoint(1,x,y,z);
	gs[i]->Draw("line:same");
	break;}
  // printf ("\tu=%f (%f %f %f) \t zG=%f \t surface=%f \n",u,x,y,z,zG,z-zG);
  }
  if (u==1) {printf ("nope \n"); nok++;}
  }
 }
printf ("ok=%d nok=%d \n",ok,nok);

c4->cd(2);
hR->Draw();


}

 TGraph2D *gray=new TGraph2D();
 gray->SetPoint(0,0,0,6e-3);
 gray->SetPoint(1,100,100,6e-3);
 gray->SetPoint(2,100,100,1);
 gray->SetPoint(3,100,200,6e-3);
 gray->SetLineColor(2);
 gray->Draw("same:line");


Double_t x0Test=0;
Double_t y0Test=0;
Double_t z0Test=7/1000.;

//Double_t x1=1;
//Double_t y1=1;
//Double_t z1=1.;
Double_t t=0;

TGraph *ga=new TGraph();
TH1F *hR=new TH1F("hR","hR;R [km]",1000,0,50);
TRandom *ra=new TRandom();
Double_t x1Test,y1Test,z1Test;
TGraph2D *gs[10000];
g2->Draw("tri2z");
int ok=0;
int nok=0;

TCanvas *c4=new TCanvas("c4","c4",1000,600);
c4->Divide(2,1);
c4->cd(1);
for (int i=0; i<1000; i++){
  gs[i]=new TGraph2D();
  gs[i]->SetPoint(0,x0Test,y0Test,z0Test);

  ra->Sphere(x1Test,y1Test,z1Test,1);
  x1Test=x1Test+x0Test;
  y1Test=y1Test+y0Test;
  z1Test=z1Test+z0Test;

  if (z1>0) {
    Float_t u;
    printf ("%f %f %f  R= %f \n ",x1Test,y1Test,z1Test,sqrt(x1Test*x1Test+y1Test*y1Test+z1Test*z1Test));
    Double_t pzG=0;
    Double_t zG=0;
    Float_t umax=(10-z0)/(z1-z0);
    for (u=0; u<=umax; u=u+0.1) {
      Double_t x=(1-u)*x0Test+x1Test*u;
      Double_t y=(1-u)*y0Test+y1Test*u;
      Double_t z=(1-u)*z0Test+z1Test*u;
      // ga->SetPoint(ga->GetN(),x,y,z);
      pzG=zG;
      Float_t hillH=h2->GetBinContent(h2->FindBin(x,y));
      zG=z-hillH;
      if (pzG<0 && zG>0) { ok++; printf ("%f %f z=%f \n",zG,pzG,z); hR->Fill(sqrt(hillH*hillH+x*x+y*y));
	gs[i]->SetPoint(1,x,y,z);
	gs[i]->Draw("line:same");
	break;}
  // printf ("\tu=%f (%f %f %f) \t zG=%f \t surface=%f \n",u,x,y,z,zG,z-zG);
  }
  if (u==1) {printf ("nope \n"); nok++;}
  }
 }
printf ("ok=%d nok=%d \n",ok,nok);

c4->cd(2);
hR->Draw();






// /snap/bin/root


}
*/
