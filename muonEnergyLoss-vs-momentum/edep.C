#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TMath.h>

#include "MediumMagboltz.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "TrackHeed.hh"
#include "Plotting.hh"
#include "Random.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {
  double min = 10.;
  double loweryrange = 0.;
  double upperyrange = 0.;
  const int points = 1; //Punti di simulazione dell'impulso per decade
  int nEvents; //Numero di tracce simulate per ogni impulso
  randomEngine.Seed(1456);
  TApplication app("app", &argc, argv);
  //  plottingEngine.SetDefaultStyle();

  // Histograms
  TCanvas* c2 = new TCanvas();
  //  c2->SetLogx();
  c2->SetGrid(1,1);

  TH1::StatOverflows(kTRUE); 
  TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons",
                              50, 0, 200);
  TGraph* h1 = new TGraph(3*points);
  TGraph* h2 = new TGraph(6);
  TH1F* hEdep = new TH1F("hEdep", "Energy Loss",
                         50, 0., 100.);

  // Make a medium
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("ar", 90., "co2", 10.);
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);

  // Detector geometry
  // Gap [cm]
  const double width = 1.;
  SolidBox* box = new SolidBox(width / 2., 0., 0., width / 2., 10., 10.);
  GeometrySimple* geo = new GeometrySimple();
  geo->AddSolid(box, gas);

  // Make a component
  ComponentConstant* comp = new ComponentConstant();
  comp->SetGeometry(geo);
  comp->SetElectricField(10., 0., 0.);

  // Make a sensor
  Sensor* sensor = new Sensor();
  sensor->AddComponent(comp);

  // Track class
  double enloss;
  double momentum;
  TString label="Muon momentum [GeV/c]";
  TrackHeed* track = new TrackHeed();
  track->SetSensor(sensor);
  track->SetParticle("mu");
  for(int k=0; k<5; k++){
    nEvents = 1 + k*40;
    
  for(int j=0; j<3*points; j++){
    momentum = min*pow(10,j/10)*(1+ j%10);   //Momentum in MeV
    //Momentum viene incrementata in step logaritmici, 10 punti per decade
    
  track->SetMomentum(momentum*1.e6);


  track->EnableDebugging();
  track->DisableDeltaElectronTransport();
  for (int i = 0; i < nEvents; ++i) {
    if (i == 0) track->EnableDebugging();
    if (i == 1) track->DisableDebugging();
     
    // Initial position and direction
    enloss = 0;
    double x0 = 0., y0 = 0., z0 = 0., t0 = 0.;
    double dx0 =1., dy0 =0., dz0 = 0.,ee=0.;
    double xe = 0., ye = 0., ze = 0., te = 0.;
    double dxe =0., dye =0., dze = 0.;
    track->NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);
    // Cluster coordinates
    double xc = 0., yc = 0., zc = 0., tc = 0.;
    // Number of electrons produced in a collision
    int nc = 0;
    // Energy loss in a collision
    double ec = 0.;
    // Dummy variable (not used at present)
    double extra = 0.;
    // Total energy loss along the track
    double esum = 0.;
    // Total number of electrons produced along the track
    int nsum = 0;
    // Loop over the clusters.
    while (track->GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
      /*	while(track->GetElectron(ii,xe, ye, ze, te,ee,dxe,dye,dze))
	  {electrontotenergy+=ee;
	    printf("\n%1f",ee);
	    ii++;}*/
      esum += ec;
      nsum += nc;      
    }
          
    enloss += esum*1.e-6;     
    if(j==0)hElectrons->Fill(nsum);
    if(j==0)hEdep->Fill(esum * 1.e-3);
    
  }
      enloss =enloss/nEvents;
    //Fill with momentum [GeV/c] , EnLoss [MeV]
  /*
    h1->SetPoint(j,momentum/1000,enloss);
    if(j==0) {upperyrange= enloss;loweryrange=enloss;}
    if(enloss < loweryrange) loweryrange= enloss;*/
  if(j==0) h2->SetPoint(k,nEvents,enloss);
  }
  }
  // TCanvas* c1 = new TCanvas();
 
  /*  c1->Divide(1,2);
  c1->cd(1);
  hElectrons->GetXaxis()->SetTitle("Pion number of electrons"); 
  hElectrons->Draw();
  c1->cd(2);
  hEdep->GetXaxis()->SetTitle("Pion En loss [keV]");
  hEdep->Draw();  
  hEdep->GetXaxis()->SetRange(0,100);
  c1->SaveAs("diff_histos.root");*/

  /* h1->GetXaxis()->SetTitle(label);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->GetXaxis()->SetTitleOffset(1.3);
  h1->GetYaxis()->SetTitle("Energy Loss [MeV/cm]");
  h1->SetMinimum(loweryrange*0.8);
  h1->SetMaximum(upperyrange*1.2);
  h1->GetXaxis()->SetLimits(.009,10.09);
  c2->cd();
  h1->SetLineColor(2);
  h1->SetLineWidth(4);
  h1->SetMarkerColor(4);
  h1->SetMarkerSize(1.5);
  h1->SetMarkerStyle(21);
  TString  title = "Energy loss averaged over "+std::to_string(nEvents)+" tracks";
  h1->SetTitle(title);

  h1->Draw("APL");*/
  h2->SetLineColor(2);
  h2->SetLineWidth(4);
  h2->SetMarkerColor(4);
  h2->SetMarkerSize(1.5);
  h2->SetMarkerStyle(21);
  h2->Draw("APL");
  c2->Update();
  c2->Modified();
  TString g ="EnLoss-vs-tracks";
  c2->SaveAs(g+".root");
  c2->SaveAs(g+".pdf");
  c2->Close();
  return 0;
  app.Run(kTRUE); 

  }
