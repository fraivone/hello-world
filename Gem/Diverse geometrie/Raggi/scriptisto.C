#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <TLegend.h>
#include <TLatex.h>
#include <sys/types.h>
#include <dirent.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <THStack.h>
#include <sstream>
#include <TH1F.h>
#include <TGeoManager.h>
#include <TPie.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TText.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>

#include "ComponentAnsys123.hh"
#include "ViewField.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "ViewFEMesh.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();
  
  TH2F* totalelectrons = new TH2F("totalelectrons","",8,20,100,8,20,60);
  TH2F* lostelectrons =  new TH2F("lostelectrons","",8,20,100,8,20,60);
  TH2F* fractionup =  new TH2F("fractionup","",8,20,100,8,20,60);
  TH2F* fractionplastic  =  new TH2F("fractionplastic","",8,20,100,8,20,60);
  TH2F* fractionlow =  new TH2F("fractionlow","",8,20,100,10,10,60);
  

  const bool debug = true;

  // Dimensions of the GEM
  const double pitch = 0.014;
  const double kapton = 50.e-4;
  const double metal = 5.e-4;
  const double outdia = 70.e-4;
  const double middia = 50.e-4;
  // Voltages on GEM
  const double drft =-200;
  const double induction =700;
  const double upmetal=-100;
  const double lowmetal=200;
  
  int od=0;
  int md=0;  
  int lost;
  int prod;
  double hold;
  double lower;
  double upper;
  std::stringstream ss;
  const char* PATH = ".";
  char buff[FILENAME_MAX];
  char buff2[FILENAME_MAX];
  getcwd(buff,sizeof(buff));
  std::string a(buff);
  std::string help;
  std::string search = "/od"; //stringa da <cercare 
    DIR *dir = opendir(PATH);

    struct dirent *entry = readdir(dir);

    while (entry != NULL)
      {
	getcwd(buff,sizeof(buff));
	a = buff;
	a +="/";
        if (entry->d_type == DT_DIR)	  
	  { a += entry->d_name;
	    if(a.find(search)!=std::string::npos)
	      {
		od=0;
		md=0;
		sscanf(entry->d_name,"%s %s",buff,buff2);
		od = std::atoi(&buff[2]);
		md = std::atoi(&buff2[2]);
		a += "/log.txt";		
		std::ifstream f ( a );
		a = "";
		while(getline( f, a))
		  {
	    	    if(a.find("prodotti")!=std::string::npos)
		      {
			lost = 0;
			prod=0;
			ss<<a;
			//std::cout<<"Sto leggendo:\t"<<ss.rdbuf()<<"\n";
			ss.ignore(26,'\n');
			ss>>prod;
			ss.ignore(25,'z');
			ss>>lost;
			totalelectrons->Fill(od,md,prod-lost);
			lostelectrons->Fill(od,md,lost);
			ss.str("");
			ss.clear();
			continue;
			
		      }

		    
		    if(a.find("upper metal:")!=std::string::npos)
		      {
			upper = 0.;
			ss<<a;
			ss.ignore(16,'z');
			ss>>upper>>buff;
			fractionup->Fill(od,md,upper*prod/100);
			// std::cout << "A VGEM: "<<voltages<<" trovo lowmetal: "<< hold <<"\n";
			ss.str("");
			ss.clear();
			continue;
			
		      }
		    
		    if(a.find("plastic:")!=std::string::npos)
		      {
			hold=0;
			ss<<a;
			ss.ignore(16,'z');
			// std::cout<<voltages<<ss.rdbuf()<<"\n";
			ss>>hold>>buff;
			fractionplastic->Fill(od,md,(hold)*prod/100);			
			// std::cout<<voltages<<" trovo plastic: "<< hold;
			ss.str("");
			ss.clear();
			continue;
			
		      }
		    if(a.find("lower")!=std::string::npos)
		      {
			lower = 0.;
			ss<<a;
			ss.ignore(16,'z');
			// std::cout<<voltages<<ss.rdbuf()<<"\n";
			ss>>lower>>buff;
			std::cout<<"Coordinate: "<<od<<md<<"Prodotti: "<< lost<<"\tpercentuale di lost: "<< lower*prod/100<<"\n";
			fractionlow->Fill(od,md,lower*prod/100);
			ss.str("");
			ss.clear();
			continue;			
		      }
		    
		  }
	      }
	  }
	

        entry = readdir(dir);
    }
    closedir(dir);
    
    lostelectrons->SetFillColor(21);
    lostelectrons->SetBarWidth(1);
    lostelectrons->SetBarOffset(0.);
    lostelectrons->SetXTitle("");
    lostelectrons->SetYTitle("");
    lostelectrons->SetStats(kFALSE);
    fractionlow->SetFillColor(30);
    fractionlow->SetBarWidth(1);
    fractionlow->SetBarOffset(0.);
    fractionlow->SetXTitle("");
    fractionlow->SetYTitle("");
    fractionlow->SetStats(kFALSE);

    fractionplastic->SetFillColor(46);
    fractionplastic->SetBarWidth(1);
    fractionplastic->SetBarOffset(0.);
    fractionplastic->SetXTitle("");
    fractionplastic->SetYTitle("");
    fractionplastic->SetStats(kFALSE);    
    

    
    fractionup->SetFillColor(41);
    fractionup->SetBarWidth(1);
    fractionup->SetBarOffset(0.);
    fractionup->SetXTitle("");
    fractionup->SetYTitle("");
    fractionup->SetStats(kFALSE);   
    totalelectrons->SetFillColor(4);
    totalelectrons->SetBarWidth(1);
    totalelectrons->SetBarOffset(0.);
    totalelectrons->SetXTitle("Outer diameter [#mum]");
    totalelectrons->SetYTitle("Mid diameter [#mum]");
    totalelectrons->SetStats(kFALSE);
    TCanvas* c4 = new TCanvas("c4", "totalelectrons",1600, 800);
    //totalelectrons->Draw("lego4");
    //lostelectrons->Draw("SAME hist b");
    // fractionlow->Draw("SAME lego2");
    // fractionplastic->Draw("SAME lego");
    // fractionup->Draw("SAME lego2");
    THStack* stack= new THStack("stack", "Stacked histograms");
    stack->Add(fractionplastic);
    stack->Add(fractionup);
    stack->Add(fractionlow);
    stack->Add(totalelectrons);
    stack->Draw("lego4");
    stack->GetYaxis()->SetTitle("Middiam [#mum]");
    stack->GetYaxis()->SetTitleOffset(1.8);
    stack->GetXaxis()->SetTitle("Outdiam [#mum]");
    stack->GetXaxis()->SetTitleOffset(1.8);
    gPad->Modified(); 
    
    auto legend = new TLegend(0.15,0.75,0.43,0.95);
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
    legend->AddEntry(totalelectrons,"Transfered electrons","f");
    legend->AddEntry(fractionlow,"Lost on lowmetal","f");
    legend->AddEntry(fractionup,"Lost on upmetal","f");
    legend->AddEntry(fractionplastic,"Lost on plastic","f");
    
    legend->SetFillColorAlpha(kWhite,0.5);
    legend->Draw();
    char c[200];
    char b[200];
    sprintf(c,"GEM dimension [#mum ]   -   kapton = %.0f,  hole diam = %.0f-%.0f,  pitch = %.0f",kapton*10000, middia*10000,outdia*10000,pitch*10000);
    // TLatex *text = new TLatex(50,50, 30, c);
    // sprintf(b,"#splitline{GEM field: %.0f kV/cm}{Drift field: %.0f kV/cm}",(lowmetal-upmetal)/5,(upmetal-drft)/100);   //Fattore 100 ingloba spazio di drift e conversione V->kV
    // TLatex *text2 = new TLatex(50,50 ,50, b);
    // text->SetTextSize(0.032);
    // text2->SetTextSize(0.039);
    // text->Draw();
    // text2->Draw();
    c4->Modified();
    c4->Update();
    

    
    app.Run(kTRUE);
    
}
    
