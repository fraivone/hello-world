#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <TLegend.h>
#include <TLatex.h>

#include <TApplication.h>
#include <TCanvas.h>
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
#include "TrackHeed.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "ViewFEMesh.hh"

using namespace Garfield;


int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();
  gROOT->SetBatch(kTRUE);



  

  // Dimensions of the GEM
  double pitch = 0.;
  double kapton = 0.;
  double metal = 0.;
  double outdia = 0.;
  double middia = 0.;
  double driftwidth = 0.;
  double induct =0.;
  double rim = 0.;

  //Voltage of the gem
  double drft =0.;
  double induction =0.;
  double upmetal=0.;
  double lowmetal=0.;


  // Histograms
  int nBinsGain = 100;
  double gmin =   0.;
  double gmax = 100.;
  char oldname[1024];
  char newname[1024];
  char holdname[1024];
  gStyle->SetOptTitle(1);
  TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons",
                              nBinsGain, gmin, gmax);
  TH1F* hIons = new TH1F("hIons", "Number of ions",
                         nBinsGain, gmin, gmax);

  int nBinsChrg = 100;
  TH1F* hChrgE = new TH1F("hChrgE", "Position of the electrons on plastic",
                          nBinsChrg, -0.5e4 * kapton, 0.5e4 * kapton);
  TH1F* hChrgI = new TH1F("hChrgI", "Position of the ions on plastic",
                          nBinsChrg, -0.5e4 * kapton, 0.5e4 * kapton);


  const char* PATH = ".";
  DIR *dir = opendir(PATH);
  struct dirent *entry = readdir(dir);
  if(!getcwd(oldname,sizeof(oldname)))

    {
      std::cout<<"Error while reading working director. Script will end right now\n";
      return 0; 
    }
  

  std::string a;
  std::string b;
  std::stringstream ss;
  const char* geometry[12];
  geometry[0] = "pitch =";
  geometry[1] = "kapton =";
  geometry[2] = "metal";
  geometry[3] = "outdia =";
  geometry[4] = "middia =";
  geometry[5] = "drift =";
  geometry[6] = "induct =";
  geometry[7] = "rim = ";
  geometry[8] = "ASEL, S, LOC, Z, drift";
  geometry[9] = "ASEL, S, LOC, Z, induct";
  geometry[10] = "! lower metal";
  geometry[11] = "! upper metal";

  
  while (entry != NULL)
    {
      if(!getcwd(newname,sizeof(newname)))
	{
	  std::cout<<"Error while reading working director. Script will end right now\n";
	  return 0; 
	}
      if (entry->d_type == DT_DIR)
	{
	  a += entry->d_name;
	  if(a.find("d-")!=std::string::npos)
	    {
	      a.clear();
	      a = oldname;
	      a += "/";
	      a += entry->d_name;
	      a += "/";

	      int t=0;
	      std::ifstream f(a+"ansys_gem_script");
	      while(getline(f,b))
		{
		  if(t>=12) break;
		  if(b.find(geometry[t])!=std::string::npos)
		    {
		      ss<<b;
		      ss.ignore(40,'=');
		      switch (t)
			{
			case 0: 
			  ss>>pitch;
			  pitch = pitch/10;
			  break;
     
			case 1: 
			  ss>>kapton;
			  kapton = kapton/10;
			  break;

			case 2:
			  ss>>metal;
			  outdia = metal/10;
			  break;

		    
			case 3:
			  ss>>outdia;
			  outdia = outdia/10;
			  break;
		    
			case 4:
			  ss>>middia;
			  middia = middia/10;
			  break;

			case 5:
			  ss>>driftwidth;
			  driftwidth = driftwidth/10;
			  break;
		    
			case 6:
			  ss>>induct;
			  induct = abs(induct);
			  induct = induct/10;
			  break;			 
		    
			case 7:
			  ss>>rim;
			  rim = rim/10;
			  break;


		    
			case 8:
			  if(getline(f,b))
			    {
			      ss.str("");
			      ss.clear();
			      ss<<b;
			      ss.ignore(16, 'z');
			      ss>>drft;
			      std::cout<<"Ho letto drift: "<<drft<<"\n";			
			    }
			  else
			    std::cout<<"Errore lettura linea\n";
			  break;

		    
			case 9:
			  if(getline(f,b))
			    {
			      ss.str("");
			      ss.clear();
			      ss<<b;
			      ss.ignore(16, 'z');
			      ss>>induction;
			      std::cout<<"Ho letto induction: "<<induction<<"\n";
			    }
			  else
			    std::cout<<"Errore lettura linea\n";
			  break;

		    
			case 10:
			  if(getline(f,b))
			    {
			      ss.str("");
			      ss.clear();
			      ss<<b;
			      ss.ignore(16, 'z');
			      ss>>lowmetal;
			      std::cout<<"Ho letto lowmetal: "<<lowmetal<<"\n";			
			    }
			  else
			    std::cout<<"Errore lettura linea\n";
			  break;
		    
			case 11:
			  if(getline(f,b))
			    {
			      ss.str("");
			      ss.clear();
			      ss<<b;
			      ss.ignore(16, 'z');
			      ss>>upmetal;
			      std::cout<<"Ho letto upmetal: "<<upmetal<<"\n";			
			    }
			  else
			    std::cout<<"Errore lettura linea\n";
			  break;
		    
		    
		    
		    
			default: 
			  break;
			}

		      ss.str("");
		      ss.clear();
		      t++;
		      continue;
		    }
		}
	
	








	      // Load the field map.
	      ComponentAnsys123* fm = new ComponentAnsys123();
	      std::string efile = a+"ELIST.lis";
	      std::string nfile = a+"NLIST.lis";
	      std::string mfile = a+"MPLIST.lis";
	      std::string sfile = a+"PRNSOL.lis";
	      fm->Initialise(efile, nfile, mfile, sfile, "mm");
	      fm->EnableMirrorPeriodicityX();
	      fm->EnableMirrorPeriodicityY();
	      fm->PrintRange();


	      //Creo questo file con i risultati della simulazione
	      strcat (newname, "/");
	      strcat (newname, entry->d_name);
	      strcat (newname, "/log.txt");
	      FILE *ptr=fopen(newname,"w");
	      fprintf(ptr,"Il campo elettrico è stato simulato da ANSYS in base a questi valori\n\n");
	      fprintf(ptr,"-Drift plane\t\t%.1f V\n",drft);
	      fprintf(ptr,"-Upper metal\t\t%.0f V\n",upmetal);
	      fprintf(ptr,"-Lower metal\t\t %.0fV\n",lowmetal);
	      fprintf(ptr,"-Induction plane\t %.0fV\n",induction);
	      fprintf(ptr,"__________________________________________________________\n\n");
	      fprintf(ptr,"-Pitch\t\t\t%.0f micron\n",pitch*10000);
	      fprintf(ptr,"-Kapton\t\t\t%.0f micron\n",kapton*10000);
	      fprintf(ptr,"-Metal\t\t\t%.0f micron\n",metal*10000);
	      fprintf(ptr,"-Outdiam\t\t%.0f micron\n",outdia*10000);
	      fprintf(ptr,"-Middiam\t\t%.0f micron\n",middia*10000);
	      fprintf(ptr,"-Driftwidth\t\t%.0f micron\n",driftwidth*10000);
	      fprintf(ptr,"-Inductwidth\t\t%.0f micron\n",induct*10000);
	      fprintf(ptr,"-Rim\t\t\t%.0f micron\n",rim*10000);



	      //Se vuoi plottare il contourplot del campo
	      const bool plotField = false;
	      if (plotField) {
		ViewField* fieldView = new ViewField();
		fieldView->SetComponent(fm);
		fieldView->SetPlane(0., -1., 0., 0., 0., 0.);
		fieldView->SetArea(-pitch / 2., -0.02, pitch / 2., 0.02);
		fieldView->SetVoltageRange(-160., 160.);
		TCanvas* cF = new TCanvas();
		fieldView->SetCanvas(cF);
		fieldView->PlotContour();
	      }

	      // Setup the gas.
	      MediumMagboltz* gas = new MediumMagboltz();
	      gas->SetComposition("ar", 70., "co2", 30.);
	      gas->SetTemperature(293.15);
	      gas->SetPressure(760.);
	      gas->EnableDebugging();
	      gas->Initialise();
	      gas->DisableDebugging();


	      // Set the Penning transfer efficiency.
	      const double rPenning = 0.57;
	      const double lambdaPenning = 0.;
	      gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");



	      // Load the ion mobilities.
	      gas->LoadIonMobility("/home/francesco/Garfield/Data/IonMobility_Ar+_Ar.txt");

	      // Associate the gas with the corresponding field map material.
	      const int nMaterials = fm->GetNumberOfMaterials();
	      for (int i = 0; i < nMaterials; ++i) {
		const double eps = fm->GetPermittivity(i);
		if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);
	      }
	      //Stampa le informazioni del materiale
	      fm->PrintMaterials();

	      // Create the sensor.
	      Sensor* sensor = new Sensor();
	      sensor->AddComponent(fm);
	      sensor->SetArea(-10 * pitch, -10 * pitch, -0.1,
			      10 * pitch,  10 * pitch,  0.1);



	      //aval è per gli elettroni
	      AvalancheMicroscopic* aval = new AvalancheMicroscopic();
	      aval->SetSensor(sensor);
	      //drift per gli ioni
	      AvalancheMC* drift = new AvalancheMC();
	      drift->SetSensor(sensor);
	      drift->SetDistanceSteps(2.e-4);

	      const bool plotDrift = true;
	      ViewDrift* driftView = new ViewDrift();
	      if (plotDrift) {
		driftView->SetArea(-10 * pitch, -10 * pitch, -0.1,
				   10 * pitch,  10 * pitch,  0.1);
		// Plot every 10 collisions (in microscopic tracking).
		aval->SetCollisionSteps(1);
		aval->EnablePlotting(driftView);
		drift->EnablePlotting(driftView);
	      }

	      //Variabili per il log finale
	      double sumIonsTotal = 0.;
	      double sumIonsDrift = 0.;
	      double sumIonsPlastic = 0.;
	      double sumElectronsTotal = 0.;
	      double sumElectronsPlastic = 0.;
	      double sumElectronsUpperMetal = 0.;
	      double sumElectronsLowerMetal = 0.;
	      double sumElectronsTransfer = 0.;
	      double sumElectronsOther = 0.;




	      //Inizio della simulazione vera e propria
	      const int nEvents = 1;
	      const bool debug = true;
	      // Track class
	      TrackHeed* track = new TrackHeed();
	      track->SetSensor(sensor);
	      track->SetParticle("mu");
	      track->SetMomentum(1.e9);
	      track->DisableDeltaElectronTransport();

	      for (int i = nEvents; i--;) {
		if (debug || i % 10 == 0) std::cout << i << "/" << nEvents << "\n";
		// Randomize the initial position.
		//La direzione giusta è verso il basso per come ho simulato la gem e il gas in ansys
		//Inoltre sensor
		track->NewTrack(-pitch +RndmUniform()*pitch,-sqrt(3)*pitch/2 + sqrt(3) * RndmUniform()*pitch/2, 0.0249, 0., 0., 0.,-1);
		double x0 = 0.;
		double y0 = 0.;
		double z0 = 0.;
		double t0 = 0.;
		double e0 = 0.;
		double dx0 = 0.;
		double dy0 = 0.;
		double dz0 = 0.;
		int ne = 0, ni = 0;
		//Variables to save the cluster position
		double xc = 0., yc = 0., zc = 0., tc = 0.;
		// Number of electrons produced in a collision
		int nc = 0;
		// Energy loss in a collision
		double ec = 0.;
		// Dummy variable (not used at present)
		double extra = 0.;

		double xe1, ye1, ze1, te1, e1;
		double xe2, ye2, ze2, te2, e2;
		double xi1, yi1, zi1, ti1;
		double xi2, yi2, zi2, ti2;
		int status;
		//Loop su tutti i cluster prodotti dalla particella
		while(track->GetCluster(xc,yc,zc,tc,nc,ec,extra)){
		  //Loop su tutti gli elettroni nel cluster in questione
		  std::cout<<"\nSto trattando il cluster prodotto al tempo: "<<tc<<"\n__________________________________________________________\n";
		  for(int k = nc; k--;)
		    { track->GetElectron(k,x0,y0,z0,t0,ec,dx0,dy0,dz0);
		      std::cout<<"L'elettrone "<<k+1<<" aveva z: "<<z0<<", Energia :"<<std::setprecision(2)<<ec<<"\n";

		      //Lancio la moltiplicazione del singolo elettrone
		      //Il campo elettrio e gli scattering determineranno il numero e la traitteoria di tutti i secondari
		      aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
		      aval->GetAvalancheSize(ne, ni);
		      const int np = aval->GetNumberOfElectronEndpoints();
		      //Riempio questi istogrammi con il numero di elettroni/ioni prodotti per questo primario
		      hElectrons->Fill(ne);
		      hIons->Fill(ni);

		      //Adesso GetElectronEndPoint mi restituisce la posizione iniziale e finale di ogni elettrone originato dalla valanga del primario
		      //Questo ciclo for looppa su tutti i secondari della valanga
		      for (int j = np; j--;) {
			aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1,
						  xe2, ye2, ze2, te2, e2, status);
			sumElectronsTotal += 1.;
			if (ze2 > -kapton / 2. && ze2 < kapton / 2.) {
			  hChrgE->Fill(ze2 * 1.e4);
			  sumElectronsPlastic += 1.;
			} else if (ze2 >= kapton / 2. && ze2 <= kapton  / 2. + metal) {
			  sumElectronsUpperMetal += 1.;
			} else if (ze2 <= -kapton / 2. && ze2 >= -kapton / 2. - metal) {
			  sumElectronsLowerMetal += 1.;
			} else if (ze2 < -kapton / 2. - metal) {
			  sumElectronsTransfer += 1.;
			} else {
			  sumElectronsOther += 1.;
			}
			//Faccio la stessa cosa per gli ioni, stavolta è un avalancheMC
			drift->DriftIon(xe1, ye1, ze1, te1);
			drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1,
					      xi2, yi2, zi2, ti2, status);
			if (zi1 < 0.01) {
			  sumIonsTotal += 1.;
			  if (zi2 > 0.01) sumIonsDrift += 1.;
			}
			if (zi2 > -kapton / 2. && zi2 < kapton / 2.) {
			  hChrgI->Fill(zi2 * 1.e4);
			  sumIonsPlastic += 1.;
			}
		      }


		      //Fine loop su cluster in questione
		    }

		  //Fine loop sui cluster prodotti dalla particella
		  std::cout<<"**********************************************************\n**********************************************************\n\n";
		}

		//Qui finisce il loop sugli eventi
	      }



	      //Stampo alcune statistiche a video
	      std::cout<<"Elettroni totali prodotti "<< sumElectronsTotal <<"\t Elettroni totali persi: " << sumElectronsPlastic+sumElectronsUpperMetal+sumElectronsLowerMetal+sumElectronsOther<< "\n\n";
	      fprintf(ptr,"__________________________________________________________\n\n");
	      fprintf(ptr,"Elettroni totali prodotti %.0f\t Elettroni totali persi: %.0f\n\n",sumElectronsTotal,sumElectronsPlastic+sumElectronsUpperMetal+sumElectronsLowerMetal+sumElectronsOther);
	      double fFeedback = 0.;
	      if (sumIonsTotal > 0.) fFeedback = sumIonsDrift / sumIonsTotal;
	      std::cout << "Fraction of ions drifting back: " << fFeedback << "\n";
	      fprintf(ptr, "Fraction of ions drifting back: %.2f %c\n",fFeedback,'%');

	      //Stampo il numero medio di ioni e elettroni prodotto per primario
	      const double neMean = hElectrons->GetMean();
	      std::cout << "Mean number of electrons: " << neMean << "\n";
	      const double niMean = hIons->GetMean();
	      std::cout << "Mean number of ions: " << niMean << "\n";

	      std::cout << "Mean number of electrons on plastic: "
			<< sumElectronsPlastic / nEvents << "\n";
	      std::cout << "Mean number of ions on plastic: "
			<< sumIonsPlastic / nEvents << "\n";

	      std::cout<< "Electron endpoints:\n";

	      //Stampo le varie porzioni di elettroni
	      const double fUpperMetal = sumElectronsUpperMetal / sumElectronsTotal;
	      const double fPlastic = sumElectronsPlastic / sumElectronsTotal;
	      const double fLowerMetal = sumElectronsLowerMetal / sumElectronsTotal;
	      const double fTransfer = sumElectronsTransfer / sumElectronsTotal;
	      const double fOther = sumElectronsOther / sumElectronsTotal;
	      std::cout << "    upper metal: " << fUpperMetal * 100. << "%\n";
	      fprintf(ptr,"    upper metal: %.2f %c\n",fUpperMetal*100,'%');
	      std::cout << "    plastic:     " << fPlastic * 100. << "%\n";
	      fprintf(ptr,"    plastic:     %.2f %c\n",fPlastic*100,'%');
	      std::cout << "    lower metal: " << fLowerMetal * 100. << "%\n";
	      fprintf(ptr,"    lower metal: %.2f %c\n",fLowerMetal*100,'%');
	      std::cout << "    transfer:    " << fTransfer * 100. << "%\n";
	      fprintf(ptr,"    transfer:    %.2f %c\n",fTransfer*100,'%');
	      std::cout << "    other:       " << fOther * 100. << "%\n";
	      fprintf(ptr,"    other:        %.2f %c\n",fOther*100,'%');


	      //PLOT PART


	      //Questo plot è pesante
	      //lo lascio tutto commentato se ti serve usalo

	      //const bool plotGeo = false;
	      // if (plotGeo && plotDrift) {
	      // 	TCanvas* cD = new TCanvas("cD","GEOCanvas",800,600);
	      // 	// Build the geometry in Root.
	      // 	TGeoManager* geoman = new TGeoManager("world", "geometry");
	      // 	TGeoMaterial* matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
	      // 	TGeoMedium* medVacuum = new TGeoMedium("Vacuum", 1, matVacuum);
	      // 	TGeoMaterial* matKapton = new TGeoMaterial("Kapton", 12, 6, 1.42);
	      // 	TGeoMedium* medKapton = new TGeoMedium("Kapton", 2, matKapton);
	      // 	TGeoMaterial* matCopper = new TGeoMaterial("Copper", 63, 29, 8.94);
	      // 	TGeoMedium* medCopper = new TGeoMedium("Copper", 3, matCopper);
	      // 	TGeoVolume* volTop = geoman->MakeBox("TOP",
	      // 					     medVacuum, pitch, pitch, 0.02);
	      // 	volTop->SetVisibility(0);
	      // 	TGeoBBox* shpKapton = new TGeoBBox("K", pitch / 2.,
	      // 					   pitch / 2.,
	      // 					   kapton / 2.);
	      // 	TGeoPcon* shpHole = new TGeoPcon("H", 0., 360., 3);
	      // 	shpHole->DefineSection(0, -kapton / 2., 0., outdia / 2.);
	      // 	shpHole->DefineSection(1,           0., 0., middia / 2.);
	      // 	shpHole->DefineSection(2,  kapton / 2., 0., outdia / 2.);

	      // 	TGeoCompositeShape* shpGem = new TGeoCompositeShape("G", "K - H");
	      // 	TGeoVolume* volKapton = new TGeoVolume("Kapton", shpGem, medKapton);
	      // 	volKapton->SetLineColor(kGreen);
	      // 	volKapton->SetTransparency(50);

	      // 	TGeoBBox* shpMetal = new TGeoBBox("M", pitch / 2.,
	      // 					  pitch / 2.,
	      // 					  metal / 2.);
	      // 	TGeoTube* shpTube = new TGeoTube("T", 0., outdia / 2., metal / 2.);
	      // 	TGeoCompositeShape* shpElectrode = new TGeoCompositeShape("E", "M - T");
	      // 	TGeoVolume* volElectrode = new TGeoVolume("Electrode",
	      // 						  shpElectrode, medCopper);
	      // 	volElectrode->SetLineColor(kBlue);
	      // 	volElectrode->SetTransparency(50);

	      // 	TGeoVolumeAssembly* volGem = new TGeoVolumeAssembly("Gem");
	      // 	const double shift =  0.5 * (metal + kapton);
	      // 	volGem->AddNode(volKapton, 1);
	      // 	volGem->AddNode(volElectrode, 2, new TGeoTranslation(0., 0.,  shift));
	      // 	volGem->AddNode(volElectrode, 3, new TGeoTranslation(0., 0., -shift));

	      // 	volTop->AddNode(volGem, 1);
	      // 	volTop->AddNode(volGem, 2, new TGeoTranslation(-pitch, 0., 0.));
	      // 	volTop->AddNode(volGem, 3, new TGeoTranslation(+pitch, 0., 0.));
	      // 	volTop->AddNode(volGem, 4,
	      // 			new TGeoTranslation(-pitch / 2., sqrt(3) * pitch / 2., 0.));
	      // 	volTop->AddNode(volGem, 5,
	      // 			new TGeoTranslation(+pitch / 2., sqrt(3) * pitch / 2., 0.));
	      // 	volTop->AddNode(volGem, 6,
	      // 			new TGeoTranslation(-pitch / 2., -sqrt(3) * pitch / 2., 0.));
	      // 	volTop->AddNode(volGem, 7,
	      // 			new TGeoTranslation(+pitch / 2., -sqrt(3) * pitch / 2., 0.));

	      // 	geoman->SetVerboseLevel(0);
	      // 	cD->cd();
	      // 	geoman->SetTopVolume(volTop);
	      // 	geoman->CloseGeometry();
	      // 	geoman->CheckOverlaps(0.1e-4);
	      // 	geoman->SetNmeshPoints(100000);
	      // 	geoman->GetTopVolume()->Draw("ogl");
	      // 	if (1) {
	      // 	  driftView->SetCanvas(cD);
	      // 	  driftView->Plot();
	      // 	}
	      // 	cD->Close();

	      // }



	      //Plot degli istogrammi con carica e luogo di arrivo
	      const bool plotHistogram = false;
	      if (plotHistogram) {
		TCanvas* cH = new TCanvas("cH", "Histograms", 800, 700);
		cH->Divide(2, 2);
		cH->cd(1);
		hElectrons->Draw();
		cH->cd(2);
		hIons->Draw();
		cH->cd(3);
		hChrgE->Draw();
		cH->cd(4);
		hChrgI->Draw();
	      }


	      //In questa  parte rinomino la cartella con i dati della geometria e voltaggio

	      
	      int result;
	      // sprintf(newname,"/d%.0f md%.0f k%.0f p%.0f",drft*10000,*10000,kapton*10000,pitch*10000);
	      // if(getcwd(oldname,sizeof(oldname))==NULL)
	      //   perror("getcwd() error");
	      // else
	      //   { result = strlen(oldname);
	      //     oldname[result] = '/';
	      //     oldname[result+1]= '\0';
	      //     result= std::rename( oldname , newname );
	      //     if ( result == 0 )
	      // 	puts ( "File successfully renamed" );
	      //     else
	      // 	perror( "Error renaming file" );
	      //   }

	      strcpy(newname , oldname);
	      sprintf(holdname, "/d%.0f um%.0f lm%.0f i%.0f",drft,upmetal,lowmetal,induction);
	      strcat (newname, holdname);
	      puts(newname);
	      strcpy(holdname , oldname);
	      strcat(holdname , "/");
	      strcat(holdname , entry->d_name);
	      
	      result = strlen(holdname);
	      holdname[result]= '\0';
	      puts(holdname);
	      result= std::rename(holdname  ,newname );
	      if ( result == 0 )
		puts ( "File successfully renamed" );
	      else
	      	perror( "Error renaming file" );
	    
	      

	      
	      fclose(ptr);

	      memset(newname, 0, sizeof(newname));
	      if(!getcwd(newname,sizeof(newname)))
		{
		  std::cout<<"Error while reading working director. Script will end right now\n";
		  return 0; 		  
		}
	      strcat (newname, "/");
	      strcat (newname, entry->d_name);
	      strcat (newname, "/Front 20 mu 1 GeV.pdf");


	      //PLOT MESH UTILE, BIDIMENSIONALE
	      TCanvas* c87 = new TCanvas("c87","mesh",800,600);
	      ViewFEMesh* meshView = new ViewFEMesh();
	      meshView->SetCanvas(c87);
	      meshView->SetComponent(fm);
	      meshView->SetPlane(0.,-1, 0.,0.,0.,0.03);
	      meshView->SetFillMesh(true);
	      meshView->SetColor(1, 5);
	      meshView->SetViewDrift(driftView);
	      meshView->SetArea(-0.02, -0.02, -0.07, 0.02, 0.02, 0.07);
	      meshView->Plot();
	      c87->SaveAs(newname);
	      c87->Clear();
	      c87->Close();

	      memset(newname, 0, sizeof(newname));
	      if(!getcwd(newname,sizeof(newname)))
		{
		  std::cout<<"Error while reading working director. Script will end right now\n";
		  return 0;
		}
	      strcat (newname, "/");
	      strcat (newname, entry->d_name);
	      strcat (newname, "/Top 20 mu 1 GeV.pdf");

	      TCanvas* c88 = new TCanvas("c88","mesh2",800,600);
	      ViewFEMesh* meshView2 = new ViewFEMesh();
	      meshView2->SetCanvas(c88);
	      meshView2->SetComponent(fm);
	      meshView2->SetPlane(0,0,-1,0,0,0);
	      meshView2->SetFillMesh(true);
	      meshView2->SetColor(1, 5);
	      meshView2->SetViewDrift(driftView);
	      meshView2->SetArea(-0.02, -0.02, -0.07, 0.02, 0.02, 0.07);
	      meshView2->Plot();
	      c88->SaveAs(newname);
	      c88->Clear();
	      c88->Close();


	    }

	}

      a.clear();
      memset(newname, 0, sizeof(newname));
      entry= readdir(dir);
    }

  closedir(dir);


  std::cout<<"\nBye...";
  app.Run(kTRUE);

}
