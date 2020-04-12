
#include <TTree.h>
#include <TRandom.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <map>
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TFile.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

using namespace std;

int makeAmptroot(string run_id)
{
	// PDG Database
	TDatabasePDG *db = new TDatabasePDG();
	TParticlePDG *p_info = new TParticlePDG();

	// Constants
	const int proton_pid = 2212;
	float p_min = 0.15;
	float pt_min = 0.3;
	float pt_max = 2.5;
	float eta_max = 1.0;

	//input file variables:
	Int_t    evn,tn,nov,npp,npt,nesp,ninesp,nest,ninest,pid; // for ampt.dat
	Float_t    px,py,pz,mass,x,y,z,t,b;                                //for ampt.dat


	//Booked variables in tree:
	const Int_t  mul  = 90000;
	Int_t  event=0, refmult2, refmult3, pmult;
	Float_t  imp, event_plane_ref2, event_plane_ref3; // event variables
	Int_t   pid_array[mul];        //particle variables
	Float_t px_array[mul], py_array[mul], pz_array[mul];


	//------------define a root file and tree :------------------------------
	TFile *fampt = new TFile(("data_" + run_id + ".root").data(),"RECREATE");
	if(!fampt){
		cout << "Output file cannot be opened" << endl;
		exit(0);
	}

	TTree *tr = new TTree("tree","AMPT Proton Data");



	//Define event branches:-------------------------------------------
	tr->Branch("event", &event, "event/I");
	tr->Branch("pmult", &pmult, "pmult/I");
	tr->Branch("refmult2",  &refmult2, "refmult2/I");           // multiplicity = tracks
	tr->Branch("refmult3",  &refmult3, "refmult3/I");
	tr->Branch("event_plane_ref2",  &event_plane_ref2, "event_plane_ref2/F");
	tr->Branch("event_plane_ref3",  &event_plane_ref3, "event_plane_ref3/F");
	tr->Branch("imp",&imp, "imp/F");


	//particle branches:
	tr->Branch("pid", &pid_array, "pid[pmult]/I");
	tr->Branch("px",  &px_array,  "px[pmult]/F");
	tr->Branch("py",  &py_array,  "py[pmult]/F");
	tr->Branch("pz",  &pz_array,  "pz[pmult]/F");

	//**************************************************************************************

	cout<<"making .root file from .dat file..."<<endl;

	ifstream infile;//[300000];
	char *fname   = new char[100];
	sprintf(fname,"ana/ampt.dat");
	infile.open(fname);

	while(infile)
	{//2 event loop // Please see ampt.dat file and npart-xy.dat to comment proper line:
		infile>>evn>>tn>>nov>>b>>npp>>npt>>nesp>>ninesp>>nest>>ninest;

		if(infile.eof())  break;

		event += 1;
		imp = b;

		//****************************ampt.dat particle loop**************
		vector<float> p_px, p_py, p_pz;
		vector<int> p_pid;
		float Qx_ref2 = 0, Qy_ref2 = 0, Qx_ref3 = 0, Qy_ref3 = 0;
		refmult2 = 0; refmult3 = 0;

		for(int j=0;j<nov;j++)                          //particle loop
		{
			infile>>pid>>px>>py>>pz>>mass>>x>>y>>z>>t;

			p_info = db->GetParticle((int)pid);
			if(fabs((int)p_info->Charge()) != 3) continue;

			TVector3 p_mom(px, py, pz);
			if(p_mom.Mag() < p_min) continue;

			double eta = p_mom.PseudoRapidity();
			if(fabs(eta > eta_max)) continue;

			double pt = p_mom.Perp();
			double phi = p_mom.Phi();

			// Refmult and event plane logic

			// ref2
			if(fabs(eta) > 0.5 && fabs(eta) < eta_max) refmult2++;
			if(fabs(eta) < eta_max && pt > 0.2 && pt < 2.0) {
				if(fabs(eta) > 0.5 || fabs(pid) != proton_pid){
					Qx_ref2 += cos(2*phi); Qy_ref2 += sin(2*phi);
				}
			}
			// ref3
			if(fabs(eta) < eta_max && pid != proton_pid) {
				refmult3++;
				if(pt > 0.2 && pt < 2.0) {
					Qx_ref3 += cos(2*phi); Qy_ref3 += sin(2*phi);
				}
			}

			// Proton selection logic
			if(fabs(pid) == proton_pid) {
				if(p_mom.Perp() >= pt_min && p_mom.Perp() <= pt_max) {
					p_px.push_back(px);
					p_py.push_back(py);
					p_pz.push_back(pz);
					p_pid.push_back(pid);
				}
			}
		}

		TVector2 Q_ref2(Qx_ref2, Qy_ref2); TVector2 Q_ref3(Qx_ref3, Qy_ref3);
		event_plane_ref2 = 0.5 * Q_ref2.Phi(); event_plane_ref3 = 0.5 * Q_ref3.Phi();
		pmult = (int)p_pid.size();

		for(int j=0; j<pmult; j++) {      //proton loop
			px_array[j] = p_px[j];
			py_array[j] = p_py[j];
			pz_array[j] = p_pz[j];
			pid_array[j] = p_pid[j];
		}

		//***************************ampt.dat particle loops end**************

		tr->Fill();


	}//2 event loop end


	fampt->cd();
	tr->Write();
	fampt->Close();
	cout<<" Tree written succesfully "<<endl;

	return 0;
}//main end


