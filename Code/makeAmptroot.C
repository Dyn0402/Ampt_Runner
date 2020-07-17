
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

int makeAmptroot_all(string run_id)
{
	// PDG Database
	TDatabasePDG *db = new TDatabasePDG();
	TParticlePDG *p_info = new TParticlePDG();

	// Constants
	const int proton_pid = 2212;
	float p_min = 0.1;  // STAR minimum momentum acceptance
	float pt_min = 0.01;  // Issues arise in eta if pt is zero
	float qvec_pt_min = 0.2;
	float qvec_pt_max = 2.0;
	float ref_eta_max = 0.5;
	float ref2_eta_min = 0.5;
	float ref2_eta_max = 1.0;
	float ref3_eta_max = 1.0;
	float eta_max = 1.0;

	// Input file variables
	Int_t evn, tn, nov, npp, npt, nesp, ninesp, nest, ninest, pid;  // for ampt.dat event
	Float_t px, py, pz, mass, x, y, z, t;  // for ampt.dat track


	// Output tree variables
	const Int_t mul = 90000;
	Int_t event=0, refmult, refmult2, refmult3, pmult;  // event variables
	Float_t imp, qx, qy;  // event variables
	Int_t pid_array[mul];  // particle variables
	Float_t px_array[mul], py_array[mul], pz_array[mul];  // particle variables


	//------------define a root file and tree :------------------------------
	TFile *fampt = new TFile(("data_" + run_id + ".root").data(),"RECREATE");
	if(!fampt){
		cout << "Output file cannot be opened" << endl;
		exit(0);
	}

	TTree *tr = new TTree("tree","AMPT Data");


	//Define event branches:-------------------------------------------
	tr->Branch("event",    &event,     "event/I");
	tr->Branch("pmult",    &pmult,     "pmult/I");
	tr->Branch("refmult",  &refmult,   "refmult/I");
	tr->Branch("refmult2", &refmult2,  "refmult2/I");
	tr->Branch("refmult3", &refmult3,  "refmult3/I");
	tr->Branch("qx",       &qx,        "qx/F");
	tr->Branch("qy",       &qy,        "qy/F");
	tr->Branch("imp",      &imp,       "imp/F");

	tr->Branch("npp",      &npp,       "npp/I");
	tr->Branch("npt",      &npt,       "npt/I");
	tr->Branch("nesp",     &nesp,      "nesp/I");
	tr->Branch("ninesp",   &ninesp,    "ninesp/I");
	tr->Branch("nest",     &nest,      "nest/I");
	tr->Branch("ninest",   &ninest,    "ninest/I");


	//particle branches:
	tr->Branch("pid",      &pid_array, "pid[pmult]/I");
	tr->Branch("px",       &px_array,  "px[pmult]/F");
	tr->Branch("py",       &py_array,  "py[pmult]/F");
	tr->Branch("pz",       &pz_array,  "pz[pmult]/F");

	//**************************************************************************************

	cout<<"making .root file from .dat file..."<<endl;

	ifstream infile;
	char *fname   = new char[100];
	sprintf(fname,"ana/ampt.dat");
	infile.open(fname);

	while(infile)
	{// event loop
		infile>>evn>>tn>>nov>>imp>>npp>>npt>>nesp>>ninesp>>nest>>ninest;

		if(infile.eof())  break;

		event += 1;

		//****************************ampt.dat particle loop**************
		qx = 0.; qy = 0.;
		refmult = 0; refmult2 = 0; refmult3 = 0; pmult = 0;

		for(int j=0;j<nov;j++)                          //particle loop
		{
			infile>>pid>>px>>py>>pz>>mass>>x>>y>>z>>t;

			px_array[j] = px;
			py_array[j] = py;
			pz_array[j] = pz;
			pid_array[j] = pid;
			pmult++;

			p_info = db->GetParticle((int)pid);
			if(!p_info) { cout << "pid: " << pid << " not in TDatabasePDG" << endl; continue; }
			else if(p_info->Mass() != (double) mass) { cout << "pid: " << pid << " with mass " << mass << "  expected mass: " << p_info->Mass() << endl; }
			if(fabs((int)p_info->Charge()) == 0) continue;

			TVector3 p_mom(px, py, pz);
			if(p_mom.Mag() < p_min) continue; // Questionable

			double pt = p_mom.Perp();
			if(pt < pt_min) continue;  // Avoid bad pseudorapidity warning, should be out of eta cut anyway.

			double eta = p_mom.PseudoRapidity();

			double phi = p_mom.Phi();

			// Refmult and event plane logic

			// ref
			if(fabs(eta) < ref_eta_max) refmult++;

			// ref2
			if(fabs(eta) > ref2_eta_min && fabs(eta) < ref2_eta_max && fabs((int)p_info->Charge()) == 3) refmult2++;

			// ref3
			if(fabs(eta) < ref3_eta_max && fabs(pid) != proton_pid && fabs((int)p_info->Charge()) == 3) refmult3++;

			// event plane q-vector
			if(pt > qvec_pt_min && pt < qvec_pt_max && fabs(eta) < ref3_eta_max) {
				qx += cos(2*phi); qy += sin(2*phi);
			}
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


