#include <vector>

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
	float p_min = 0.1;  // STAR minimum momentum acceptance
	float pt_min = 0.01;  // Issues arise in eta if pt is zero
	float qvec_pt_min = 0.2;
	float qvec_pt_max = 2.0;
	float ref_eta_max = 0.5;
	float ref2_eta_min = 0.5;
	float ref2_eta_max = 1.0;
	float ref3_eta_max = 1.0;
	float eta_max = 1.0;
	float mass_qa_percent = 1.0;  // % difference in mass to output warning
	int buffer_size = 5000000;
	int split_level = 1;

	// Input file variables
	int evn, tn, nov, npp, npt, nesp, ninesp, nest, ninest, pid;  // for ampt.dat event
	float px, py, pz, mass, x, y, z, t;  // for ampt.dat track
	float empty;  // For last entry on event line


	// Output tree variables
	int event=0, refmult, refmult2, refmult3;  // event variables
	float imp, qx, qy;  // event variables
	vector<int> pid_vec;  // track variables
	vector<float> px_vec;  // Apparently need to be declared on separate lines
	vector<float> py_vec;
	vector<float> pz_vec;  // track variables


	//------------define a root file and tree :------------------------------
	TFile *fampt = new TFile(("data_" + run_id + ".root").data(),"RECREATE");
	if(!fampt){
		cout << "Output file cannot be opened" << endl;
		exit(0);
	}

	TTree *tr = new TTree("tree","AMPT Data");


	//Define event branches:-------------------------------------------
	tr->Branch("event",    &event,     "event/I");
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
	tr->Branch("pid",      &pid_vec, buffer_size, split_level);
	tr->Branch("px",       &px_vec,  buffer_size, split_level);
	tr->Branch("py",       &py_vec,  buffer_size, split_level);
	tr->Branch("pz",       &pz_vec,  buffer_size, split_level);

	//**************************************************************************************

	cout<<"making .root file from .dat file..."<<endl;

	ifstream infile;
	char *fname   = new char[100];
	sprintf(fname,"ana/ampt.dat");
	infile.open(fname);

	while(infile)
	{// event loop
		infile>>evn>>tn>>nov>>imp>>npp>>npt>>nesp>>ninesp>>nest>>ninest>>empty;

		cout << "New event: " << evn << endl;

		if (infile.eof()) { cout << "end of file true. event number " << evn << endl; }
		else { cout << "end of file false. event number " << evn << endl; }

		if(infile.eof())  break;

		event += 1;

		//****************************ampt.dat particle loop**************
		qx = 0.; qy = 0.;
		refmult = 0; refmult2 = 0; refmult3 = 0;
		pid_vec.clear(); px_vec.clear(); py_vec.clear(); pz_vec.clear();

		for(int j=0;j<nov;j++)                          //particle loop
		{
			infile>>pid>>px>>py>>pz>>mass>>x>>y>>z>>t;

			cout << "New particle: " << pid << " " << px << " " << py << " " << pz << endl;

			px_vec.push_back(px);
			py_vec.push_back(py);
			pz_vec.push_back(pz);
			pid_vec.push_back(pid);

			p_info = db->GetParticle((int)pid);
			if(!p_info) { cout << "pid: " << pid << " not in TDatabasePDG" << endl; continue; }
			else if(p_info->Mass() * 1+mass_qa_percent/100 < mass || p_info->Mass() * 1-mass_qa_percent/100 > mass) {
				cout << "pid: " << pid << " with mass " << mass << "  expected mass: " << p_info->Mass() << endl;
			}
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


