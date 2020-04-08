
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

using namespace std;

int makeAmptroot()
{
	//Constants
	const int proton_pid = 2212;
	float p_min = 0.15;
	float pt_min = 0.3;
	float pt_max = 2.5;
	float eta_max = 1.0;

	//input file variables:
	Int_t    evn,tn,nov,npp,npt,nesp,ninesp,nest,ninest,pid,counter=0; // for ampt.dat
	Float_t    px,py,pz,mass,x,y,z,t,b;                                //for ampt.dat


	//Booked variables in tree:
	const Int_t  mul  = 90000;
	Int_t  event=0, pmult;
	Float_t  imp; // event variables
	Int_t   pid_array[mul];        //particle variables
	Float_t px_array[mul], py_array[mul], pz_array[mul];


	//------------define a root file and tree :------------------------------
	TFile *fampt = new TFile("ana/test.root","RECREATE");
	if(!fampt){
		cout << "Output file cannot be opened" << endl;
		exit(0);
	}

	TTree *tr = new TTree("tree","AMPT Proton Data");



	//Define event branches:-------------------------------------------
	tr->Branch("event", &event,"event/I");
	tr->Branch("pmult",  &pmult, "pmult/I");           // multiplicity = tracks
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
		for(int j=0;j<nov;j++)                          //particle loop
		{
			infile>>pid>>px>>py>>pz>>mass>>x>>y>>z>>t;
			if(fabs(pid) == proton_pid) {
				TVector3 p_mom(px, py, pz);
				if(p_mom.Mag() >= p_min && p_mom.Perp() >= pt_min && p_mom.Perp() <= pt_max && fabs(p_mom.PseudoRapidity()) <= eta_max) {
					p_px.push_back(px);
					p_py.push_back(py);
					p_pz.push_back(pz);
					p_pid.push_back(pid);
				}
			}
		}

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


