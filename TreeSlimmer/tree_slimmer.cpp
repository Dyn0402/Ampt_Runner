
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

using namespace std;


void slim_tree(string file);

struct ampt_tree_branches {
	int event;
	int refmult;
	int refmult2;
	int refmult3;

	float imp;
	float qx, qy;

	int npp, npt, nesp, ninesp, nest, ninest;

	vector<int>* pid = 0;
	vector<float>* px = 0;
	vector<float>* py = 0;
	vector<float>* pz = 0;

	TBranch* branch_pid = 0;
	TBranch* branch_px = 0;
	TBranch* branch_py = 0;
	TBranch* branch_pz = 0;
};


void set_ampt_tree_branches(TTree* tree, ampt_tree_branches& branches) {
	tree->SetBranchAddress("event", &branches.event);

	tree->SetBranchAddress("refmult", &branches.refmult);
	tree->SetBranchAddress("refmult2", &branches.refmult2);
	tree->SetBranchAddress("refmult3", &branches.refmult3);

	tree->SetBranchAddress("imp", &branches.imp);
	tree->SetBranchAddress("qx", &branches.qx);
	tree->SetBranchAddress("qy", &branches.qy);

	tree->SetBranchAddress("npp", &branches.npp);
	tree->SetBranchAddress("npt", &branches.npt);
	tree->SetBranchAddress("nesp", &branches.nesp);
	tree->SetBranchAddress("ninesp", &branches.ninesp);
	tree->SetBranchAddress("nest", &branches.nest);
	tree->SetBranchAddress("ninest", &branches.ninest);

	tree->SetBranchAddress("pid", &branches.pid, &branches.branch_pid);
	tree->SetBranchAddress("px", &branches.px, &branches.branch_px);
	tree->SetBranchAddress("py", &branches.py, &branches.branch_py);
	tree->SetBranchAddress("pz", &branches.pz, &branches.branch_pz);
}


void tree_slimmer(string input_file_list) {
	ifstream in_file(input_file_list);
	if (in_file.is_open()) {
		//stringstream ss(input_file_list);
		string file;

		if (input_file_list != NULL) {
			while (getline(in_file, file)) {
				cout << "Slim " << file << endl;
				slim_tree(file);
			}
		}
	}
}


void slim_tree(string file) {
	TFile *f_in = new TFile(file.data(), "READ");
	TTree *tree_in = (TTree*)f_in->Get("tree");
	
	ampt_tree_branches branches;
	set_ampt_tree_branches(tree_in, branches);

	string f_out_name = file.substr(0, file.rfind('.')) + "_protons.root";
	TFile* f_out = new TFile(f_out_name.data(), "RECREATE");
	TTree* tree_out = new TTree("tree", "AMPT Data");

	cout << "Slimming " << file << " to " << f_out_name << endl;

	const int proton_pid = 2212;
	float pt_min = 0.3;  // GeV
	float pt_max = 2.2;  // GeV
	float eta_max = 2.1; 

	int buffer_size = 5000000;
	int split_level = 1;

	vector<int> pid_vec;  // track variables
	vector<float> px_vec;  // Apparently need to be declared on separate lines
	vector<float> py_vec;
	vector<float> pz_vec;  // track variables

	//Define event branches:-------------------------------------------
	tree_out->Branch("event", &branches.event, "event/I");
	tree_out->Branch("refmult", &branches.refmult, "refmult/I");
	tree_out->Branch("refmult2", &branches.refmult2, "refmult2/I");
	tree_out->Branch("refmult3", &branches.refmult3, "refmult3/I");
	tree_out->Branch("qx", &branches.qx, "qx/F");
	tree_out->Branch("qy", &branches.qy, "qy/F");
	tree_out->Branch("imp", &branches.imp, "imp/F");

	tree_out->Branch("npp", &branches.npp, "npp/I");
	tree_out->Branch("npt", &branches.npt, "npt/I");
	tree_out->Branch("nesp", &branches.nesp, "nesp/I");
	tree_out->Branch("ninesp", &branches.ninesp, "ninesp/I");
	tree_out->Branch("nest", &branches.nest, "nest/I");
	tree_out->Branch("ninest", &branches.ninest, "ninest/I");

	//particle branches:
	tree_out->Branch("pid", &pid_vec, buffer_size, split_level);
	tree_out->Branch("px", &px_vec, buffer_size, split_level);
	tree_out->Branch("py", &py_vec, buffer_size, split_level);
	tree_out->Branch("pz", &pz_vec, buffer_size, split_level);

	int event_index = 0;
	while (tree_in->GetEvent(event_index)) {
		for (int particle_index = 0; particle_index < (int)branches.pid->size(); particle_index++) {
			TVector3 p_mom(branches.px->at(particle_index), branches.py->at(particle_index), branches.pz->at(particle_index));
			if (p_mom.Perp() < pt_min || p_mom.Perp() > pt_max) continue;
			if (fabs(p_mom.PseudoRapidity()) > eta_max) continue;
			px_vec.push_back(branches.px->at(particle_index));
			py_vec.push_back(branches.py->at(particle_index));
			pz_vec.push_back(branches.pz->at(particle_index));
			pid_vec.push_back(branches.pid->at(particle_index));
		}
		tree_out->Fill();
	}

	f_out->cd();
	tree_out->Write();
	f_out->Close();

	cout << file << " slimmed successfully to " << f_out_name << endl;

	tree_in->ResetBranchAddresses();
	f_in->Close();
}