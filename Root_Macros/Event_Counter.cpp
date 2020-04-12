#include <iostream>
#include <vector>
#include <dirent.h>
#include <sys/stat.h>
#include <fstream>
#include <string>
#include <sstream>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TPad.h"

using namespace std;


vector<string> get_files_in_dir(string, string, string);
vector<string> split_string_by_char(string, char);


int Event_Counter() {
	int target_events = 2000;
	string path = "/media/dylan/SSD_Storage/Research/Trees_Ampt/19GeV/";

	vector<string> in_files = get_files_in_dir(path, "root", "path");
	TH1I *num_events = new TH1I("Num_Events", "Number of AMPT Events", target_events / 100 + 1, -0.5, target_events + 0.5);
	int num_quit = 0, total_events = 0;
	for(string file:in_files) {
		TFile *f = new TFile(file.data(), "READ");
		TTree *tree = (TTree*)f->Get("tree");
		num_events->Fill(tree->GetEntries());
		total_events += tree->GetEntries();
		if(tree->GetEntries() < target_events) { cout << "Entries: " << tree->GetEntries() << endl; num_quit++; }
		f->Close();
		delete f;
	}
	cout << endl << "Number died early: " << num_quit << endl;
	cout << "Total number of events: " << total_events << endl;
	num_events->Draw();
	return 0;
}


vector<string> split_string_by_char(string str, char del) {
	vector<string> split;
	stringstream stream(str);
	string hold;
	while(getline(stream, hold, del)) {
		split.push_back(hold);
	}

	return(split);
}


//Return name (not path) of all files in dir_path with extension ext.
vector<string> get_files_in_dir(string dir_path, string ext, string out) {
	vector<string> files;
	struct stat info;

	if( stat(dir_path.data(), &info) !=0 ) {
		cout << "Can't access path: " << dir_path << " Skipping following read/write." << endl;
	} else {
		DIR* files_dir = opendir(dir_path.data());
		struct dirent* dp;
		while((dp=readdir(files_dir)) != NULL) {
			string file = dp->d_name;
			vector<string> fields = split_string_by_char(file, '.');
			if(fields.back() == ext) {
				files.push_back(file);
			}
		}
		if(out == "path") {
			for(unsigned i=0; i<files.size(); i++) {
				files[i] = dir_path + files[i];
			}
		}
	}

	return(files);
}
