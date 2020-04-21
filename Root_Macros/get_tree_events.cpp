
#include <iostream>
#include <dirent.h>
#include <sys/stat.h>
#include <fstream>
#include <string>
#include <sstream>
#include <map>


#include "TFile.h"
#include "TTree.h"

using namespace std;


vector<string> get_files_in_dir(string, string, string);
vector<string> split_string_by_char(string, char);


int get_tree_events(string path) {

	vector<string> in_files = get_files_in_dir(path, "root", "path");
	for(string file:in_files) {
		TFile *f = new TFile(file.data(), "READ");
		TTree *tree = (TTree*)f->Get("tree");
		cout << "Number of events in tree: " << file << " " << tree->GetEntries() << endl;
		f->Close();
		delete f;
	}

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
