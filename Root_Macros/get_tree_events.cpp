
#include <iostream>

#include "TFile.h"
#include "TTree.h"

using namespace std;

int get_tree_events(string path) {

	TFile *file = new TFile(path.data(), "READ");
	TTree *tree = (TTree*)file->Get("tree");
	cout << "Number of events in tree: " << tree->GetEntries() << endl;
	file->Close();
	delete file;

	return 0;
}
