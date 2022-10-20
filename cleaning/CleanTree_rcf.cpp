#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TMath.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TList.h"
#include "TString.h"
#include "TRandom.h" 
#include "TRandom3.h" 
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include <algorithm>
#include <stdlib.h>
#include "TComplex.h"
#include "TVector2.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include <string>
#include "TEntryList.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"

using namespace std;

const float etamax = 1;
const int max_ident = 2;

bool IdenticalTrack(float px1, float py1, float pz1, float pid1, float px2, float py2, float pz2, float pid2);

void CleanTree_rcf(const Char_t *inFile = "placeholder.list", const TString newdir = "./", const char *ext=".root")
{       
    int nfile = 0;
    ifstream fin(inFile);
    string line;
    TString fullname, dirname, fname;
    while(getline(fin, line)) //loop wiill run till end of file
    {
        cout << line << endl;
        size_t found = line.find_last_of("/");
        dirname = (line.substr(0,found)).c_str();
        fname   = (line.substr(found+1)).c_str();
        fullname = line.c_str();
        nfile++;         

        TFile *f = new TFile(fullname);
        TTree *oldtree, *newtree;
        gDirectory->GetObject("tree", oldtree);

        std::vector<float> *px_vec = nullptr;
        std::vector<float> *py_vec = nullptr;
        std::vector<float> *pz_vec = nullptr;
        std::vector<int>   *pid_vec= nullptr;
        oldtree->SetBranchAddress("px", &px_vec);
        oldtree->SetBranchAddress("py", &py_vec);
        oldtree->SetBranchAddress("pz", &pz_vec);
        oldtree->SetBranchAddress("pid",&pid_vec);

        bool isCorrupted = false;
        std::vector<int> BadEventList;
        /******* checking the file *******/
        // event loop
        for (int i = 0; i < oldtree->GetEntries(); ++i)
        {
            oldtree->GetEntry(i);

            int num_ident = 0;
            float px1, py1, pz1, px2, py2, pz2;
            int   pid1, pid2;
            float eta1, eta2, theta1, theta2, pt1, pt2;

            // particle loop
            for (int j = 0; j < px_vec->size(); ++j)
            {
                pid1   = pid_vec->at(j);
                px1    = px_vec->at(j);
                py1    = py_vec->at(j);
                pz1    = pz_vec->at(j);
                pt1    = sqrt(px1*px1 + py1*py1);
                theta1 = atan2(pt1,pz1);
                eta1   = -log(tan(theta1/2.));

                // avoid some particles
                if (isnan(eta1)) continue;
                if (fabs(eta1) > etamax) continue;
                if (abs(pid1) == 111 || abs(pid1) == 313) continue;

                // cout << "PID || " << pid1 << " || px || " << px1 << " || py || " << py1 <<" || pz || " << pz1 << " || eta || " << eta1 << endl;

                for (int k = j+1; k < px_vec->size(); ++k)
                {   
                    pid2   = pid_vec->at(k);
                    px2    = px_vec->at(k);
                    py2    = py_vec->at(k);
                    pz2    = pz_vec->at(k);
                    pt2    = sqrt(px2*px2 + py2*py2);
                    theta2 = atan2(pt2,pz2);
                    eta2   = -log(tan(theta2/2.));

                    if (isnan(eta2)) continue;
                    if (fabs(eta2) > etamax) continue;
                    if (abs(pid2) == 111 || abs(pid2) == 313) continue;

                    if (IdenticalTrack(px1, py1, pz1, pid1, px2, py2, pz2, pid2)) num_ident += 1;
                }
            }
            if (num_ident >= max_ident) {
                cout << "Event # " << i << " corrupt, num identical: " << num_ident << endl;
                isCorrupted = true;
                BadEventList.push_back(i);
            }
        }

        cout << setw(25) << fname << " || Is Corrupted? || " << setw(4) << (bool)isCorrupted << " || # of Bad Events || " << setw(4) << BadEventList.size() << " || # of Total Events || " << oldtree->GetEntries() << endl; 
        /******* updating the file *******/
        if (isCorrupted) 
        {
            // match branches in old files
            int event=0, refmult, refmult2, refmult3;
            float imp, qx, qy;
            int npp, npt, nesp, ninesp, nest, ninest;
            int split_level = 1;
            int buffer_size = 5000000;

            oldtree->SetBranchAddress("event",    &event);
            oldtree->SetBranchAddress("refmult",  &refmult);
            oldtree->SetBranchAddress("refmult2", &refmult2);
            oldtree->SetBranchAddress("refmult3", &refmult3);
            oldtree->SetBranchAddress("qx",       &qx);
            oldtree->SetBranchAddress("qy",       &qy);
            oldtree->SetBranchAddress("imp",      &imp);
            oldtree->SetBranchAddress("npp",      &npp);
            oldtree->SetBranchAddress("npt",      &npt);
            oldtree->SetBranchAddress("nesp",     &nesp);
            oldtree->SetBranchAddress("ninesp",   &ninesp);
            oldtree->SetBranchAddress("nest",     &nest);
            oldtree->SetBranchAddress("ninest",   &ninest);

            //------------define a root file and tree :------------------------------
            TFile *fnew = new TFile(newdir+"fix_"+fname,"RECREATE");
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

            for (int i = 0; i < oldtree->GetEntries(); i++)
            {
                oldtree->GetEntry(i);
                if (std::find(BadEventList.begin(), BadEventList.end(), i) == BadEventList.end()) // if not bad event
                {
                    tr->Fill();
                }
            }

            fnew->cd();
            tr->Write();
            fnew->Close();    
        }
        
        f->Close();
        delete f; 

        // delete good input files from ./INPUTFILES/
        if (!isCorrupted) gSystem->Unlink(fullname);
    }
    fin.close();    
    
}

bool IdenticalTrack(float px1, float py1, float pz1, float pid1, float px2, float py2, float pz2, float pid2)
{
    if (px1 == px2 && py1 == py2 && pz1 == pz2 && pid1 == pid2) return true;
    return false;
}
