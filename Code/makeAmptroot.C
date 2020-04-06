
#include<TTree.h>
#include<TRandom.h>
#include<iostream>
#include<math.h>
#include<stdio.h>
#include<fstream>
#include<map>
#include"TH2F.h"
#include"TCanvas.h"
#include"TLorentzVector.h"
#include"TROOT.h"
#include "TFile.h"

int makeAmptroot()
{
//input file variables:
  Int_t    evn,tn,nov,npp,npt,nesp,ninesp,nest,ninest,pid,counter=0; // for ampt.dat
  Float_t    px,py,pz,mass,X,Y,Z,t,b;                                //for ampt.dat
//  Int_t    evn2, it, na, nb, nab, sr, stat, tpp,id1,id2;
  Float_t    nx, ny, zz, theta_p , phi_p , theta_t , phi_t , imp, psi;
  psi = 0;

 


//Booked variables in tree: 
  const Int_t  mul  = 90000;
  const Int_t  nucl = 600;
  Int_t  Event=0, Na, Nb ,Nab, Mult, Npartp, Npartt, Nesp, Ninesp, Nest, Ninest;
  Float_t  Imp,  Psi; // event variables
  Float_t   Theta_p , Phi_p , Theta_t , Phi_t; // event variables
  Int_t  Stat[nucl],PID1[nucl], PID2[nucl], PID[mul];        //particle variables
  Float_t  Nx[nucl], Ny[nucl], Nz[nucl], Px[mul], Py[mul], Pz[mul];
  Float_t   XX[mul],  YY[mul],  ZZ[mul], TT[mul], Mass[mul]; //particle variables
  
  Char_t outfile[100];

//------------define a root file and tree :------------------------------
  TFile *fampt = new TFile("ana/test.root","RECREATE");
  if(!fampt){
    cout << "Output file cannot be opened" << endl;
    exit(0);
  }


  // TFile *fampt = new TFile("AuAu_7GeV_SM_ntmax150_3mb_0_15fm_file0.root","recreate");
  TTree *tr    = new TTree("tr","Reconst ntuple");



//Define event branches:-------------------------------------------
 tr->Branch("Event", &Event,"Event/I");
 tr->Branch("Mult",  &Mult, "Mult/I");           // multiplicity = tracks
 tr->Branch("Npartp",&Npartp,"Npartp/I");
 tr->Branch("Npartt",&Npartt,"Npartt/I");
 tr->Branch("Nesp",  &Nesp, "Nesp/I");
 tr->Branch("Ninesp",&Ninesp, "Ninesp/I");
 tr->Branch("Nest",  &Nest, "Nest/I");
 tr->Branch("Ninest",&Ninest, "Ninest/I");
 tr->Branch("Imp",&Imp, "Imp/F");
// tr->Branch("Na", &Na, "Na/I");       //Na = Projectile mass no.
// tr->Branch("Nb", &Nb, "Nb/I");       //Nb = Target mass no.
// tr->Branch("Nab",&Nab, "Nab/I");
// tr->Branch("Psi",&Psi, "Psi/F");



 //--for U+U only:
 /* tr->Branch("Theta_p", &Theta_p,"Theta_p/F");
 tr->Branch("Phi_p",    &Phi_p,    "Phi_p/F");
 tr->Branch("Theta_t", &Theta_t, "Theta_t/F");
 tr->Branch("Phi_t",    &Phi_t,     "Phi_t/F");

*/


//particle branches:
//  tr->Branch("Nx",  &Nx,  "Nx[Nab]/F");   // Nab = na+nb; kore dio niche;
//  tr->Branch("Ny",  &Ny,  "Ny[Nab]/F");
//  tr->Branch("Nz",  &Nz,  "Nz[Nab]/F");
//  tr->Branch("Stat",&Stat,"Stat[Nab]/I");
//  tr->Branch("PID1",&PID1,"PID1[Nab]/I");
//  tr->Branch("PID2",&PID2,"PID2[Nab]/I");

  tr->Branch("PID", &PID, "PID[Mult]/I");
  tr->Branch("Px",  &Px,  "Px[Mult]/F");
  tr->Branch("Py",  &Py,  "Py[Mult]/F");
  tr->Branch("Pz",  &Pz,  "Pz[Mult]/F");
  tr->Branch("Mass",&Mass,"Mass[Mult]/F");
  tr->Branch("XX",  &XX,  "XX[Mult]/F");
  tr->Branch("YY",  &YY,  "YY[Mult]/F");
  tr->Branch("ZZ",  &ZZ,  "ZZ[Mult]/F");
  tr->Branch("TT",  &TT,  "TT[Mult]/F");

//**************************************************************************************

cout<<"making .root file from .dat file..."<<endl;


ifstream infile;//[300000]; 

char *fname   = new char[100];

 int mull,signal,signal2,oldpid,oldstat;
 double oldnx,oldny,oldnz;
 double oldxx,oldyy,oldzz,oldpx,oldpy,oldpz;


       sprintf(fname,"ana/ampt.dat");
       infile.open(fname); 

while(infile)   
    {//2 event loop // Please see ampt.dat file and npart-xy.dat to comment proper line:
      infile>>evn>>tn>>nov>>b>>npp>>npt>>nesp>>ninesp>>nest>>ninest;//>>psi;    // for newer Ampt version


      if(infile.eof())  break;

 
        Event   = Event+1;  
        Mult    = nov; 
        Imp     = b;
        Npartp  = npp;
        Npartt  = npt;         //tpp = npp+npt;
        Nesp    = nesp;
        Ninesp  = ninesp;
        Nest    = nest;
        Ninest  = ninest ;
	Psi     = psi;

//for U+U only: //if you use for AuAu, these values will be set zero 
	/*      Theta_p = theta_p; 
        Phi_p    = phi_p;
        Theta_t = theta_t;
        Phi_t    = phi_t;
*/
	mull = 0; //counter for produced particle.
	oldxx = 0;
	oldyy = 0;
	oldzz = 0;
	oldpx = 0;
	oldpy = 0;
	oldpz = 0;
	signal = 0;
	oldpid = 0;
 //****************************ampt.dat particle loop************** 
 for(int j=0;j<nov;j++)                          //particle loop
    {
            infile>>pid>>px>>py>>pz>>mass>>X>>Y>>Z>>t;
            Px[j] = px;    
            Py[j] = py;
            Pz[j] = pz;
          Mass[j] = mass;
           PID[j] = pid;
 	    XX[j] = X;
	    YY[j] = Y;
	    ZZ[j] = Z;
	    TT[j] = t;

	    if(oldpid==pid && oldxx==X && oldyy==Y && oldzz==Z && 
               oldpx==px && oldpy==py && oldpz==pz){
               signal++;
	       //break;
              } 

	    mull++;
	    oldxx = X;
	    oldyy = Y;
	    oldzz = Z;
	    oldpx = px;
	    oldpy = py;
	    oldpz = pz;
	    oldpid = pid;
    }
 
//***************************ampt.dat particle loop ends**************
 
 oldnx = 0; oldny = 0; oldnz = 0; oldstat=0;
 signal2 = 0; 
 
tr->Fill();      


}//2 event loop end

//}//1 folder loop end

fampt->cd();
tr->Write();
fampt->Close();
cout<<" Tree written succesfully "<<endl;

return 0;
}//main end


