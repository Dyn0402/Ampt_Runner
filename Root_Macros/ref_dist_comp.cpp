
#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraph.h"

using namespace std;

float get_chi_diff(float);
int plot_ref_dist_comp();
void minimize_chi_b();

int ref_dist_comp() {
	plot_ref_dist_comp();
//	minimize_chi_b();
	return 0;
}

void minimize_chi_b() {
	float b_cut_min = 14.0;
	float b_cut_max = 16.5;
	float b_cut_step = 0.1;

	vector<float> b_cuts;
	for(float b = b_cut_min; b <= b_cut_max; b += b_cut_step) {
		b_cuts.push_back(b);
	}

	vector<float> chi_diff;
	for(float b:b_cuts) {
		cout << "Calculating chi^2 for b = " << b << endl;
		chi_diff.push_back(get_chi_diff(b));
	}

	TGraph *graph = new TGraph((int)b_cuts.size(), b_cuts.data(), chi_diff.data());
	graph->Draw();
}


int plot_ref_dist_comp() {
	string data_path = "/media/dylan/SSD_Storage/Research/Data_Old_Ref3/eta050/19GeV/QA_19GeV.root";
	string ampt_path = "/media/dylan/SSD_Storage/Research/Trees_Ampt/19.root";
	float b_cut = 14.75;

	TFile *data_file = new TFile(data_path.data(), "READ");
	TH1I *int_data_dist = (TH1I*)data_file->Get("pre_reftwo_eta050_19");
	TH1D* data_dist = new TH1D("data_dist", "BES1 Distribution", int_data_dist->GetNbinsX(), int_data_dist->GetXaxis()->GetXmin(), int_data_dist->GetXaxis()->GetXmax());
	for(int i =1; i<=int_data_dist->GetNbinsX(); i++) { data_dist->SetBinContent(i, int_data_dist->GetBinContent(i)); }
	data_dist->SetLineColor(kBlack);

	TFile *ampt_file = new TFile(ampt_path.data(), "READ");
	TTree *ampt_tree = (TTree*)ampt_file->Get("tree");
	cout << data_dist->GetXaxis()->GetXmin() << " " << data_dist->GetXaxis()->GetXmax() << endl;
	TH1D *ampt_dist = new TH1D("ampt_dist", "Ampt Distribution", data_dist->GetNbinsX(), data_dist->GetXaxis()->GetXmin(), data_dist->GetXaxis()->GetXmax());
	ampt_dist->SetLineColor(kBlue);
	TLeaf *dist_leaf = ampt_tree->GetLeaf("refmult3");
	TLeaf *b_leaf = ampt_tree->GetLeaf("imp");
	int event_index = 0;
	while(ampt_tree->GetEntry(event_index)) {
		if(b_leaf->GetValue() < b_cut) {
			ampt_dist->Fill(dist_leaf->GetValue());
		}
		event_index++;
	}

//	for(int i=0; i<=802; i++) { cout << "Bin " << i << ": " << data_dist->GetBinContent(i) << endl; }

	auto leg1 = new TLegend();
	leg1->AddEntry(data_dist, "BES1 Distribution", "l");
	leg1->AddEntry(ampt_dist, "Ampt Distribution", "l");
	TCanvas *can1 = new TCanvas("Raw Canvas");
	data_dist->Draw();
	ampt_dist->Draw("SAMES");
	leg1->Draw();

	TH1D *norm_data_dist = (TH1D*) data_dist->Clone("norm_data_dist");
	TH1D *norm_ampt_dist = (TH1D*) ampt_dist->Clone("norm_ampt_dist");
	norm_data_dist->Scale(1.0/data_dist->Integral());
	norm_ampt_dist->Scale(1.0/ampt_dist->Integral());
	auto leg2 = new TLegend();

//	for(int i=0; i<=802; i++) { cout << "Bin " << i << ": " << norm_data_dist->GetBinContent(i) << endl; }

	leg2->AddEntry(norm_data_dist, "BES1 Distribution", "l");
	leg2->AddEntry(norm_ampt_dist, "Ampt Distribution", "l");
	TCanvas *can2 = new TCanvas("Norm Canvas");
	norm_data_dist->Draw("HIST");
	norm_ampt_dist->Draw("HISTSAMES");
	leg2->Draw();

	TH1D *ratio_dist = new TH1D("ratio_dist", "Ampt / BES1 Ratio (Normed Values)", data_dist->GetNbinsX(), data_dist->GetXaxis()->GetXmin(), data_dist->GetXaxis()->GetXmax());
	for(int i =1; i<=int_data_dist->GetNbinsX(); i++) {
		if(norm_data_dist->GetBinContent(i) > 0) {
			ratio_dist->SetBinContent(i, norm_ampt_dist->GetBinContent(i) / norm_data_dist->GetBinContent(i));
		} else {
			ratio_dist->SetBinContent(i, 1.0);
		}
	}
	TCanvas *can3 = new TCanvas("Ratio Canvas");
	ratio_dist->Draw();

	return 0;
}


float get_chi_diff(float b_cut) {
	string data_path = "/media/dylan/SSD_Storage/Research/Data_Old_Ref3/eta050/7GeV/QA_7GeV.root";
	string ampt_path = "/media/dylan/SSD_Storage/Research/Trees_Ampt/7.1.root";

	TFile *data_file = new TFile(data_path.data(), "READ");
	TH1I *int_data_dist = (TH1I*)data_file->Get("pre_reftwo_eta050_7");
	TH1D* data_dist = new TH1D("data_dist", "BES1 Distribution", int_data_dist->GetNbinsX(), int_data_dist->GetXaxis()->GetXmin(), int_data_dist->GetXaxis()->GetXmax());
	for(int i =1; i<=int_data_dist->GetNbinsX(); i++) { data_dist->SetBinContent(i, int_data_dist->GetBinContent(i)); }

	TFile *ampt_file = new TFile(ampt_path.data(), "READ");
	TTree *ampt_tree = (TTree*)ampt_file->Get("tree");
	TH1D *ampt_dist = new TH1D("ampt_dist", "Ampt Distribution", data_dist->GetNbinsX(), data_dist->GetXaxis()->GetXmin(), data_dist->GetXaxis()->GetXmax());
	TLeaf *dist_leaf = ampt_tree->GetLeaf("refmult3");
	TLeaf *b_leaf = ampt_tree->GetLeaf("imp");
	int event_index = 0;
	while(ampt_tree->GetEntry(event_index)) {
		if(b_leaf->GetValue() < b_cut) {
			ampt_dist->Fill(dist_leaf->GetValue());
		}
		event_index++;
	}

	TH1D *norm_data_dist = (TH1D*) data_dist->Clone("norm_data_dist");
	TH1D *norm_ampt_dist = (TH1D*) ampt_dist->Clone("norm_ampt_dist");
	norm_data_dist->Scale(1.0/data_dist->Integral());
	norm_ampt_dist->Scale(1.0/ampt_dist->Integral());

	float chi2 = 0;
	for(int i = 1; i <= norm_data_dist->GetNbinsX(); i++) {
		chi2 += pow(norm_data_dist->GetBinContent(i) - norm_ampt_dist ->GetBinContent(i), 2);
	}

	delete int_data_dist;
	delete data_dist;
	delete ampt_dist;
	delete norm_data_dist;
	delete norm_ampt_dist;

	data_file->Close();
	ampt_file->Close();
	delete data_file;
	delete ampt_file;

	return chi2;
}
