
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
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLine.h"

using namespace std;

float get_chi_diff(float b, int energy);
int plot_ref_dist_comp(int energy);
void minimize_chi_b(vector<int> energy);

int ref_dist_comp() {
//	plot_ref_dist_comp(energy);
	minimize_chi_b({7, 11, 19, 27, 39, 62});
	return 0;
}

void minimize_chi_b(vector<int> energy_vec) {
	float acc = 0.01;
	map<int, int> energy_color {{7, 1}, {11, 2}, {19, 4}, {27, 6}, {39, 8}, {62, 9}};
	map<int, float> b_min {{7, 14}, {11, 14}, {19, 13.5}, {27, 13.5}, {39, 13}, {62, 13}};
	map<int, float> opt;

//	vector<float> b_cuts;
//	for(float b = b_cut_min; b <= b_cut_max; b += b_cut_step) {
//		b_cuts.push_back(b);
//	}

	auto *mg = new TMultiGraph();
	auto *leg = new TLegend();
	float chi_max = 0;

	for(int energy:energy_vec) {
		float b_cut_step = 0.5;
		vector<float> chi_diff;
		vector<float> b;
		map<float, float> b_chi_map;
		b.push_back(b_min[energy]);
		cout << "Calculating chi^2 for " + to_string(energy) + "GeV, b = " << b.back() << endl;
		chi_diff.push_back(get_chi_diff(b.back(), energy));
		b_chi_map[b.back()] = chi_diff.back();

		while(fabs(b_cut_step) >= acc) {
			b.push_back(b.back() + b_cut_step);
			cout << "Calculating chi^2 for " + to_string(energy) + "GeV, b = " << b.back() << endl;
			chi_diff.push_back(get_chi_diff(b.back(), energy));
			b_chi_map[b.back()] = chi_diff.back();
			while(chi_diff.back() < chi_diff[(int)chi_diff.size() - 2]) {
				b.push_back(b.back() + b_cut_step);
				cout << "Calculating chi^2 for " + to_string(energy) + "GeV, b = " << b.back() << endl;
				chi_diff.push_back(get_chi_diff(b.back(), energy));
				b_chi_map[b.back()] = chi_diff.back();
			}
			b_cut_step /= -2;
		}
		opt[energy] = b[(int)b.size() - 2];
		if(chi_diff[(int)chi_diff.size() - 2] > chi_max) { chi_max = chi_diff[(int)chi_diff.size() - 2]; }

		chi_diff.clear(); b.clear();
		for(pair<float, float> b_chi:b_chi_map) {
			b.push_back(b_chi.first);
			chi_diff.push_back(b_chi.second);
		}

		TGraph *graph = new TGraph((int)b.size(), b.data(), chi_diff.data());
		graph->SetMarkerColor(energy_color[energy]);
		graph->SetLineColor(energy_color[energy]);
		graph->SetLineWidth(3);
		graph->SetMarkerStyle(8);
		mg->Add(graph, "l");
		leg->AddEntry(graph, (to_string(energy) + "GeV").data(), "l");
	}

	mg->GetXaxis()->SetTitle("Impact Parameter b (fm)");
	mg->GetYaxis()->SetTitle("Sum( (BES1 - AMPT)^2 )");
	mg->SetTitle("BES1 to AMPT Ref3 distribution Chi^2 vs b cutoff (max)");
	mg->Draw("a");
	leg->Draw();

	for(pair<int, float> opt_pars:opt) {
		cout << "Energy " << opt_pars.first << "GeV optimal b = " << opt_pars.second << "fm" << endl;
		auto *line = new TLine(opt_pars.second, 0.0, opt_pars.second, chi_max);
		line->SetLineColor(energy_color[opt_pars.first]);
		line->Draw();
	}
}


int plot_ref_dist_comp(int energy) {
	string data_path = "/media/dylan/SSD_Storage/Research/Data_Old_Ref3/eta050/" + to_string(energy) + "GeV/QA_" + to_string(energy) + "GeV.root";
	string ampt_path = "/media/dylan/SSD_Storage/Research/Trees_Ampt/" + to_string(energy) + ".root";
	float b_cut = 20.0;

	TFile *data_file = new TFile(data_path.data(), "READ");
	TH1I *int_data_dist = (TH1I*)data_file->Get(("pre_reftwo_eta050_" + to_string(energy)).data());
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


float get_chi_diff(float b_cut, int energy) {
	string data_path = "/media/dylan/SSD_Storage/Research/Data_Old_Ref3/eta050/" + to_string(energy) + "GeV/QA_" + to_string(energy) + "GeV.root";
	string ampt_path = "/media/dylan/SSD_Storage/Research/Trees_Ampt/" + to_string(energy) + ".root";

	TFile *data_file = new TFile(data_path.data(), "READ");
	TH1I *int_data_dist = (TH1I*)data_file->Get(("pre_reftwo_eta050_" + to_string(energy)).data());
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
