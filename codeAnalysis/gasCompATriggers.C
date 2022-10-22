#include "helpers/fitFunctions.C"
#include "helpers/getSPE.C"

int gasCompATriggers(){

	// 3V, medium
	TFile * fm3 = TFile::Open("../procData/outProcStatic2_20210723_3.root", "READ");
	TTree * tm3 = (TTree*) fm3->Get("tt");

	// 3V, mega
	TFile * fM3 = TFile::Open("../procData/outProcStatic2_20210723_8.root", "READ");
	TTree * tM3 = (TTree*) fM3->Get("tt");

	// 4V, medium
	TFile * fm4 = TFile::Open("../procData/50koutProcStatic2_20210723_4.root", "READ");
	TTree * tm4 = (TTree*) fm4->Get("tt");

	// 4V, mega
	TFile * fM4 = TFile::Open("../procData/outProcStatic2_20210723_9.root", "READ");
	TTree * tM4 = (TTree*) fM4->Get("tt");

	// 5V, medium
	TFile * fm5 = TFile::Open("../procData/50koutProcStatic2_20210723_5.root", "READ");
	TTree * tm5 = (TTree*) fm5->Get("tt");

	// 5V, mega
	TFile * fM5 = TFile::Open("../procData/50koutProcStatic2_20210723_10.root", "READ");
	TTree * tM5 = (TTree*) fM5->Get("tt");
	
	// tM4->Draw("blineV>>h(100, -0.1, 0.1");
	// tM4->Draw("tsEnd>>h(500, 1500, 2000");
	// return 0;

	///////////// building histograms /////////////
	int nqBins = 300;
	int rLowQ = -10;
	int rHiQ = 600;

	///////////////// runA, med trigger /////////////////
	TH1 * hm3 = new TH1F("hm3", "hm3", nqBins, rLowQ, rHiQ);
	tm3->Draw("qs50>>hm3", "", "goff");

	TH1 * hm3c = new TH1F("hm3c", "hm3 cuts", nqBins, rLowQ, rHiQ);
	tm3->Draw("qs50>>hm3c", "blineV < 0.044 && tsEnd < 1945", "goff");

	TH1 * hM3 = new TH1F("hM3", "hM3", nqBins, rLowQ, rHiQ);
	tM3->Draw("qs50>>hM3", "", "goff");

	TH1 * hM3c = new TH1F("hM3c", "hM3c", nqBins, rLowQ, rHiQ);
	tM3->Draw("qs50>>hM3c", "blineV < 0.044 && tsEnd < 1945", "goff");

	TH1 * hm4 = new TH1F("hm4", "hm4", nqBins, rLowQ, rHiQ);
	tm4->Draw("qs50>>hm4", "", "goff");

	TH1 * hm4c = new TH1F("hm4c", "hm4c", nqBins, rLowQ, rHiQ);
	tm4->Draw("qs50>>hm4c", "blineV < 0.044 && tsEnd < 1945", "goff");

	TH1 * hM4 = new TH1F("hM4", "hM4", nqBins, rLowQ, rHiQ);
	tM4->Draw("qs50>>hM4", "", "goff");

	TH1 * hM4c = new TH1F("hM4c", "hM4c", nqBins, rLowQ, rHiQ);
	tM4->Draw("qs50>>hM4c", "blineV < 0.044 && tsEnd < 1945", "goff");

	TH1 * hm5 = new TH1F("hm5", "hm5", nqBins, rLowQ, rHiQ);
	tm5->Draw("qs50>>hm5", "", "goff");

	TH1 * hm5c = new TH1F("hm5c", "hm5c", nqBins, rLowQ, rHiQ);
	tm5->Draw("qs50>>hm5c", "blineV < 0.044 && tsEnd < 1945", "goff");

	TH1 * hM5 = new TH1F("hM5", "hM5", nqBins, rLowQ, rHiQ);
	tM5->Draw("qs50>>hM5", "", "goff");

	TH1 * hM5c = new TH1F("hM5c", "hM5c", nqBins, rLowQ, rHiQ);
	tM5->Draw("qs50>>hM5c", "blineV < 0.044 && tsEnd < 1945", "goff");

	TCanvas * c1 = new TCanvas();
	TLegend * leg1 = new TLegend(0.85, 0.85, "Comparing data", "ARC NDC");
 
	hm3->SetLineColor(1);
	hm3->Draw("hist");

	hm3c->SetLineColor(1);
	hm3c->SetLineStyle(9);
	hm3c->SetLineWidth(3);
	hm3c->Draw("hist same");

	hM3->SetLineColor(2);
	hM3->Draw("hist same");

	hM3c->SetLineColor(2);
	hM3c->SetLineStyle(9);
	hM3c->SetLineWidth(3);
	hM3c->Draw("hist same");

	hm4->SetLineColor(3);
	hm4->Draw("hist same");

	hm4c->SetLineColor(3);
	hm4c->SetLineStyle(9);
	hm4c->SetLineWidth(3);
	hm4c->Draw("hist same");

	hM4->SetLineColor(4);
	hM4->Draw("hist same");

	hM4c->SetLineColor(4);
	hM4c->SetLineStyle(9);
	hM4c->SetLineWidth(3);
	hM4c->Draw("hist same");

	hm5->SetLineColor(6);
	hm5->Draw("hist same");

	hm5c->SetLineColor(6);
	hm5c->SetLineStyle(9);
	hm5c->SetLineWidth(3);
	hm5c->Draw("hist same");

	hM5->SetLineColor(7);
	hM5->Draw("hist same");

	hM5c->SetLineColor(7);
	hM5c->SetLineStyle(9);
	hM5c->SetLineWidth(3);
	hM5c->Draw("hist same");

	leg1->AddEntry(hm3, "medium 3V", "l");
	leg1->AddEntry(hm3c, "cuts medium 3V", "l");

	leg1->AddEntry(hM3, "mega 3V", "l");
	
	leg1->AddEntry(hm4, "medium 4V", "l");

	leg1->AddEntry(hM4, "mega 4V", "l");

	leg1->AddEntry(hm5, "medium 5V", "l");

	leg1->AddEntry(hM5, "mega 5V", "l");

	leg1->Draw("same");

	return 0;
}