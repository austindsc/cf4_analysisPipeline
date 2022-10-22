/*
load in all 'med' data files

for the 1.25cm run, remove events with bad baseline or end of pulse too late
*/
#include "helpers/fitFunctions.C"
#include "helpers/getSPE.C"

int gasComp2Medium(){

	// 3V data
	TFile * fAmed = TFile::Open("../procData/outProcStatic2_20210723_3.root", "READ");
	TTree * tAmed = (TTree*) fAmed->Get("tt");
	tAmed->Show(0);

	TFile * fBmed = TFile::Open("../procData/outProcStatic2_20210803_5.root", "READ");
	TTree * tBmed = (TTree*) fBmed->Get("tt");

	TFile * fCmed = TFile::Open("../procData/outProcStatic2_20210812_0.root", "READ");
	TTree * tCmed = (TTree*) fCmed->Get("tt");

	///////////// building histograms /////////////
	int nqBins = 300;
	int rLowQ = -10;
	int rHiQ = 600;

	// cut values (require this)
	float cutBline = 0.045; //blineV
	float cutTend = 1945.;	//tsEnd
	string scBline = "blineV < " + to_string(cutBline);
	string scTend = "tsEnd < " + to_string(cutTend) + " && tsEnd > 20";
	string scBoth = scBline + " && " + scTend;

	///////////////// runA, med trigger /////////////////
	TH1 * hAmed_qsRaw = new TH1F("hAmed_qsRaw", "runA, medTrig, smart Q", nqBins, rLowQ, rHiQ);
	tAmed->Draw("qs50>>hAmed_qsRaw", "", "goff");

	TH1 * hAmed_qsBoth = new TH1F("hAmed_qsBoth", "runA, medTrig, smart Q", nqBins, rLowQ, rHiQ);
	tAmed->Draw("qs50>>hAmed_qsBoth", scBoth.c_str(), "goff");

	///////////////// runB, med trigger /////////////////
	TH1 * hBmed_tEnd = new TH1F("hBmed_tend", "runB, uncut pulse tail", 1000, 0, 2000);
	tBmed->Draw("tsEnd>>hBmed_tend", "", "goff");

	TH1 * hBmed_qsRaw = new TH1F("hBmed_qsRaw", "runB, medTrig, smartQ raw", nqBins, rLowQ, rHiQ);
	TH1 * hBmed_qsBline = (TH1*)hBmed_qsRaw->Clone("hBmed_qsBline");
	TH1 * hBmed_qsTend = (TH1*)hBmed_qsRaw->Clone("hBmed_qsTend");
	TH1 * hBmed_qsBoth = (TH1*)hBmed_qsRaw->Clone("hBmed_qsBoth");

	tBmed->Draw("qs50>>hBmed_qsRaw", "", "goff");
	tBmed->Draw("qs50>>hBmed_qsBline", scBline.c_str(), "goff");
	tBmed->Draw("qs50>>hBmed_qsTend", scTend.c_str(), "goff");
	tBmed->Draw("qs50>>hBmed_qsBoth", scBoth.c_str(), "goff");

	///////////////// runC, med trigger /////////////////
	TH1 * hCmed_qsRaw = new TH1F("hCmed_qsRaw", "runC, medTrig, smart Q", nqBins, rLowQ, rHiQ);
	tCmed->Draw("qs50>>hCmed_qsRaw", "", "goff");

	TH1 * hCmed_qsBoth = new TH1F("hCmed_qsBoth", "runC, medTrig, smart Q", nqBins, rLowQ, rHiQ);
	tCmed->Draw("qs50>>hCmed_qsBoth", scBoth.c_str(), "goff");

	////////////// fitting //////////////
	float scaleFitBack = 3.;
	float scaleFitForward = 0.5;
	float startFit2 = 120.;

	//---- fitting A ----//
	// TF1 * fAg1 = new TF1("fAg1", "gaus", 175, 230);
	TF1 * fAg1 = new TF1("fAg1", "gaus", 80, 120);
	hAmed_qsBoth->Fit("fAg1", "+RMN0", "goff");

	//---- fitting B ----//
	// TF1 * fBg1 = new TF1("fBg1", "gaus", 250, 600);	
	TF1 * fBg1 = new TF1("fBg1", "gaus", 150, 400);	
	hBmed_qsRaw->Fit("fBg1", "+RMN0", "goff");

	//---- fitting C ----//
	TF1 * fCg1 = new TF1("fCg1", "gaus", 130, 230);
	hCmed_qsBoth->Fit("fCg1", "+RMN0", "goff");

	// how to improve this strategy
	/*
	1. fit start with exponential, include exponential in all gaussian fits
	2. repeat similar analysis on medium-trigger data (embrace weird shape)

	*/

	//////// plotting fit  R2 data ////////

	float xr[3] = { 2.5, 1.9, 1.25 }; //opposite order
	Double_t er = 2 * 0.2; // multiply by 1/R^2
	Double_t exr[3] = { er, er, er }; //opposite order

	float yLight[3] = {(float)fAg1->GetParameter(1), (float)fCg1->GetParameter(1), (float)fBg1->GetParameter(1)};
	
	///// -- building TGraph -- \\\\\

	gStyle->SetOptFit(1);
	TGraph * tgR = new TGraph(3, xr, yLight);
	tgR->SetTitle("Light vs Distance");

	// actual correct function
	TF1 * fLYr_accurate = new TF1("fLYr_accurate", MFF::fLY_omegaSquare, 1., 4, 2);
	Double_t sipmSide = 0.6*2+0.05; //in cm
	fLYr_accurate->SetParameter(0, sipmSide);
	fLYr_accurate->SetParLimits(0, sipmSide, sipmSide);
	tgR->Fit("fLYr_accurate", "+M0", "", 2.2, 4);

	TCanvas * cR = new TCanvas();
	tgR->SetMarkerStyle(20);
	tgR->Draw("AP");

	fLYr_accurate->SetLineColor(2);
	fLYr_accurate->Draw("same");
	cR->Update();

	///////////////// plotting /////////////////

	//// runA, B, C, comparing cuts ////
	TCanvas * c1 = new TCanvas();

	TLegend * leg1 = new TLegend(0.85, 0.85, "Comparing cuts", "ARC NDC");

	hAmed_qsRaw->SetLineColor(1);
	hAmed_qsRaw->Draw("hist");

	hAmed_qsBoth->SetLineColor(1);
	hAmed_qsBoth->SetLineStyle(9);
	hAmed_qsBoth->SetLineWidth(3);
	hAmed_qsBoth->Draw("hist same");

	hBmed_qsRaw->SetLineColor(3);
	hBmed_qsRaw->Draw("hist same");

	hBmed_qsBoth->SetLineColor(3);
	hBmed_qsBoth->SetLineStyle(9);
	hBmed_qsBoth->SetLineWidth(3);
	hBmed_qsBoth->Draw("hist same");

	hCmed_qsRaw->SetLineColor(4);
	hCmed_qsRaw->Draw("hist same");

	hCmed_qsBoth->SetLineColor(4);
	hCmed_qsBoth->SetLineStyle(9);
	hCmed_qsBoth->SetLineWidth(3);
	hCmed_qsBoth->Draw("hist same");

	// ---- drawing fits ---- //
	fAg1->SetLineWidth(3);
	fAg1->Draw("same");

	fBg1->SetLineWidth(3);
	fBg1->Draw("same");

	fCg1->SetLineWidth(3);
	fCg1->Draw("same");

	/// ---- formatting ---- ///
	leg1->AddEntry(hAmed_qsRaw, "2.5cm, raw", "l");
	leg1->AddEntry(hAmed_qsBoth, "2.5cm, cuts", "l");
	leg1->AddEntry(hCmed_qsRaw, "1.9cm, raw", "l");
	// leg1->AddEntry(hCmed_qsBoth, "1.9cm, cuts", "l");
	leg1->AddEntry(hBmed_qsRaw, "1.25cm, raw", "l");
	// leg1->AddEntry(hBmed_qsBoth, "1.25cm, cuts", "l");
	leg1->Draw("same");

	hAmed_qsRaw->SetTitle("Gas - Medium Trigger - alpha spectra");
	hAmed_qsRaw->GetXaxis()->SetTitle("Raw Charge (V x ns)");
	hAmed_qsRaw->GetYaxis()->SetTitle("Counts");

	
	return 0;
}