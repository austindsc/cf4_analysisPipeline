#include "helpers/getSPE.C"
#include "helpers/fitFunctions.C"

// TH1* ozGetHist(TTree * tt, string variable, string sTh1, string cuts){
	
// 	tt->Draw(sTh1+">>", cuts, "goff");
// 	return (TH1*) gDirectory->Get
// }

/*
load in all 'mega' data files

for the 1.25cm run, remove events with bad baseline or end of pulse too late
*/

/*
1. load in data files
1. get SPE values from other files
3. 

???? WHERE DOES +1 IN SCALING COME FROM (in picking bin) ????


*/


int liqAlphaYield1(){

	TFile * files5V[3] = {
				TFile::Open("$CF4HOME/procData/liqMedFull_OPS2_20210726_8.root", "READ"),
				TFile::Open("$CF4HOME/procData/liqMedFull_OPS2_20210804_5.root", "READ"),
				TFile::Open("$CF4HOME/procData/liqMedFull_OPS2_20210813_5.root", "READ")
	};

	TTree * trees5V[3];
	for(int irun = 0; irun < 3; irun++){
		trees5V[irun] = (TTree*) files5V[irun]->Get("tt");
		cout << "entries = " << trees5V[irun]->GetEntries() << endl;
	}

	//////////// getting SPE values
	speStuff spe5v[3];
	for(int irun = 0; irun < 3; irun++){
		TTree * t = trees5V[irun];
		bool verbose = 1;
		if (verbose == 1) new TCanvas();
		speStuff speInfo = getSPEval(t, 0.1, verbose);
		spe5v[irun] = speInfo;
		cout << speInfo.spe << endl;
	}

	////////////// histogram building
	int peMax = 30;
	float nbins = 300;
	float nStart = -0.5;
	string var = "qr200/";

	/// run A
	TH1* hAqr200 = new TH1F("hAqr200", "run A, qr200", nbins, nStart, peMax-0.5);
	string sA = var + to_string(spe5v[0].spe) + ">> hAqr200";
	trees5V[0]->Draw(sA.c_str(), "", "goff");

	TH1* hiAqr200 = new TH1F("hiAqr200", "run A, qr200", peMax, nStart, peMax-0.5);
	sA = var + to_string(spe5v[0].spe) + ">> hiAqr200";
	trees5V[0]->Draw(sA.c_str(), "", "goff");

	/// run B
	TH1* hBqr200 = new TH1F("hBqr200", "run B, qr200", nbins, nStart, peMax-0.5);
	string sB = var + to_string(spe5v[1].spe) + ">> hBqr200";
	trees5V[1]->Draw(sB.c_str(), "", "goff");

	TH1* hiBqr200 = new TH1F("hiBqr200", "run B, qr200", peMax, nStart, peMax-0.5);
	sB = var + to_string(spe5v[1].spe) + ">> hiBqr200";
	trees5V[1]->Draw(sB.c_str(), "", "goff");

	/// run C
	TH1* hCqr200 = new TH1F("hCqr200", "run C, qr200", nbins, nStart, peMax-0.5);
	string sC = var + to_string(spe5v[2].spe) + ">> hCqr200";
	trees5V[2]->Draw(sC.c_str(), "", "goff");

	TH1* hiCqr200 = new TH1F("hiCqr200", "run C, qr200", peMax, nStart, peMax-0.5);
	sC = var + to_string(spe5v[2].spe) + ">> hiCqr200";
	trees5V[2]->Draw(sC.c_str(), "", "goff");

	TH1* hqr[6] = {hAqr200, hBqr200, hCqr200, hiAqr200, hiBqr200, hiCqr200};

	///////// scaling and manipulating /////////
	Double_t scaleFromPE = 5; // match areas from 3 onwards (3PE bin)
	Double_t endFit = 10;

	for(int irun = 0; irun < 6; irun++){
		TH1 * h = hqr[irun];
		Double_t scaleVal = h->Integral(h->GetXaxis()->FindBin(scaleFromPE), h->GetXaxis()->FindBin(endFit));
		h->Scale(1/scaleVal);
	}

	// float testHist1 = 0;
	// for(int i = scaleFromPE; i <= endFit; i++){
	// 	testHist1 += hiAqr200->GetBinContent(hiAqr200->GetXaxis()->FindBin(i));
	// }
	// cout << "test hist1 val = " << testHist1 << endl;

	///////// PDF work and Fitting /////////

	int pMade = 25;
	float eff = 0.1;
	Double_t sipmSide = 6*2 + 0.5; // in mm
	Double_t dA = 25; // in mm
	Double_t dB = 12.5; // in mm
	Double_t dC = 19.0; // in mm

	TF1 * fA = new TF1("fA", MFF::myBinomial, 0, 30, 2);
	fA->SetParameters(eff, pMade);

	// TF1 * f1 = new TF1("f1", MFF::myBinomial, 0, 30, 2);
	// f1->SetParameters(eff*2, pMade);

	// TF1 * f2 = new TF1("f2", MFF::myBinomial, 0, 30, 2);
	// f2->SetParameters(eff, pMade*2);

	// float bFactor = MFF::omegaSquare(&dB, &sipmSide)/MFF::omegaSquare(&dA, &sipmSide);
	// TF1 * fB = new TF1("fB", MFF::myBinomial, 0, 30, 2);
	// fB->SetParameters(eff, bFactor * pMade);

	///// poisson fitting
	TF1 * fpoisA = new TF1("fpoisA", MFF::scaledPoisson,-0.5, 30, 3);
	fpoisA->SetParameter(0, 2.5);
	fpoisA->SetParameter(1, scaleFromPE);
	fpoisA->SetParLimits(1, scaleFromPE, scaleFromPE);
	fpoisA->SetParameter(2, endFit);
	fpoisA->SetParLimits(2, endFit, endFit);
	hiAqr200->Fit("fpoisA", "+M0", "", scaleFromPE, endFit);

	// float testFunc = 0;
	// float testHist = 0;
	// for(int i = scaleFromPE; i <= endFit; i++){
	// 	testFunc += fpoisA->Eval(i);
	// 	testHist += hiAqr200->GetBinContent(hiAqr200->GetXaxis()->FindBin(i));
	// }
	// cout << "test func val = " << testFunc << endl;
	// cout << "test hist val = " << testHist << endl;

	TF1 * fpoisB = new TF1("fpoisB", MFF::scaledPoisson,-0.5, 30, 3);
	fpoisB->SetParameter(0, 4.5);
	fpoisB->SetParameter(1, scaleFromPE);
	fpoisB->SetParLimits(1, scaleFromPE, scaleFromPE);
	fpoisB->SetParameter(2, endFit);
	fpoisB->SetParLimits(2, endFit, endFit);
	hiBqr200->Fit("fpoisB", "+M0", "", scaleFromPE, endFit);

	TF1 * fpoisC = new TF1("fpoisC", MFF::scaledPoisson,-0.5, 30, 3);
	fpoisC->SetParameter(0, 3);
	fpoisC->SetParameter(1, scaleFromPE);
	fpoisC->SetParLimits(1, scaleFromPE, scaleFromPE);
	fpoisC->SetParameter(2, endFit);
	fpoisC->SetParLimits(2, endFit, endFit);
	hiCqr200->Fit("fpoisC", "+M0", "", scaleFromPE, endFit);

	///// R2 fitting /////

	Double_t push = 0.; // in mm
	Double_t xr[3] = { dA + push, dC + push, dB + push }; //opposite order {A, C, B}
	Double_t testD = 0;
	Double_t er = 0.2; // multiply by 1/R^2
	Double_t exr[3] = { er, er, er }; //opposite order
	Double_t yMu[3] = { fpoisA->GetParameter(0), fpoisC->GetParameter(0), fpoisB->GetParameter(0) };
	Double_t eyMu[3]= { fpoisA->GetParError(0), fpoisC->GetParError(0), fpoisB->GetParError(0) };
	TGraphErrors * tgR = new TGraphErrors(3, xr, yMu, exr, eyMu);

	gStyle->SetOptFit(1111);
	// actual correct function
	TF1 * fLYr_omega = new TF1("fLYr_omega", MFF::fLY_omegaSquare, 0, 60, 2);
	fLYr_omega->SetParameter(0, sipmSide);
	fLYr_omega->SetParLimits(0, sipmSide, sipmSide);
	tgR->Fit("fLYr_omega", "+M0", "", 0, 60);

	// actual correct function
	TF1 * fLYr_omega2 = new TF1("fLYr_omega2", MFF::fLY_omegaSquare2, 0, 60, 3);
	fLYr_omega2->SetParameter(0, sipmSide);
	fLYr_omega2->SetParLimits(0, sipmSide, sipmSide);
	fLYr_omega2->SetParameter(1, 300.);
	fLYr_omega2->SetParameter(2, 40.);
	tgR->Fit("fLYr_omega2", "+M0", "same", 0, 60);

	// with attenuation
	// TF1 * fLYr_atten = new TF1("fLYr_atten", MFF::fLY_omegaAttenuation, 10, 60, 3);
	// fLYr_atten->SetParameter(0, sipmSide);
	// fLYr_atten->SetParLimits(0, sipmSide, sipmSide);

	// fLYr_atten->SetParameter(1, 4 * fLYr_omega->GetParameter(1)); // LY factor
	// fLYr_atten->SetParameter(2, 20); // attenuation k (1/mm)
	// tgR->Fit("fLYr_atten", "+MN0", "", 0, 40);

	TCanvas * cR = new TCanvas();
	tgR->SetMarkerStyle(20);
	tgR->Draw("AP");
	fLYr_omega->SetLineColor(2);
	fLYr_omega->Draw("same");

	fLYr_omega2->SetLineColor(3);
	fLYr_omega2->Draw("same");

	tgR->GetXaxis()->SetTitle("True Distance (mm)");
	tgR->GetYaxis()->SetTitle("Poisson fit average PE");
	tgR->SetTitle("LCF4 - alpha - LY vs Distance");

	cR->Update();

	// fLYr_atten->Draw("same");

	////////// drawing //////////////
	TCanvas * c2 = new TCanvas();
	TLegend * leg1 = new TLegend();
	gStyle->SetOptStat(0);

	hAqr200->GetXaxis()->SetTitle("PE");
	hAqr200->SetTitle("Liquid-alpha, comparing Runs");
	hAqr200->SetLineColor(1);
	hAqr200->Draw("hist");

	hBqr200->SetLineColor(3);
	hBqr200->Draw("hist same");

	hCqr200->SetLineColor(4);
	hCqr200->Draw("hist same");

	leg1->SetTextSize(0.04);
	leg1->SetHeader("Source - SiPM Distance");
	leg1->AddEntry(hAqr200, "25 mm", "L");
	leg1->AddEntry(hCqr200, "19 mm", "L");
	leg1->AddEntry(hBqr200, "12.5 mm", "L");
	leg1->Draw("same");

	c2->SetLogy();

	///// drawing /////
	TCanvas * c1 = new TCanvas();

	hiAqr200->SetFillColorAlpha(1, 0.3);
	hiAqr200->Draw("hist");
	// fpoisA->SetLineColor(1);
	fpoisA->Draw("same");
	c1->SetLogy();

	TCanvas * c1n2 = new TCanvas();
	hiBqr200->SetFillColorAlpha(3, 1);
	hiBqr200->Draw("hist");
	fpoisB->Draw("same");

	c1n2->SetLogy();
	c1n2->Update();

	TCanvas * c1n3 = new TCanvas();
	hiCqr200->SetFillColorAlpha(4, 1);
	hiCqr200->Draw("hist same");
	fpoisC->Draw("same");

	c1n3->SetLogy();


	// fA->SetLineStyle(9);
	// fA->SetLineColor(1);
	// fA->Draw("same");

	// fB->SetLineColor(4);
	// fB->SetLineStyle(9);
	// fB->Draw("same");

	// f1->SetLineColor(2);
	// f1->Draw("same");

	// f2->SetLineColor(3);
	// f2->Draw("same");


	// fpoisB->SetLineColor(3);
	// fpoisC->SetLineColor(4);

	return 0;
}