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
		bool verbose = 0;
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
	float scaleFromPE = 2.6; // match areas from 2.6PE onwards (3PE bin)

	for(int irun = 0; irun < 6; irun++){
		TH1 * h = hqr[irun];
		float scaleVal = h->Integral(h->GetXaxis()->FindBin(scaleFromPE), h->GetNbinsX());
		h->Scale(1/scaleVal);
		// h->Scale(1/h->GetBinWidth(10));
	}

	///////// PDF work (all manual fitting) /////////

	auto c1 = new TCanvas();
	int pMade = 25;
	float eff = 0.1;
	Double_t sipmSide = 6*2 + 0.5; // in mm
	Double_t dA = 25; // in mm
	Double_t dB = 12.5; // in mm

	TF1 * fA = new TF1("fA", MFF::myBinomial, 0, 30, 2);
	fA->SetParameters(eff, pMade);

	// TF1 * f1 = new TF1("f1", MFF::myBinomial, 0, 30, 2);
	// f1->SetParameters(eff*2, pMade);

	// TF1 * f2 = new TF1("f2", MFF::myBinomial, 0, 30, 2);
	// f2->SetParameters(eff, pMade*2);

	float bFactor = MFF::omegaSquare(&dB, &sipmSide)/MFF::omegaSquare(&dA, &sipmSide);
	TF1 * fB = new TF1("fB", MFF::myBinomial, 0, 30, 2);
	fB->SetParameters(eff, bFactor * pMade);

	///// poisson fitting
	TF1 * fpoisA = new TF1("fpoisA", MFF::scaledPoisson,-0.5, 30, 2);
	fpoisA->SetParameter(0, 2.5);
	fpoisA->SetParameter(1, scaleFromPE);
	fpoisA->SetParLimits(1, scaleFromPE, scaleFromPE);
	hiAqr200->Fit("fpoisA", "+M0", "", scaleFromPE, 15);

	TF1 * fpoisB = new TF1("fpoisB", MFF::scaledPoisson,-0.5, 30, 2);
	fpoisB->SetParameter(0, 3.);
	fpoisB->SetParameter(1, scaleFromPE);
	fpoisB->SetParLimits(1, scaleFromPE, scaleFromPE);
	hiBqr200->Fit("fpoisB", "+M0", "", scaleFromPE, 15);

	TF1 * fpoisC = new TF1("fpoisC", MFF::scaledPoisson,-0.5, 30, 2);
	fpoisC->SetParameter(0, 4.);
	fpoisC->SetParameter(1, scaleFromPE);
	fpoisC->SetParLimits(1, scaleFromPE, scaleFromPE);
	hiCqr200->Fit("fpoisC", "+M0", "", scaleFromPE, 20);

	////////// drawing //////////////
	hiAqr200->SetLineColor(1);
	hiAqr200->Draw("hist");

	// fA->SetLineStyle(9);
	// fA->SetLineColor(1);
	// fA->Draw("same");

	hiBqr200->SetLineColor(3);
	hiBqr200->Draw("hist same");

	hiCqr200->SetLineColor(4);
	hiCqr200->Draw("hist same");

	// fB->SetLineColor(4);
	// fB->SetLineStyle(9);
	// fB->Draw("same");

	// f1->SetLineColor(2);
	// f1->Draw("same");

	// f2->SetLineColor(3);
	// f2->Draw("same");

	fpoisA->Draw("same");
	fpoisB->Draw("same");
	fpoisC->Draw("same");

	c1->SetLogy();
	
	// return 0;

	/////// drawing ///////
	TCanvas * c2 = new TCanvas();
	TLegend * leg1 = new TLegend(0.85, 0.85, "Liquid spectra", "ARC NDC");

	hAqr200->GetXaxis()->SetTitle("PE");
	hAqr200->SetTitle("Liquid-alpha, comparing Runs");
	hAqr200->SetLineColor(1);
	hAqr200->Draw("hist");

	hBqr200->SetLineColor(3);
	hBqr200->Draw("hist same");

	hCqr200->SetLineColor(4);
	hCqr200->Draw("hist same");

	c2->SetLogy();
	return 0;
}