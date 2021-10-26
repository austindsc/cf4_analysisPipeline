
// TH1* ozGetHist(TTree * tt, string variable, string sTh1, string cuts){
	
// 	tt->Draw(sTh1+">>", cuts, "goff");
// 	return (TH1*) gDirectory->Get
// }

/*
load in all 'med' data files

for the 1.25cm run, remove events with bad baseline or end of pulse too late
*/


Double_t omedSquare(Double_t *x, Double_t *par){

	float d = x[0];
	float a = par[0]; // side length
	Double_t top = a*a;
	Double_t bottom = 2*d * TMath::Sqrt(4*d*d + 2*a*a);

	return 4*TMath::ATan(top/bottom) /4 /TMath::Pi();
}

Double_t fLY_omedSquare(Double_t *x, Double_t *par){

	float d = x[0];
	float a = par[0]; // side length
	Double_t top = a*a;
	Double_t bottom = 2*d * TMath::Sqrt(4*d*d + 2*a*a);

	Double_t Rfrac = 4*TMath::ATan(top/bottom) /4 /TMath::Pi();

	return par[1] * Rfrac;
}

Double_t omedSimpleSquare(Double_t *x, Double_t *par){

	Double_t d = x[0];
	Double_t a = par[0]; // sipm side length

	Double_t sipmArea = a*a;
	Double_t circleArea = 4 * TMath::Pi() * d * d;

	return sipmArea / circleArea;
}

int gasCompMedium(){
	// TFile * fmed = TFile::Open("../procData/outProcStatic1_20210812_1.root", "READ");
	// TTree * tmed = (TTree*) fmed->Get("tt");
	// tmed->Show(0);

	TFile * fAmed = TFile::Open("../procData/outProcStatic2_20210723_3.root", "READ");
	TTree * tAmed = (TTree*) fAmed->Get("tt");
	tAmed->Show(0);

	TFile * fBmed = TFile::Open("../procData/outProcStatic2_20210803_5.root", "READ");
	TTree * tBmed = (TTree*) fBmed->Get("tt");

	TFile * fCmed = TFile::Open("../procData/outProcStatic2_20210812_0.root", "READ");
	TTree * tCmed = (TTree*) fCmed->Get("tt");

	///////////// building histograms /////////////
	int nqBins = 200;
	int rLowQ = -10;
	int rHiQ = 600;

	// cut values (require this)
	float cutBline = 0.045; //blineV
	float cutTend = 1945.;	//tsEnd
	string scBline = "blineV < " + to_string(cutBline);
	string scTend = "tsEnd < " + to_string(cutTend);
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
	TF1 * fAg1 = new TF1("fAg1", "gaus", 100, 175);
	hAmed_qsBoth->Fit("fAg1", "+RMN0", "goff");

	// TF1 * fAg2 = new TF1("fAg2", "gaus", fAg1->GetParameter(1) - scaleFitBack*fAg1->GetParameter(2), fAg1->GetParameter(1)+ scaleFitForward*fAg1->GetParameter(2));
	TF1 * fAg2 = new TF1("fAg2", "gaus", 60, fAg1->GetParameter(1) - scaleFitForward*fAg1->GetParameter(2));
	hAmed_qsBoth->Fit("fAg2", "+RMN0", "goff");

	//---- fitting B ----//
	// TF1 * fBg1 = new TF1("fBg1", "gaus", 250, 600);	
	TF1 * fBg1 = new TF1("fBg1", "gaus", 150, 400);	
	hBmed_qsRaw->Fit("fBg1", "+RMN0", "goff");

	TF1 * fBg2 = new TF1("fBg2", "gaus", 80, 300);	
	hBmed_qsBoth->Fit("fBg2", "+RMN0", "goff");

	//---- fitting C ----//
	TF1 * fCg1 = new TF1("fCg1", "gaus", 150, 250);
	hCmed_qsBoth->Fit("fCg1", "+RMN0", "goff");

	// TF1 * fCg2 = new TF1("fCg2", "gaus", fCg1->GetParameter(1) - scaleFitBack*fCg1->GetParameter(2), fCg1->GetParameter(1)+ scaleFitForward*fCg1->GetParameter(2));
	TF1 * fCg2 = new TF1("fCg2", "gaus", 80, fCg1->GetParameter(1) - scaleFitForward*fCg1->GetParameter(2));
	hCmed_qsBoth->Fit("fCg2", "+RMN0", "goff");

	// how to improve this strategy
	/*
	1. fit start with exponential, include exponential in all gaussian fits
	2. repeat similar analysis on medium-trigger data (embrace weird shape)

	*/

	//////// plotting fit  R2 data ////////

	float scaleA = fAg1->GetParameter(1)/fAg2->GetParameter(1);

	float scaleC = fCg1->GetParameter(1)/fCg2->GetParameter(1);

	float testScaleB = fBg1->GetParameter(1)/fBg2->GetParameter(1);

	float scaleB = (scaleA + scaleC)/2;

	float runBlight = fBg2->GetParameter(1) * scaleB;
	cout << "estimated light = " << runBlight << endl;

	cout << "scaleA = " << scaleA << endl;
	cout << "scaleC = " << scaleC << endl;
	cout << "testScaleB = " << testScaleB << endl;


	// plot linearly, x = 1/R^2
	// plot pure, x = R
	// when x = R^2, funtionality is 1/x

	float xr[3] = { 2.5, 1.9, 1.25 }; //opposite order
	float xr2[3] = { 2.5*2.5, 1.9*1.9, 1.25*1.25 }; //opposite order
	float xr2lin[3] = { 1/(2.5*2.5), 1/(1.9*1.9), 1/(1.25*1.25) }; //opposite order
	Double_t er = 2 * 0.2; // multiply by 1/R^2
	Double_t exr2[3] = { er/(2.5*2.5), er/(1.9*1.9), er/(1.25*1.25) }; //opposite order

	//float yLightScaled[3] = {(float)fAg1->GetParameter(1), (float)fCg1->GetParameter(1), runBlight};
	float yLightScaled[3] = {(float)fAg1->GetParameter(1), (float)fCg1->GetParameter(1), (float)fBg1->GetParameter(1)};

	///// -- building TGraphs

	////// 1/x^2
	TGraph * tgR = new TGraph(3, xr, yLightScaled);
	tgR->SetTitle("1/x^2");
	TF1 * fLYr = new TF1("fLYr", "[0]/(x*x)", 1., 10);
	tgR->Fit("fLYr", "+RMN0");

	// correct function?
	// TF1 * fLYr_good1 = new TF1("fLYr_good", "[0] * 1/(1 + x*x/[1])", 0, 4);
	// tgR->Fit("fLYr_good", "+RMN0");

	// actual correct function
	TF1 * fLYr_accurate = new TF1("fLYr_accurate", fLY_omedSquare, 1., 4, 2);
	Double_t sipmSide = 0.6*2+0.05; //in cm
	fLYr_accurate->SetParameter(0, sipmSide);
	fLYr_accurate->SetParLimits(0, sipmSide, sipmSide);
	tgR->Fit("fLYr_accurate", "+RMN0");

	TCanvas * cR = new TCanvas();
	tgR->SetMarkerStyle(20);
	tgR->Draw("AP");
	fLYr->Draw("same");

	// fLYr_good1->SetLineColor(4);
	// fLYr_good1->Draw("same");

	fLYr_accurate->SetLineColor(3);
	fLYr_accurate->Draw("same");
	cR->Update();


	////// 1/x
	TGraph * tgR2 = new TGraph(3, xr2, yLightScaled);
	tgR2->SetTitle("1/x");
	TF1 * fLYr2 = new TF1("fLYr2", "[0]/x", 2., 10);
	tgR2->Fit("fLYr2", "+RMNO");

	TCanvas * cR2 = new TCanvas();
	tgR2->SetMarkerStyle(20);
	tgR2->Draw("AP");
	fLYr2->Draw("same");

	cR2->Update();

	////// x
	TGraph * tgRlin = new TGraph(3, xr2lin, yLightScaled);
	tgRlin->SetTitle("x");
	TF1 * fLYr2lin = new TF1("fLYr2lin", "[0] * x", 0, 0.6);
	tgRlin->Fit("fLYr2lin", "+RMN0");

	TCanvas * cRlin = new TCanvas();

	tgRlin->SetMarkerStyle(20);
	tgRlin->Draw("AP");
	fLYr2lin->Draw("same");

	cRlin->Update();

	///////////////// plotting /////////////////

	// TCanvas *cg = new TCanvas();
	// tgR2->SetMarkerStyle(20);
	// tgR2->Draw("AP");

	// cg->Update();


	// TCanvas *cg2 = new TCanvas();
	// tgR_rawlight->SetMarkerStyle(20);
	// tgR_rawlight->Draw("AP");
	// fLYr2->Draw("same");


	// cg2->Update();

	//// runA, B, C, comparing cuts ////
	TCanvas * c1 = new TCanvas();

	TLegend * leg1 = new TLegend(0.85, 0.85, "Comparing cuts", "ARC NDC");

	hAmed_qsRaw->SetLineColor(1);
	hAmed_qsRaw->Draw("hist");

	hAmed_qsBoth->SetLineColor(1);
	hAmed_qsBoth->SetLineStyle(9);
	hAmed_qsBoth->SetLineWidth(3);
	hAmed_qsBoth->Draw("hist same");

	hCmed_qsRaw->SetLineColor(4);
	hCmed_qsRaw->Draw("hist same");

	hCmed_qsBoth->SetLineColor(4);
	hCmed_qsBoth->SetLineStyle(9);
	hCmed_qsBoth->SetLineWidth(3);
	hCmed_qsBoth->Draw("hist same");

	hBmed_qsRaw->SetLineColor(3);
	hBmed_qsRaw->Draw("hist same");

	hBmed_qsBoth->SetLineStyle(9);
	hBmed_qsBoth->SetLineWidth(3);
	hBmed_qsBoth->SetLineColor(3);
	hBmed_qsBoth->Draw("hist same");

	// ---- drawing fits ---- //
	fAg1->SetLineWidth(3);
	fAg1->Draw("same");

	fAg2->SetLineColor(6);
	fAg2->SetLineWidth(3);
	fAg2->Draw("same");

	fBg1->SetLineWidth(3);
	fBg1->Draw("same");

	fBg2->SetLineWidth(3);
	fBg2->SetLineColor(6);
	fBg2->Draw("same");

	fCg1->SetLineWidth(3);
	fCg1->Draw("same");

	fCg2->SetLineWidth(3);
	fCg2->SetLineColor(6);
	fCg2->Draw("same");

	/// ---- formatting ---- ///
	leg1->AddEntry(hAmed_qsRaw, "2.5cm, raw", "l");
	leg1->AddEntry(hAmed_qsBoth, "2.5cm, cuts", "l");
	leg1->AddEntry(hCmed_qsRaw, "1.9cm, raw", "l");
	// leg1->AddEntry(hCmed_qsBoth, "1.9cm, cuts", "l");
	leg1->AddEntry(hBmed_qsRaw, "1.25cm, raw", "l");
	// leg1->AddEntry(hBmed_qsBoth, "1.25cm, cuts", "l");
	leg1->Draw("same");

	hAmed_qsRaw->SetTitle("Gas - Mega Trigger - alpha spectra");
	hAmed_qsRaw->GetXaxis()->SetTitle("Raw Charge (V x ns)");
	hAmed_qsRaw->GetYaxis()->SetTitle("Counts");

	c1->Update();

	//// runB qs cuts ////
	/*
	TCanvas * c2 = new TCanvas();

	hBmed_qsRaw->SetLineColor(1);
	hBmed_qsRaw->Draw("hist");

	hBmed_qsBline->SetLineColor(2);
	hBmed_qsBline->Draw("hist same");

	hBmed_qsTend->SetLineColor(3);
	hBmed_qsTend->Draw("hist same");

	hBmed_qsBoth->SetLineColor(4);
	hBmed_qsBoth->Draw("hist same");

	hBmed_qsRaw->SetTitle("1.25cm Run, Mega trig, QS, comparing cuts");
	hBmed_qsRaw->GetXaxis()->SetTitle("raw charge (V x ns)");

	TLegend * leg2 = new TLegend(0.85, 0.85, "Comparing cuts", "ARC NDC");
	leg2->AddEntry(hBmed_qsRaw, "no cuts", "l");
	leg2->AddEntry(hBmed_qsBline, "baseline cut", "l");
	leg2->AddEntry(hBmed_qsTend, "time cut (end of pulse)", "l");
	leg2->AddEntry(hBmed_qsBoth, "both cuts", "l");
	leg2->Draw("same");
	c2->Update();
	*/
	
	//// all runs, cuts ////






	return 0;
}