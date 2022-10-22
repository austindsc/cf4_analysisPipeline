
// TH1* ozGetHist(TTree * tt, string variable, string sTh1, string cuts){
	
// 	tt->Draw(sTh1+">>", cuts, "goff");
// 	return (TH1*) gDirectory->Get
// }

/*
load in all 'mega' data files

for the 1.25cm run, remove events with bad baseline or end of pulse too late
*/

int gasCompMegaCuts(){

	// 3V data
	TFile * fAmega = TFile::Open("../procData/outProcStatic2_20210723_8.root", "READ");
	TTree * tAmega = (TTree*) fAmega->Get("tt");
	tAmega->Show(0);

	TFile * fBmega = TFile::Open("../procData/outProcStatic2_20210803_9.root", "READ");
	TTree * tBmega = (TTree*) fBmega->Get("tt");

	TFile * fCmega = TFile::Open("../procData/outProcStatic2_20210812_4.root", "READ");
	TTree * tCmega = (TTree*) fCmega->Get("tt");

	///////////// building histograms /////////////
	int nqBins = 300;
	int rLowQ = -10;
	int rHiQ = 600;

	// cut values (require this)
	float cutBline = 0.045; //blineV
	float cutTend = 1945.;	//tsEnd
	string scBline = "blineV < " + to_string(cutBline);
	string scTend = "tsEnd < " + to_string(cutTend);
	string scBoth = scBline + " && " + scTend;

	///////////////// runA, mega trigger /////////////////
	TH1 * hAmega_qS = new TH1F("hAmega_qS", "runA, megaTrig, smart Q", nqBins, rLowQ, rHiQ);
	tAmega->Draw("qs50>>hAmega_qS", "blineV < 0.045", "goff");

	///////////////// runB, mega trigger /////////////////
	TH1 * hBmega_tEnd = new TH1F("hBmega_tend", "runB, uncut pulse tail", 1000, 0, 2000);
	tBmega->Draw("tsEnd>>hBmega_tend", "", "goff");

	TH1 * hBmega_qsRaw = new TH1F("hBmega_qsRaw", "runB, megaTrig, smartQ raw", nqBins, rLowQ, rHiQ);
	TH1 * hBmega_qsBline = (TH1*)hBmega_qsRaw->Clone("hBmega_qsBline");
	TH1 * hBmega_qsTend = (TH1*)hBmega_qsRaw->Clone("hBmega_qsTend");
	TH1 * hBmega_qsBoth = (TH1*)hBmega_qsRaw->Clone("hBmega_qsBoth");

	tBmega->Draw("qs50>>hBmega_qsRaw", "", "goff");
	tBmega->Draw("qs50>>hBmega_qsBline", scBline.c_str(), "goff");
	tBmega->Draw("qs50>>hBmega_qsTend", scTend.c_str(), "goff");
	tBmega->Draw("qs50>>hBmega_qsBoth", scBoth.c_str(), "goff");


	///////////////// plotting /////////////////

	TLegend * leg1 = new TLegend(0.85, 0.85, "Comparing cuts", "ARC NDC");

	//// runB qS cuts ////
	hBmega_qsRaw->SetLineColor(1);
	hBmega_qsRaw->Draw("hist");

	hBmega_qsBline->SetLineColor(2);
	hBmega_qsBline->Draw("hist same");

	hBmega_qsTend->SetLineColor(3);
	hBmega_qsTend->Draw("hist same");

	hBmega_qsBoth->SetLineColor(4);
	hBmega_qsBoth->Draw("hist same");

	leg1->AddEntry(hBmega_qsRaw, "no cuts", "l");
	leg1->AddEntry(hBmega_qsBline, "baseline cut", "l");
	leg1->AddEntry(hBmega_qsTend, "time end", "l");
	leg1->AddEntry(hBmega_qsBoth, "both cuts", "l");

	return 0;
}