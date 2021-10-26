
TH1* ozGetHist(TFile * f, string name){

	TFile * fmed = TFile::Open("../procData/outProcStatic1_20210812_1.root", "READ");
	return NULL;
}

int gasComp1(){
	TFile * fmed = TFile::Open("../procData/outProcStatic1_20210812_1.root", "READ");
	TTree * tmed = (TTree*) fmed->Get("tt");
	tmed->Show(0);

	TFile * fmega = TFile::Open("../procData/outProcStatic2_20210812_5.root", "READ");
	TTree * tmega = (TTree*) fmega->Get("tt");
	cout << tmega->GetEntries() << endl;

	///////////////// medium trigger /////////////////
	tmed->Draw("qr800>>hmed_qr800(500, -10, 600)", "blineV < 0.045", "goff");
	TH1 * hmed_qr800 = (TH1*) gDirectory->Get("hmed_qr800");

	tmed->Draw("qs50>>hmed_qs(500, -10, 600)", "blineV < 0.045 && tsEnd < 1945", "goff");
	TH1 * hmed_qs = (TH1*) gDirectory->Get("hmed_qs");

	tmed->Draw("qr200>>hmed_qr200(500, -10, 600)", "blineV < 0.045", "goff");
	TH1 * hmed_qr200 = (TH1*) gDirectory->Get("hmed_qr200");

	tmed->Draw("hbline_min>>hmed_h(200, -0.1, 1)", "", "goff");
	TH1 * hmed_h = (TH1*) gDirectory->Get("hmed_h");

	tmed->Draw("blineV>>hmedB(500, -0.25, .25)", "", "goff");
	TH1 * hmedB = (TH1*) gDirectory->Get("hmedB");

	///////////////// mega trigger /////////////////
	tmega->Draw("qr200>>hmega_qr200(500, -10, 600)", "blineV < 0.045", "goff");
	TH1 * hmega_qr200 = (TH1*) gDirectory->Get("hmega_qr200");

	tmega->Draw("qs50>>hmega_qs(500, -10, 600)", "blineV < 0.045 && tsEnd < 1945", "goff");
	TH1 * hmega_qs = (TH1*) gDirectory->Get("hmega_qs");

	tmega->Draw("qr800>>hmega_qr800(500, -10, 600)", "blineV < 0.045", "goff");
	TH1 * hmega_qr800 = (TH1*) gDirectory->Get("hmega_qr800");

	tmega->Draw("hbline_min>>hmega_h(200, -0.1, 1)", "", "goff");
	TH1 * hmega_h = (TH1*) gDirectory->Get("hmega_h");

	tmega->Draw("blineV>>hmegaB(500, -0.25, .25)", "", "goff");
	TH1 * hmegaB = (TH1*) gDirectory->Get("hmegaB");


	///////////////// plotting /////////////////
	TCanvas * c1 = new TCanvas();
	TLegend * leg1 = new TLegend(0.85, 0.85, "Trigger and Q", "ARC NDC");

	hmega_qs->SetLineColor(1);
	hmega_qs->Draw("hist");

	// hmega_qr800->SetLineColor(2);
	// hmega_qr800->Draw("hist same");

	hmed_qs->SetLineColor(3);
	hmed_qs->Draw("hist same");

	// hmed_qr800->SetLineColor(4);
	// hmed_qr800->Draw("hist same");

	leg1->AddEntry(hmega_qs, "mega trigger, QS, cuts applied", "l");
	// leg1->AddEntry(hmega_qr800, "mega trigger, QR800", "l");
	leg1->AddEntry(hmed_qs, "medium trig, QS, cuts applied", "l");
	// leg1->AddEntry(hmed_qr800, "medium trig, QR800", "l");
	leg1->Draw("same");

	hmega_qs->SetTitle("Gas - 1.9cm - comparing trigger type - alpha spectra");
	hmega_qs->GetXaxis()->SetTitle("Raw Charge (V x ns)");
	hmega_qs->GetYaxis()->SetTitle("Counts");

	c1->Update();

	// new TCanvas();

	// hmed_qr200->SetLineColor(1);
	// hmed_qr200->Draw("hist");

	// hmega_qr200->SetLineColor(2);
	// hmega_qr200->Draw("hist same");

	new TCanvas();
	hmega_h->SetLineColor(1);
	hmega_h->Draw("hist");

	hmed_h->SetLineColor(2);
	hmed_h->Draw("hist same");

	new TCanvas();
	TLegend * leg2 = new TLegend(0.85, 0.85, "Trigger", "ARC NDC");

	hmedB->SetLineColor(1);
	hmedB->Draw("hist");

	hmegaB->SetLineColor(2);
	hmegaB->Draw("hist same");

	leg2->AddEntry(hmedB, "medium", "l");
	leg2->AddEntry(hmegaB, "mega", "l");
	leg2->Draw("same");	

return 0;
}