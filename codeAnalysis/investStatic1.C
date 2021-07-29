#include <stdio.h>
using namespace std;
using namespace TMath;

/*
1. load in 3 voltage files

2. get their SPE

3. make Q vs V plot, fit, get VBR



*/

struct speStuff{
  speStuff(float a, float b){
    spe = a;
    espe = b;
  }

  float spe;
  float espe;
  float V;
  float eV;
};

speStuff getSPEval(TTree * tt, bool verbose = false){
  const char * dOpt = "goff";
  const char * fOpt = "QRN0M";
  if (verbose == true) { dOpt = ""; fOpt = "+RM same"; }

  Double_t gausWidth = 0.36;
  TH1F * histQ = new TH1F("histQ", "histQ", 200, 0, 10);
  tt->Draw(string("qr300 >> histQ").c_str(), "", dOpt);
  
  TSpectrum * s = new TSpectrum(4);
  s->Search(histQ, 2., dOpt, 0.05);
  Double_t * xPeaks = s->GetPositionX();
  int n = 4;
  sort(xPeaks, xPeaks+n);

  Double_t * xP1 = xPeaks;
  // dealing with noise
  if (xPeaks[0] < 0.5){
    xP1 = xPeaks + 1;
  }

  TF1 * fgaus = new TF1("fgaus", "gaus", *xP1 - gausWidth, *xP1 + gausWidth);
  fgaus->SetParameters(histQ->GetMaximum(), *xP1, gausWidth/2);
  histQ->Fit("fgaus", fOpt);
  Double_t q1PE = fgaus->GetParameter(1);

  TF1 * fgaus2 = new TF1("fgaus2", "gaus", *(xP1+1) - gausWidth, *(xP1+1) + gausWidth);
  fgaus2->SetParameters(histQ->GetMaximum(), *(xP1+1), gausWidth/2);
  histQ->Fit("fgaus2", fOpt);
  Double_t q2PE = fgaus2->GetParameter(1);

  TF1 * fgaus3 = new TF1("fgaus3", "gaus", *(xP1+2) - gausWidth, *(xP1+2) + gausWidth);
  fgaus3->SetParameters(histQ->GetMaximum(), *(xP1+2), gausWidth/2);
  histQ->Fit("fgaus3", fOpt);
  Double_t q3PE = fgaus3->GetParameter(1);

  printf("chi2 = %0.2f   QSPE = %0.2f\n", fgaus->GetChisquare()/fgaus->GetNDF(), q1PE);

  float SPE = q3PE - q2PE;

  float eSPE = Sqrt(Power(fgaus3->GetParError(1), 2) + Power(fgaus2->GetParError(1), 2));

  speStuff val1(SPE, eSPE);

  if(verbose == true){
    printf("%0.2f\n", q1PE);
    printf("%0.2f\n", q2PE);
    printf("%0.2f\n", q3PE);
    gStyle->SetOptFit(1);
  } else {
    delete fgaus;
    delete histQ;
  }
  return val1;
}

int investStatic1(){
  TFile * f4 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210726_0.root");
  TTree * tt4 = (TTree*)f4->Get("tt");

  TFile * f5 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210726_1.root");
  TTree * tt5 = (TTree*)f5->Get("tt");

  TFile * f6 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210726_2.root");
  TTree * tt6 = (TTree*)f6->Get("tt");

  speStuff v4 = getSPEval(tt4);
  speStuff v5 = getSPEval(tt5);
  speStuff v6 = getSPEval(tt6);

  // liquid data voltages
  v4.V = 48.25;
  v4.eV = 0.01;
  v5.V = 49.25;
  v5.eV = 0.01;
  v6.V = 50.26;
  v6.eV = 0.01;

  Double_t aQr[3] = { v4.spe, v5.spe, v6.spe };
  Double_t aeQr[3] = { v4.espe, v5.espe, v6.espe };
  Double_t aV[3] = { v4.V, v5.V, v6.V };
  Double_t aeV[3] = { v4.eV, v5.eV, v6.eV };

  TGraphErrors * tgVbr_qr = new TGraphErrors(3, aV, aQr, aeV, aeQr);

  TF1 * fLine = new TF1("fLine", "[1]*(x - [0])", 0, 55);
  TF1 * flineVbr_qr = new TF1("flineVbr_qr", "fLine", 0, 55);
  flineVbr_qr->SetParameters(48., 0.25);
  tgVbr_qr->Fit("flineVbr_qr", "MR", "");

  /*
  gROOT->Reset();
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1);
  auto c3 = new TCanvas();

  tgVbr_qr->SetTitle("SPE charge vs Voltage");
  tgVbr_qr->GetXaxis()->SetTitle("Applied voltage");
  tgVbr_qr->GetYaxis()->SetTitle("SPE charge");

  tgVbr_qr->GetXaxis()->SetLimits(40, 51);
  tgVbr_qr->GetXaxis()->SetRangeUser(40, 51);

  tgVbr_qr->GetYaxis()->SetLimits(0, 1.5);
  tgVbr_qr->GetYaxis()->SetRangeUser(0, 1.5);

  tgVbr_qr->Draw("APE");
  //flineVbr_qr->Draw("same");

  c3->Update();

  auto c2 = new TCanvas();
  hqr5->SetTitle("Comparing charge variables");
  hqr5->GetXaxis()->SetTitle("charge (Vns)");
  hqr5->Draw("hist");

  hqs5->SetLineColor(3);
  hqs5->Draw("hist same");

  hqf5->SetLineColor(2);
  hqf5->Draw("hist same");

  auto leg2 = new TLegend();
  leg2->AddEntry(hqr5, "fixed range Q", "L");
  leg2->AddEntry(hqs5, "'smart' range Q", "L");
  leg2->AddEntry(hqf5, "fit Q", "L");
  leg2->Draw("same");

  c2->SetGrid(1);
  c2->Update();
  */

  //////// large trigger stuff ////////
  // load large trigger data, gas
  TFile * fga4 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210723_0.root");
  TTree * ttga4 = (TTree*)fga4->Get("tt");

  TFile * fga5 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210723_1.root");
  TTree * ttga5 = (TTree*)fga5->Get("tt");

  TFile * fga6 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210723_2.root");
  TTree * ttga6 = (TTree*)fga6->Get("tt");

  auto c3 = new TCanvas();
  TH1F * hga4 = new TH1F("hga4", "hga4", 150, 0, 500);
  hga4->SetTitle("Gas alpha scintillation");
  hga4->GetXaxis()->SetTitle("NPE (with liquid SPE values)");
  ttga4->Draw(("qs50 * " + to_string(1/v4.spe) + " >> hga4").c_str(), "", "");

  TF1 * fitga4 = new TF1("fitga4", "gaus", 120, 200);
  hga4->Fit("fitga4", "+MR", "");
  cout << fitga4->GetParameter(1)/fitga4->GetParameter(2) << endl;

  TH1F * hga5 = new TH1F("hga5", "hga5", 150, 0, 500);
  hga5->SetLineColor(3);
  ttga5->Draw(("qs50 * " + to_string(1/v5.spe) + " >> hga5").c_str(), "", "same");
  TF1 * fitga5 = new TF1("fitga5", "gaus", 140, 220);
  hga5->Fit("fitga5", "+MR", "");
  cout << fitga5->GetParameter(1)/fitga5->GetParameter(2) << endl;

  TH1F * hga6 = new TH1F("hga6", "hga6", 150, 0, 500);
  hga6->SetLineColor(6);
  ttga6->Draw(("qs50 * " + to_string(1/v6.spe) + " >> hga6").c_str(), "", "same");
  TF1 * fitga6 = new TF1("fitga6", "gaus", 160, 260);
  hga6->Fit("fitga6", "+MR", "");
  cout << fitga6->GetParameter(1)/fitga6->GetParameter(2) << endl;

  auto leg3 = new TLegend();
  leg3->AddEntry(hga4, "4V", "L");
  leg3->AddEntry(hga5, "5V", "L");
  leg3->AddEntry(hga6, "6V", "L");
  leg3->Draw("same");
  c3->Update();

  // load large trigger data, liquid
  TFile * fla4 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210726_7.root");
  TTree * ttla4 = (TTree*)fla4->Get("tt");

  TFile * fla5 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210726_8.root");
  TTree * ttla5 = (TTree*)fla5->Get("tt");

  TFile * flag5 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210727_2.root");
  TTree * ttlag5 = (TTree*)flag5->Get("tt");

  auto c4 = new TCanvas();
  TH1F * hla4 = new TH1F("hla4", "hla4", 150, 0, 50);
  hla4->SetTitle("Liquid alpha scintillation");
  hla4->GetXaxis()->SetTitle("NPE");
  ttla4->Draw(("qs50 * " + to_string(1/v4.spe) + " >> hla4").c_str(), "", "");

  TH1F * hla5 = new TH1F("hla5", "hla5", 150, 0, 50);
  hla5->SetLineColor(3);
  ttla5->Draw(("qs50 * " + to_string(1/v5.spe) + " >> hla5").c_str(), "", "same");

  TH1F * hlag5 = new TH1F("hlag5", "hlag5", 150, 0, 100);
  hlag5->SetLineColor(6);
  ttlag5->Draw(("qs50 * " + to_string(1/v5.spe) + " >> hlag5").c_str(), "", "same");

  auto leg4 = new TLegend();
  leg4->AddEntry(hla4, "4V, alpha", "L");
  leg4->AddEntry(hla5, "5V, alpha", "L");
  leg4->AddEntry(hga6, "5V, gamma (Na22)", "L");
  leg4->Draw("same");

  c4->SetLogy();
  c4->Update();

  /// drawing ///


	return 0;
}