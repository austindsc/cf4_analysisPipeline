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

int investStatic1Gas(){
  TFile * f4 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210723_0.root");
  TTree * tt4 = (TTree*)f4->Get("tt");

  TFile * f5 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210723_1.root");
  TTree * tt5 = (TTree*)f5->Get("tt");

  TFile * f6 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210723_2.root");
  TTree * tt6 = (TTree*)f6->Get("tt");

  speStuff v4 = getSPEval(tt4, true);
  speStuff v5 = getSPEval(tt5, true);
  speStuff v6 = getSPEval(tt6, true);

  // liquid data voltages
  v4.V = 48.25;
  v4.eV = 0.01;
  v5.V = 49.25;
  v5.eV = 0.01;
  v6.V = 50.25;
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

  // can't get SPE from gas data 
  // 1. digitizer is 1V
  // will assume actual SiPM temp also unchanged (no pressure reading either)

  //////// large trigger stuff ////////
  // load large trigger data
  TFile * fga4 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210723_0.root");
  TTree * ttga4 = (TTree*)fga5->Get("tt");

  TFile * fga5 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210723_1.root");
  TTree * ttga4 = (TTree*)fga5->Get("tt");

  TFile * fga6 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210723_2.root");
  TTree * ttga6 = (TTree*)fga6->Get("tt");

  /// drawing ///
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

  /*
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

	return 0;
}