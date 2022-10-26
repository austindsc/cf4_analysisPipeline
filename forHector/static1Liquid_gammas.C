#include <stdio.h>
using namespace std;
using namespace TMath;

/*
1. load in 3 voltage files

2. get their SPE

3. make Q vs V plot, fit, get VBR


run in command line with "root -l nameOfMacro.C"

*/

// a custom object to store SPE information ()
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

// get SPE information from a ttree object
speStuff getSPEval(TTree * tt, bool verbose = false){
  const char * dOpt = "goff"; //don't draw
  const char * fOpt = "QRN0M"; //Quiet, use Range of function, No graphic, 0 - don't draw, M - improve fit
  if (verbose == true) { dOpt = ""; fOpt = "+RM same"; }

  TH1F * histQ = new TH1F("histQ", "histQ", 200, 0, 10);
  tt->Draw(string("qr300 >> histQ").c_str(), "", dOpt); // putting data into histQ TH1
  
  TSpectrum * s = new TSpectrum(4); // used to find the peaks in histQ
  s->Search(histQ, 2., dOpt, 0.05);
  Double_t * xPeaks = s->GetPositionX();
  int n = 4;
  sort(xPeaks, xPeaks+n); // sorting the peaks 

  Double_t * xP1 = xPeaks;
  if (xPeaks[0] < 0.5){ // dealing with noise, skip any peaks below 0.5PE
    xP1 = xPeaks + 1;
  }

  // defining gaussians to fit each peak with
  Double_t gausWidth = 0.36; //empirical value to help with fitting
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

  float SPE = q3PE - q2PE; // get one SPE value
  float eSPE = Sqrt(Power(fgaus3->GetParError(1), 2) + Power(fgaus2->GetParError(1), 2)); // error on SPE

  speStuff val1(SPE, eSPE); // object holding SPE information

  if(verbose == true){
    printf("%0.2f\n", q1PE);
    printf("%0.2f\n", q2PE);
    printf("%0.2f\n", q3PE);
    gStyle->SetOptFit(1);
  } else {
    delete fgaus;
    delete histQ;
  }
  return val1; // pass speStuff object
}

int static1Liquid_gammas(){

  ////////////////////// getting breakdown voltages //////////////////////
  TFile * f4 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210726_0.root");
  TTree * tt4 = (TTree*)f4->Get("tt"); //open file and get TTree

  TFile * f5 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210726_1.root");
  TTree * tt5 = (TTree*)f5->Get("tt");

  TFile * f6 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210726_2.root");
  TTree * tt6 = (TTree*)f6->Get("tt");

  speStuff v4 = getSPEval(tt4, true); //set to true, add 'return 0;' at next line to see SPE distribution
  speStuff v5 = getSPEval(tt5);
  speStuff v6 = getSPEval(tt6);

  // liquid data voltages, hardcoded for now
  v4.V = 48.25;
  v4.eV = 0.01;
  v5.V = 49.25;
  v5.eV = 0.01;
  v6.V = 50.26;
  v6.eV = 0.01;

  // makinng arrays to plot TGraphErrors
  Double_t aQr[3] = { v4.spe, v5.spe, v6.spe };
  Double_t aeQr[3] = { v4.espe, v5.espe, v6.espe };
  Double_t aV[3] = { v4.V, v5.V, v6.V };
  Double_t aeV[3] = { v4.eV, v5.eV, v6.eV };

  TGraphErrors * tgVbr_qr = new TGraphErrors(3, aV, aQr, aeV, aeQr);

  // making breakdown voltage plot
  TF1 * fLine = new TF1("fLine", "[1]*(x - [0])", 0, 55);
  TF1 * flineVbr_qr = new TF1("flineVbr_qr", "fLine", 0, 55);
  flineVbr_qr->SetParameters(48., 0.25);
  tgVbr_qr->Fit("flineVbr_qr", "QMR", "");

  gROOT->Reset();
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1);
  
  /*
  // drawing breakdown voltage graph
  auto c3 = new TCanvas();

  tgVbr_qr->SetTitle("SPE charge vs Voltage");
  tgVbr_qr->GetXaxis()->SetTitle("Applied voltage");
  tgVbr_qr->GetYaxis()->SetTitle("SPE charge");

  tgVbr_qr->GetXaxis()->SetLimits(40, 51);
  tgVbr_qr->GetXaxis()->SetRangeUser(40, 51);

  tgVbr_qr->GetYaxis()->SetLimits(0, 1.5);
  tgVbr_qr->GetYaxis()->SetRangeUser(0, 1.5);

  tgVbr_qr->Draw("APE");
  flineVbr_qr->Draw("same");

  c3->Update();
  */

  //// done getting breakdown voltages ////


  //////// large trigger stuff ////////

  // load large trigger data, liquid
  TFile * fla4 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210726_7.root");
  TTree * ttla4 = (TTree*)fla4->Get("tt");

  TFile * fla5 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210726_8.root");
  TTree * ttla5 = (TTree*)fla5->Get("tt");

  TFile * flag5 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210729_4.root");
  TTree * ttlag5 = (TTree*)flag5->Get("tt");

  TFile * flag5Co = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210729_6.root");
  TTree * ttlag5Co = (TTree*)flag5Co->Get("tt");

  TH1F * hlaN22 = new TH1F("hlaN22", "hlaN22", 150, 0, 1.); // "hist liquid alpha Na22 --> hlaN22"
  ttlag5->Draw("hbline_min:qs50 >> hlaN22", "", ""); // put into TH1 named hlaN22

  auto c5 = new TCanvas();
  TH1F * hla4 = new TH1F("hla4", "hla4", 150, 0, 300);
  hla4->SetTitle("Liquid alpha scintillation");
  hla4->GetXaxis()->SetTitle("NPE");
  ttla4->Draw(("qs50 * " + to_string(1/v4.spe) + " >> hla4").c_str(), "", ""); // draw each TTree's data and normalize by SPE value

  TH1F * hla5 = new TH1F("hla5", "hla5", 150, 0, 300);
  hla5->SetLineColor(3);
  ttla5->Draw(("qs50 * " + to_string(1/v5.spe) + " >> hla5").c_str(), "", "same");

  TH1F * hlag5 = new TH1F("hlag5", "hlag5", 150, 0, 300);
  hlag5->SetLineColor(6);
  ttlag5->Draw(("qs50 * " + to_string(1/v5.spe) + " >> hlag5").c_str(), "", "same");

  TH1F * hlag5Co = new TH1F("hlag5Co", "hlag5Co", 150, 0, 300);
  hlag5Co->SetLineColor(7);
  ttlag5Co->Draw(("qs50 * " + to_string(1/v5.spe) + " >> hlag5Co").c_str(), "", "same");

  auto leg5 = new TLegend();
  leg5->AddEntry(hla4, "4V, alpha", "L");
  leg5->AddEntry(hla5, "5V, alpha", "L");
  leg5->AddEntry(hlag5, "5V, gamma (Na22)", "L");
  leg5->AddEntry(hlag5Co, "5V, gamma (Co57)", "L");
  leg5->Draw("same");

  c5->SetLogy();
  c5->Update();


  /// drawing ///

	return 0;
}