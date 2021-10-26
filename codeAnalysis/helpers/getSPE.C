#include <stdio.h>
using namespace std;
using namespace TMath;

/*
1. load in 3 voltage files

2. get their SPE

3. make Q vs V plot, fit, get VBR



*/

struct speStuff{
  speStuff(){
    this->spe = -1;
    this->espe = -1;
  }

  speStuff(float a, float b){
    spe = a;
    espe = b;
  }

  float spe;
  float espe;
  float V;
  float eV;
};

speStuff getSPEval(TTree * tt, float width = 0.2, bool verbose = false){
  const char * dOpt = "goff";
  const char * fOpt = "QRN0M";
  if (verbose == true) { dOpt = ""; fOpt = "+RM same"; }

  Double_t gausWidth = 0.2;
  float findPeak1Factor = 1.4;

  TH1F * histQ = new TH1F("histQ", "histQ", 200, 0, 10);
  histQ->GetXaxis()->SetTitle("charge (V x ns)");
  histQ->SetTitle("raw charge histograms");
  tt->Draw(string("qr200 >> histQ").c_str(), "", dOpt);
  
  TSpectrum * s = new TSpectrum(4);
  s->Search(histQ, 2., dOpt, 0.05);
  Double_t * xPeaks = s->GetPositionX();
  Double_t * yPeaks = s->GetPositionY();

  int n = 4;
  sort(xPeaks, xPeaks+n);

  Double_t * xP1 = xPeaks;
  // dealing with noise
  // if (xPeaks[0] < 0.3){
  //   xP1 = xPeaks + 1;
  // }
  if(histQ->GetBinContent(histQ->GetXaxis()->FindBin(xPeaks[0])) * findPeak1Factor  < histQ->GetMaximum()
    || xPeaks[0] < 0.4){
    // cout << *xP1 << endl;
    xP1 = xPeaks + 1;
  }
  if(*xP1 < 0.4) xP1 = xPeaks + 2;

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

  float SPE = q2PE - q1PE;

  float eSPE = Sqrt(Power(fgaus2->GetParError(1), 2) + Power(fgaus->GetParError(1), 2));
  
  printf("chi2 = %0.2f   QSPE = %0.2f\n", fgaus->GetChisquare()/fgaus->GetNDF(), SPE);

  speStuff val1(SPE, eSPE);

  if(verbose == true){
    printf("qpe 1 = %0.2f\n", q1PE);
    printf("qpe 2 = %0.2f\n", q2PE);
    printf("qpe 3 = %0.2f\n", q3PE);
    gStyle->SetOptFit(1);
  } else {
    delete fgaus;
    delete histQ;
  }
  return val1;
}