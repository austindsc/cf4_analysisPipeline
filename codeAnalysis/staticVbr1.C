#include <stdio.h>
using namespace std;

/*
1. load in 3 voltage files

2. 



*/

Double_t getSPEval(int chan, TTree * tt, bool verbose = false){
  Double_t gausWidth = 40;
  TH1F * histQ = new TH1F("histQ", "histQ", 500, 0, 10000);
  tt->Draw(("ch" + to_string(chan) + ".pulse.charge >> histQ").c_str(), ("ch" + to_string(chan) + ".pulse.charge > 100").c_str(), "goff");
  
  TSpectrum * s = new TSpectrum(10);
  s->Search(histQ, 2., "goff", 0.05);
  Double_t * xPeaks = s->GetPositionX();
  TF1 * fgaus = new TF1("fgaus", "gaus", *xPeaks - gausWidth, *xPeaks + gausWidth);
  fgaus->SetParameters(histQ->GetMaximum(), *xPeaks, gausWidth/2);
  histQ->Fit("fgaus", "QRN0");
  Double_t qSPE = fgaus->GetParameter(1);

  printf("ch = %i  chi2 = %0.2f   QSPE = %0.2f\n", chan, fgaus->GetChisquare()/fgaus->GetNDF(), qSPE);
  if(verbose == true){
    gStyle->SetOptFit(1);
    histQ->Draw("hist");
    fgaus->Draw("same");
  } else {
    delete fgaus;
    delete histQ;
  }
  return qSPE;
}

int staticVbr1(){
	cout << "hi yo" << endl;

	return 0;
}