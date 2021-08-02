 

#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorT.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TFile.h"
#include "TMath.h"

#include "TApplication.h"

#include <stdio.h>
using namespace std;

// R__LOAD_LIBRARY(gsl);
// R__LOAD_LIBRARY(libgsl);
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>

/*
looking at spectra in Gas

! no pure Cs spectrum in gas!
*/

// Constants
Double_t fine = 1./137.;
Double_t me = 511.0;
Double_t Z = 83.;
Double_t A = 210.;
Double_t c = TMath::C();

Double_t qEnd = 1162.; // keV


// scaled full energy
Double_t getW(Double_t *x){
    Double_t xx = x[0];
    return xx/me + 1.;
}

// *x is kinetic energy in keV
// momentum in units of me * c
Double_t getP(Double_t *x){
    Double_t xx = x[0];
    //printf("xx energy = %f\n", xx);
    return TMath::Sqrt(TMath::Power(xx + me, 2) - me*me)/me;
}

// Nuclear Transitions
Double_t getCW(Double_t *x){
    Double_t W = getW(x);
    return 1. - 0.4523 * W + 0.056 * W*W;
}

//  scaled energy over momentum
Double_t getY(Double_t *x){ //is this right?
    return Z*fine * getW(x)/getP(x);
}
Double_t getY(Double_t *x, Double_t *par){
    return getY(x);
}

Double_t getFermi(Double_t *x){
    Double_t xx = x[0];

    Double_t gamma = TMath::Sqrt(1 - fine*fine * Z*Z);
    Double_t R = 0.0028*TMath::Power(A, 1/3) + 0.0063*TMath::Power(A, 2/3) - 0.017/A;
    Double_t Y = getY(x);
    Double_t p = getP(x);

    Double_t A = 4*TMath::Power(2*p*R, 2*(gamma - 1));

    Double_t B = TMath::Exp(TMath::Pi() * Y);
    
    gsl_sf_result Clog, Carg;
    gsl_sf_lngamma_complex_e(gamma, Y, &Clog, &Carg);
    Double_t C = TMath::Exp(Clog.val);
    
    Double_t D = TMath::Power(TMath::Gamma(1. + 2.*gamma), 2);

    Double_t val = A * B * C / D;

    // cout << "K = "  << *x << endl;
    // cout << "W = " << getW(x) << endl;
    // cout << "P = " << p << endl;
    // cout << "Y = " << Y << endl;

    // printf("A = %f\n", A);
    // printf("B = %f\n", B);
    // printf("C = %f\n", C);
    // printf("D = %f\n", D);
    // printf("Fermi = %f\n", val);
    
    return val;
}

Double_t getFermi(Double_t *x, Double_t *par){
    return getFermi(x);
}


// no  fermi function here
Double_t getSpectrum1(Double_t *x, Double_t *par){
    Double_t W = getW(x);
    Double_t Wmax = getW(&qEnd);
    return  getP(x) * W * TMath::Power(W - Wmax, 2) * getCW(x) + par[0];
}

// no  fermi function here
Double_t getSpectrumFull(Double_t *x, Double_t *par){
    Double_t W = getW(x);
    Double_t Wmax = getW(&qEnd);
    return  getP(x) * W * TMath::Power(W - Wmax, 2) * getCW(x) * getFermi(x);
}


void drawSpectrum(){

    Double_t en1 = 10;
    Double_t en2 = qEnd;


    // TF1 * fspec1 = new TF1("fspec1", getSpectrum1, 1, qEnd, 1);

    // TF1 * fFermi = new TF1("fFermi", getFermi, 1, qEnd, 0);
    // TF1 * fY = new TF1("fY", getY, 1, qEnd, 0);
    TF1 * fspecFull = new TF1("fspecFull", getSpectrumFull, 0, qEnd, 1);


    auto cTest = new TCanvas("cTest", "cTest", 400, 400);

    fspecFull->Draw("L*");
    //fspec1->Draw("L same");
    cTest->SetLogy();
    cTest->Update();

    return;
}

void data_ana(){

    // drawSpectrum();

    cout << "second one" << endl;
    // return;

    string path = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/proc2/";

    TH1 *hBkg, *hPo, *hCsPo;
    TH1 *kBkg, *kPo, *kCsPo;
    TH2 *k2_bkg, *k2_Po, *k2_CsPo;
    TFile * f;
    TVectorD* v;

    /////////// loading data 1 ///////////
    string nPo = "gasBig_Po_53_1123.root";
    string nCsPo = "gasBig_CsPo_53_1123.root";
    string nBkg = "gasBig_bkg_53_1203.root";

    f = TFile::Open((path + nPo).c_str(), "READ");
    kPo = (TH1*) f->Get( "hF_cuts" )->Clone(); // Po charge cut
    k2_Po = (TH2*) f->Get( "h2_hq" )->Clone();

    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t lifetimePo = (*v)[0];

    f = TFile::Open((path + nCsPo).c_str(), "READ");
    kCsPo = (TH1*) f->Get("hF_cuts")->Clone(); // CsPo charge cut
    k2_CsPo = (TH2*) f->Get( "h2_hq_cuts" )->Clone();

    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t lifetimeCsPo = (*v)[0];

    f = TFile::Open((path + nBkg).c_str(), "READ");
    kBkg = (TH1*) f->Get("hF_cuts")->Clone(); // CsPo charge cut
    hBkg = (TH1*) f->Get("hF_all")->Clone(); // CsPo charge cut

    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t lifetimeBkg = (*v)[0];

    cout << kPo->GetEntries() << endl;
    cout << kCsPo->GetEntries() << endl;
    cout << kBkg->GetEntries() << "\n" << endl;

    int rebin = 10;

    kPo->Rebin(rebin);
    kCsPo->Rebin(rebin);
    kBkg->Rebin(rebin);
    hBkg->Rebin(rebin);

    kPo->Scale(1/lifetimePo);
    kCsPo->Scale(1/lifetimeCsPo);
    kBkg->Scale(1/lifetimeBkg);
    hBkg->Scale(1/lifetimeBkg);

    // printf("%f\n", lifetimePo);
    // printf("%f\n", lifetimeCsPo);
    // printf("%f\n", lifetimeBkg);

    // counting peak size
    float firstBin = 0.7;

    float a_raw_Bkg = hBkg->Integral( hBkg->GetXaxis()->FindBin(firstBin), hBkg->GetNbinsX() );
    float aBkg = kBkg->Integral( kBkg->GetXaxis()->FindBin(firstBin), kBkg->GetNbinsX() );
    float aPo_N = kPo->Integral(kPo->GetXaxis()->FindBin(10.), kPo->GetNbinsX());
    float aPo_sub = kPo->Integral(kPo->GetXaxis()->FindBin(firstBin), kPo->GetNbinsX()) - aBkg;
    float aCsPo_sub = kCsPo->Integral(kCsPo->GetXaxis()->FindBin(firstBin), kCsPo->GetNbinsX()) - aBkg;
    float aCs = aCsPo_sub - aPo_sub;

    cout << "uncut bkg rate (Hz) = " << a_raw_Bkg << endl;
    cout << "bkg rate (Hz) = " << aBkg << endl;
    cout << "Integral Peak Activity(Hz) = " << aPo_N << endl;
    cout << "sub, Po Activity (Hz) = " << aPo_sub << endl;
    cout << "sub, Po+Cs Activity (Hz) = " << aCsPo_sub << endl;
    cout << "sub sub, Cs Activity (Hz) = " << aCs << endl;

    // issue is this depends on DN (V_over) not changing (temperature)

    // temperatures are different by 3C
    float Vov_extra = 1* 3/20;
    float gainExtra = (53-48.5)*(1+Vov_extra);

    // getting Cs rate
    float uCu_to_Bq = 37000.;
    float tau_cs = 30.08; // years
    TF1* fExp = new TF1("fExp", "[0] * TMath::Exp(-x/[1])");
    fExp -> SetParameters(uCu_to_Bq * 11.10, 30.08);

    cout << fExp->Eval(45.4) << endl;; // 45 years later

    // fitting
    TF1* FGaus = (TF1*) gROOT->FindObjectAny("FGaus");
    if(!FGaus) TF1 * gausPo = new TF1("gausPo","[0]/sqrt(2*TMath::Pi())/[2]*exp(-(x-[1])*(x-[1])/2/[2]/[2])");

    TF1 * fgausPo = new TF1("fgausPo", "gaus");
    fgausPo->SetParameters(0.1, 80., 40);
    kPo->Fit("fgausPo",  "",  "", 60, 100);



    // look at peak shape vs 'leftover' rate, see if fraction is comparable to solid angle
    // (maybe long tail comes from 'shadowing')

    // Cs must contribute a couple PE

    // k1 peak energy 5.3MeV

    // k3 peak energy 0.66 MeV

    // getting the Polonium rate

    //////////// c2 ////////////
    // TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

    // auto leg1 = new TLegend();

    // // Po charge cut
    // kPo->SetLineColor(1); // black
    // // kPo->Draw("hist");

    // // kPo->GetXaxis()->SetTitle("qPE (nVs)");
    // // kPo->GetYaxis()->SetTitle("rate (Hz)");

    // // // CsPo charge cut
    // // kCsPo->SetLineColor(2); // red
    // // kCsPo->Draw("hist same");

    // // kBkg->SetLineColor(3); // green
    // // kBkg->Draw("hist same");

    // fgausPo->Draw();


    // leg1->AddEntry(kBkg, "bkg - 1203", "L");
    // leg1->AddEntry(kPo, "Polonium - 1123", "L");
    // leg1->AddEntry(kCsPo, "Cs + Po - 1123", "L");
    // leg1->Draw("same");

    // c1->SetGrid();
    // c1->SetLogy();
    // c1->Update();

    return;
}


int main(int argc, char** argv){

    TApplication theApp("App",&argc, argv);
    
    data_ana();

    drawSpectrum();

    theApp.Run();

    return 0;
    }