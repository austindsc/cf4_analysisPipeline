

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
#include "TPad.h"

#include "TApplication.h"

#include <stdio.h>
using namespace std;

// R__LOAD_LIBRARY(gsl);
// R__LOAD_LIBRARY(libgsl);

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
    string nPo = "gasBig_Po_51_1123.root";
    string nCsPo = "gasBig_CsPo_51_1123.root";
    string nBkg = "gasBig_bkg_51_1203.root";

    // need bkg data from same day as gasBig_Po (1123)
    // for more correct 1PE estimation (same temp?)

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

    /////////// manipulating ///////////
    // make a copy, for a finer 1pe fit
    TH1* kPo_spe = (TH1*) kPo->Clone();
    TH1* kBkg_spe = (TH1*) kBkg->Clone();

    int rebinSpe = 1;
    kBkg_spe->Rebin(rebinSpe); // xx test
    kPo_spe->Rebin(rebinSpe);

    kPo_spe->Scale(1/kPo_spe->Integral(0, kPo_spe->GetXaxis()->FindBin(15)));
    kBkg_spe->Scale(1/kBkg_spe->Integral(0, kBkg_spe->GetXaxis()->FindBin(15)));


    int rebin = 20;

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

    /////////// FITTING ///////////
    float rNPE_low = 0.5;
    float rNPE_high = 5.;
    float r1PE_high = 2.5;

    TF1* FGaus = (TF1*) gROOT->FindObjectAny("FGaus");
    if(!FGaus) TF1 * gausPo = new TF1("gausPo","[0]/sqrt(2*TMath::Pi())/[2]*exp(-(x-[1])*(x-[1])/2/[2]/[2])");

    // try to fit the bkg spectrum
    // 1 pe fit
    TF1 * fgBkg_1PE = new TF1("fgBkg_1PE", "gaus", rNPE_low, r1PE_high);
    fgBkg_1PE->SetParameters(kBkg_spe->GetMaximum(), 2, 1);
    kBkg_spe->Fit("fgBkg_1PE", "MR0", "");
    float speBkg_Height = fgBkg_1PE->GetParameter(0);
    float speBkg = fgBkg_1PE->GetParameter(1);
    float speBkg_width = fgBkg_1PE->GetParameter(2);

    //cout << "XX YOO speBkg = " << speBkg << endl;

    // fitting NPE together
    TF1 * fgBkg_NPE = new TF1("fgBkg_NPE", "gaus(0) + gaus(3) + gaus(6)", rNPE_low, rNPE_high);
    fgBkg_NPE->SetParNames("1height", "1PE", "1sigma",  "2 height", "2PE", "2sigma", "AP height", "AP mean", "AP sigma");

    // 1PE peak
    fgBkg_NPE->SetParameter(0, speBkg_Height);
    //fgBkg_NPE->SetParLimits(0, speBkg_Height, speBkg_Height);
    fgBkg_NPE->SetParameter(1, speBkg);
    fgBkg_NPE->SetParLimits(1, 0.8*speBkg, 1.1* speBkg);
    fgBkg_NPE->SetParameter(2, speBkg_width);
    //fgBkg_NPE->SetParLimits(2, speBkg_width, speBkg_width);

    // 2PE peak
    fgBkg_NPE->SetParameter(3, speBkg_Height/8);
    //fgBkg_NPE->SetParLimits(3, speBkg_Height/15, speBkg_Height/15);
    fgBkg_NPE->SetParameter(4, 2 * speBkg);
    fgBkg_NPE->SetParLimits(4, 1.5* speBkg, 2.5*speBkg);
    fgBkg_NPE->SetParameter(5, speBkg_width/2);
    
    // 1PE afterpulses
    fgBkg_NPE->SetParameter(6, speBkg_Height/8);
    fgBkg_NPE->SetParLimits(6, speBkg_Height/20, speBkg_Height/2);
    fgBkg_NPE->SetParameter(7, 1.2 * speBkg);
    fgBkg_NPE->SetParLimits(7, 1.* speBkg, 1.8*speBkg);
    fgBkg_NPE->SetParameter(8, speBkg_width);

    kBkg_spe->Fit("fgBkg_NPE", "MR0", "");

    cout << "1PE FIT = " << fgBkg_1PE->GetChisquare()/fgBkg_1PE->GetNDF() << endl;
    cout << "NPE FIT = " << fgBkg_NPE->GetChisquare()/fgBkg_NPE->GetNDF() << endl;

    ////////////////// fitting Po SPE data //////////////////

    // Po 1PE
    TF1 * fgPo_1PE = new TF1("fgPo_1PE", "gaus", 0.5, r1PE_high);
    fgPo_1PE->SetParameters(kPo_spe->GetMaximum(), 2, 0.25);
    kPo_spe->Fit("fgPo_1PE", "MR0", "");
    float spePo = fgPo_1PE->GetParameter(1);

    // NPE fit
    TF1 * fgPo_NPE = (TF1*) fgBkg_NPE->Clone();
    fgPo_NPE->SetRange(rNPE_low, rNPE_high);
    fgPo_NPE->SetName("fgPo_NPE");
    kPo_spe->Fit("fgPo_NPE", "MR0", "");

    cout << "1PE FIT = " << fgPo_1PE->GetChisquare()/fgPo_1PE->GetNDF() << endl;
    cout << "NPE FIT = " << fgPo_NPE->GetChisquare()/fgPo_NPE->GetNDF() << endl;

    // Po alpha
    TF1 * fgausPo = new TF1("fgausPo", "gaus", 60, 120);
    fgausPo->SetParameters(0.1, 90., 20);
    kPo->Fit("fgausPo",  "NMR0",  "");

    ////////////////// getting extra functions for plotting //////////////////
    TF1 * fgBkg_ex1 = new TF1("fgBkg_ex1", "gaus", 1, 5.);
    fgBkg_ex1->SetParameters(fgBkg_NPE->GetParameter(0), fgBkg_NPE->GetParameter(1),  fgBkg_NPE->GetParameter(2));
    TF1 * fgBkg_ex2 = new TF1("fgBkg_ex2", "gaus", 1, 7.);
    fgBkg_ex2->SetParameters(fgBkg_NPE->GetParameter(3), fgBkg_NPE->GetParameter(4),  fgBkg_NPE->GetParameter(5));
    TF1 * fgBkg_ex1AP = new TF1("fgBkg_ex1AP", "gaus", 1, 5.);
    fgBkg_ex1AP->SetParameters(fgBkg_NPE->GetParameter(6), fgBkg_NPE->GetParameter(7),  fgBkg_NPE->GetParameter(8));

    TF1 * fgPo_ex1 = new TF1("fgPo_ex1", "gaus", 1, 5.);
    fgPo_ex1->SetParameters(fgPo_NPE->GetParameter(0), fgPo_NPE->GetParameter(1),  fgPo_NPE->GetParameter(2));
    TF1 * fgPo_ex2 = new TF1("fgPo_ex2", "gaus", 1, 7.);
    fgPo_ex2->SetParameters(fgPo_NPE->GetParameter(3), fgPo_NPE->GetParameter(4),  fgPo_NPE->GetParameter(5));
    TF1 * fgPo_ex1AP = new TF1("fgPo_ex1AP", "gaus", 1, 5.);
    fgPo_ex1AP->SetParameters(fgPo_NPE->GetParameter(6), fgPo_NPE->GetParameter(7),  fgPo_NPE->GetParameter(8));

    ////////// SPE plots ////////////
    TCanvas *c2 = new TCanvas("c2", "c2", 850, 600);
    c2->Divide(1, 2);
    TVirtualPad *c2A, *c2B;
    c2A = c2->GetPad(1);
    c2B = c2->GetPad(2);
    c2A->SetBottomMargin(0.15);
    c2B->SetBottomMargin(0.15);

    c2->cd(1);
    kBkg_spe->SetName("full fit");
    kBkg_spe->SetTitle("Charge histogram - Bkg, 53V");
    kBkg_spe->GetXaxis()->SetTitle("charge (nsV)");
    kBkg_spe->GetXaxis()->SetTitleSize(0.06);
    kBkg_spe->GetYaxis()->SetTitle("normalized counts");

    kBkg_spe->GetXaxis()->SetRangeUser(-.5, 7.);
    kBkg_spe->SetFillColorAlpha(41, 0.75);
    kBkg_spe->Draw("hist");

    fgBkg_1PE->SetLineColor(1);
    fgBkg_1PE->Draw("same");

    fgBkg_NPE->SetLineColor(2);
    fgBkg_NPE->SetLineWidth(4);
    fgBkg_NPE->Draw("same");

    fgBkg_ex1->SetLineColor(3);
    fgBkg_ex1->Draw("same");
    fgBkg_ex2->SetLineColor(4);
    fgBkg_ex2->Draw("same");
    fgBkg_ex1AP->SetLineColor(6);
    fgBkg_ex1AP->Draw("same");

    auto leg2A = new TLegend();
    leg2A->AddEntry(fgBkg_1PE, "1PE only", "L");
    leg2A->AddEntry(fgBkg_NPE, "full: A+B+C", "L");
    leg2A->AddEntry(fgBkg_ex1, "A: 1PE", "L");
    leg2A->AddEntry(fgBkg_ex2, "B: 2PE", "L");
    leg2A->AddEntry(fgBkg_ex1AP, "C: 1PE CDA (delayed)", "L");
    leg2A->Draw("same");

    gStyle->SetOptFit();

    c2->cd(2);
    kPo_spe->SetName("full fit");
    kPo_spe->SetTitle("Charge histogram - Pb210, 53V");
    kPo_spe->GetXaxis()->SetTitle("charge (nsV)");
    kPo_spe->GetXaxis()->SetTitleSize(0.06);
    kPo_spe->GetYaxis()->SetTitle("normalized counts");

    kPo_spe->GetXaxis()->SetRangeUser(-.5, 7.);
    kPo_spe->SetFillColorAlpha(41, 0.75);
    kPo_spe->Draw("hist same");

    fgPo_1PE->SetLineColor(1);
    fgPo_1PE->Draw("same");

    fgPo_NPE->SetLineColor(2);
    fgPo_NPE->SetLineWidth(4);
    fgPo_NPE->Draw("same");

    fgPo_ex1->SetLineColor(3);
    fgPo_ex1->Draw("same");
    fgPo_ex2->SetLineColor(4);
    fgPo_ex2->Draw("same");
    fgPo_ex1AP->SetLineColor(6);
    fgPo_ex1AP->Draw("same");

    auto leg2B = new TLegend();
    leg2B->AddEntry(fgPo_1PE, "1PE only", "L");
    leg2B->AddEntry(fgPo_NPE, "full: A+B+C", "L");
    leg2B->AddEntry(fgPo_ex1, "A: 1PE", "L");
    leg2B->AddEntry(fgPo_ex2, "B: 2PE", "L");
    leg2B->AddEntry(fgPo_ex1AP, "C: 1PE CDA (delayed)", "L");
    leg2B->Draw("same");

    kPo_spe->GetXaxis()->SetTitleOffset(0.);
    kBkg_spe->GetXaxis()->SetTitleOffset(0.);

    c2A->Update();

    c2A->SetGrid();
    c2B->SetGrid();
    c2->Update();

/*
    ///// Polonium alpha stuff /////
    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    c3->Divide(1, 2);
    TVirtualPad *c3A, *c3B;
    c3A = c3->GetPad(1);
    c3->cd(1);

    // Po charge cut
    kPo->SetLineColor(1); // black
    kPo->Draw("hist");

    kPo->GetXaxis()->SetTitle("qPE (nVs)");
    kPo->GetYaxis()->SetTitle("rate (Hz)");

    // // CsPo charge cut
    kCsPo->SetLineColor(2); // red
    kCsPo->Draw("hist same");

    kBkg->SetLineColor(3); // green
    kBkg->SetLineWidth(2);
    kBkg->Draw("hist same");

    // functions
    fgausPo->SetLineWidth(5);
    fgausPo->Draw("same");

    auto leg2 = new TLegend();
    leg2->AddEntry(kBkg, "bkg - 1203", "L");
    leg2->AddEntry(kPo, "Polonium - 1123", "L");
    leg2->AddEntry(kCsPo, "Cs + Po - 1123", "L");
    leg2->Draw("same");

    c3A->SetGrid();
    c3A->SetLogy();
*/

    return;
}


int compGas2_betaFit(){
    
    data_ana();

    //drawSpectrum();

    return 0;
    }