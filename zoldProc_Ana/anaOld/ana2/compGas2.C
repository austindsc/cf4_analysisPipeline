// Builds a graph with errors, displays it and saves it as
// image. First, include some header files
// (not necessary for Cling)

#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"

/*
looking at spectra in Gas

! no pure Cs spectrum in gas!
*/


void compGas2(){

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

    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t lifetimeBkg = (*v)[0];

    cout << kPo->GetEntries() << endl;
    cout << kCsPo->GetEntries() << endl;
    cout << kBkg->GetEntries() << "\n" << endl;

    int rebin = 5;

    kPo->Rebin(rebin);
    kCsPo->Rebin(rebin);
    kBkg->Rebin(rebin);

    kPo->Scale(1/lifetimePo);
    kCsPo->Scale(1/lifetimeCsPo);
    kBkg->Scale(1/lifetimeBkg);

    // printf("%f\n", lifetimePo);
    // printf("%f\n", lifetimeCsPo);
    // printf("%f\n", lifetimeBkg);

    // counting peak size
    float firstBin = 2;

    float aBkg = kBkg->Integral( kBkg->GetXaxis()->FindBin(firstBin), kBkg->GetNbinsX() );

    float aPo_N = kPo->Integral(kPo->GetXaxis()->FindBin(10.), kPo->GetNbinsX());

    float aPo_sub = kPo->Integral(kPo->GetXaxis()->FindBin(firstBin), kPo->GetNbinsX()) - aBkg;
    float aCsPo_sub = kCsPo->Integral(kCsPo->GetXaxis()->FindBin(firstBin), kCsPo->GetNbinsX()) - aBkg;
    float aCs = aCsPo_sub - aPo_sub;

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
    if(!FGaus) TF1 * FGaus = new TF1("FGaus","[0]/sqrt(2*TMath::Pi())/[2]*exp(-(x-[1])*(x-[1])/2/[2]/[2])");

    // look at peak shape vs 'leftover' rate, see if fraction is comparable to solid angle
    // (maybe long tail comes from 'shadowing')

    // Cs must contribute a couple PE

    // k1 peak energy 5.3MeV

    // k3 peak energy 0.66 MeV

    // getting the Polonium rate

    //////////// c2 ////////////
    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 800);

    auto leg1 = new TLegend();

    // Po charge cut
    kPo->SetLineColor(1); // black
    kPo->Draw("hist");

    kPo->GetXaxis()->SetTitle("qPE (nVs)");
    kPo->GetYaxis()->SetTitle("rate (Hz)");

    // CsPo charge cut
    kCsPo->SetLineColor(2); // red
    kCsPo->Draw("hist same");

    kBkg->SetLineColor(3); // green
    kBkg->Draw("hist same");

    leg1->AddEntry(kBkg, "bkg - 1203", "L");
    leg1->AddEntry(kPo, "Polonium - 1123", "L");
    leg1->AddEntry(kCsPo, "Cs + Po - 1123", "L");
    leg1->Draw("same");

    c1->SetGrid();
    c1->SetLogy();
    c1->Update();

}











int main(){
    compGas2();
    }