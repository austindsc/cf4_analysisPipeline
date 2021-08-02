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

void compLiq2(){


    string path = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/proc2/";

    TH1 *hPo, *hBkg, *hCs, *hCsPo;
    TVectorD* v;
    TFile *f;

    /////////// loading data 1 ///////////
    // 51V Po
    string nbkg = "liqBig_bkg_53_1205.root";
    string nPo = "liqBig_Po_53_1117.root";
    string nCs = "liqBig_Cs_53_1205.root";
    string nCsPo = "liqBig_CsPo_53_1117.root";

    f = TFile::Open((path + nbkg).c_str(), "READ");
    hBkg = (TH1*) f->Get( "hF_cuts" )->Clone(); // Po charge cut
    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t ltBkg = (*v)[0];

    f = TFile::Open((path + nPo).c_str(), "READ");
    hPo = (TH1*) f->Get( "hF_cuts" )->Clone(); // Po charge cut
    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t ltPo = (*v)[0];

    f = TFile::Open((path + nCs).c_str(), "READ");
    hCs = (TH1*) f->Get( "hF_cuts" )->Clone(); // Po charge cut
    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t ltCs = (*v)[0];

    f = TFile::Open((path + nCsPo).c_str(), "READ");
    hCsPo = (TH1*) f->Get( "hF_cuts" )->Clone(); // Po charge cut
    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t ltCsPo = (*v)[0];
    
    // cout << hBkg->GetEntries() << endl;
    // cout << hCs->GetEntries() << endl;
    // cout << hPo->GetEntries() << endl;
    // cout << hCsPo->GetEntries() << endl;

    int rebin = 10;
    hBkg->Rebin(rebin);
    hPo->Rebin(rebin);
    hCs->Rebin(rebin);
    hCsPo->Rebin(rebin);

    hBkg->Scale(1/ltBkg); // y-scale is in Hz
    hPo->Scale(1/ltPo);   // y-scale is Hz
    hCs->Scale(1/ltCs);
    hCsPo->Scale(1/ltCsPo);

    float firstBin=1.;

    float aBkg = hBkg->Integral(hBkg->GetXaxis()->FindBin(firstBin), 1000);

    float aPo = hPo->Integral(hPo->GetXaxis()->FindBin(firstBin), 1000) - aBkg;
    float aCs = hCs->Integral(hCs->GetXaxis()->FindBin(firstBin), 1000) - aBkg;
    float aCsPo = hCsPo->Integral(hCsPo->GetXaxis()->FindBin(firstBin), 1000) - aBkg;

    cout << "Po Activity (Hz) = " << aPo << endl;
    cout << "Cs Activity (Hz) = " << aCs << endl;
    cout << "CsPo Activity (Hz) = " << aCsPo << endl;


    // // 2d histograms
    // h2 = (TH1*) fc_51->Get( "h2_hq" )->Clone();
    // k2 = (TH1*) fc_51->Get( "h2_hq_cuts" )->Clone();

    // // counting peak size
    // float NPo_51 = k1->Integral(k1->GetXaxis()->FindBin(10.), k1->GetNbinsX());

    // cout << NPo_51 << endl;

    // cout << "Peak Activity (Hz) = " << NPo_51/lifetime51Po << endl;

    // // getting Cs rate
    // float uCu_to_Bq = 37000.;
    // float tau_cs = 30.08; // years
    // TF1* fExp = new TF1("fExp", "[0] * TMath::Exp(-x/[1])");
    // fExp -> SetParameters(uCu_to_Bq * 11.10, 30.08);

    // cout << fExp->Eval(45.4) << endl;; // 45 years later


    // fitting
    // TF1* FGaus = (TF1*) gROOT->FindObjectAny("FGaus");
    // if(!FGaus) TF1 * FGaus = new TF1("FGaus","[0]/sqrt(2*TMath::Pi())/[2]*exp(-(x-[1])*(x-[1])/2/[2]/[2])");



    // k1 peak energy 5.3MeV

    // k3 peak energy 0.66 MeV

    // getting the Polonium rate

    // lifetime corrections? pure Cs is 1 month later


    // hq_CsPo_51->Scale(1/hq_CsPo_51->Integral(hq_CsPo_51->GetXaxis()->FindBin(scaleLow), hq_CsPo_51->GetXaxis()->FindBin(scaleHigh)));
    // hq_CsPo_53->Scale(1/hq_CsPo_53->Integral(hq_CsPo_53->GetXaxis()->FindBin(scaleLow), hq_CsPo_53->GetXaxis()->FindBin(scaleHigh)));
    
    //////////// c1 ////////////
    TCanvas * c1 = new TCanvas("c1", "c1", 900, 700);

    // first canvas //
    c1->cd(1);

    auto leg1 = new TLegend();

    hBkg->GetYaxis()->SetRangeUser(0.01, 500);
    c1->SetLogy();

    hBkg->SetLineColor(1);
    hBkg->Draw("hist");

    hBkg->GetXaxis()->SetTitle("qPE (nVs)");
    hBkg->GetYaxis()->SetTitle("rate (Hz/bin)");

    hPo->SetLineColor(2);
    hPo->Draw("hist same");

    hCs->SetLineColor(3);
    hCs->Draw("hist same");

    hCsPo->SetLineColor(4);
    hCsPo->Draw("hist same");

    leg1->AddEntry(hBkg, "bkg - 1205", "L");
    leg1->AddEntry(hPo, "Po - 1117", "L");
    leg1->AddEntry(hCs, "Cs - 1205", "L");
    leg1->AddEntry(hCsPo, "Cs + Po - 1117", "L");
    leg1->Draw("same");

    c1->Update();

}











int main(){
    compLiq2();
    }