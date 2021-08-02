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

void compLiq1(){


    string path = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/proc1/";
    
    string gas = "gas/";
    string liq = "liquid/big/";

    string nPo51 = "liqBig_51_Po_1117_x20.root";
    string nPo53 = "liqBig_53_Po_1117_x20.root";
    //string nPo53 = "gasBig_55_Po_1123_x20.root";

    string nBkg51 = "liqBig_51_bkg_1205_x20.root";
    string nBkg53 = "liqBig_53_bkg_1205_x20.root";

    string nCsPo51 = "liqBig_51_Cs_1117_x20.root";
    string nCsPo53 = "liqBig_53_Cs_1117_x20.root";

    string nCs51 = "liqBig_51_CsnoPo_1205_x20.root";
    string nCs53 = "liqBig_53_CsnoPo_1205_x20.root";

    //string nCsPo53 = "gasBig_55_Cs_1123_x20.root";

    // data data
    TFile* f_Po_51 = TFile::Open((path + liq + nPo51).c_str(), "READ");
    TH1* hq_Po_51 = (TH1*) ((TH1*) f_Po_51->Get("hF_cuts")->Clone() );
    TFile* f_Po_53 = TFile::Open((path + liq + nPo53).c_str(), "READ");
    TH1* hq_Po_53 = (TH1*) ((TH1*) f_Po_53->Get("hF_cuts")->Clone() );

    TFile* f_bkg_51 = TFile::Open((path + liq + nBkg51).c_str(), "READ");
    TH1* hq_bkg_51 = (TH1*) ((TH1*) f_bkg_51->Get("hF_cuts")->Clone() );
    TFile* f_bkg_53 = TFile::Open((path + liq + nBkg53).c_str(), "READ");
    TH1* hq_bkg_53 = (TH1*) ((TH1*) f_bkg_53->Get("hF_cuts")->Clone() );

    TFile* f_Cs_51 = TFile::Open((path + liq + nCs51).c_str(), "READ");
    TH1* hq_Cs_51 = (TH1*) ((TH1*) f_Cs_51->Get("hF_cuts")->Clone() );
    TFile* f_Cs_53 = TFile::Open((path + liq + nCs53).c_str(), "READ");
    TH1* hq_Cs_53 = (TH1*) ((TH1*) f_Cs_53->Get("hF_cuts")->Clone() );

    TFile* f_CsPo_51 = TFile::Open((path + liq + nCsPo51).c_str(), "READ");
    TH1* hq_CsPo_51 = (TH1*) ((TH1*) f_CsPo_51->Get("hF_cuts")->Clone() );
    TFile* f_CsPo_53 = TFile::Open((path + liq + nCsPo53).c_str(), "READ");
    TH1* hq_CsPo_53 = (TH1*) ((TH1*) f_CsPo_53->Get("hF_cuts")->Clone() );

    int rebin = 40;
    hq_Po_51->Rebin(rebin);
    hq_Po_53->Rebin(rebin);

    hq_bkg_51->Rebin(rebin);
    hq_bkg_53->Rebin(rebin);

    hq_Cs_51->Rebin(rebin);
    hq_Cs_53->Rebin(rebin);

    hq_CsPo_51->Rebin(rebin);
    hq_CsPo_53->Rebin(rebin);

    TH1* allHists[8] = {hq_Po_51, hq_Po_53, hq_bkg_51, hq_bkg_53, hq_Cs_51, hq_Cs_53, hq_CsPo_51, hq_CsPo_53};
    for (int i = 0; i < 8; ++i){
        allHists[i]->Scale(1/allHists[i]->GetEntries());
    }

    // // scale from 10 onwards
    // float scaleLow = 10;
    // float scaleHigh = 300;

    // hq_Po_51->Scale(1/hq_Po_51->Integral(hq_Po_51->GetXaxis()->FindBin(scaleLow), hq_Po_51->GetXaxis()->FindBin(scaleHigh)));
    // hq_Po_53->Scale(1/hq_Po_53->Integral(hq_Po_53->GetXaxis()->FindBin(scaleLow), hq_Po_53->GetXaxis()->FindBin(scaleHigh)));
    
    // hq_CsPo_51->Scale(1/hq_CsPo_51->Integral(hq_CsPo_51->GetXaxis()->FindBin(scaleLow), hq_CsPo_51->GetXaxis()->FindBin(scaleHigh)));
    // hq_CsPo_53->Scale(1/hq_CsPo_53->Integral(hq_CsPo_53->GetXaxis()->FindBin(scaleLow), hq_CsPo_53->GetXaxis()->FindBin(scaleHigh)));
    

    //////////// c1 ////////////
    TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
    c1->Divide(1, 2);

    // first canvas //
    c1->cd(1);
    auto leg1 = new TLegend();

    hq_Po_51->SetLineColor(2);
    hq_Po_51->Draw("hist");

    hq_bkg_51->SetLineColor(6);
    hq_bkg_51->Draw("hist same");
    
    hq_Cs_51->SetLineColor(3);
    hq_Cs_51->Draw("hist same");

    hq_CsPo_51->SetLineColor(4);
    hq_CsPo_51->Draw("hist same");

    leg1->AddEntry(hq_bkg_51, "bkg", "L");
    leg1->AddEntry(hq_Po_51, "Po", "L");
    leg1->AddEntry(hq_Cs_51, "Cs", "L");
    leg1->AddEntry(hq_CsPo_51, "Cs + Po", "L");

    leg1->Draw("same");
    c1->cd(1)->SetLogy();
    //hq_Po_51->GetYaxis()->SetRangeUser(1, max(hq_Po_51->GetMaximum(), hq_bkg_51->GetMaximum()));
    //hq_Po_51->GetYaxis()->SetRangeUser(0.0001, 1);

    // second canvas //
    c1->cd(2);
    auto leg2 = new TLegend();

    hq_Po_53->SetLineColor(2);
    hq_Po_53->Draw("hist");

    hq_bkg_53->SetLineColor(6);
    hq_bkg_53->Draw("hist same");
    
    hq_Cs_53->SetLineColor(3);
    hq_Cs_53->Draw("hist same");

    hq_CsPo_53->SetLineColor(4);
    hq_CsPo_53->Draw("hist same");

    leg2->AddEntry(hq_bkg_53, "bkg", "L");
    leg2->AddEntry(hq_Po_53, "Po", "L");
    leg2->AddEntry(hq_Cs_53, "Cs", "L");
    leg2->AddEntry(hq_CsPo_53, "Cs + Po", "L");

    //hq_Po_53->GetYaxis()->SetRangeUser(0.0001, 1);

    //hq_Po_53->GetYaxis()->SetRangeUser(1, max(hq_Po_53->GetMaximum(), hq_bkg_53->GetMaximum()));
    leg2->Draw("same");

    c1->cd(2)->SetLogy();
    c1->SetGrid();

    c1->Update();
}

int main(){
    compLiq1();
    }