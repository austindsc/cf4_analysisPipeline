// Builds a graph with errors, displays it and saves it as
// image. First, include some header files
// (not necessary for Cling)

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
#include "TPaveText.h"

#include "TApplication.h"

#include <stdio.h>
using namespace std;

/*
1. get Po events in gas
2. 
*/

void getEff(){

    // new processed data
    string path = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/proc2/";    
    string n53 = "gasBig_Po_53_1123.root";
    string nL_CsPo_53 = "liqBig_CsPo_53_1117.root";
    string nL_Po_53 = "liqBig_Po_53_1117.root";

    TH1 *hq53, *hqL_Po_53, *hqL_CsPo_53;
    TVectorD* v;
    TFile* f;

    f = TFile::Open((path + n53).c_str(), "READ");
    hq53 = (TH1*) ((TH1*) f->Get("hF_cuts")->Clone() );
    hq53->SetName(n53.c_str());
    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t ltGPo = (*v)[0];

    f = TFile::Open((path + nL_CsPo_53).c_str(), "READ");
    hqL_CsPo_53 = (TH1*) ((TH1*) f->Get("hF_cuts")->Clone() );
    hqL_CsPo_53->SetName(nL_CsPo_53.c_str());
    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t ltLCsPo = (*v)[0];

    f = TFile::Open((path + nL_Po_53).c_str(), "READ");
    hqL_Po_53 = (TH1*) ((TH1*) f->Get("hF_cuts")->Clone() );
    hqL_Po_53->SetName(nL_Po_53.c_str());
    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t ltLPo = (*v)[0];

    // spe analysis
    TH1 * hq53_spe = (TH1*) hq53->Clone();
    hq53_spe->SetName("hq53_spe");

    TF1 * fg53_spe = new TF1("fg53_spe", "gaus", 0.6, 2.5);
    hq53_spe->Fit("fg53_spe", "MR0", "");
    float qspe = fg53_spe->GetParameter(1);
    
    TH1 * hqL_Cs_spe = (TH1*) hqL_CsPo_53->Clone();
    hqL_Cs_spe->SetName("hqL_Cs_spe");

    TF1 * fgCs_spe = new TF1("fgCs_spe", "gaus", 0.6, 2.5);
    hqL_Cs_spe->Fit("fgCs_spe", "MR0", "");
    float qspeCs = fgCs_spe->GetParameter(1);
    
    cout << "HELLO CS spe = " << qspeCs << endl;

    // manipulatinng
    int rebin = 40;
    hq53->Rebin(rebin);
    hqL_Po_53->Rebin(rebin);
    hqL_CsPo_53->Rebin(rebin);

    hq53->Scale(1/ltGPo);
    hqL_Po_53->Scale(1/ltLPo);
    hqL_CsPo_53->Scale(1/ltLCsPo);

    TF1 * fg53 = new TF1("fg53", "gaus", 60, 130);
    hq53->Fit("fg53", "MR0", "");
    float Ndet = fg53->GetParameter(1)/qspe;

    // estimate total  yield //
    float NexM = 5.3*3000;
    float NexL = 5.3*225;

    float effM = Ndet/NexM;
    float effL = Ndet/NexL;

    cout << "efficiency M = " << effM << endl;
    cout << "efficiency L = " << effL << endl;

    TH1 * liqSub = (TH1*) hqL_CsPo_53->Clone();
    liqSub->SetName("hq_Cs_sub_Po");

    liqSub->Add(hqL_Po_53, -1.);

    //////////// c1 ////////////
    TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
    c1->cd();

    gStyle->SetOptFit(1);

    hq53->SetTitle("Gas CF4 - Pb210 Source");
    hq53->GetXaxis()->SetTitle("Charge (nVs)");

    hq53->SetLineColor(4);
    hq53->SetLineWidth(3);
    hq53->Draw("hist");

    fg53->SetLineWidth(3);
    fg53->Draw("same");

    c1->SetLogy();
    c1->SetGrid();
    c1->Update();

        //////////// c1 ////////////
    // TCanvas * c2 = new TCanvas("c2", "c2", 600, 600);
    // c2->cd();

    // gStyle->SetOptFit(1);

    // hq53_spe->SetTitle("Gas CF4 - Pb210 Source");
    // hq53_spe->GetXaxis()->SetTitle("Charge (nVs)");

    // hq53_spe->Draw("hist");
    // fg53_spe->Draw("same");

    // c2->SetLogy();
    // c2->SetGrid();
    // c2->Update();

    TCanvas * c3 = new TCanvas("c3", "c3", 600, 600);
    c3->cd();

    gStyle->SetOptFit(1);

    hqL_Po_53->SetTitle("53V data, Cs and Pb sources");
    hqL_Po_53->GetXaxis()->SetTitle("charge (nVs)");

    hqL_Po_53->SetLineColor(2);
    hqL_Po_53->Draw("hist");
    hqL_CsPo_53->Draw("hist same");

    auto leg3 = new TLegend();
    leg3->AddEntry(hqL_Po_53, "Pb source", "L");
    leg3->AddEntry(hqL_CsPo_53, "Pb + Cs source", "L");
    leg3->Draw("same");

    // liqSub->SetLineColor(3);
    // liqSub->SetLineWidth(3);
    // liqSub->Draw("hist same");

    c3->SetLogy();
    c3->SetGrid();
    c3->Update();
}

int ly1_rough(){
    getEff();

    return 0;
    }