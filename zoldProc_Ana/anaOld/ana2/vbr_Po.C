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

void vbr_Po(){


    // string path = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/proc1/";
    
    // string fold = "gas/";
    // string n51 = "gasBig_51_Po_1123_x20.root";
    // string n53 = "gasBig_53_Po_1123_x20.root";
    // string n55 = "gasBig_55_Po_1123_x20.root";

    // new processed data
    string path = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/proc2/";    
    string fold = "";
    string n51 = "gasBig_Po_51_1123.root";
    string n53 = "gasBig_Po_53_1123.root";
    string n55 = "gasBig_Po_55_1123.root";

    // data data
    TFile* f_Po_51 = TFile::Open((path + fold + n51).c_str(), "READ");
    TH1* hq51 = (TH1*) ((TH1*) f_Po_51->Get("hF_cuts")->Clone() );
    hq51->SetName(n51.c_str());

    TFile* f_Po_53 = TFile::Open((path + fold + n53).c_str(), "READ");
    TH1* hq53 = (TH1*) ((TH1*) f_Po_53->Get("hF_cuts")->Clone() );
    hq53->SetName(n53.c_str());

    TFile* f_Po_55 = TFile::Open((path + fold + n55).c_str(), "READ");
    TH1* hq55 = (TH1*) ((TH1*) f_Po_55->Get("hF_cuts")->Clone() );
    hq55->SetName(n55.c_str());

    // rebinnned items
    TH1 *hqr51, *hqr53, *hqr55;
    hqr51 = (TH1*) hq51->Clone();
    hqr53 = (TH1*) hq53->Clone();
    hqr55 = (TH1*) hq55->Clone();

    int rebin = 40;
    hqr51->Rebin(rebin);
    hqr53->Rebin(rebin);
    hqr55->Rebin(rebin);

    // charges in array
    Double_t qx[3] = {51., 53., 55.};
    Double_t qex[3] = {0.1, 0.1, 0.1};
    Double_t qy[3];
    Double_t qey[3];

    // gaussian fitting
    TF1 * fg51 = new TF1("fg51", "gaus", 30, 70);
    TFitResultPtr r51 = hqr51->Fit("fg51", "SMR0", "");
    qy[0] = *(r51->GetParams()+1);
    qey[0] = *(r51->GetErrors()+1);

    TF1 * fg53 = new TF1("fg53", "gaus", 60, 130);
    TFitResultPtr r53 = hqr53->Fit("fg53", "SMR0", "");
    qy[1] = *(r53->GetParams()+1);
    qey[1] = *(r53->GetErrors()+1);

    TF1 * fg55 = new TF1("fg55", "gaus", 100, 200);
    TFitResultPtr r55 = hqr55->Fit("fg55", "SMR0", "");
    qy[2] = *(r55->GetParams()+1);
    qey[2] = *(r55->GetErrors()+1);

    // cout << *r51->GetParams() << endl;
    // cout << *(r51->GetParams()+1) << endl;
    // cout << *(r51->GetErrors()+1) << endl;

    TGraphErrors * gr_q1 = new TGraphErrors(3, qx, qy, qex, qey);
    TGraphErrors * gr_q1Clone = (TGraphErrors *) gr_q1->Clone();

    TF1 * fLine = new TF1("fLine", "[1]*(x - [0])", 0, 60);

    // numbers are from fitting VUV4 paper data with Ae^((V - Vbr)/b)
    TF1 * fLAP = new TF1("fLAP", "[1]*(x - [0])*(1 + 0.0438*(TMath::Exp((x - [0])/1.937) -1))", 0, 60);
    TFitResultPtr rfLine = gr_q1->Fit("fLine", "SMR0", "");

    fLAP->SetParameters(fLine->GetParameter(0), fLine->GetParameter(1));
    TFitResultPtr rfLAP = gr_q1Clone->Fit("fLAP", "SMR0", "");

    cout << *(rfLine->GetParams()) << endl;
    cout << *(rfLAP->GetParams()) << endl;

    //////////// c1 ////////////
    TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
    c1->cd();
    //c1->Divide(1, 2);

    // first canvas //
    /*
    c1->cd(1);

    hq51->SetLineColor(2);
    hq51->Draw("hist");
    
    hq53->SetLineColor(3);
    hq53->Draw("hist same");

    hq55->SetLineColor(4);
    hq55->Draw("hist same");

    hq51->GetYaxis()->SetRangeUser(1, hq55->GetMaximum());

    // second canvas //
    c1->cd(2);
    */
    gStyle->SetOptFit(1);

    hqr51->SetTitle("Gas CF4 - Pb210 Source");
    hqr51->GetXaxis()->SetTitle("Charge (nVs)");
    hqr51->SetLineColor(1);
    hqr51->SetLineWidth(3);
    hqr51->Draw("hist");

    // TPaveStats * st2 = (TPaveStats*)c1->GetPrimitive("stats");
    // st2->SetX1NDC(0.6);
    // st2->SetX2NDC(0.9);
    // st2->SetY1NDC(0.6);
    // st2->SetY2NDC(0.9);

    hqr53->SetLineColor(3);
    hqr53->SetLineWidth(3);
    hqr53->Draw("hist sames");

    hqr55->SetLineColor(4);
    hqr55->SetLineWidth(3);
    hqr55->Draw("hist sames");

    fg51->Draw("same");
    fg53->Draw("same");
    fg55->Draw("same");

    auto leg2 = new TLegend();
    leg2->AddEntry(hqr51, "51V", "L");
    leg2->AddEntry(hqr53, "53V", "L");
    leg2->AddEntry(hqr55, "55V", "L");
    leg2->Draw("same");

    c1->SetLogy();
    c1->SetGrid();
    c1->Update();


    //////////// c2 ////////////
    TCanvas * c2 = new TCanvas();
    gStyle->SetOptFit(1);

    gr_q1->SetTitle("Po peak vs voltage");
    gr_q1->GetXaxis()->SetLimits(40, 58);
    gr_q1->GetYaxis()->SetRangeUser(0., 150.);

    gr_q1->GetXaxis()->SetTitle("requested Voltage (V)");
    gr_q1->GetYaxis()->SetTitle("alpha peak charge (nVs)");

    gr_q1->Draw("APE*");
    fLine->Draw("same");
    gr_q1Clone->Draw("sames");
    fLAP->SetLineColor(3);
    fLAP->Draw("same");

    c2->SetGrid();
    c2->Update();
}

int main(){
    vbr_Po();
    }