

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
struct myFitData{

    float voltage;
    string source;
    TF1* f1PE;
    TF1* fNPE;
    TF1* fH1PE; // height gaussian

    float spe1, espe1; // spe from 1gaus fit
    float speN1, espeN1; //spe, first term in 3gaus fit
    float speN2, espeN2; //2pe - 1pe from 3gaus fit
    float tpe, etpe; // 2pe peak, 3 gaus fit

    float hpe1, ehpe1; // 1PE heights

    myFitData();
    myFitData(float V, string src);

    void SetQ1PEvars();
    void SetH1PEvars();
};

myFitData::myFitData(){}
myFitData::myFitData(float V, string src){
    voltage = V;
    source = src;
}

void myFitData::SetQ1PEvars(){
    spe1 = f1PE->GetParameter(1);
    espe1 = f1PE->GetParError(1);
    speN1 = fNPE->GetParameter(1);
    espeN1 = fNPE->GetParError(1);
    tpe = fNPE->GetParameter(4);
    etpe = fNPE->GetParError(4);
    speN2  = tpe - speN1;
    espeN2 = TMath::Sqrt(etpe*etpe + espeN1*espeN1);
}

void myFitData::SetH1PEvars(){
    hpe1 = fH1PE->GetParameter(1);
    ehpe1= fH1PE->GetParError(1);
}

///////////////////// done with data holding object /////////////////////

void getRates(string fileName, bool plot, myFitData * fitR, float limLow, float limHigh){
    const char * fitOption = "MR0";
    string path = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/proc2/";

    TH1 *hQ, *kQ;
    TFile * f;
    TVectorD* v;

    /////////// loading data 1 ///////////
    string nPo = fileName;

    f = TFile::Open((path + nPo).c_str(), "READ");
    hQ = (TH1*) f->Get( "hF_all" )->Clone(); // Po charge cut
    kQ = (TH1*) f->Get( "hF_cuts" )->Clone(); // Po charge cut

    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t lifetime = (*v)[0];

    int rebin = 1;
    if (rebin > 1){
        kQ->Rebin(rebin);
    }

    // kQ->Scale(1/lifetime);
    // hQ->Scale(1/lifetime);

    float area = 4*5.85*5.95;

    cout << "rate 1, uncut (Hz/mm^2) = " << hQ->GetEntries()/lifetime/area << endl;
    cout << "rate 1, cuts (Hz/mm^2) = " << kQ->GetEntries()/lifetime/area << endl;

    // counting peak size
    // float a_raw_Bkg = hQ->Integral( hBkg->GetXaxis()->FindBin(firstBin), hBkg->GetNbinsX() );
    
    if (plot == true){
        TCanvas *c = new TCanvas(("Q: " + fileName).c_str(), ("Q spe" + fileName).c_str(), 800, 600);

        hQ->SetLineColor(3);
        hQ->Draw("hist");
        kQ->Draw("hist same");

        c->Update();
    }
}

void speAna_h(string fileName,  bool plot, myFitData * fitR){
    
    const char * fitOption = "MR0";
    string path = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/proc2/";

    TH1 *hH, *kH;
    TH2 *hQH, *kQH;
    TFile * f;
    TVectorD* v;

    /////////// loading data 1 ///////////
    string nPo = fileName;

    f = TFile::Open((path + nPo).c_str(), "READ");
    kH = (TH1*) f->Get( "hH_cuts" )->Clone(); // Po charge cut
    kQH = (TH2*) f->Get( "h2_hq" )->Clone();

    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t lifetimePo = (*v)[0];

    /////////// manipulating ///////////
    int rebin = 1;
    if (kH->GetEntries() < 10000) rebin = 2;
    kH->Rebin(rebin);
    //kH->Scale(1/kH->Integral(0, kH->GetXaxis()->FindBin(0.05))); //normalize

    /////////// fitting ///////////
    float maxX = kH->GetBinCenter(kH->GetMaximumBin());
    float max = kH->GetMaximum();
    float rangeG1 = 0.003;
    TF1 * fg_H1PE = new TF1("fg_H1PE", "gaus(0)", maxX-1.5*rangeG1, maxX+rangeG1);
    fg_H1PE->SetParameters(max, maxX, rangeG1);

    cout << "height = " << max << endl;
    cout << "mean = " << maxX << endl;
    kH->Fit("fg_H1PE", fitOption, "");


    if (plot == true){
        TCanvas *c = new TCanvas(("H spe: " + fileName).c_str(), ("H spe" + fileName).c_str(), 800, 600);

        kH->SetTitle(("H hist: " + fileName).c_str());
        kH->GetXaxis()->SetRangeUser(-0.01, 0.05);
        kH->Draw("hist");
        fg_H1PE->Draw("same");
        gStyle->SetOptFit();
        //kQH->Draw("colz");

        c->SetGrid();
        c->Update();
    }

    fitR->fH1PE = fg_H1PE;
    fitR->SetH1PEvars();
}

void speAna_q(string fileName,  bool plot, myFitData * fitR){

    const char * fitOption = "MR0Q";

    string path = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4_coding/proc2/";

    TH1 *hBkg, *hPo, *hCsPo;
    TH1 *kBkg, *kPo, *kCsPo;
    TH2 *k2_bkg, *k2_Po, *k2_CsPo;
    TFile * f;
    TVectorD* v;

    /////////// loading data 1 ///////////
    string nPo = fileName;

    // need bkg data from same day as gasBig_Po (1123)
    // for more correct 1PE estimation (same temp?)

    f = TFile::Open((path + nPo).c_str(), "READ");
    kPo = (TH1*) f->Get( "hF_cuts" )->Clone(); // Po charge cut
    k2_Po = (TH2*) f->Get( "h2_hq" )->Clone();

    v = (TVectorD* ) f->Get("liveTime")->Clone();
    Double_t lifetimePo = (*v)[0];

    /////////// manipulating ///////////
    // make a copy, for a finer 1pe fit
    TH1* kPo_spe = (TH1*) kPo->Clone();

    int rebinSpe = 1;
    kPo_spe->Rebin(rebinSpe);
    kPo_spe->Scale(1/kPo_spe->Integral(0, kPo_spe->GetXaxis()->FindBin(15))); //normalize

    int rebin = 20;
    kPo->Rebin(rebin);
    kPo->Scale(1/lifetimePo);

    // counting peak size
    float firstBin = 0.7;

    // issue is this depends on DN (V_over) not changing (temperature)

    ////////////////// fitting small SPE data //////////////////
    float rNPE_low = 0.5;
    float rNPE_high = 5.;
    float r1PE_high = 3.;
    float rNPE_factor = 1.;

    if(kPo_spe->GetEntries() > 20000){
        rNPE_factor = 0.75;
    }

        // 1gaus 1PE fit
    TF1 * fgPo_1PE = new TF1("fgPo_1PE", "gaus", 0.5, kPo_spe->GetBinCenter(kPo_spe->GetMaximumBin()) + 0.5 );
    fgPo_1PE->SetParameters(kPo_spe->GetMaximum(), 2, 0.75);
    kPo_spe->Fit("fgPo_1PE", fitOption, "");
    float spe_Height = fgPo_1PE->GetParameter(0);
    float spe = fgPo_1PE->GetParameter(1);
    float spe_width = fgPo_1PE->GetParameter(2);


        ////// NPE fit
    TF1 * fgPo_NPE = new TF1("fgPo_NPE", "gaus(0) + gaus(3) + gaus(6)", rNPE_low, 2*(spe + rNPE_factor * spe_width));
    fgPo_NPE->SetParNames("1height", "1PE", "1sigma",  "2 height", "2PE", "2sigma", "AP height", "AP mean", "AP sigma");

        // 1PE peak
    fgPo_NPE->SetParameter(0, 0.9 * spe_Height);
    fgPo_NPE->SetParLimits(0, 0.6 * spe_Height, spe_Height);
    fgPo_NPE->SetParameter(1, 0.9 * spe);
    fgPo_NPE->SetParLimits(1, 0.8 * spe, 1. * spe);
    fgPo_NPE->SetParameter(2, spe_width);
    //fgPo_1PE->SetParLimits(2, speBkg_width, speBkg_width);

        // 2PE peak
    fgPo_NPE->SetParameter(3, spe_Height/15);
    fgPo_NPE->SetParLimits(3, 0, spe_Height/2);
    fgPo_NPE->SetParameter(4, 2.2 * spe);
    fgPo_NPE->SetParLimits(4, 1.6* spe, 2.4*spe);
    fgPo_NPE->SetParameter(5, spe_width);
    
        // 1PE afterpulses
    fgPo_NPE->SetParameter(6, spe_Height/15);
    fgPo_NPE->SetParLimits(6, spe_Height/20, spe_Height/2);
    fgPo_NPE->SetParameter(7, 1.2 * spe);
    fgPo_NPE->SetParLimits(7, 1.* spe, 1.5*spe);
    fgPo_NPE->SetParameter(8, spe_width);

    kPo_spe->Fit("fgPo_NPE", fitOption, "");

    //cout << "1PE FIT = " << fgPo_1PE->GetChisquare()/fgPo_1PE->GetNDF() << endl;
    //cout << "NPE FIT = " << fgPo_NPE->GetChisquare()/fgPo_NPE->GetNDF() << endl;

    ////////////////// getting extra functions for plotting //////////////////
    TF1 * fgPo_ex1 = new TF1("fgPo_ex1", "gaus", 1, 5.);
    fgPo_ex1->SetParameters(fgPo_NPE->GetParameter(0), fgPo_NPE->GetParameter(1),  fgPo_NPE->GetParameter(2));
    TF1 * fgPo_ex2 = new TF1("fgPo_ex2", "gaus", 1, 7.);
    fgPo_ex2->SetParameters(fgPo_NPE->GetParameter(3), fgPo_NPE->GetParameter(4),  fgPo_NPE->GetParameter(5));
    TF1 * fgPo_ex1AP = new TF1("fgPo_ex1AP", "gaus", 1, 5.);
    fgPo_ex1AP->SetParameters(fgPo_NPE->GetParameter(6), fgPo_NPE->GetParameter(7),  fgPo_NPE->GetParameter(8));

    // TF1 * fgPo_ex2AP = new TF1("fgPo_ex2AP", "gaus", 1, 8.);
    // fgPo_ex2AP->SetParameters(fgPo_NPE->GetParameter(9), fgPo_NPE->GetParameter(10),  fgPo_NPE->GetParameter(11));


    if(plot == true){
        ////////// SPE plots ////////////
        TCanvas *c = new TCanvas(("Q spe: " + fileName).c_str(), ("Q spe" + fileName).c_str(), 800, 600);
        
        kPo_spe->SetName("full fit");
        kPo_spe->SetTitle(("Q hist: " + fileName).c_str());
        kPo_spe->GetXaxis()->SetTitle("charge (nsV)");
        kPo_spe->GetXaxis()->SetTitleSize(0.06);
        kPo_spe->GetYaxis()->SetTitle("normalized counts");

        kPo_spe->GetXaxis()->SetRangeUser(-.5, 10.);
        kPo_spe->SetFillColorAlpha(41, 0.75);
        kPo_spe->Draw("hist");

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

        // fgPo_ex2AP->SetLineColor(7);
        // fgPo_ex2AP->Draw("same");

        auto leg2B = new TLegend();
        leg2B->AddEntry(fgPo_1PE, "1PE only", "L");
        leg2B->AddEntry(fgPo_NPE, "full: A+B+C", "L");
        leg2B->AddEntry(fgPo_ex1, "A: 1PE", "L");
        leg2B->AddEntry(fgPo_ex2, "B: 2PE", "L");
        leg2B->AddEntry(fgPo_ex1AP, "C: 1PE CDA (delayed)", "L");
        leg2B->Draw("same");

        gStyle->SetOptFit();

        c->SetBottomMargin(0.15);

        kPo_spe->GetXaxis()->SetTitleOffset(0.);

        c->SetGrid();
        c->SetGrid();
        //c->SetLogy();
        c->Update();
    } // plotting if statement

    fitR->f1PE = fgPo_1PE;
    fitR->fNPE = fgPo_NPE;
    fitR->SetQ1PEvars();
}

// have three SPE values for each voltage.
// f1PE->GetParameter(1)
// fNPE->GetParameter(1)
// fNPE-> par2 - par1

// par2 is very sensitive to the range (esp at 51V)

int compGas2_speVbr(){
    
    myFitData *fake1, *fake2, *fake3;
    bool view = true;
    myFitData * d51 = new myFitData(51., "Pb");
    myFitData * d53 = new myFitData(53., "Pb");
    myFitData * d55 = new myFitData(55., "Pb");

    getRates("gasBig_bkg_51_1203.root", true,  fake1, 0.,  10.);
    getRates("gasBig_bkg_53_1203.root", true,  fake2, 0.,  10.);
    getRates("gasBig_bkg_55_1203.root", true,  fake3, 0.,  10.);
/*
    ////// VBr from Height //////
    TF1 * fLine = new TF1("fLine", "[1]*(x - [0])", 0, 60);
    Double_t V[3] = {51., 53., 55.};
    Double_t eV[3] = {0.1, 0.1, 0.1};

    speAna_h("gasBig_Po_51_1123.root",  view, d51);
    d51->fH1PE->SetName("fg_Po_H1pe_51V"); // MUST CHANGE NAME

    speAna_h("gasBig_Po_53_1123.root",  view, d53);
    d53->fH1PE->SetName("fg_Po_H1pe_53V"); // MUST CHANGE NAME

    speAna_h("gasBig_Po_55_1123.root",  view, d55);
    d55->fH1PE->SetName("fg_Po_H1pe_55V"); // MUST CHANGE NAME
    
    float vShift = 0.00;
    Double_t h1PE[3] = {d51->hpe1 - vShift, d53->hpe1 - vShift, d55->hpe1 -vShift};
    Double_t eh1PE[3] = {d51->ehpe1, d53->ehpe1, d55->ehpe1};

    TGraphErrors * gVbr_H1pe = new TGraphErrors(3, V, h1PE, eV, eh1PE);
    TF1 * flineH1PE = new TF1("flineH1PE", "fLine", 0, 60);
    flineH1PE->SetParameters(48., d53->hpe1/8.);
    gVbr_H1pe->Fit("flineH1PE", "MR0", "");

    // plotting //
    auto c3 = new TCanvas("c3", "c3", 600, 600);
    gVbr_H1pe->GetYaxis()->SetRangeUser(0., 0.008);
    gVbr_H1pe->GetXaxis()->SetLimits(30., 58.);

    gVbr_H1pe->SetMarkerColor(2);
    gVbr_H1pe->SetMarkerSize(1.5);
    gVbr_H1pe->SetMarkerStyle();
    gVbr_H1pe->Draw("AP");

    flineH1PE->SetLineColor(2);
    flineH1PE->Draw("same");

    gStyle->SetOptFit();
    c3->SetGrid();
    c3->Update();

    ////////// getting Q VBR data //////////
    //speAna_q("gasBig_Po_51_1123.root",  view, d51);
    speAna_q("gasBig_Po_51_1123.root",  view, d51);
    d51->f1PE->SetName("fg_Po_npe_51V"); // MUST CHANGE NAME

    speAna_q("gasBig_Po_53_1123.root",  view, d53);
    d53->f1PE->SetName("fg_Po_npe_53V"); // MUST CHANGE NAME

    speAna_q("gasBig_Po_55_1123.root",  view, d55);
    d55->f1PE->SetName("fg_Po_npe_55V"); // MUST CHANGE NAME

    ////// building Vbr graphs //////
    Double_t q1PE[3] = {d51->spe1, d53->spe1, d55->spe1};
    Double_t qN1PE[3] = {d51->speN1, d53->speN1, d55->speN1};
    Double_t qN2PE[3] = {d51->speN2, d53->speN2, d55->speN2};

    Double_t eq1PE[3] = {d51->espe1, d53->espe1, d55->espe1};
    Double_t eqN1PE[3] = {d51->espeN1, d53->espeN1, d55->espeN1};
    Double_t eqN2PE[3] = {d51->espeN2, d53->espeN2, d55->espeN2};

    // breakdownn voltage graphs
    TGraphErrors * gVbr_1PE = new TGraphErrors(3, V, q1PE, eV, eq1PE); //  1PE
    TGraphErrors * gVbr_N1PE = new TGraphErrors(3, V, qN1PE, eV, eqN1PE); // 1PE from NPE
    TGraphErrors * gVbr_N2PE = new TGraphErrors(3, V, qN2PE, eV, eqN2PE); // 2PE-1PE from NPE

    ////////// fitting vbr graphs //////////

    TF1 * fline1PE = new TF1("fline1PE", "fLine", 0, 60);
    fline1PE->SetParameters(48., 0.5);
    gVbr_1PE->Fit("fline1PE", "MR0", "");

    TF1 * flineN1PE = new TF1("flineN1PE", "fLine", 0, 60);
    flineN1PE->SetParameters(48., 0.5);
    gVbr_N1PE->Fit("flineN1PE", "MR0", "");

    TF1 * flineN2PE = new TF1("flineN2PE", "fLine", 0, 60);
    flineN2PE->SetParameters(48., 0.5);
    gVbr_N2PE->Fit("flineN2PE", "MR0", "");

    ////////// plotting //////////
    auto c2 = new TCanvas("c2", "c2",  600, 600);

    gVbr_1PE->SetTitle("SPE charge vs V");
    gVbr_1PE->GetYaxis()->SetRangeUser(0., 3.);
    gVbr_1PE->GetXaxis()->SetLimits(30., 58.);

    gVbr_1PE->GetYaxis()->SetTitle("SPE charge");
    gVbr_1PE->GetXaxis()->SetTitle("Requested Voltage (V)");

    gVbr_1PE->SetMarkerColor(2);
    gVbr_1PE->SetMarkerSize(1.5);
    gVbr_1PE->SetMarkerStyle(21);
    gVbr_1PE->Draw("AP");

    fline1PE->SetLineColor(2);
    fline1PE->Draw("same");

    gVbr_N1PE->SetMarkerColor(3);
    gVbr_N1PE->SetMarkerSize(1.5);
    gVbr_N1PE->SetMarkerStyle(43);
    gVbr_N1PE->Draw("P same");

    flineN1PE->SetLineColor(3);
    flineN1PE->Draw("same");

    gVbr_N2PE->SetMarkerColor(4);
    gVbr_N2PE->SetMarkerSize(1.5);
    gVbr_N2PE->SetMarkerStyle(47);
    gVbr_N2PE->Draw("P same");

    flineN2PE->SetLineColor(4);
    flineN2PE->Draw("same");

    auto leg2 = new TLegend();
    leg2->AddEntry(gVbr_1PE, "1gaus fit", "P");
    leg2->AddEntry(gVbr_N1PE, "full fit, 1PE", "P");
    leg2->AddEntry(gVbr_N2PE, "full fit, 2PE-1PE", "P");
    leg2->Draw("same");

    c2->SetGrid();
    c2->Update();
*/
    return 0;
    }