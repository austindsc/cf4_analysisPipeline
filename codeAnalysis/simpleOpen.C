#include <stdio.h>
using namespace std;
using namespace TMath;

/*
run in command line with "root -l simpleOpen.C"

*/

void drawStraightFromTtree(TTree * t){
  
  auto c5 = new TCanvas("c5", "c5", 600, 400); // canvas to draw in
  TH1F * hla4 = new TH1F("hla4", "hla4", 150, 0, 300);

  hla4->SetTitle("Liquid alpha scintillation");
  hla4->GetXaxis()->SetTitle("raw charge");

  t->Draw("qs50 >> hla4", "", ""); // draw each TTree's data in hla4

  auto leg5 = new TLegend(); // legend
  leg5->AddEntry(hla4, "4V, alpha", "L");
  leg5->Draw("same"); // draw on same canvas

  c5->SetLogy();
  c5->Update();
}

void drawWithEachEntry(TTree *t){

  float qSmart50;
  t->SetBranchAddress("qs50", &qSmart50); // tell the TTree to load its data into the qSmart50 variable

  TH1 * hExampleQ = new TH1F("hEx", "example histogram title", 200, 0, 30); //histogram to fill

  for (int i = 0; i < t->GetEntries(); ++i){ //loop on ttree entries (rows of data)
    t->GetEntry(i); // load next data into variables


    hExampleQ->Fill(qSmart50/100); //fill histogram bin and do some math (divide valuee by 100)
  }

  TCanvas * c1 = new TCanvas("c1", "c1 canvas name", 600, 400);
  hExampleQ->Draw("hist");
  c1->Update();
}

int simpleOpen(){


  // load large trigger data, liquid
  TFile * fla4 = new TFile("/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/outProcStatic1_20210726_7.root");
  TTree * ttla4 = (TTree*)fla4->Get("tt");

  cout << "showing ttree entry 4" << endl;
  ttla4->Show(4);

  drawStraightFromTtree(ttla4);

  drawWithEachEntry(ttla4);

	return 0;
}