#include "helpers/getSPE.C"
#include "helpers/fitFunctions.C"

// TH1* ozGetHist(TTree * tt, string variable, string sTh1, string cuts){
	
// 	tt->Draw(sTh1+">>", cuts, "goff");
// 	return (TH1*) gDirectory->Get
// }

/*
load in all 'mega' data files

for the 1.25cm run, remove events with bad baseline or end of pulse too late
*/

/*
1. load in data files
1. get SPE values from other files
3. 


*/


int liqVBRspe(){

	TFile * files[3][3] = {
						{TFile::Open("$CF4HOME/procData/liqSPE1outProcStatic2_20210726_0.root", "READ"),
						TFile::Open("$CF4HOME/procData/liqSPE1outProcStatic2_20210726_1.root", "READ"),
						TFile::Open("$CF4HOME/procData/liqSPE1outProcStatic2_20210726_2.root", "READ")},

						{TFile::Open("$CF4HOME/procData/liqSPE1outProcStatic2_20210804_0.root", "READ"),
						TFile::Open("$CF4HOME/procData/liqSPE1outProcStatic2_20210804_1.root", "READ"),
						TFile::Open("$CF4HOME/procData/liqSPE1outProcStatic2_20210804_2.root", "READ")},

						{TFile::Open("$CF4HOME/procData/liqSPE1outProcStatic2_20210813_0.root", "READ"),
						TFile::Open("$CF4HOME/procData/liqSPE1outProcStatic2_20210813_1.root", "READ"),
						TFile::Open("$CF4HOME/procData/liqSPE1outProcStatic2_20210813_2.root", "READ")}
					};

	TTree * trees[3][3];
	for(int irun = 0; irun < 3; irun++){
		for(int ivolt = 0; ivolt < 3; ivolt++){
			// cout << files[irun][ivolt]->GetName() << endl;
			trees[irun][ivolt] = (TTree*) files[irun][ivolt]->Get("tt");
			// cout << trees[irun][ivolt]->GetEntries() << endl;
		}
	}

	speStuff speVals[3][3];
	Double_t voltages[3] = {48.25, 49.25, 50.25};
	Double_t errV = 0.05;
	for(int irun = 0; irun < 3; irun++){
		for(int ivolt = 0; ivolt < 3; ivolt++){
			TTree * t = trees[irun][ivolt];
			int verbose = 0;

			if(irun == 0) {new TCanvas(); verbose = 1;}

			cout << "\nrun " << irun << "    voltage " << ivolt << endl;
			speStuff speInfo = getSPEval(t, 0.15, verbose);
			speInfo.V = voltages[ivolt];
			speInfo.eV = errV;

			speVals[irun][ivolt] = speInfo;
		}
	}

	TF1 * flinearvbr = new TF1("flinearvbr", "[0]*(x - [1])");
	TF1 * fits[3];

	TGraphErrors * tgRuns[3];
	for (int irun = 0; irun < 3; ++irun){
		TGraphErrors * tg = new TGraphErrors(3);
		tgRuns[irun] = tg;
		for(int ivolt = 0; ivolt < 3; ivolt++){
			speStuff spe = speVals[irun][ivolt];
			
			tg->SetPoint(ivolt, spe.V, spe.spe);
			tg->SetPointError(ivolt, spe.eV, spe.espe);
		}

		//fitting tgraphs
		string fname = "f" + to_string(irun);
		fits[irun] = new TF1(fname.c_str(), "flinearvbr", 0, 60);
		tg->Fit(fname.c_str(), "+RM0Q", "goff");
	}

	TCanvas * cg = new TCanvas();
	TLegend * l1 = new TLegend();

	tgRuns[0]->GetYaxis()->SetRangeUser(0, 1.25);
	tgRuns[0]->GetXaxis()->SetLimits(42., 51.);
	tgRuns[0]->GetXaxis()->SetTitle("Applied Voltage");
	tgRuns[0]->GetYaxis()->SetTitle("SPE charge (V x ns)");
	tgRuns[0]->SetTitle("SPE runs: Applied voltage vs SPE charge");

	tgRuns[0]->SetMarkerColor(1);
	tgRuns[0]->Draw("A*");
	fits[0]->Draw("same");

	tgRuns[1]->SetMarkerColor(3);
	tgRuns[1]->Draw("* same");
	fits[1]->Draw("same");

	tgRuns[2]->SetMarkerColor(4);
	tgRuns[2]->Draw("* same");
	fits[2]->Draw("same");

	l1->AddEntry(tgRuns[0], "Run A - 1.5cm", "P");
	l1->AddEntry(tgRuns[1], "Run B - 1.25cm", "P");
	l1->AddEntry(tgRuns[2], "Run C - 1.9cm", "P");
	l1->Draw("same");

	cg->Update();

	return 0;
}