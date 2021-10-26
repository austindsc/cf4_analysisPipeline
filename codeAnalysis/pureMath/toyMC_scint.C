
/* Toy Monte Carlo for building charge spectra

all for 1 event{
	1. generate Number of photons
		gaussian with mean E * LY

		ignore shadowing

	2. propagate photons in liquid CF4
		treat as transparent (for now)

	3. geometric factor
		use MFF functions

	4. is it detected (binomial)
		treat as 1 average PDE value

	5. charge from SPE gaussian distribution

	6. correlated pulses within few 100ns
}

*/

#include "../helpers/fitFunctions.C"

namespace MMC{
	
	int getNPE(TRandom * r, float E, float LY, float spread){
		// spread as percentage of mean
		// Energy in MeV
		// Light Yield in photons/MeV

		float mean = E * LY;
		float sig = mean * spread;
		return (int) r->Gaus(mean, sig);
	}

	bool bTry(TRandom * r, float eff){ //eff = 1, return get true, bernouli trial
		if(r->Rndm() <= eff) return true;
		else return false;
	}

	float getSipmCharge1(TRandom * r, float sig){
		return r->Gaus(1, sig);
	}

}

int toyMC_scint(){

	TRandom3 * rr = new TRandom3(0);

	int nbins = 300;
	int peMax = 30;


	// building histograms
	TH1 * hNPEraw = new TH1F("hNPEraw", "raw PE spectrum", peMax, -0.5, peMax-0.5);

	TH1 * hNPEsipm = new TH1F("hNPEsipm", "raw PE spectrum", nbins, 0, peMax);

	TH1 * hNPExtSipm = new TH1F("hNPExtsipm", "raw PE spectrum", nbins, 0, peMax);

	TH1 * hNPEfull = new TH1F("hNPEfull", "raw PE spectrum", nbins, 0, peMax);


	// Light production parameters
	float E = 5.4;
	float LY1 = 400;
	float LYsig = 0.1;

	// detection parameters
	Double_t sipmSide = 6*2 + 0.5; // in mm
	Double_t d = 25; // in mm

	float angleFactor = MFF::omegaSquare(&d, &sipmSide);
	float transportEff = 1.0;
	float sipmEff = 0.1; // sipm efficiency

	// sipm internal parameters
	float speValue = 1.0; // duh!
	float speSigma = 0.2;

	float sipmNoise = 0.1; // correlated avalanche

	// trigger and DAQ scheme
	float trigPE = 1.5;
	float trigSig = 0.2;

	// Running simulation
	int Nev = 100000; // total number of events

	for(int i = 0; i < Nev; i++){ // event loop

		int npe = MMC::getNPE(rr, E, LY1, LYsig); // number of photons
		int nDetRaw = 0;
		float nDetSipm = 0;
		float nDetSipmNoise = 0;

		// photon loop
		for(int iphoton = 0; iphoton < npe; iphoton++){

			if (! MMC::bTry(rr, angleFactor)) continue; // sipm solid angle
			if (! MMC::bTry(rr, transportEff)) continue; // sipm transport

			if (! MMC::bTry(rr, sipmEff)) continue; // sipm efficiency

			// photon detected!
			nDetRaw++;

			// add gaussian sipm response
			Double_t newCharge = MMC::getSipmCharge1(rr, speSigma);

			nDetSipm += newCharge;
			nDetSipmNoise += newCharge;

			// sipm crosstalk/noise
			while(MMC::bTry(rr, sipmNoise)){
				nDetSipmNoise += MMC::getSipmCharge1(rr, speSigma);
			}
		}

		if(nDetRaw == 0) continue;
		hNPEraw->Fill(nDetRaw);
		hNPEsipm->Fill(nDetSipm);
		hNPExtSipm->Fill(nDetSipmNoise);

		if(nDetSipm >= rr->Gaus(trigPE, trigSig)) hNPEfull->Fill(nDetSipm);
	}

	// scaling
	hNPEraw->Scale(hNPEsipm->GetBinWidth(5)/hNPEraw->GetBinWidth(5));
	hNPEraw->Scale(1/hNPEraw->GetEntries());
	hNPEsipm->Scale(1/hNPEsipm->GetEntries());
	hNPEfull->Scale(1/hNPEfull->GetEntries());
	hNPExtSipm->Scale(1/hNPExtSipm->GetEntries());

	hNPEsipm->SetLineColor(2);
	hNPEsipm->Draw("hist");

	hNPEfull->SetLineColor(3);
	hNPEfull->Draw("hist same");

	hNPExtSipm->SetLineColor(4);
	hNPExtSipm->Draw("hist same");

	hNPEraw->SetLineColor(1);
	hNPEraw->Draw("hist same");

	return 0;
}