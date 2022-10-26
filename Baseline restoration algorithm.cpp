//Baseline restoration algorithm
void BaselineRestoration(TH1D* hCorr, Double_t tau = 3300,Double_t t0=0.0)
{

  Double_t binsize = hCorr->GetBinWidth(0); // Bin width in nano-seconds
  vector<Double_t> vShift;
  vShift.push_back(0.0);// Initial condition
  Int_t binstart = 1;
  int n = 1;

  //Recursive shift
  for(Int_t iB=binstart;iB<hCorr->GetNbinsX();iB++){
    Double_t input = hCorr->GetBinContent(iB);//VExp(n)
    Double_t correction = input-vShift[n-1]*exp(-binsize/tau);//Vinp(n)
    Double_t shift = vShift[n-1]*exp(-binsize/tau)-correction*binsize/tau; // Vshift(n)
    vShift.push_back(shift);
    hCorr->SetBinContent(iB,correction);
    n++;
  }
  vShift.clear();

	Int_t nbins = hCorr->GetNbinsX();
  //Calculate DC offset at end of pulse
  Double_t offset = 0.0;
  Int_t bstart = nbins-7;
  Int_t bend = nbins-2;
  for(int iB =bstart;iB<bend;iB++)
  {
    offset+=hCorr->GetBinContent(iB);
  }

  offset/=Double_t(bend-bstart);

  //T0 of first pulse
  bstart = hCorr->FindBin(t0)+1;

	//Shift down post-pulse baseline
  for(int iB=bstart;iB<nbins;iB++)
  {
    Double_t bincontent = hCorr->GetBinContent(iB);
    hCorr->SetBinContent(iB,bincontent-offset);
  }

}