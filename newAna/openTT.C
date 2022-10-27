

void openTT(){

  _file0->ls();

  Float_t bline[3];

  tt->SetBranchAddress("blineV", bline);

  tt->GetEntry(0);

  tt->Show(0);

  cout << bline[0] << endl;

  tt->Draw("blineV[0]");
  tt->Draw("blineV[1]", "", "hist same");

}
