#include <stdio.h>
using namespace std;
using namespace TMath;

Double_t omegaSquare(Double_t *x, Double_t *par){

	float d = x[0];
	float a = par[0]; // side length
	Double_t top = a*a;
	Double_t bottom = 2*d * TMath::Sqrt(4*d*d + 2*a*a);

	return 4*TMath::ATan(top/bottom) /4 /TMath::Pi();
}

Double_t omegaSimpleSquare(Double_t *x, Double_t *par){

	Double_t d = x[0];
	Double_t a = par[0]; // sipm side length

	Double_t sipmArea = a*a;
	Double_t circleArea = 4 * TMath::Pi() * d * d;

	return sipmArea / circleArea;
}

void solidAngle(){

	cout << "YOOOO " << endl;

	float sipmSide = 6*2 + 0.5; // in mm
	
	TF1 * fSquare = new TF1("fSquare", omegaSquare, 0, 40, 1);
	fSquare->SetParameter(0, sipmSide);

	TF1 * fApprox = new TF1("fApprox", omegaSimpleSquare, 0, 40, 1);
	fApprox->SetParameter(0, sipmSide);

	fSquare->SetLineColor(4);
	fSquare->GetXaxis()->SetTitle("SiPM distance (mm)");
	fSquare->GetYaxis()->SetTitle("Solid Angle Fraction");
	fSquare->SetTitle("Square solid angle functions");
	fSquare->Draw();
	fApprox->Draw("same");

	return;
}