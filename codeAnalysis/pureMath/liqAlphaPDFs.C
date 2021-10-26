#include "../helpers/fitFunctions.C"

void liqAlphaPDFs(){

	TF1 * f = new TF1("f", "ROOT::Math::binomial_pdf(x, 0.1, 30)", 0, 40);
	
	f->Draw();

	return;
}