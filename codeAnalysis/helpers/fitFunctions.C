
/*
- my fitting functions

*/

namespace MFF{

	/// pdfs ///
	Double_t myBinomial(Double_t *x, Double_t *par){

		int k = (int) x[0];
		Double_t p = par[0];
		int n = (int) par[1];

		return ROOT::Math::binomial_pdf(k, p, n);
	}

	Double_t scaledPoisson(Double_t *x, Double_t *par){

		int n = (int) (x[0] + 0.5);
		Double_t mu = par[0];
		int frontLimit = (int) par[1];
		Double_t scaling = ROOT::Math::poisson_cdf_c(frontLimit, mu);

		return ROOT::Math::poisson_pdf(n, mu)/scaling;
	}

	////////////// solid angle stuff //////////////
	Double_t omegaSquare(Double_t *x, Double_t *par){

		float d = x[0]; // source distance
		float a = par[0]; // side length
		Double_t top = a*a;
		Double_t bottom = 2*d * TMath::Sqrt(4*d*d + 2*a*a);

		return 4*TMath::ATan(top/bottom) /4 /TMath::Pi();
	}

	Double_t fLY_omegaSquare(Double_t *x, Double_t *par){

		float d = x[0]; // source distance

		float a = par[0]; // side length
		float scale = par[1]; //scaling factor

		Double_t top = a*a;
		Double_t bottom = 2*d * TMath::Sqrt(4*d*d + 2*a*a);

		Double_t Rfrac = 4*TMath::ATan(top/bottom) /4 /TMath::Pi();

		return scale * Rfrac;
	}

	Double_t omegaSimpleSquare(Double_t *x, Double_t *par){

		Double_t d = x[0]; //source distance
		Double_t a = par[0]; // sipm side length

		Double_t sipmArea = a*a;
		Double_t circleArea = 4 * TMath::Pi() * d * d;

		return sipmArea / circleArea;
	}

}