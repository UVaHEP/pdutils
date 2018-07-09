#include "TMath.h"
#include "TF1.h"
#include "TH1.h"



// "simple" function to fit multiple peaks in pulse height distribution
// par[0] : # of peaks to fit 0=noise only, 1=noise+1pe, ....
// par[1] : noise peak normalization
// par[2] : noise peak mean
// par[3] : noise peak width
// par[4] : enf
// par[5] : gain
// par[6] : np1 normalization
// par[7] : np2 normalization
// ...
Double_t fcn(Double_t *xp, Double_t *par){
  double x=xp[0];
  int npeFit=par[0];
  double noise = par[1]*TMath::Gaus(x,par[2],par[3]);
  double val=noise;
  double enf=par[4];
  double M=par[5];
  for (int npe=1; npe<=npeFit; npe++){
    double mu=par[2]+M*npe;
    double sig=TMath::Sqrt(par[3]*par[3]+npe*enf*enf);
    val+=par[5+npe]*TMath::Gaus(x,mu,sig);
  }
  return val;
}


TF1 *tf_npefcn = new TF1("npefcn",fcn,0,10,20);


