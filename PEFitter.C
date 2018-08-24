#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include <iostream>



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
  double M=par[5];
  double s0=par[3]*M;    // noise/enf are given as a fraction of 1pe signal
  double enf=par[4]*M;
  double noise = par[1]*TMath::Gaus(x,par[2],s0);
  double val=noise;
  for (int npe=1; npe<=npeFit; npe++){
    double mu=par[2]+M*npe;
    double sig=TMath::Sqrt(s0*s0+npe*enf*enf);
    val+=par[5+npe]*TMath::Gaus(x,mu,sig);
  }
  return val;
}


// function x-range is set later in fitting code
TF1 *tf_npefcn = new TF1("npefcn",fcn,0,10,15);


