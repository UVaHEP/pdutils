from ROOT import *
from math import sqrt
from math import log
from math import pi


# light source
# assume Poisson distribution of photons based on mean mu
# and additional Gaussian variation of mu on each pulse
class LightSource():
    def __init__(self,muGamma=20, pulseJitter=0.10):
        self.muGamma=muGamma   # average photons per pulse
        self.pulseJitter=pulseJitter # fluctuation in pulse intensity
        self.trand=TRandom(0)
        self.count=0
        self.pulsehist=TH1F()
        self.pulsehist.SetTitle("Photons per pulse")
        self.pulsehist.SetName("nPhotons")

    def Print(self):
        print "Light source parameters"
        print "<photons/pulse>",self.muGamma
        print "SD of pulse intensity:",self.pulseJitter
        
    def Pulse(self):         # pulse intensity in # of photons
        if (self.count==0):
            xmin=0
            xmax=int(self.muGamma+sqrt(self.muGamma)*4)
            self.pulsehist.SetBins(xmax,xmin,xmax)
        nPhot=self.muGamma*(1+self.trand.Gaus()*self.pulseJitter)
        nPhot=self.trand.Poisson(nPhot)
        self.pulsehist.Fill(nPhot)
        return nPhot
        
# very simple SIPM model
class SIPM():
    def __init__(self,pde=0.2,noise=0.1,enf=0.1,dcr=2,M=1,pulseWid=25):
        # device properties
        self.pde=pde
        self.noise=noise        # electronic noise relative to 1PE peak
        self.enf=enf            # excess noise factor (eg uniformity)
        self.dcr=dcr            # in MHz
        self.M=M                # Gain = counts/1pe
        self.pulseWid=pulseWid  # in ns
        self.trand=TRandom(0)

    def Print(self):
        print "---"
        print "SIPM input parameters"
        print "pde:",self.pde
        print "gain",self.M
        print "noise, enf",self.noise,self.enf
        print "---"

    def NoiseSample(self):
        return self.trand.Gaus()*self.noise*self.M

    # calculate signal contribution from DCR
    def DarkSample(self,gate):
        # dark pulses expected in integration window
        extendedWindow=self.pulseWid*1e-9*(gate+1)
        reducedWindow=self.pulseWid*1e-9*(gate-1)
        pulseAny=self.dcr*1e6 * extendedWindow
        pulseFull=self.dcr*1e6 * reducedWindow
        r=self.trand.Uniform()
        darkSig=0
        if r<pulseFull: darkSig=1
        elif r<pulseAny: darkSig=self.trand.Uniform()
        darkSig=darkSig*self.M*(1+self.trand.Gaus()*self.enf)
        return darkSig

    # simulate charge collection from pulses in gate window, including effect of random overlaps fronm DCR 
    # can also add afterpulsing, but how to model?  Probably needs recharge time + falling expo probabiity
    def SimPhD(self,pulser,gate,ntrials=20000):  # get pulse height distribution
        if gate<1: print "Warning gate < pulse width"
        muGamma=pulser.muGamma
        xmin=-3*self.noise
        xmax=self.pde*pulser.muGamma+4*sqrt(self.pde*pulser.muGamma)*self.M
        hPhD = TH1F("hPhD","Pulse height dist",300,xmin,xmax)
        for nt in range(ntrials):
            signal=self.NoiseSample()+self.DarkSample(gate) # noise+dark counts
            nPhot=pulser.Pulse()      # ignore gate assume in time, wid>=1.0
            for np in range(nPhot):
                if self.trand.Uniform()<self.pde:
                    signal=signal+self.M*(1+self.trand.Gaus()*self.enf)
            hPhD.Fill(signal)
        return hPhD

    # gate describes integration window relative to pulse width
    # eg gate = 1, means gate = pulse width
    # very simple dark count model
    def SimDarkPhD(self, gate, ntrials=10000):
        if gate<1: print "Warning gate < pulse width"
        xmin=-3*self.noise
        xmax=self.M*4
        hDark = TH1F("hDark","Dark pulse height dist",300,xmin,xmax)
        for nt in range(ntrials):
            signal=self.NoiseSample() + self.DarkSample(gate)
            hDark.Fill(signal)
        return hDark


class PeakFitter():
    def __init__(self):
        code="""
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
        val+=TMath::Gaus(x,mu,sig);
        }
        return val;
        };
        """
        gInterpreter.ProcessLine(code)

    
class PhDAnalyzier():
    def __init__(self,hPhD,hDark=None):
        self.hPhD=hPhD
        self.ts=TSpectrum()
        self.hPhD0=hDark    # dark count pulse height distributions
        self.xylist=[]
        
    # use this method for the dark pulse height distribution     
    def Fit0Peak(self):
        maxbin=self.hPhD0.GetMaximumBin()
        max=self.hPhD0.GetBinContent(maxbin)
        mu=self.hPhD0.GetBinCenter(maxbin)
        for ibin in range(maxbin-1,0,-1):
            y=self.hPhD0.GetBinContent(ibin)
            if y<max/2:
                sig=mu-self.hPhD0.GetBinCenter(ibin)
                break
        xmin=self.hPhD0.GetBinCenter(1)
        xmax=self.hPhD0.GetBinCenter(maxbin)+sig*1.5 # 1.5 is a hack!
        self.hPhD0.Fit("gaus","","",xmin,xmax)
        
    def FindPeaks(self):
        hfft=hPhD.FFT(0,"RE")
        hfft.SetBinContent(1,0)  # suppress DC component
        max=hfft.GetMaximumBin()
        if max>hfft.GetNbinsX()/2: max = hfft.GetNbinsX()-max
        self.peakWid=(hPhD.GetNbinsX()/max)/4  # est. peak distance / 4 in Nbins
        self.npeaks=self.ts.Search(hPhD,self.peakWid)
        print "found",self.npeaks,"peaks"
        xvals=self.ts.GetPositionX()
        yvals=self.ts.GetPositionY()
        for i in range(self.npeaks):
            self.xylist.append([xvals[i],yvals[i]])
        self.xylist.sort()        
        return self.npeaks

    # warning must call Fit0Peak and FindPeaks first
    def FitPhD(self):
        parnames=["nfit","a0","mu0","sig0","enf","gain"]
        gROOT.ProcessLine(".L PEFitter.C+")
        npefcn.SetRange(self.hPhD.GetXaxis().GetXmin(),self.hPhD.GetXaxis().GetXmax())
        npefcn.SetNpx(self.hPhD.GetNbinsX())
        npePeaks=self.npeaks-1
        fcn0=self.hPhD0.GetFunction("gaus")
        npefcn.FixParameter(0,npePeaks) # of peaks to fit 0=noise only, 1=noise+1pe, ....
        npefcn.SetParameter(1,self.xylist[0][1]) # noise peak normalization
        npefcn.SetParameter(2,fcn0.GetParameter(1)) # noise peak mean
        npefcn.SetParameter(3,fcn0.GetParameter(2)) # noise peak width
        npefcn.SetParameter(4,fcn0.GetParameter(2)) # enf -- starting guess
        npefcn.SetParameter(5,self.xylist[1][0]-self.xylist[0][0]) # gain
        npefcn.Print()
        for i in range(npePeaks):
            npefcn.SetParameter(6+i,self.xylist[1+i][1])
            parnames.append("a"+str(i+1))
        for i in range(len(parnames)): npefcn.SetParName(i,parnames[i])
        for i in range(len(parnames),npefcn.GetNpar()): npefcn.FixParameter(i,0)
        self.hPhD.Fit("npefcn")
        
    
    def CalcNpe(self):  # calculate average # of detected photons per pulse
        npeaks=self.FindPeaks()
        zeropeak=(self.xylist[0])[0]
        xmin=zeropeak-self.peakWid*self.hPhD.GetBinWidth(1)*1.5 # 1.5 is a hack!
        xmax=zeropeak+self.peakWid*self.hPhD.GetBinWidth(1)*1.5
        self.hPhD.Fit("gaus","","",xmin,xmax)
        fcn=self.hPhD.GetFunction("gaus")
        A=fcn.GetParameter(0)
        mu=fcn.GetParameter(1)
        self.noise=fcn.GetParameter(2)
        nPed=A/self.hPhD.GetBinWidth(1) * sqrt(2*pi) * self.noise
        nDarkPed=1
        nDarkTot=1
        if self.hPhD0:
            self.Fit0Peak()
            fcn=self.hPhD0.GetFunction("gaus")
            A=fcn.GetParameter(0)
            mu=fcn.GetParameter(1)
            noise=fcn.GetParameter(2)
            nDarkPed=A/self.hPhD0.GetBinWidth(1) * sqrt(2*pi) * noise
            nDarkTot=self.hPhD0.GetEntries() 
        self.npe = -log(nPed/self.hPhD.GetEntries()) + log(nDarkPed/nDarkTot)
        return self.npe

    def GetPDE(self,pulser):
        return self.npe/pulser.muGamma

    def GetNoise(self):
        return npefcn.GetParameter("sig0")
    def GetENF(self):
        return npefcn.GetParameter("enf")
    def GetGain(self):
        return npefcn.GetParameter("gain")

    

if __name__ == "__main__":
    pulser=LightSource()
    s=SIPM()
    s.dcr=2
    gate=1.5  # gate width relative to pulse width
    hLight=s.SimPhD(pulser,gate)
    hDark=s.SimDarkPhD(gate)
    ana=PhDAnalyzier(hLight,hDark)
    npe=ana.CalcNpe()
    ana.FitPhD()


    screenY=TGClient.Instance().GetDisplayHeight()
    c1=TCanvas("results","results",int(screenY*.75),int(screenY*.75))
    c1.Divide(2,2)
    c1.cd(1)
    ana.hPhD.Draw()
    c1.cd(2)
    pulser.pulsehist.Draw()
    c1.cd(3)
    if ana.hPhD0: ana.hPhD0.Draw()
    else: hDark.Draw()
    c1.cd(4)
    #ana.FitPhD()
    npefcn.Draw()


    s.Print()

    print "Calculated PDE =",'{0:.1f}%'.format(ana.GetPDE(pulser)*100)
    print "Calculated noise =",'{0:.2f}'.format(ana.GetNoise())
    print "Calculated ENF =",'{0:.2f}'.format(ana.GetENF())
    print "Calculated Gain =",'{0:.2f}'.format(ana.GetGain())
    
    raw_input("Press Enter to continue...")