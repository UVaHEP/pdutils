from ROOT import *
from pdeTest import *
import sys, os


if len(sys.argv)<2:
    print "No input file given"
    sys.exit()

tf=TFile(sys.argv[1])
hLight=tf.Get("hpulses1")
hDark=tf.Get("hpulses0")

hLight.Draw()
hDark.Draw("same")

raw_input("Press Enter to continue...")


ana=PhDAnalyzier(hLight.Clone(),hDark.Clone())
npe=ana.CalcNpe()

ana.FitPhD() # do a nice fit to the peaks

screenY=TGClient.Instance().GetDisplayHeight()
c1=TCanvas("results",os.path.basename(sys.argv[1]),int(screenY*.75),int(screenY*.75))
c1.Divide(1,2)
c1.cd(1)
ana.hPhD.Draw()
ana.fcn0.Draw("same")
#ana.fcn0.Print()
c1.cd(2)
ana.hPhD0.Draw()


print "Mean NPE detectected",npe

raw_input("Press Enter to continue...")
