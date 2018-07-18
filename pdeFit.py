from ROOT import *
from pdeTest import *
import sys


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
ana.FitPhD()

screenY=TGClient.Instance().GetDisplayHeight()
c1=TCanvas("results","results",int(screenY*.75),int(screenY*.35))
c1.Divide(1,2)
c1.cd(1)
ana.hPhD.Draw()
c1.cd(2)
ana.hPhD0.Draw()


raw_input("Press Enter to continue...")
