

import numpy as np
from array import array
import ROOT
from ROOT import TFile, TTree

defOpen = "/Users/austindestcroix/Desktop/school_holders/projects/sbc/CF4/procData/npv1_202203101333.root"

file = ROOT.TFile.Open(defOpen, "READ")

file.ls()

tt = file.Get("tt")

tt.Show(0)

bline = np.empty((3), dtype = "float32")
tt.SetBranchAddress("blineV", bline)

tt.GetEntry(0)
print(bline)

tt.Draw("blineV[0]")
# c1.Update()

# ROOT.gApplication.Run()