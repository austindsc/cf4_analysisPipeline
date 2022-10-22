import ROOT
from ROOT import TFile, TObject, TVectorD


hFile = TFile("testOutput.root", 'RECREATE', 'testWrite.py output file')

hFile.cd()

num = TVectorD(1)
num[0] = 200.

num.Write("blinePre")

hFile.Close()
