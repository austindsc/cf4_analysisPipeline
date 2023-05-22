
import ROOT as r
import numpy as np

class A(object):
  def __init__(self):
    self.b = 0
    
  def myFun(self, x, par):
    return par[0] * (x[0]-par[1])

a = A()
afunc = a.myFun
f1 = r.TF1('f1', afunc, 0, 10, 2)
f1.SetParameters(1, 0)
f1.Eval(0)

# Define the user-defined function
def falone(x, p):
  return p[0]*(x[0] - p[1])


# Create a TF1 object from the user-defined function
f2 = r.TF1("f2", falone, 0, 10, 2)
f2.SetParameters(1, 0)
f2.Eval(69)

# # Set the initial parameters for the function
# my_tf1.SetParameters(1, 0.1)

# Draw the function on a canvas
# canvas = r.TCanvas("canvas")
# my_tf1.Draw()
# canvas.Draw()