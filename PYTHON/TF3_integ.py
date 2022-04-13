import ROOT

func = ROOT.TF3("func", "1*x+2*y+3*z", 0, 10, 0, 10, 0, 10)

print("Integral : ", func.Integral(0, 1, 0, 1, 0, 2) )
