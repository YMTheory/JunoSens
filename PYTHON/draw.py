import numpy as np
import matplotlib.pyplot as plt
import ROOT

f1 = ROOT.TFile("PredEvisSpec_NO.root", "read")
hNO = f1.Get("hPredEvisSpec")

f2 = ROOT.TFile("PredEvisSpec_IO.root", "read")
hIO = f2.Get("hPredEvisSpec")


cont1, cont2, cent = [], [], []
for i in range(340):
    cent.append(hNO.GetBinCenter(i+1))
    cont1.append(hNO.GetBinContent(i+1))
    cont2.append(hIO.GetBinContent(i+1))



fig, ax = plt.subplots()
ax.plot(cent, cont1, label="Noraml Ordering ")
ax.plot(cent, cont2, label="Inverted Ordering")

ax.legend()
plt.show()




