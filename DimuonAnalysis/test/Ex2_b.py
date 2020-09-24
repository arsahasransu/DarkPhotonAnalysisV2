import ROOT as r

file = r.TFile("scout_cmsdas.root","READ")
h_dimuM = file.Get("dimuonMass")

h_dimuM.GetXaxis().SetRange(50, 65)
h_dimuM.Fit("gaus","","",2.8,3.3)

h_dimuM.SetMaximum(28000)

c1 = r.TCanvas()
h_dimuM.Draw()
c1.SaveAs("JpsimassFit.pdf")
