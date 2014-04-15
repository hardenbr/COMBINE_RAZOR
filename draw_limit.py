from  optparse  import OptionParser
import ROOT as rt
import array, os
import random, sys
import rootlogon

rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptTitle(0)
rootlogon.style()

input_file = sys.argv[1]

file = rt.TFile(input_file,"READ")
exp = file.Get("exp_graph")
exp_0 = file.Get("exp_0_graph")
exp_p = file.Get("exp_p_graph")
exp_m = file.Get("exp_m_graph")
obs = file.Get("obs_graph")
obs_p = file.Get("obs_p_graph")
obs_m = file.Get("obs_m_graph")

exp_hist = file.Get("exp_theory")
#exp_hist.GetXaxis().SetRangeUser(450,2050)
canvas = rt.TCanvas("limit_plot","limit", 1000,900)

high_pad =  rt.TPad("top_pad","top_pad",0,0.8,1,1)
high_pad.SetTopMargin(0)
high_pad.SetBottomMargin(.4)
high_pad.Draw()

high_pad.cd()
lumi_text= rt.TPaveText(.15,.95,.95,0.0125,"NDC")
lumi_text.SetTextFont(42)
lumi_text.SetTextSize(0.172)
lumi_text.SetFillColor(0)
lumi_text.SetTextAlign(12)
lumi_text.AddText("GGM Bino-like #tilde{#chi}^{0}: pp#rightarrow (#tilde{g}#tilde{g})/(#tilde{q}#tilde{q}) #rightarrow (#tilde{\chi}^{0}#rightarrow #tilde{G}#gamma) (#tilde{\chi}^{0}#rightarrow #tilde{G}#gamma)+jet(s)")
lumi_text.AddText("CMS Preliminary #sqrt{s} = 8 TeV, L = 19.7 fb^{-1}")
lumi_text.SetShadowColor(0)
lumi_text.SetLineWidth(3)
lumi_text.Draw()

#switch bad to the full canvas
canvas.cd()

low_pad = rt.TPad("low_bad", "low_pad", 0, 0, 1, 0.8);

low_pad.SetTopMargin(0)
low_pad.SetBottomMargin(.2)
low_pad.Draw()
low_pad.cd()

low_pad.SetGridx(1)
low_pad.SetGridy(1)

exp_hist.GetXaxis().SetTitleSize(.065)
exp_hist.GetYaxis().SetTitleSize(.065)
exp_hist.GetXaxis().SetLabelSize(.05)
exp_hist.GetYaxis().SetLabelSize(.05)
exp_hist.GetXaxis().SetTitle("m_{#tilde{q}} [GeV]")
exp_hist.GetYaxis().SetTitle("m_{#tilde{g}} [GeV]")
exp_hist.Draw("colz")

#low_pad.SetLogz(1)

exp.Draw("F")
exp_0.Draw("L")
exp_m.Draw("L")
exp_p.Draw("L")
#obs.Draw("L")
#obs_p.Draw("L")
#obs_m.Draw("L")

leg =  rt.TLegend(.58,.7,.88,.95,"Razor #gamma#gamma + #geq 1 jet: Exp. #sigma Excluded (fb)")
leg.SetTextSize(.03)
leg.AddEntry(exp_0, "Exp. Exclusion","l")
leg.AddEntry(exp, "Exp. Exclusion #pm 1 #sigma","f")
#leg.AddEntry(obs, "Obs. Exclusion","l")
#leg.AddEntry(obs_p, "Obs. Exclusion #pm 1 #sigma_{theory}","l")
leg.SetFillColor(0)
leg.SetLineWidth(4)
leg.Draw()

#low_pad.BuildLegend()
#lumi_text.Draw("same")

raw_input("RAW INPUT:")
