import ROOT,math,os
ROOT.gROOT.LoadMacro('./RooCMSShape.cc+')

from math import sqrt
from ROOT import RooCMSShape,TCanvas, TPad
import CMSTDRStyle
CMSTDRStyle.setTDRStyle().cd()
import CMSstyle
from array import array

def eleSF(ismc,filename):

#  filein=ROOT.TFile("./Pt10To20Etam0p0Top0p8.root","READ")
  filein=ROOT.TFile(filename,"READ")
  
  dy_pass=ROOT.TH1D()
  dy_fail=ROOT.TH1D()
  EG_pass=ROOT.TH1D()
  EG_fail=ROOT.TH1D()
  filein.GetObject("TnP_mass_DYpass",dy_pass)
  filein.GetObject("TnP_mass_DYfail",dy_fail)
  filein.GetObject("TnP_mass_EGpass",EG_pass)
  filein.GetObject("TnP_mass_EGfail",EG_fail)
  
  DY_pass_error=ROOT.Double(0.)
  DY_fail_error=ROOT.Double(0.)
  EG_pass_error=ROOT.Double(0.)
  EG_fail_error=ROOT.Double(0.)
  DY_pass_total=dy_pass.IntegralAndError(1,60,DY_pass_error)
  DY_fail_total=dy_fail.IntegralAndError(1,60,DY_fail_error)
  EG_pass_total=EG_pass.IntegralAndError(1,60,EG_pass_error)
  EG_fail_total=EG_fail.IntegralAndError(1,60,EG_fail_error)
  
  x = ROOT.RooRealVar("x", "x", 60, 120)
  x.setRange("FitRange", 61, 119)
  
  # make RootDataHist 
  GenPass = ROOT.RooDataHist("GenPass","GenPass",ROOT.RooArgList(x),dy_pass)
  GenFail = ROOT.RooDataHist("GenFail","GenFail",ROOT.RooArgList(x),dy_fail)
  DataPass = ROOT.RooDataHist("DataPass","DataPass",ROOT.RooArgList(x),EG_pass)
  DataFail = ROOT.RooDataHist("DataFail","DataFail",ROOT.RooArgList(x),EG_fail)
  
  ZPassShape = ROOT.RooHistPdf("ZPassShape","ZPassShape",ROOT.RooArgSet(x), GenPass)
  ZFailShape = ROOT.RooHistPdf("ZFailShape","ZFailShape",ROOT.RooArgSet(x), GenFail)
  
  # gauss smearing of DY shape
  meanPass = ROOT.RooRealVar("meanP", "mean of Reg", -0.0, -5.0, 5.0)
  sigmaPass = ROOT.RooRealVar("sigmaP", "width of Reg", 0.9,0.5,5.0)
  gauss_Pass = ROOT.RooGaussian("gaussP", "gaussian Reg", x, meanPass, sigmaPass)

  meanFail = ROOT.RooRealVar("meanF", "mean of Reg", -0.0, -5.0, 5.0)
  sigmaFail = ROOT.RooRealVar("sigmaF", "width of Reg", 0.9,0.5,5.0)
  gauss_Fail = ROOT.RooGaussian("gaussF", "gaussian Reg", x, meanFail, sigmaFail)
  
  sig_pass = ROOT.RooFFTConvPdf("sigP", "signal shape", x, ZPassShape, gauss_Pass)
  sig_fail = ROOT.RooFFTConvPdf("sigF", "signal shape", x, ZFailShape, gauss_Fail)
  
  # parameter of RooCMSShape, pdf of background
  acmsP = ROOT.RooRealVar("acmsP", "acms", 60.,50.,80.)
  betaP = ROOT.RooRealVar("betaP", "beta", 0.05,0.01,0.08)
  gammaP = ROOT.RooRealVar("gammaP", "gamma", 0.1, -2, 2)
  peakP = ROOT.RooRealVar("peakP", "peak", 90.0)

  bkgP = RooCMSShape("bkgP", "bkg shape", x, acmsP, betaP, gammaP, peakP)

  acmsF = ROOT.RooRealVar("acmsF", "acms", 60.,50.,80.)
  betaF = ROOT.RooRealVar("betaF", "beta", 0.05,0.01,0.08)
  gammaF = ROOT.RooRealVar("gammaF", "gamma", 0.1, -2, 2)
  peakF = ROOT.RooRealVar("peakF", "peak", 90.0)

  bkgF = RooCMSShape("bkgF", "bkg shape", x, acmsF, betaF, gammaF, peakF)
  
  
  nSigP = ROOT.RooRealVar("nSigP","nSigP",0.9*EG_pass_total,0.5*EG_pass_total,1.5*EG_pass_total)
  nBkgP = ROOT.RooRealVar("nBkgP","nBkgP",0.1*EG_pass_total,0.,1.5*EG_pass_total)
  
  nSigF = ROOT.RooRealVar("nSigF","nSigF",0.9*EG_fail_total,0.5*EG_fail_total,1.5*EG_fail_total)
  nBkgF = ROOT.RooRealVar("nBkgF","nBkgF",0.1*EG_fail_total,0.,1.5*EG_fail_total)
  
  modelP=ROOT.RooAddPdf("modelP","modelP", ROOT.RooArgList(sig_pass,bkgP), ROOT.RooArgList(nSigP,nBkgP))
  modelF=ROOT.RooAddPdf("modelF","modelF", ROOT.RooArgList(sig_fail,bkgF), ROOT.RooArgList(nSigF,nBkgF))
  
  Pframe = x.frame(ROOT.RooFit.Title("passing probe"))
  Fframe = x.frame(ROOT.RooFit.Title("failing probe"))
  
  if ismc:
    GenPass.plotOn(Pframe)
    rPass = sig_pass.fitTo(GenPass, ROOT.RooFit.Range("FitRange"), ROOT.RooFit.Save())
    sig_pass.plotOn(Pframe, ROOT.RooFit.LineColor(ROOT.kRed))
    GenFail.plotOn(Fframe)
    rFail = sig_fail.fitTo(GenFail, ROOT.RooFit.Range("FitRange"), ROOT.RooFit.Save())
    sig_fail.plotOn(Fframe, ROOT.RooFit.LineColor(ROOT.kRed))
    nTot=DY_pass_total+DY_fail_total
    eff=DY_pass_total/nTot
    e_eff = 1./(nTot*nTot) * sqrt( DY_pass_total*DY_pass_total* DY_fail_error*DY_fail_error + DY_fail_total*DY_fail_total * DY_pass_error*DY_pass_error )


  if not ismc:
    DataPass.plotOn(Pframe)
    rPass = modelP.fitTo(DataPass, ROOT.RooFit.Range("FitRange"), ROOT.RooFit.Save())
    modelP.plotOn(Pframe, ROOT.RooFit.LineColor(ROOT.kRed))
    modelP.plotOn(Pframe, ROOT.RooFit.Components("bkgP"),ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.LineStyle(ROOT.kDashed))
    DataFail.plotOn(Fframe)
    rFail = modelF.fitTo(DataFail, ROOT.RooFit.Range("FitRange"), ROOT.RooFit.Save())
    modelF.plotOn(Fframe, ROOT.RooFit.LineColor(ROOT.kRed))
    modelF.plotOn(Fframe, ROOT.RooFit.Components("bkgF"),ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.LineStyle(ROOT.kDashed))
    nTot=nSigP.getVal()+nSigF.getVal()
    eff=nSigP.getVal()/nTot
    e_eff = 1./(nTot*nTot)*sqrt(nSigP.getVal()*nSigP.getVal()*nSigF.getError()*nSigF.getError() + nSigF.getVal()*nSigF.getVal() * nSigP.getError()*nSigP.getError() )


  c1=TCanvas("TnP","TnP",1200,600)
  c1.Divide(3,1)
  c1.cd(2)
  Pframe.Draw()
  c1.cd(3)
  Fframe.Draw()

  # Add text results
  text1 = ROOT.TPaveText(0.05,0.75,0.3,0.95)
  text1.SetFillColor(0)
  text1.SetBorderSize(0)
  text1.SetTextAlign(12)
  if ismc:
    text1.AddText('* MC Fit status:')
    text1.AddText('passing: '+str(rPass.status())+', '+'failing: '+str(rFail.status()))
    text1.AddText('* Eff = '+str('%1.4f'%eff)+' #pm '+str('%1.4f'%e_eff))
    text1.SetTextSize(0.08)
  else:
    text1.AddText('* Data Fit status:')
    text1.AddText('passing: '+str(rPass.status())+', '+'failing: '+str(rFail.status()))
    text1.AddText('* Eff = '+str('%1.4f'%eff)+' #pm '+str('%1.4f'%e_eff))
    text1.SetTextSize(0.08)

  text2 = ROOT.TPaveText(0.05,0.05,0.3,0.72)
  text2.SetFillColor(0)
  text2.SetBorderSize(0)
  text2.SetTextAlign(12)
  if ismc:
    text2.AddText('  --- parameters ')
    text2.AddText('-meanP = '+str('%1.3f'%meanPass.getVal())+' #pm '+str('%1.3f'%meanPass.getError()))
    text2.AddText('-sigmaP = '+str('%1.3f'%sigmaPass.getVal())+' #pm '+str('%1.3f'%sigmaPass.getError()))
    text2.AddText('-nSigP = '+str('%1.3f'%DY_pass_total)+' #pm '+str('%1.3f'%DY_pass_error))
    text2.AddText('-meanF = '+str('%1.3f'%meanFail.getVal())+' #pm '+str('%1.3f'%meanFail.getError()))
    text2.AddText('-sigmaF = '+str('%1.3f'%sigmaFail.getVal())+' #pm '+str('%1.3f'%sigmaFail.getError()))
    text2.AddText('-nSigF = '+str('%1.3f'%DY_fail_total)+' #pm '+str('%1.3f'%DY_fail_error))
    text2.SetTextSize(0.06)
  else:
    text2.AddText('  --- parameters ')
    text2.AddText('-nBkgP = '+str('%1.3f'%nBkgP.getVal())+' #pm '+str('%1.3f'%nBkgP.getError()))
    text2.AddText('-nSigP = '+str('%1.3f'%nSigP.getVal())+' #pm '+str('%1.3f'%nSigP.getError()))
    text2.AddText('-meanP = '+str('%1.3f'%meanPass.getVal())+' #pm '+str('%1.3f'%meanPass.getError()))
    text2.AddText('-sigmaP = '+str('%1.3f'%sigmaPass.getVal())+' #pm '+str('%1.3f'%sigmaPass.getError()))
    text2.AddText('-acmsP = '+str('%1.3f'%acmsP.getVal())+' #pm '+str('%1.3f'%acmsP.getError()))
    text2.AddText('-betaP = '+str('%1.3f'%betaP.getVal())+' #pm '+str('%1.3f'%betaP.getError()))
    text2.AddText('-gammaP = '+str('%1.3f'%gammaP.getVal())+' #pm '+str('%1.3f'%gammaP.getError()))
    text2.AddText('-nBkgF = '+str('%1.3f'%nBkgF.getVal())+' #pm '+str('%1.3f'%nBkgF.getError()))
    text2.AddText('-nSigF = '+str('%1.3f'%nSigF.getVal())+' #pm '+str('%1.3f'%nSigF.getError()))
    text2.AddText('-meanF = '+str('%1.3f'%meanFail.getVal())+' #pm '+str('%1.3f'%meanFail.getError()))
    text2.AddText('-sigmaF = '+str('%1.3f'%sigmaFail.getVal())+' #pm '+str('%1.3f'%sigmaFail.getError()))
    text2.AddText('-acmsF = '+str('%1.3f'%acmsF.getVal())+' #pm '+str('%1.3f'%acmsF.getError()))
    text2.AddText('-betaF = '+str('%1.3f'%betaF.getVal())+' #pm '+str('%1.3f'%betaF.getError()))
    text2.AddText('-gammaF = '+str('%1.3f'%gammaF.getVal())+' #pm '+str('%1.3f'%gammaF.getError()))
    text2.SetTextSize(0.05)

  c1.cd(1)
  text1.Draw()
  text2.Draw()
  if ismc:
    c1.SaveAs("mc_"+filename+".png")
  else:
    c1.SaveAs("data_"+filename+".png")

  return eff, e_eff

if __name__ == "__main__":

  tdptbin=array('d',[10,20,35,50,100,200,500])
  tdptbin_plain=array('d',[1,2,3,4,5,6,7])
  tdptbinname=['10~20','20~35','35~50','50~100','100~200','200~500']
  tdetabin=array('d',[0.0,0.8,1.4442,1.566,2.0,2.5])

  h2_SF = ROOT.TH2D('EleIDSF', 'EleIDSF', 6, tdptbin_plain, 5, tdetabin)
  h2_SF.Sumw2()
  h2_SF.SetStats(0)
  h2_SF.GetXaxis().SetTitle('Electron P_{T} [GeV]')
  h2_SF.GetYaxis().SetTitle('Electron #||{#eta}')
  h2_SF.SetTitle('')
  for ib in range(1,7):
    h2_SF.GetXaxis().SetBinLabel(ib,tdptbinname[ib-1])

  h2_data = ROOT.TH2D('EleIDDataEff', 'EleIDDataEff', 6, tdptbin, 5, tdetabin)
  h2_data.Sumw2()
  h2_data.SetStats(0)
  h2_data.GetXaxis().SetTitle('Electron P_{T} [GeV]')
  h2_data.GetYaxis().SetTitle('Electron #||{#eta}')
  h2_data.SetTitle('')

  h2_mc = ROOT.TH2D('EleIDMCEff', 'EleIDMCEff', 6, tdptbin, 5, tdetabin)
  h2_mc.Sumw2()
  h2_mc.SetStats(0)
  h2_mc.GetXaxis().SetTitle('Electron P_{T} [GeV]')
  h2_mc.GetYaxis().SetTitle('Electron #||{#eta}')
  h2_mc.SetTitle('')

  ptbinnames=['Pt10To20','Pt20To35','Pt35To50','Pt50To100','Pt100To200','Pt200To500']
  etabinnames=['Etam0p0Top0p8','Etap0p8Top1p4442','Etap1p4442Top1p566','Etap1p566Top2p0','Etap2p0Top2p5']

  for files in os.listdir('./'):
    if not (files.startswith('Pt') and files.endswith('.root')):continue
    eff,eff_err = eleSF(0,files)
    if ptbinnames[0] in files:
      if etabinnames[0] in files:
        h2_data.SetBinContent(1,1,eff)
        h2_data.SetBinError(1,1,eff_err)
      if etabinnames[1] in files:
        h2_data.SetBinContent(1,2,eff)
        h2_data.SetBinError(1,2,eff_err)
      if etabinnames[2] in files:
        h2_data.SetBinContent(1,3,eff)
        h2_data.SetBinError(1,3,eff_err)
      if etabinnames[3] in files:
        h2_data.SetBinContent(1,4,eff)
        h2_data.SetBinError(1,4,eff_err)
      if etabinnames[4] in files:
        h2_data.SetBinContent(1,5,eff)
        h2_data.SetBinError(1,5,eff_err)
    if ptbinnames[1] in files:
      if etabinnames[0] in files:
        h2_data.SetBinContent(2,1,eff)
        h2_data.SetBinError(2,1,eff_err)
      if etabinnames[1] in files:
        h2_data.SetBinContent(2,2,eff)
        h2_data.SetBinError(2,2,eff_err)
      if etabinnames[2] in files:
        h2_data.SetBinContent(2,3,eff)
        h2_data.SetBinError(2,3,eff_err)
      if etabinnames[3] in files:
        h2_data.SetBinContent(2,4,eff)
        h2_data.SetBinError(2,4,eff_err)
      if etabinnames[4] in files:
        h2_data.SetBinContent(2,5,eff)
        h2_data.SetBinError(2,5,eff_err)
    if ptbinnames[2] in files:
      if etabinnames[0] in files:
        h2_data.SetBinContent(3,1,eff)
        h2_data.SetBinError(3,1,eff_err)
      if etabinnames[1] in files:
        h2_data.SetBinContent(3,2,eff)
        h2_data.SetBinError(3,2,eff_err)
      if etabinnames[2] in files:
        h2_data.SetBinContent(3,3,eff)
        h2_data.SetBinError(3,3,eff_err)
      if etabinnames[3] in files:
        h2_data.SetBinContent(3,4,eff)
        h2_data.SetBinError(3,4,eff_err)
      if etabinnames[4] in files:
        h2_data.SetBinContent(3,5,eff)
        h2_data.SetBinError(3,5,eff_err)
    if ptbinnames[3] in files:
      if etabinnames[0] in files:
        h2_data.SetBinContent(4,1,eff)
        h2_data.SetBinError(4,1,eff_err)
      if etabinnames[1] in files:
        h2_data.SetBinContent(4,2,eff)
        h2_data.SetBinError(4,2,eff_err)
      if etabinnames[2] in files:
        h2_data.SetBinContent(4,3,eff)
        h2_data.SetBinError(4,3,eff_err)
      if etabinnames[3] in files:
        h2_data.SetBinContent(4,4,eff)
        h2_data.SetBinError(4,4,eff_err)
      if etabinnames[4] in files:
        h2_data.SetBinContent(4,5,eff)
        h2_data.SetBinError(4,5,eff_err)
    if ptbinnames[4] in files:
      if etabinnames[0] in files:
        h2_data.SetBinContent(5,1,eff)
        h2_data.SetBinError(5,1,eff_err)
      if etabinnames[1] in files:
        h2_data.SetBinContent(5,2,eff)
        h2_data.SetBinError(5,2,eff_err)
      if etabinnames[2] in files:
        h2_data.SetBinContent(5,3,eff)
        h2_data.SetBinError(5,3,eff_err)
      if etabinnames[3] in files:
        h2_data.SetBinContent(5,4,eff)
        h2_data.SetBinError(5,4,eff_err)
      if etabinnames[4] in files:
        h2_data.SetBinContent(5,5,eff)
        h2_data.SetBinError(5,5,eff_err)
    if ptbinnames[5] in files:
      if etabinnames[0] in files:
        h2_data.SetBinContent(6,1,eff)
        h2_data.SetBinError(6,1,eff_err)
      if etabinnames[1] in files:
        h2_data.SetBinContent(6,2,eff)
        h2_data.SetBinError(6,2,eff_err)
      if etabinnames[2] in files:
        h2_data.SetBinContent(6,3,eff)
        h2_data.SetBinError(6,3,eff_err)
      if etabinnames[3] in files:
        h2_data.SetBinContent(6,4,eff)
        h2_data.SetBinError(6,4,eff_err)
      if etabinnames[4] in files:
        h2_data.SetBinContent(6,5,eff)
        h2_data.SetBinError(6,5,eff_err)

    eff_mc, eff_err_mc = eleSF(1,files)
    if ptbinnames[0] in files:
      if etabinnames[0] in files:
        h2_mc.SetBinContent(1,1,eff_mc)
        h2_mc.SetBinError(1,1,eff_err_mc)
      if etabinnames[1] in files:
        h2_mc.SetBinContent(1,2,eff_mc)
        h2_mc.SetBinError(1,2,eff_err_mc)
      if etabinnames[2] in files:
        h2_mc.SetBinContent(1,3,eff_mc)
        h2_mc.SetBinError(1,3,eff_err_mc)
      if etabinnames[3] in files:
        h2_mc.SetBinContent(1,4,eff_mc)
        h2_mc.SetBinError(1,4,eff_err_mc)
      if etabinnames[4] in files:
        h2_mc.SetBinContent(1,5,eff_mc)
        h2_mc.SetBinError(1,5,eff_err_mc)
    if ptbinnames[1] in files:
      if etabinnames[0] in files:
        h2_mc.SetBinContent(2,1,eff_mc)
        h2_mc.SetBinError(2,1,eff_err_mc)
      if etabinnames[1] in files:
        h2_mc.SetBinContent(2,2,eff_mc)
        h2_mc.SetBinError(2,2,eff_err_mc)
      if etabinnames[2] in files:
        h2_mc.SetBinContent(2,3,eff_mc)
        h2_mc.SetBinError(2,3,eff_err_mc)
      if etabinnames[3] in files:
        h2_mc.SetBinContent(2,4,eff_mc)
        h2_mc.SetBinError(2,4,eff_err_mc)
      if etabinnames[4] in files:
        h2_mc.SetBinContent(2,5,eff_mc)
        h2_mc.SetBinError(2,5,eff_err_mc)
    if ptbinnames[2] in files:
      if etabinnames[0] in files:
        h2_mc.SetBinContent(3,1,eff_mc)
        h2_mc.SetBinError(3,1,eff_err_mc)
      if etabinnames[1] in files:
        h2_mc.SetBinContent(3,2,eff_mc)
        h2_mc.SetBinError(3,2,eff_err_mc)
      if etabinnames[2] in files:
        h2_mc.SetBinContent(3,3,eff_mc)
        h2_mc.SetBinError(3,3,eff_err_mc)
      if etabinnames[3] in files:
        h2_mc.SetBinContent(3,4,eff_mc)
        h2_mc.SetBinError(3,4,eff_err_mc)
      if etabinnames[4] in files:
        h2_mc.SetBinContent(3,5,eff_mc)
        h2_mc.SetBinError(3,5,eff_err_mc)
    if ptbinnames[3] in files:
      if etabinnames[0] in files:
        h2_mc.SetBinContent(4,1,eff_mc)
        h2_mc.SetBinError(4,1,eff_err_mc)
      if etabinnames[1] in files:
        h2_mc.SetBinContent(4,2,eff_mc)
        h2_mc.SetBinError(4,2,eff_err_mc)
      if etabinnames[2] in files:
        h2_mc.SetBinContent(4,3,eff_mc)
        h2_mc.SetBinError(4,3,eff_err_mc)
      if etabinnames[3] in files:
        h2_mc.SetBinContent(4,4,eff_mc)
        h2_mc.SetBinError(4,4,eff_err_mc)
      if etabinnames[4] in files:
        h2_mc.SetBinContent(4,5,eff_mc)
        h2_mc.SetBinError(4,5,eff_err_mc)
    if ptbinnames[4] in files:
      if etabinnames[0] in files:
        h2_mc.SetBinContent(5,1,eff_mc)
        h2_mc.SetBinError(5,1,eff_err_mc)
      if etabinnames[1] in files:
        h2_mc.SetBinContent(5,2,eff_mc)
        h2_mc.SetBinError(5,2,eff_err_mc)
      if etabinnames[2] in files:
        h2_mc.SetBinContent(5,3,eff_mc)
        h2_mc.SetBinError(5,3,eff_err_mc)
      if etabinnames[3] in files:
        h2_mc.SetBinContent(5,4,eff_mc)
        h2_mc.SetBinError(5,4,eff_err_mc)
      if etabinnames[4] in files:
        h2_mc.SetBinContent(5,5,eff_mc)
        h2_mc.SetBinError(5,5,eff_err_mc)
    if ptbinnames[5] in files:
      if etabinnames[0] in files:
        h2_mc.SetBinContent(6,1,eff_mc)
        h2_mc.SetBinError(6,1,eff_err_mc)
      if etabinnames[1] in files:
        h2_mc.SetBinContent(6,2,eff_mc)
        h2_mc.SetBinError(6,2,eff_err_mc)
      if etabinnames[2] in files:
        h2_mc.SetBinContent(6,3,eff_mc)
        h2_mc.SetBinError(6,3,eff_err_mc)
      if etabinnames[3] in files:
        h2_mc.SetBinContent(6,4,eff_mc)
        h2_mc.SetBinError(6,4,eff_err_mc)
      if etabinnames[4] in files:
        h2_mc.SetBinContent(6,5,eff_mc)
        h2_mc.SetBinError(6,5,eff_err_mc)

  h2_data.Divide(h2_mc)
  for ix in range(h2_data.GetNbinsX()):
    for iy in range(h2_data.GetNbinsY()):
      h2_SF.SetBinContent(ix+1,iy+1,h2_data.GetBinContent(ix+1,iy+1))
      h2_SF.SetBinError(ix+1,iy+1,h2_data.GetBinError(ix+1,iy+1))

  c1 = TCanvas()
  pad1 = TPad()
  pad1.Draw()
  pad1.cd()
  h2_data.Draw('COLZ TEXT E')
  CMSstyle.SetStyle(pad1)
  pad1.SetRightMargin(0.15)
  c1.SetGridx(False);
  c1.SetGridy(False);
  c1.SaveAs('./SF.png')
  c1.SaveAs('./SF.pdf')
  pad1.Close()

  c2 = TCanvas()
  pad2 = TPad()
  pad2.Draw()
  pad2.cd()
  h2_SF.Draw('COLZ TEXT E')
  CMSstyle.SetStyle(pad2)
  pad2.SetRightMargin(0.15)
  c2.SetGridx(False);
  c2.SetGridy(False);
  c2.SaveAs('./SF_plainX.png')
  c2.SaveAs('./SF_plainX.pdf')
  pad2.Close()

  fout = ROOT.TFile('output.root','recreate')
  fout.cd()
  h2_data.Write()
  fout.Close()
