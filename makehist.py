import ROOT
import time

path="/eos/user/m/melu/Electron_SF/"

SingleEG_names = ROOT.std.vector('string')()
for f in ["Tree_SingleElectron_Run2017B_0.root","Tree_SingleElectron_Run2017C_0.root","Tree_SingleElectron_Run2017D_0.root","Tree_SingleElectron_Run2017E_0.root","Tree_SingleElectron_Run2017E_1.root","Tree_SingleElectron_Run2017F_0.root","Tree_SingleElectron_Run2017F_1.root"]:
  SingleEG_names.push_back(path+f)

DY_name = ROOT.std.vector('string')()
for f in ["Tree_DYJetsToLL_M50_Electron2017_0.root","Tree_DYJetsToLL_M50_Electron2017_1.root","Tree_DYJetsToLL_M50_Electron2017_2.root"]:
  DY_name.push_back(path+f)

def get_mcEventnumber(filename):
  print 'opening file ', filename
  nevent_temp=0
  for i in range(0,len(filename)):
    ftemp=ROOT.TFile.Open(filename[i])
    htemp=ftemp.Get('Events')
    nevent_temp=htemp.GetEntriesFast()
  return nevent_temp

def makehist(ptbin, etabin):


  lumi = 41480.

  DY_xs = 6077.22
  DY_ev = get_mcEventnumber(DY_name)

  filters_tag="Tag_pt > 35 && abs(Tag_eta) < 2.17 && Tag_mvaFall17V2Iso_WP90"
  filters_other="TnP_mass>60 && TnP_mass<120 && nElectron==2 && (Tag_charge*Probe_charge==-1) && abs(Probe_eta)<2.5 && Probe_pt>10"
  filters_probe=" (Probe_sip3d<8 && abs(Probe_dxy)<0.05 && abs(Probe_dz)<0.1 && Probe_miniPFRelIso_all<0.4 && Probe_mvaFall17V2noIso_WPL && Probe_lostHits<2 && Probe_convVeto && Probe_hoe<0.1 && Probe_eInvMinusPInv>-0.04) && ((abs(Probe_deltaEtaSC+Probe_eta)<1.479 && Probe_sieie<0.011) || (abs(Probe_deltaEtaSC+Probe_eta)>1.479 && Probe_sieie<0.03)) && Probe_mvaTTH>0.25 && Probe_tightCharge==2"
  antifilters_probe="!((Probe_sip3d<8 && abs(Probe_dxy)<0.05 && abs(Probe_dz)<0.1 && Probe_miniPFRelIso_all<0.4 && Probe_mvaFall17V2noIso_WPL && Probe_lostHits<2 && Probe_convVeto && Probe_hoe<0.1 && Probe_eInvMinusPInv>-0.04) && ((abs(Probe_deltaEtaSC+Probe_eta)<1.479 && Probe_sieie<0.011) || (abs(Probe_deltaEtaSC+Probe_eta)>1.479 && Probe_sieie<0.03)) && Probe_mvaTTH>0.25 && Probe_tightCharge==2)"

  filters_pass=filters_tag+" && " + filters_other+ " && " +filters_probe
  filters_fail=filters_tag+" && " + filters_other+ " && " +antifilters_probe

  eta_bin=[0.0, 0.8, 1.4442, 1.566, 2.0, 2.5]
  etabin_names=['m0p0','p0p8','p1p4442','p1p566','p2p0','p2p5']
  pt_bin=[10,20,35,50,100,200,500]
  ptbin_names=['10','20','35','50','100','200','500']

  outputname='Pt'+ptbin_names[ptbin]+'To'+ptbin_names[ptbin+1]+'Eta'+etabin_names[etabin]+'To'+etabin_names[etabin+1]+'.root'
  fileout=ROOT.TFile(outputname,"RECREATE")


  filters_pass_final=filters_pass+ "&& Probe_pt>" + str(pt_bin[ptbin]) + " && Probe_pt<" +str(pt_bin[ptbin+1])+"&& abs(Probe_eta)<" +str(eta_bin[etabin+1]) + " && abs(Probe_eta)>" + str(eta_bin[etabin])
  filters_fail_final=filters_fail+ "&& Probe_pt>" + str(pt_bin[ptbin]) + " && Probe_pt<" +str(pt_bin[ptbin+1])+"&& abs(Probe_eta)<" +str(eta_bin[etabin+1]) + " && abs(Probe_eta)>" + str(eta_bin[etabin])
  title_pass_temp="Pass_eta"+etabin_names[etabin]+"To"+etabin_names[etabin+1]+"pt"+ptbin_names[ptbin]+"To"+ptbin_names[ptbin+1]
  title_fail_temp="Fail_eta"+etabin_names[etabin]+"To"+etabin_names[etabin+1]+"pt"+ptbin_names[ptbin]+"To"+ptbin_names[ptbin+1]

  df_DY_tree_pass = ROOT.RDataFrame("Events",DY_name)
  df_DY_pass = df_DY_tree_pass.Filter(filters_pass_final)
  df_DY_pass_histo = df_DY_pass.Histo1D(("TnP_mass_DYpass",title_pass_temp ,60,60,120), "TnP_mass",'puWeight')

  df_DY_tree_fail = ROOT.RDataFrame("Events",DY_name)
  df_DY_fail = df_DY_tree_fail.Filter(filters_fail_final)
  df_DY_fail_histo = df_DY_fail.Histo1D(("TnP_mass_DYfail",title_fail_temp ,60,60,120), "TnP_mass",'puWeight')

  df_SingleEle_tree_pass = ROOT.RDataFrame("Events", SingleEG_names)
  df_SingleEle_pass = df_SingleEle_tree_pass.Filter(filters_pass_final)
  df_SingleEle_pass_histo = df_SingleEle_pass.Histo1D(("TnP_mass_EGpass",title_pass_temp,60,60,120), "TnP_mass")

  df_SingleEle_tree_fail = ROOT.RDataFrame("Events", SingleEG_names)
  df_SingleEle_fail = df_SingleEle_tree_fail.Filter(filters_fail_final)
  df_SingleEle_fail_histo = df_SingleEle_fail.Histo1D(("TnP_mass_EGfail",title_fail_temp,60,60,120), "TnP_mass")

  df_DY_pass_histo.Draw()
  df_DY_fail_histo.Draw()
  df_SingleEle_pass_histo.Draw()
  df_SingleEle_fail_histo.Draw()

  fileout.cd()
  df_DY_pass_histo.Write()
  df_DY_fail_histo.Write()
  df_SingleEle_pass_histo.Write()
  df_SingleEle_fail_histo.Write()

  fileout.Close()

if __name__ == "__main__":
  start = time.time()
  start1 = time.clock()
  for ptbin in range(6):
    for etabin in range(5):
      makehist(ptbin,etabin)
  end = time.time()
  end1 = time.clock()
  print "wall time:", end-start
  print "process time:", end1-start1

