import sys
sys.path.append("/uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/python_analysis/analysisTools")
from analysisImports import *
setLibParams()

cfgDir = "/uscms/home/sbrightt/nobackup/iDM/iDMe_analysis/CMSSW_10_6_26/src/iDMe/python_analysis/configs/"
plt.style.use(cfgDir+'plt_settings.mplstyle')

def genMatch_electrons(file):
    print("******************************************")
    print("Analyzing "+file)
    print("------------------------------------------")
    
    importlib.reload(aTools)
    events = aTools.loadNano(file,"outT")
    nev = len(events)
    
    genE_msk = events.GenPart.ID == 11
    genP_msk = events.GenPart.ID == -11

    genE = ak.flatten(events.GenPart[genE_msk])
    genP = ak.flatten(events.GenPart[genP_msk])

    genE_pt = genE.pt
    genE_eta = genE.eta
    genE_phi = genE.phi

    genP_pt = genP.pt
    genP_eta = genP.eta
    genP_phi = genP.phi

    es = events.Electron
    les = events.LptElectron

    e_pt = es.pt
    e_eta = es.eta
    e_phi = es.phi

    le_pt = les.pt
    le_eta = les.eta
    le_phi = les.phi

    genE_pos = ak.zip({"x":genE_eta,"y":genE_phi},with_name="TwoVector")
    genP_pos = ak.zip({"x":genP_eta,"y":genP_phi},with_name="TwoVector")
    e_pos = ak.zip({"x":e_eta,"y":e_phi},with_name="TwoVector")
    le_pos = ak.zip({"x":le_eta,"y":le_phi},with_name="TwoVector")
    
    dRe_e = ak.fill_none(ak.min((e_pos - genE_pos).r,axis=1),999)
    dRe_le = ak.fill_none(ak.min((le_pos - genE_pos).r,axis=1),999)

    dRp_e = ak.fill_none(ak.min((e_pos - genP_pos).r,axis=1),999)
    dRp_le = ak.fill_none(ak.min((le_pos - genP_pos).r,axis=1),999)
    
    e_matchReg = dRe_e < 0.1
    e_matchLow = dRe_le < 0.1

    p_matchReg = dRp_e < 0.1
    p_matchLow = dRp_le < 0.1

    e_match = e_matchReg | e_matchLow
    p_match = p_matchReg | p_matchLow

    regreg = e_matchReg & p_matchReg
    lowreg = e_matchLow & p_matchReg
    reglow = e_matchReg & p_matchLow
    lowlow = e_matchLow & p_matchLow

    bothMatch = e_match & p_match

    print("Electron matched to regular in {0} events ({1:.2f}%)".format(ak.count_nonzero(e_matchReg),100*ak.count_nonzero(e_matchReg)/nev))
    print("Electron matched to low-pT in {0} events ({1:.2f}%)".format(ak.count_nonzero(e_matchLow),100*ak.count_nonzero(e_matchLow)/nev))
    print("Electron matched to either in {0} events ({1:.2f}%)".format(ak.count_nonzero(e_match),100*ak.count_nonzero(e_match)/nev))
    print("---------------------------------------------------------")
    print("Positron matched to regular in {0} events ({1:.2f}%)".format(ak.count_nonzero(p_matchReg),100*ak.count_nonzero(p_matchReg)/nev))
    print("Positron matched to low-pT in {0} events ({1:.2f}%)".format(ak.count_nonzero(p_matchLow),100*ak.count_nonzero(p_matchLow)/nev))
    print("Positron matched to either in {0} events ({1:.2f}%)".format(ak.count_nonzero(p_match),100*ak.count_nonzero(p_match)/nev))
    print("---------------------------------------------------------")
    print("e = reg, p = reg in {0} events ({1:.2f}%)".format(ak.count_nonzero(regreg),100*ak.count_nonzero(regreg)/nev))
    print("e = low, p = reg in {0} events ({1:.2f}%)".format(ak.count_nonzero(lowreg),100*ak.count_nonzero(lowreg)/nev))
    print("e = reg, p = low in {0} events ({1:.2f}%)".format(ak.count_nonzero(reglow),100*ak.count_nonzero(reglow)/nev))
    print("e = low, p = low in {0} events ({1:.2f}%)".format(ak.count_nonzero(lowlow),100*ak.count_nonzero(lowlow)/nev))
    print("---------------------------------------------------------")
    print("e & p both matched to something in {0} events ({1:.2f}%)".format(ak.count_nonzero(bothMatch),100*ak.count_nonzero(bothMatch)/nev))