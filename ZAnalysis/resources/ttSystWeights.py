from DataFormats.FWLite import Handle, Runs
import ROOT

lheruninfo=Handle('LHERunInfoProduct')

runs=Runs('root://eoscms//eos/cms/store/group/phys_top/gkrintir/TopHI/HINPbPbAutumn18DR_skims/WJetsToLNu_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/B02D2CBB-AFBB-B741-B810-6D00B71B2AAD.root')
#runs=Runs('root://cms-xrd-global.cern.ch//store/himc/HINPbPbAutumn18DR/TT_TuneCP5_HydjetDrumMB_5p02TeV-powheg-pythia8/AODSIM/mva98_103X_upgrade2018_realistic_HI_v11-v1/270000/EF779097-4F67-024E-A13E-5EE36B3E5101.root')
for r in runs:
    r.getByLabel('externalLHEProducer',lheruninfo)
    it=lheruninfo.product().headers_begin()
    while it!=lheruninfo.product().headers_end():
        lines=it.lines()
        allowPrint=False
        wgtCtr=0
        for i in xrange(0,lines.size()):
            linestr=lines.at(i)
            if '<weightgroup' in linestr : allowPrint=True
            if '</weightgroup' in linestr : allowPrint=False
            if not allowPrint : continue
            if 'weightgroup' in linestr :
                print '*'*50
                print linestr
                print '*'*50
            else:
                if not 'weight' in linestr : continue
                print wgtCtr,linestr
                wgtCtr+=1
        it.next()
