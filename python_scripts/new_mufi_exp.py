import matplotlib.pyplot as plt
import ROOT,os,subprocess
ROOT.gStyle.SetOptStat(0)
import atexit
import time
import ctypes
from array import array
import numpy as np
import ana_kit as ana
# from multiprocessing import Pool
from scipy.optimize import curve_fit
import json
import pandas as pd
import uproot

def linearFunc(x,intercept,slope):
    y = intercept + slope * x
    return y
# for fixing a root bug,  will be solved in the forthcoming 6.26 release.
ROOT.gInterpreter.Declare("""
#include "MuFilterHit.h"
#include "AbsMeasurement.h"
#include "TrackPoint.h"

void fixRoot(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value, bool mask) {
   std::map<int,float> m = aHit.GetAllSignals(false);
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRootT(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value, bool mask) {
   std::map<int,float> m = aHit.GetAllTimes(false);
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRoot(MuFilterHit& aHit, std::vector<TString>& key,std::vector<float>& value, bool mask) {
   std::map<TString, float> m = aHit.SumOfSignals();
   std::map<TString, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}

void fixRoot(std::vector<genfit::TrackPoint*>& points, std::vector<int>& d,std::vector<int>& k, bool mask) {
      for(std::size_t i = 0; i < points.size(); ++i) {
        genfit::AbsMeasurement*  m = points[i]->getRawMeasurement();
        d.push_back( m->getDetId() );
        k.push_back( int(m->getHitId()/1000) );
    }
}
""")
def pyExit():
       print("Make suicide until solution found for freezing")
       os.system('kill '+str(os.getpid()))
atexit.register(pyExit)

Tkey  = ROOT.std.vector('TString')()
Ikey   = ROOT.std.vector('int')()
Value = ROOT.std.vector('float')()

def map2Dict(aHit,T='GetAllSignals',mask=False):
     if T=="SumOfSignals":
        key = Tkey
     elif T=="GetAllSignals" or T=="GetAllTimes":
        key = Ikey
     else: 
           print('use case not known',T)
           1/0
     key.clear()
     Value.clear()
     if T=="GetAllTimes": ROOT.fixRootT(aHit,key,Value,mask)
     else:                         ROOT.fixRoot(aHit,key,Value,mask)
     theDict = {}
     for k in range(key.size()):
         if T=="SumOfSignals": theDict[key[k].Data()] = Value[k]
         else: theDict[key[k]] = Value[k]
     return theDict

import rootUtils as ut
import shipunit as u
h={}
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,required=True)
parser.add_argument("-p", "--path", dest="path", help="run number",required=False,default="")
parser.add_argument("-f", "--inputFile", dest="fname", help="file name for MC", type=str,default=None,required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
parser.add_argument("-b", "--heartBeat", dest="heartBeat", help="heart beat", default=10000,type=int)
parser.add_argument("-c", "--command", dest="command", help="command", default="")
parser.add_argument("-n", "--nEvents", dest="nEvents", help="number of events", default=-1,type=int)
parser.add_argument("-q", "--nStart", dest="Nstart", help="first event", default=0,type=int)
parser.add_argument("-t", "--trackType", dest="trackType", help="DS or Scifi", default="DS")
parser.add_argument("--remakeScifiClusters", dest="remakeScifiClusters", help="remake ScifiClusters", default=False)
options = parser.parse_args()

# runNr   = str(options.runNumber).zfill(6)
runNr = str(options.runNumber)
# /eos/experiment/sndlhc/convertedData/commissioning/TI18/
# python new_mufi_exp.py -r 6069 -p /eos/experiment/sndlhc/convertedData/physics/2023/ -g geofile_sndlhc_TI18_V2_2023.root
# runNr   = "00" + str(options.runNumber)
print(runNr)
path     = options.path
partitions = 0
if options.runNumber > 0: 
 if path.find('eos')>0:
    path     = os.environ['EOSSHIP']+options.path
    dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+options.path,shell=True) )
# single file, before Feb'22
    data = "sndsw_raw_"+runNr+".root"
    if  dirlist.find(data)<0:
# check for partitions
       dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+options.path+"run_"+runNr,shell=True) )
       while 1>0:
        data = "raw-"+ str(partitions).zfill(4)
        if dirlist.find(data)>0:
            partitions+=1
        else: break
 else:
# check for partitions
       data = "sndsw_raw_"+runNr+".root"
       dirlist = os.listdir(options.path)
       if  not data in dirlist:
          dirlist  = os.listdir(options.path+"run_"+runNr)
          for x in dirlist:
            data = "raw-"+ str(partitions).zfill(4)
            if x.find(data)>0:
               partitions+=1

import SndlhcGeo
if (options.geoFile).find('../')<0: geo = SndlhcGeo.GeoInterface(path+options.geoFile)
else:                                         geo = SndlhcGeo.GeoInterface(options.geoFile[3:])
MuFilter = geo.modules['MuFilter']
Scifi       = geo.modules['Scifi']
nav = ROOT.gGeoManager.GetCurrentNavigator()

A,B,locA,locB,globA,globB    = ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3(),ROOT.TVector3()
latex = ROOT.TLatex()


# MuFilter mapping of planes and bars 
systemAndPlanes  = {1:2,2:5,3:7}
systemAndBars     = {1:7,2:10,3:60}
def systemAndOrientation(s,plane):
   if s==1 or s==2: return "horizontal"
   if plane%2==1 or plane == 6: return "vertical"
   return "horizontal"

systemAndChannels     = {1:[8,0],2:[6,2],3:[1,0]}
sdict                     = {1:'Veto',2:'US',3:'DS'}

freq      = 160.316E6
TDC2ns = 1E9/freq

# some helper functions

def linear_function(x, y):
    return a * x + b * y + c
 
def getAverageZpositions():
   zPos={'MuFilter':{},'Scifi':{}}
   for s in systemAndPlanes:
       for plane in range(systemAndPlanes[s]):
          bar = 4
          p = plane
          if s==3 and (plane%2==0 or plane==7): 
             bar = 90
             p = plane//2
          elif s==3 and plane%2==1:
             bar = 30
             p = plane//2
          MuFilter.GetPosition(s*10000+p*1000+bar,A,B)
          zPos['MuFilter'][s*10+plane] = (A.Z()+B.Z())/2.
   for s in range(1,6):
      mat   = 2
      sipm = 1
      channel = 64
      for o in range(2):
          Scifi.GetPosition(channel+1000*sipm+10000*mat+100000*o+1000000*s,A,B)
          zPos['Scifi'][s*10+o] = (A.Z()+B.Z())/2.
   return zPos

def smallSiPMchannel(i):
    if i==2 or i==5 or i==10 or i==13: return True
    else: return False

def largeSiPMchannel(i):
    return not smallSiPMchannel(i)
def fit_langau(hist,o,bmin,bmax,opt=''):
   if opt == 2:
     params = {0:'Width(scale)',1:'mostProbable',2:'norm',3:'sigma',4:'N2'}
     F = ROOT.TF1('langau',langaufun,0,200,len(params))
   else:
     params = {0:'Width(scale)',1:'mostProbable',2:'norm',3:'sigma'}
     F = ROOT.TF1('langau',twoLangaufun,0,200,len(params))
   for p in params: F.SetParName(p,params[p])
   rc = hist.Fit('landau','S'+o,'',bmin,bmax)
   res = rc.Get()
   if not res: return res
   F.SetParameter(2,res.Parameter(0))
   F.SetParameter(1,res.Parameter(1))
   F.SetParameter(0,res.Parameter(2))
   F.SetParameter(3,res.Parameter(2))
   F.SetParameter(4,0)
   F.SetParLimits(0,0,100)
   F.SetParLimits(1,0,100)
   F.SetParLimits(3,0,10)
   rc = hist.Fit(F,'S'+o,'',bmin,bmax)
   res = rc.Get()
   return res

def twoLangaufun(x,par):
   N1 = langaufun(x,par)
   par2 = [par[0],par[1]*2,par[4],par[3]]
   N2 = langaufun(x,par2)
   return N1+abs(N2)

def  langaufun(x,par):
   #Fit parameters:
   #par[0]=Width (scale) parameter of Landau density
   #par[1]=Most Probable (MP, location) parameter of Landau density
   #par[2]=Total area (integral -inf to inf, normalization constant)
   #par[3]=Width (sigma) of convoluted Gaussian function
   #
   #In the Landau distribution (represented by the CERNLIB approximation),
   #the maximum is located at x=-0.22278298 with the location parameter=0.
   #This shift is corrected within this function, so that the actual
   #maximum is identical to the MP parameter.
#
      # Numeric constants
      invsq2pi = 0.3989422804014   # (2 pi)^(-1/2)
      mpshift  = -0.22278298       # Landau maximum location
#
      # Control constants
      np = 100.0      # number of convolution steps
      sc =   5.0      # convolution extends to +-sc Gaussian sigmas
#
      # Variables
      summe = 0.0
#
      # MP shift correction
      mpc = par[1] - mpshift * par[0]
#
      # Range of convolution integral
      xlow = max(0,x[0] - sc * par[3])
      xupp = x[0] + sc * par[3]
#
      step = (xupp-xlow) / np
#
      # Convolution integral of Landau and Gaussian by sum
      i=1.0
      if par[0]==0 or par[3]==0: return 9999
      while i<=np/2:
         i+=1
         xx = xlow + (i-.5) * step
         fland = ROOT.TMath.Landau(xx,mpc,par[0]) / par[0]
         summe += fland * ROOT.TMath.Gaus(x[0],xx,par[3])
#
         xx = xupp - (i-.5) * step
         fland = ROOT.TMath.Landau(xx,mpc,par[0]) / par[0]
         summe += fland * ROOT.TMath.Gaus(x[0],xx,par[3])
#
      return (par[2] * step * summe * invsq2pi / par[3])

def myPrint(tc,name,withRootFile=True):
     tc.Update()
     tc.Print(name+'-run'+str(options.runNumber)+'.png')
     tc.Print(name+'-run'+str(options.runNumber)+'.pdf')
     if withRootFile: tc.Print(name+'-run'+str(options.runNumber)+'.root')

def makeAnimation(histname,j0=1,j1=2,animated=True, findMinMax=True, lim = 50):
    ut.bookCanvas(h,'ani','',900,800,1,1)
    tc = h['ani'].cd()
    jmin,jmax = j0,j1
    if findMinMax:
       jmin,jmax = 0,0
       for j in range(j0,j1):
            tmp = histname.replace('$$$',str(j))
            if tmp in h:
                  if h[tmp].GetEntries()>lim:
                       jmin = j
                       print(j,tmp,h[tmp].GetEntries())
                       break
       for j in range(j1,j0,-1):
            tmp = histname.replace('$$$',str(j))
            if tmp in h:
                  if h[tmp].GetEntries()>lim:
                       jmax = j
                       break
    for j in range(jmin,jmax):
            tmp = histname.replace('$$$',str(j))
            h[tmp].Draw()
            tc.Update()
            stats = h[tmp].FindObject('stats')
            stats.SetOptFit(1111111)
            h[tmp].Draw()
            if animated: 
               h['ani'].Print('picAni'+str(j)+'.png')
            else:
               rc = input("hit return for next event or q for quit: ")
               if rc=='q': break
    if animated and jmax > jmin: 
           os.system("convert -delay 60 -loop 0 picAni*.png "+histname+".gif")
           os.system('rm picAni*.png')



# initialize 

zPos = getAverageZpositions()

if options.runNumber>0:
              eventChain = ROOT.TChain('rawConv')
              if partitions==0:
                   eventChain.Add(path+'sndsw_raw_'+str(options.runNumber).zfill(6)+'.root')
              else:
                   for p in range(partitions):
                       eventChain.Add(path+'run_'+runNr+'/sndsw_raw-'+str(p).zfill(4)+'.root')

else:
# for MC data and other files
            #   FNAME = options.fname
              f=ROOT.TFile.Open(options.fname)
              if f.Get('rawConv'):   eventChain = f.rawConv
              else:                        eventChain = f.cbmsim
if options.remakeScifiClusters: eventChain.SetBranchStatus("Cluster_Scifi*",0)
rc = eventChain.GetEvent(0)
run      = ROOT.FairRunAna()
ioman = ROOT.FairRootManager.Instance()
ioman.SetTreeName(eventChain.GetName())
outFile = ROOT.TMemFile('dummy','CREATE')
source = ROOT.FairFileSource(eventChain.GetCurrentFile())

if partitions>0:
    for p in range(1,partitions):
                       source.AddFile(path+'run_'+runNr+'/sndsw_raw-'+str(p).zfill(4)+'.root')
run.SetSource(source)
sink = ROOT.FairRootFileSink(outFile)
run.SetSink(sink)

houghTransform = False # under construction, not yet tested
if houghTransform:
   muon_reco_task = SndlhcMuonReco.MuonReco()
   muon_reco_task.SetName("houghTransform")
   run.AddTask(muon_reco_task)
else:
   import SndlhcTracking
   trackTask = SndlhcTracking.Tracking()
   trackTask.SetName('simpleTracking')
   run.AddTask(trackTask)

#avoiding some error messages
xrdb = ROOT.FairRuntimeDb.instance()
xrdb.getContainer("FairBaseParSet").setStatic()
xrdb.getContainer("FairGeoParSet").setStatic()

run.Init()
if partitions>0:  eventTree = ioman.GetInChain()
else:                 eventTree = ioman.GetInTree()

OT = ioman.GetSink().GetOutTree()
Reco_MuonTracks = trackTask.fittedTracks
Cluster_Scifi         = trackTask.clusScifi
# wait for user action 

def help():
    print(" following methods exist")
    print("     make and plot  hitmaps, signal distributions for MuFIlter and Scifi:")
    print("              Scifi_hitMaps(Nev) and Mufi_hitMaps(Nev)     if no number of events is given, loop over all events in file.")
    print(" ")
    print("     plot time between events and time since first event")
    print("              eventTime(Nev=-1)")
    print(" ")
    print("     MuFilter residuals, efficiency estimates with DS or Scifi tracks")
    print("              Mufi_Efficiency(Nev=-1,optionTrack='DS' or 'Scifi")
    print(" ")
    print("     analyze and plot historgams made withMufi_Efficiency ")
    print("              analyze_EfficiencyAndResiduals(readHists=False), with option True, histograms are read from file")
    print(" ")
    print("     Scifi unbiased residuals for an optional set of alignment constants")
    print("              Scifi_residuals(Nev=-1,NbinsRes=100,xmin=-500.,alignPar=False), xmin = low edge of residual histos in microns")
    print(" ")
    print("     Minimization of Scifi alignment constants")
    print("              minimizeAlignScifi(first=True,level=1,minuit=False)")
    print(" ")
    print("     first attempt to look at hadronic showers ")
    print("              USshower(Nev=-1)")
    print(" ")
    print("     make histograms with QDC distributions and fit all distributions with Langau ")
    print("              mips()")
    print("     plot the above distributions directly or reading from file ")
    print("              plotMips(readhisto=True)")
    print(" ")
    print("     beam illumination ")
    print("             scifi_beamSpot() and beamSpot() for MuFilter")
    print(" ")
    print("     rough estimate of Scifi resolution by comparing slopes  build with two stations, x and y projection")
    print("             Scifi_slopes(Nev=-1)")


def Mufi_hitMaps(Nev = options.nEvents):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 for s in systemAndPlanes:
    for l in range(systemAndPlanes[s]):
       ut.bookHist(h,'hitmult_'+str(s*10+l),'hit mult / plane '+str(s*10+l),61,-0.5,60.5)
       ut.bookHist(h,'hit_'+str(s*10+l),'channel map / plane '+str(s*10+l),160,-0.5,159.5)
       if s==3:  ut.bookHist(h,'bar_'+str(s*10+l),'hit map / plane '+str(s*10+l),60,-0.5,59.5)
       else:       ut.bookHist(h,'bar_'+str(s*10+l),'hit map / plane '+str(s*10+l),10,-0.5,9.5)
       ut.bookHist(h,'sig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,200.)
       if s==2:    ut.bookHist(h,'sigS_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,200.)
       ut.bookHist(h,'sigL_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,200.)
       ut.bookHist(h,'sigR_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,200.)
       ut.bookHist(h,'Tsig_'+str(s*10+l),'signal / plane '+str(s*10+l),200,0.0,200.)
       ut.bookHist(h,'occ_'+str(s*10+l),'channel occupancy '+str(s*10+l),100,0.0,200.)
       ut.bookHist(h,'occTag_'+str(s*10+l),'channel occupancy '+str(s*10+l),100,0.0,200.)

 ut.bookHist(h,'leftvsright_1','Veto hits in left / right',10,-0.5,9.5,10,-0.5,9.5)
 ut.bookHist(h,'leftvsright_2','US hits in left / right',10,-0.5,9.5,10,-0.5,9.5)
 ut.bookHist(h,'leftvsright_3','DS hits in left / right',2,-0.5,1.5,2,-0.5,1.5)
 ut.bookHist(h,'leftvsright_signal_1','Veto signal in left / right',100,-0.5,200.,100,-0.5,200.)
 ut.bookHist(h,'leftvsright_signal_2','US signal in left / right',100,-0.5,200.,100,-0.5,200.)
 ut.bookHist(h,'leftvsright_signal_3','DS signal in left / right',100,-0.5,200.,100,-0.5,200.)

 ut.bookHist(h,'dtime','delta event time; dt [ns]',100,0.0,1000.)
 ut.bookHist(h,'dtimeu','delta event time; dt [us]',100,0.0,1000.)
 ut.bookHist(h,'dtimem','delta event time; dt [ms]',100,0.0,1000.)

 ut.bookHist(h,'bs','beam spot',100,-100.,10.,100,0.,80.)
 ut.bookHist(h,'bsDS','beam spot',60,-0.5,59.5,60,-0.5,59.5)
 ut.bookHist(h,'slopes','track slopes',100,-0.1,0.1,100,-0.1,0.1)

 for s in range(1,4):
    ut.bookHist(h,sdict[s]+'Mult','QDCs vs nr hits',200,0.,800.,200,0.,300.)

 N=-1
 Tprev = 0
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)
 t0 =  eventTree.EventHeader.GetEventTime()/freq
 listOfHits = {1:[],2:[],3:[]}
 mult = {}
 for event in eventTree:
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    withX = False
    planes = {}
    listOfHits[1].clear()
    listOfHits[2].clear()
    listOfHits[3].clear()
    for s in systemAndPlanes:
        for l in range(systemAndPlanes[s]):   mult[s*10+l]=0
    for aHit in event.Digi_MuFilterHits:
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()
        if aHit.isVertical():     withX = True
        s = detID//10000
        l  = (detID%10000)//1000  # plane number
        bar = (detID%1000)
        if s>2:
               l=2*l
               if bar>59:
                    bar=bar-60
                    if l<6: l+=1
        mult[s*10+l]+=1
        key = s*100+l
        if not key in planes: planes[key] = {}
        sumSignal = map2Dict(aHit,'SumOfSignals')
        planes[key][bar] = [sumSignal['SumL'],sumSignal['SumR']]
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
# check left/right
        allChannels = map2Dict(aHit,'GetAllSignals')
        for c in allChannels:
           listOfHits[s].append(allChannels[c])
        if nSides==2:
           Nleft    = 0
           Nright = 0
           Sleft    = 0
           Sright = 0
           for c in allChannels:
              if  nSiPMs > c:  # left side
                    Nleft+=1
                    Sleft+=allChannels[c]
              else:
                    Nright+=1
                    Sright+=allChannels[c]
           rc = h['leftvsright_'+str(s)].Fill(Nleft,Nright)
           rc = h['leftvsright_signal_'+str(s)].Fill(Sleft,Sright)
#
        for c in allChannels:
            channel = bar*nSiPMs*nSides + c
            rc = h['hit_'+str(s)+str(l)].Fill( int(channel))
            rc = h['bar_'+str(s)+str(l)].Fill(bar)
            if s==2 and smallSiPMchannel(c) : rc  = h['sigS_'+str(s)+str(l)].Fill(allChannels[c])
            elif c<nSiPMs: rc  = h['sigL_'+str(s)+str(l)].Fill(allChannels[c])
            else             :             rc  = h['sigR_'+str(s)+str(l)].Fill(allChannels[c])
            rc  = h['sig_'+str(s)+str(l)].Fill(allChannels[c])
        allChannels.clear()
#
    for s in listOfHits:
         nhits = len(listOfHits[s])
         qcdsum = 0
         for i in range(nhits):
             rc = h[sdict[s]+'Mult'].Fill(nhits, listOfHits[s][i])
    for s in systemAndPlanes:
        for l in range(systemAndPlanes[s]):   
           rc = h['hitmult_'+str(s*10+l)].Fill(mult[s*10+l])

    maxOneBar = True
    for key in planes:
        if len(planes[key]) > 2: maxOneBar = False
    if withX and maxOneBar:  beamSpot()
    
 S = {1:[1800,800,2,1],2:[1800,1500,2,3],3:[1800,1800,2,4]}
 for s in S:
   ut.bookCanvas(h,'hitmaps' +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
   ut.bookCanvas(h,'barmaps'+str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
   ut.bookCanvas(h,'signal'    +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])
   ut.bookCanvas(h,'Tsignal'   +str(s),' ',S[s][0],S[s][1],S[s][2],S[s][3])

   for l in range(systemAndPlanes[s]):
      n = l+1
      if s==3 and n==7: n=8
      tc = h['hitmaps'+str(s)].cd(n)
      tag = str(s)+str(l)
      h['hit_'+tag].Draw()
      tc = h['barmaps'+str(s)].cd(n)
      h['bar_'+tag].Draw()
      tc = h['signal'+str(s)].cd(n)
      h['sig_'+tag].Draw()
      tc = h['Tsignal'+str(s)].cd(n)
      h['Tsig_'+tag].Draw()

 ut.bookCanvas(h,'hitmult','hit multiplicities per plane',2000,1600,4,3)
 k=1
 for s in systemAndPlanes:
     for l in range(systemAndPlanes[s]):
           tc = h['hitmult'].cd(k)
           tc.SetLogy(1)
           k+=1
           rc = h['hitmult_'+str(s*10+l)].Draw()

 ut.bookCanvas(h,'VETO',' ',1200,1800,1,2)
 for l in range(2):
    tc = h['VETO'].cd(l+1)
    hname = 'hit_'+str(1)+str(l)
    h[hname].SetStats(0)
    h[hname].Draw()
    for n in range(7):
       x = (n+1)*16-0.5
       y = h['hit_'+str(1)+str(l)].GetMaximum()
       lname = 'L'+str(n)+hname
       h[lname] = ROOT.TLine(x,0,x,y)
       h[lname].SetLineColor(ROOT.kRed)
       h[lname].SetLineStyle(9)
       h[lname].Draw('same')

 ut.bookCanvas(h,'USBars',' ',1200,900,1,1)
 colours = {0:ROOT.kOrange,1:ROOT.kRed,2:ROOT.kGreen,3:ROOT.kBlue,4:ROOT.kMagenta,5:ROOT.kCyan,
                  6:ROOT.kAzure,7:ROOT.kPink,8:ROOT.kSpring}
 for i in range(5): 
       h['bar_2'+str(i)].SetLineColor(colours[i])
       h['bar_2'+str(i)].SetLineWidth(2)
       h['bar_2'+str(i)].SetStats(0)
 h['bar_20'].Draw()
 h['bar_21'].Draw('same')
 h['bar_22'].Draw('same')
 h['bar_23'].Draw('same')
 h['bar_24'].Draw('same')
 h['lbar2']=ROOT.TLegend(0.6,0.6,0.99,0.99)
 for i in range(5): 
    h['lbar2'].AddEntry(h['bar_2'+str(i)],'plane '+str(i+1),"f")
 h['lbar2'].Draw()
 for i in range(7): 
       h['hit_3'+str(i)].SetLineColor(colours[i])
       h['hit_3'+str(i)].SetLineWidth(2)
       h['hit_3'+str(i)].SetStats(0)
 h['hit_30'].Draw()
 for i in range(1,7):
     h['hit_3'+str(i)].Draw('same')
 h['lbar3']=ROOT.TLegend(0.6,0.6,0.99,0.99)
 for i in range(7): 
    h['lbar3'].AddEntry(h['hit_3'+str(i)],'plane '+str(i+1),"f")
 h['lbar3'].Draw()

 ut.bookCanvas(h,'LR',' ',1800,900,3,2)
 h['LR'].cd(1)
 h['leftvsright_'+str(1)].Draw('textBox')
 h['LR'].cd(2)
 h['leftvsright_'+str(2)].Draw('textBox')
 h['LR'].cd(3)
 h['leftvsright_'+str(3)].Draw('textBox')
 h['LR'].cd(4)
 h['leftvsright_signal_1'].SetMaximum(h['leftvsright_signal_1'].GetBinContent(10,10))
 h['leftvsright_signal_2'].SetMaximum(h['leftvsright_signal_2'].GetBinContent(10,10))
 h['leftvsright_signal_3'].SetMaximum(h['leftvsright_signal_3'].GetBinContent(10,10))
 h['leftvsright_signal_'+str(1)].Draw('colz')
 h['LR'].cd(5)
 h['leftvsright_signal_'+str(2)].Draw('colz')
 h['LR'].cd(6)
 h['leftvsright_signal_'+str(3)].Draw('colz')

 ut.bookCanvas(h,'LRinEff',' ',1800,450,3,1)
 for s in range(1,4):
   h['lLRinEff'+str(s)]=ROOT.TLegend(0.6,0.54,0.99,0.93)
   name = 'leftvsright_signal_'+str(s)
   h[name+'0Y'] = h[name].ProjectionY(name+'0Y',1,1)
   h[name+'0X'] = h[name].ProjectionX(name+'0X',1,1)
   h[name+'1X'] = h[name].ProjectionY(name+'1Y')
   h[name+'1Y'] = h[name].ProjectionX(name+'1X')
   tc = h['LRinEff'].cd(s)
   tc.SetLogy()
   h[name+'0X'].SetStats(0)
   h[name+'0Y'].SetStats(0)
   h[name+'1X'].SetStats(0)
   h[name+'1Y'].SetStats(0)
   h[name+'0X'].SetLineColor(ROOT.kRed)
   h[name+'0Y'].SetLineColor(ROOT.kGreen)
   h[name+'1X'].SetLineColor(ROOT.kMagenta)
   h[name+'1Y'].SetLineColor(ROOT.kCyan)
   h[name+'0X'].SetMaximum(max(h[name+'1X'].GetMaximum(),h[name+'1Y'].GetMaximum()))
   h[name+'0X'].Draw()
   h[name+'0Y'].Draw('same')
   h[name+'1X'].Draw('same')
   h[name+'1Y'].Draw('same')
   # Fill(Sleft,Sright)
   h['lLRinEff'+str(s)].AddEntry(h[name+'0X'],'left with no signal right',"f")
   h['lLRinEff'+str(s)].AddEntry(h[name+'0Y'],'right with no signal left',"f")
   h['lLRinEff'+str(s)].AddEntry(h[name+'1X'],'left all',"f")
   h['lLRinEff'+str(s)].AddEntry(h[name+'1Y'],'right all',"f")
   h['lLRinEff'+str(s)].Draw()

 ut.bookCanvas(h,'signalUSVeto',' ',1200,1600,3,7)
 s = 1
 l = 1
 for plane in range(2):
     for side in ['L','R','S']:
         tc = h['signalUSVeto'].cd(l)
         l+=1
         if side=='S': continue
         h['sig'+side+'_'+str( s*10+plane)].Draw()
 s=2
 for plane in range(5):
     for side in ['L','R','S']:
         tc = h['signalUSVeto'].cd(l)
         l+=1
         h['sig'+side+'_'+str( s*10+plane)].Draw()
 ut.bookCanvas(h,'signalDS',' ',900,1600,2,7)
 s = 3
 l = 1
 for plane in range(7):
     for side in ['L','R']:
         tc = h['signalDS'].cd(l)
         l+=1
         h['sig'+side+'_'+str( s*10+plane)].Draw()


 for canvas in ['signalUSVeto','LR','USBars']:
      h[canvas].Update()
      myPrint(h[canvas],canvas)
 for canvas in ['hitmaps','barmaps','signal','Tsignal']:
      for s in range(1,4):
         h[canvas+str(s)].Update()
         myPrint(h[canvas+str(s)],canvas+sdict[s])
         
         
def smallVsLargeSiPMs(Nev=-1, offset_qdc = 8.):
 S = 2
 sipm_cut = "large_11"
 sipm_cut_type, num_of_sipms = sipm_cut.split("_") 
 conv_mech = "dir"
 for l in range(systemAndPlanes[S]):
    ut.bookHist(h,'SVSl_'+str(l),'QDC large vs small sum',200,0.,200.,200,0.,200.)
    ut.bookHist(h,'sVSl_'+str(l),'QDC large vs small average',200,0.,200.,200,0.,200.)
    for side in ['L','R']:
          for i1 in range(7):
             for i2 in range(i1+1,8):
               tag=''
               if S==2 and smallSiPMchannel(i1): tag = 's'+str(i1)
               else:                              tag = 'l'+str(i1)
               if S==2 and smallSiPMchannel(i2): tag += 's'+str(i2)
               else:                              tag += 'l'+str(i2)
               ut.bookHist(h,'cor'+tag+'_'+side+str(l),'cor_'+tag+'_'+side+str(l),200,-20.,200.,200,-20.,200.)
               for bar in range(systemAndBars[S]):
                     ut.bookHist(h,'cor'+tag+'_'+side+str(l)+str(bar),'cor_'+tag+'_'+side+str(l)+str(bar) + ";Signal [QDC];Signal [QDC]",200,-20.,200.,200,-20.,200.)

 N=-1
 if Nev < 0 : Nev = eventTree.GetEntries()
 for event in eventTree:
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    for aHit in event.Digi_MuFilterHits:
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()
        s = detID//10000
        bar = (detID%1000)
        if s!=2: continue
        l  = (detID%10000)//1000  # plane number
        sumL,sumS,SumL,SumS = 0,0,0,0
        if conv_mech == "dir":
           allChannels = {i:aHit.GetSignal(i) + offset_qdc if smallSiPMchannel(i) else aHit.GetSignal(i) for i in range(16)}
         #   for a in allChannels:
         #      if allChannels[a] < 0:
         #       print("negative!")
        else:
           allChannels = map2Dict(aHit,'GetAllSignals')
      #   allChannels = map2Dict(aHit,'GetAllSignals',mask=False)
        nS = 0
        nL = 0
        qdc_value = ana.av_qdc(aHit, sipm_cut = sipm_cut_type, sat = False, cut = int(num_of_sipms), MPVs_sipm = {i:0 for i in range(16)}, sqrt = True, offset = False, conv_mech = "any", small_count = False)
        if qdc_value == -999:
            continue
        for c in allChannels:
            if s==2 and smallSiPMchannel(c) : 
                sumS+= allChannels[c]
                nS += 1
            else:                                              
                sumL+= allChannels[c]
                nL+=1
        if nL>0: SumL=sumL/nL
        if nS>0: SumS=sumS/nS
        rc = h['sVSl_'+str(l)].Fill(SumS,SumL)
        rc = h['SVSl_'+str(l)].Fill(sumS/4.,sumL/12.)
#
        for side in ['L','R']:
          offset = 0
          if side=='R': offset = 8
          for i1 in range(offset,offset+7):
             if not i1 in allChannels: continue
             qdc1 = allChannels[i1]
             for i2 in range(i1+1,offset+8):
               if not (i2) in allChannels: continue
               if s==2 and smallSiPMchannel(i1): tag = 's'+str(i1-offset)
               else: tag = 'l'+str(i1-offset)
               if s==2 and smallSiPMchannel(i2): tag += 's'+str(i2-offset)
               else: tag += 'l'+str(i2-offset)
               qdc2 = allChannels[i2]
               rc = h['cor'+tag+'_'+side+str(l)].Fill(qdc1,qdc2)
               rc = h['cor'+tag+'_'+side+str(l)+str(bar)].Fill(qdc1,qdc2)
        allChannels.clear()

 ut.bookCanvas(h,'TSL','',1800,1400,3,2)
 ut.bookCanvas(h,'STSL','',1800,1400,3,2)
 for l in range(systemAndPlanes[S]):
     tc = h['TSL'].cd(l+1)
     tc.SetLogz(1)
     aHist = h['sVSl_'+str(l)]
     aHist.SetTitle(';small SiPM QCD av:large SiPM QCD av')
     nmax = aHist.GetBinContent(aHist.GetMaximumBin())
     aHist.SetMaximum( 0.1*nmax )
     tc = h['sVSl_'+str(l)].Draw('colz')
 myPrint(h['TSL'],"largeSiPMvsSmallSiPM")
 for l in range(systemAndPlanes[S]):
     tc = h['STSL'].cd(l+1)
     tc.SetLogz(1)
     aHist = h['SVSl_'+str(l)]
     aHist.SetTitle(';small SiPM QCD sum/2:large SiPM QCD sum/6')
     nmax = aHist.GetBinContent(aHist.GetMaximumBin())
     aHist.SetMaximum( 0.1*nmax )
     tc = h['SVSl_'+str(l)].Draw('colz')
 myPrint(h['STSL'],"SumlargeSiPMvsSmallSiPM")

 for l in range(systemAndPlanes[S]):
    for side in ['L','R']:
      for bar in range(systemAndBars[S]):
            ut.bookCanvas(h,'cor' + side+str(l)+str(bar),'',1800,1400,7,4)
      ut.bookCanvas(h,'cor'+side+str(l),'',1800,1400,7,4)
      k=1
      for i1 in range(7):
             for i2 in range(i1+1,8):
               tag=''
               if S==2 and smallSiPMchannel(i1): tag = 's'+str(i1)
               else:                              tag = 'l'+str(i1)
               if S==2 and smallSiPMchannel(i2): tag += 's'+str(i2)
               else:                              tag += 'l'+str(i2)
               tc = h['cor'+side+str(l)].cd(k)
               for bar in range(systemAndBars[S]):
                    tc = h['cor' + side+str(l)+str(bar)].cd(k)
                    if bar == 0: h['cor'+tag+'_'+side+str(l)+str(bar)].Draw('colz')
                    else: h['cor'+tag+'_'+side+str(l)+str(bar)].Draw('colzsame')
               k+=1
      for bar in range(systemAndBars[S]):
            myPrint(h['cor' + side+str(l)+str(bar)],'cor' + side+str(l)+str(bar))
      myPrint(h['cor'+side+str(l)],'QDCcor'+side+str(l))




def qdc_dist(Nev = options.nEvents, fit = "langau", title = "US_QDC_distributions"):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 bin_min = 0.
 bin_max = 50.
 hist_list = {}
 hist_list_lr = {}
 for l in range(20, 25):
    #ut.bookHist(h,'sig_'+str(l),'signal / plane '+str(l),200,0.0,50.)
    h =  ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}", 200, bin_min, bin_max)
    h_lr =  ROOT.TH2I("plane" + f"_{l}_lr", "plane" + f"_{l}_lr", 100,-0.5,200.,100,-0.5,200.)
    hist_list[l] = h
    hist_list_lr[l] = h_lr

 for l in range(30, 37):
    h_ds =  ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}", 200, bin_min, bin_max)
    hist_list[l] = h_ds


 N=-1
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)

 for event in eventTree:
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    for aHit in event.Digi_MuFilterHits:
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()

        s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
        l  = (detID%10000)//1000  # plane number
        bar = (detID%1000)

        if s != 3 : continue
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
        qdc_value = ana.av_qdc(aHit)
        if qdc_value == -1:
            continue
        print(s*10 + l)
        hist_list[s*10 + l].Fill(qdc_value)
        q_l, q_r = ana.qdc_left_right(aHit)
        hist_list_lr[s*10 + l].Fill(q_l, q_r)

 # langau fit
#  c  = ROOT.TCanvas("US QDC distribution","US QDC distribution",0,0,1000,800)
#  c.Divide(2,3)
#  for i, plane in enumerate(hist_list.keys()):
#     #print(i)
#     c.cd(i+1)
#    #  if fit == "langau":
#    #    ana.fit_langau(hist_list[plane], str(plane), bin_min, bin_max)
#    #  else:
#    #    hist_list[plane].Fit("pol5")
#     hist_list_lr[plane].Draw("COLZ")
#  c.SaveAs(title + "_" + fit + ".root")

 File = ROOT.TFile.Open(f"{options.runNumber}_run.root", "RECREATE")
 def write_dict_to_file(dict_obj, File):
    for key in dict_obj.keys():
      File.WriteObject(dict_obj[key], dict_obj[key].GetTitle())
 
 hists_all = [hist_list]
 for H in hists_all:
   write_dict_to_file(H, File)
 File.Close()

def qdc_dist_per_bar(Nev = options.nEvents, fit = "langau", title = "US_QDC_distributions"):
 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 bin_min = 0.
 bin_max = 50.
 hist_list = {}
 hist_list_lr = {}
 for l in range(20, 25):
    #ut.bookHist(h,'sig_'+str(l),'signal / plane '+str(l),200,0.0,50.)
    h =  ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}", 200, bin_min, bin_max)
    hist_list[l] = h
    hist_list_lr[l] = {}
    for bar in range(10):
      h_lr =  ROOT.TH2I("plane" + f"_{l}_{bar}_lr", "plane" + f"_{l}_{bar}_lr", 100,-0.5,200.,100,-0.5,200.)
      hist_list_lr[l][bar] = h_lr
 N=-1
 if Nev < 0 : Nev = eventTree.GetEntries()
 eventTree.GetEvent(0)

 for event in eventTree:
    N+=1
    if N%options.heartBeat == 0: print('event ',N,' ',time.ctime())
    if N>Nev: break
    for aHit in event.Digi_MuFilterHits:
        if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()

        s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
        l  = (detID%10000)//1000  # plane number
        bar = (detID%1000)

        if s != 2: continue
        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()
        qdc_value = ana.av_qdc(aHit)
        if qdc_value == -1:
            continue
        hist_list[s*10 + l].Fill(qdc_value)
        q_l, q_r = ana.qdc_left_right(aHit)
        hist_list_lr[s*10 + l][bar].Fill(q_l, q_r)

#  langau fit
#  c  = ROOT.TCanvas("US QDC distribution","US QDC distribution",0,0,1000,800)
#  c.Divide(2,3)
#  for i, plane in enumerate(hist_list.keys()):
#     #print(i)
#     c.cd(i+1)
#    #  if fit == "langau":
#    #    ana.fit_langau(hist_list[plane], str(plane), bin_min, bin_max)
#    #  else:
#    #    hist_list[plane].Fit("pol5")
#     hist_list_lr[plane].Draw("COLZ")
#  c.SaveAs(title + "_" + fit + ".root")
# langau fit
 c  = ROOT.TCanvas("US QDC distribution","US QDC distribution",0,0,1200,1200)
 c.Divide(3,4)
 for i, plane in enumerate(hist_list_lr[23].keys()):
    #print(i)
    c.cd(i+1)
   #  if fit == "langau":
   #    ana.fit_langau(hist_list[plane], str(plane), bin_min, bin_max)
   #  else:
   #    hist_list[plane].Fit("pol5")
    hist_list_lr[23][i].Draw("COLZ")
 c.SaveAs(title + "_" + fit + ".root")
 c.SaveAs(title + "_" + fit + ".pdf")

 File = ROOT.TFile.Open(f"{options.runNumber}_run.root", "RECREATE")
 def write_dict_to_file(dict_obj, File):
    for key in dict_obj.keys():
      File.WriteObject(dict_obj[key], dict_obj[key].GetTitle())
 
 #hists_all = [hist_list, hist_list_lr]
 hists_all = [hist_list_lr[key] for key in hist_list_lr.keys()]
 for H in hists_all:
   write_dict_to_file(H, File)
 File.Close()





def beamSpot():
   trackTask.ExecuteTask()
   Xbar = -10
   Ybar = -10
   for  aTrack in Reco_MuonTracks:
         state = aTrack.getFittedState()
         pos    = state.getPos()
         rc = h['bs'].Fill(pos.x(),pos.y())
         points = aTrack.getPoints()
         keys     = ROOT.std.vector('int')()
         detIDs = ROOT.std.vector('int')()
         ROOT.fixRoot(points, detIDs,keys,True)
         for k in range(keys.size()):
             #                                     m = p.getRawMeasurement()
             detID =detIDs[k] # m.getDetId()
             key = keys[k]          # m.getHitId()//1000 # for mufi
             aHit = eventTree.Digi_MuFilterHits[key]
             if aHit.GetDetectorID() != detID: continue # not a Mufi hit
             s = detID//10000
             l  = (detID%10000)//1000  # plane number
             bar = (detID%1000)
             if s>2: 
               l=2*l
               if bar>59:
                    bar=bar-60
                    if l<6: l+=1
             if s==3 and l%2==0: Ybar=bar
             if s==3 and l%2==1: Xbar=bar
             nSiPMs = aHit.GetnSiPMs()
             nSides  = aHit.GetnSides()
             for p in range(nSides):
                 c=bar*nSiPMs*nSides + p*nSiPMs
                 for i in range(nSiPMs):
                      signal = aHit.GetSignal(i+p*nSiPMs)
                      if signal > 0:
                           rc  = h['Tsig_'+str(s)+str(l)].Fill(signal)
         mom = state.getMom()
         slopeY= mom.X()/mom.Z()
         slopeX= mom.Y()/mom.Z()
         h['slopes'].Fill(slopeX,slopeY)
         if not Ybar<0 and not Xbar<0 and abs(slopeY)<0.01: rc = h['bsDS'].Fill(Xbar,Ybar)



# def OneHitPerUS(DigiHits):
#    USPlanes={}
# 	# for i, aHit in enumerate(eventTree.Digi_MuFilterHit):
#    for i, aHit in enumerate(DigiHits):
#       if not aHit.isValid(): continue
#       detID=aHit.GetDetectorID()
#       subsystem, plane, bar = parseDetID(detID)
#       if subsystem != 2: continue
#       if plane in USPlanes: 
#         continue
#       else:
#         USPlanes[plane] = 0
#       # !!!!!!!!!!!!
#       USPlanes[plane] += 1
#       # !!!!!!!!!!!1 :)
#    for plane in USPlanes:
#       if USPlanes[plane] != 1:
#          return False
#    return True


def scifi_cut(DigiHits, cut_type, hist_scifi_hit, hist_scifi_qdc, hist_us_scifi):
   N_plane_ZY = {i: 0 for i in range(1, 5)}
   N_plane_ZX = {i: 0 for i in range(1, 5)}
   scifi_num_of_hits = 0
   scifi_signal_sum = 0
   for scifiHit in DigiHits:
         scifi_num_of_hits += 1
         scifi_signal_sum += scifiHit.GetSignal(0)
         if scifiHit.isVertical(): 
               hist_scifi_qdc[scifiHit.GetStation() - 1].Fill(scifiHit.GetSignal(0))
               N_plane_ZX[scifiHit.GetStation()] += 1
         else:
               hist_scifi_qdc[scifiHit.GetStation() + 3].Fill(scifiHit.GetSignal(0))
               N_plane_ZY[scifiHit.GetStation()] += 1
   #print('number of scifi hits: ', scifi_num_of_hits, scifi_signal_sum)
   for plane in N_plane_ZY:
      if "scifiY_" + str(plane) not in hist_us_scifi:
         hist_us_scifi["scifiY_" + str(plane)] = []
      if N_plane_ZY[plane] > 0:
         hist_scifi_hit[plane - 1].Fill(N_plane_ZY[plane])
         hist_us_scifi["scifiY_" + str(plane)].append(N_plane_ZY[plane])

   for plane in N_plane_ZX:
      if "scifiX_" + str(plane) not in hist_us_scifi:
         hist_us_scifi["scifiX_" + str(plane)] = []
      if N_plane_ZX[plane] > 0:
         hist_scifi_hit[plane + 3].Fill(N_plane_ZX[plane])
         hist_us_scifi["scifiX_" + str(plane)].append(N_plane_ZX[plane])
   if cut_type == "pion":
      # print(N_plane_ZX)
      if N_plane_ZX[1] == 0:
         return False
      if N_plane_ZY[1] == 0:
         return False

   elif cut_type == "muon":
      if scifi_num_of_hits >= 25:
         return False
      if N_plane_ZX[1] == 0:
         return False
      if N_plane_ZY[1] == 0:
         return False
      return True
def beamSpot():
   trackTask.ExecuteTask()
   Xbar = -10
   Ybar = -10
   for  aTrack in Reco_MuonTracks:
         state = aTrack.getFittedState()
         pos    = state.getPos()
         rc = h['bs'].Fill(pos.x(),pos.y())
         points = aTrack.getPoints()
         keys     = ROOT.std.vector('int')()
         detIDs = ROOT.std.vector('int')()
         ROOT.fixRoot(points, detIDs,keys,True)
         for k in range(keys.size()):
             #                                     m = p.getRawMeasurement()
             detID =detIDs[k] # m.getDetId()
             key = keys[k]          # m.getHitId()//1000 # for mufi
             aHit = eventTree.Digi_MuFilterHits[key]
             if aHit.GetDetectorID() != detID: continue # not a Mufi hit
             s = detID//10000
             l  = (detID%10000)//1000  # plane number
             bar = (detID%1000)
             if s>2: 
               l=2*l
               if bar>59:
                    bar=bar-60
                    if l<6: l+=1
             if s==3 and l%2==0: Ybar=bar
             if s==3 and l%2==1: Xbar=bar
             nSiPMs = aHit.GetnSiPMs()
             nSides  = aHit.GetnSides()
             for p in range(nSides):
                 c=bar*nSiPMs*nSides + p*nSiPMs
                 for i in range(nSiPMs):
                      signal = aHit.GetSignal(i+p*nSiPMs)
                      if signal > 0:
                           rc  = h['Tsig_'+str(s)+str(l)].Fill(signal)
         mom = state.getMom()
         slopeY= mom.X()/mom.Z()
         slopeX= mom.Y()/mom.Z()
         h['slopes'].Fill(slopeX,slopeY)
         if not Ybar<0 and not Xbar<0 and abs(slopeY)<0.01: rc = h['bsDS'].Fill(Xbar,Ybar)


def find_max_bin(hist, range_min, range_max):
   max_bin = 0
   max_bin_content = 0
   for i in range(range_min, range_max+1):
      bin_content = hist.GetBinContent(i)
      if bin_content >= max_bin_content:
         max_bin_content = bin_content
         max_bin = hist.GetBinCenter(i)
      return max_bin


def draw_qdc_sipm(hist_list_sipm, label, wall_num, write_sipm = False):
   mpv_sipms = []
   end_points = []
   sipmID = []
   k_sipm = 0
   with open(f"sipm_endpoints_{wall_num}wall", "w") as f:
      for pl in hist_list_sipm:
         for br in hist_list_sipm[pl]:
            c  = ROOT.TCanvas(f"{label}. Bar {pl*10 + br}", f"{label}. Bar {pl*10 + br}",0,0,1000,1000)
            c.Divide(4,4)
            for i, Si in enumerate(hist_list_sipm[pl][br]):
               c.cd(i+1)
               if write_sipm:
                  res = ana.fit_langau(hist_list_sipm[pl][br][Si], str(Si), 0, 200)
                  sipmID.append(k_sipm)
                  mpv_sipms.append(res)
                  print(pl*10 + br, Si, k_sipm, res, file = f)
                  
               #fit_function = ROOT.TF1("fit_function_1n", "gaus", 100, 200)

               # # Fit the histogram with the Gaussian function
               #hist_list_sipm[pl][br][Si].Fit("fit_function_1n")
               #fit_parameters = fit_function.GetParameters()
               mean = 0
               sigma = 0
               n_init = 10
               last_nonzero_bin = hist_list_sipm[pl][br][Si].FindLastBinAbove(0) 
               print([hist_list_sipm[pl][br][Si].GetBinContent(last_nonzero_bin-i) for i in range(n_init)])
               list_down = [hist_list_sipm[pl][br][Si].GetBinContent(last_nonzero_bin-i) for i in range(n_init)]
               # while sum(list_down) < 10 and n_init < 15:
               list_up = [hist_list_sipm[pl][br][Si].GetBinCenter(last_nonzero_bin-i)*hist_list_sipm[pl][br][Si].GetBinContent(last_nonzero_bin-i) for i in range(10)]
               try:
                  last_bin_value = sum(list_up)/sum(list_down)
                              # Set the range for the search
               except:
                  last_bin_value = -999
               max_bin = find_max_bin(hist_list_sipm[pl][br][Si], 100, 200)
               sipmID.append(k_sipm)
               end_points.append(last_bin_value)
               print(pl*10 + br, Si, k_sipm, last_bin_value, max_bin,  mean, sigma, file = f)
               hist_list_sipm[pl][br][Si].Draw()   
               k_sipm += 1       
            c.Write()

def draw_pcb(hist_pcb, label, write_sipm = False):
   mpv_sipms = []
   sipmID = []
   k_sipm = 0
   pl = 4
   with open("sipm_mpvs", "w") as f:
      c  = ROOT.TCanvas(f"{label}. Plane {pl}. PCB", f"{label}. Plane {pl}. PCB",0,0,2000,2000)
      c.Divide(6,10)
      k = 0
      for br in range(10):
         for i, Si in enumerate(hist_pcb[br]):
            c.cd(k+1)
            k += 1
            hist_pcb[br][Si].Draw()          
      c.Write()

def draw_pcb_1(hist_pcb, label, write_sipm = False):
    # Create a canvas
    pl = 4
    canvas  = ROOT.TCanvas(f"{label}. Plane {pl}. PCB", f"{label}. Plane {pl}. PCB",0,0,2000,2000)

    # Define the number of rows and columns
    rows = 10
    columns = 6

    # Calculate the width and height of each pad
    pad_width = 1.0 / columns
    pad_height = 1.0 / rows
    large_sipms = [i for i in filter(largeSiPMchannel, range(16))]
    # Create pads and arrange them in rows and columns
    for i_row in range(rows):
        for i_col in range(columns):
            # Calculate the pad coordinates
            x1 = i_col * pad_width
            y1 = 1.0 - (i_row + 1) * pad_height
            x2 = (i_col + 1) * pad_width
            y2 = 1.0 - i_row * pad_height

            # Create a pad and draw it
            pad = ROOT.TPad("pad_{0}_{1}".format(i_row, i_col), "", x1, y1, x2, y2)
            pad.Draw()
            # Draw the histogram on the pad
            pad.cd()
            hist_pcb[i_row][large_sipms[i_col]].Draw()
    canvas.Write()

def draw_qdc_bar(hist_list_bar, label, fit_key = False):
   for pl in hist_list_bar:
      c  = ROOT.TCanvas(f"{label}. Plane {pl}", f"{label}. Bar {pl}",0,0,1000,1000)
      c.Divide(3,4)
      for i, br in enumerate(hist_list_bar[pl]):
         c.cd(i+1)
         if fit_key:
            fit_function = ROOT.TF1("fit_function", linear_function, -10, 100, -10, 100, 3)
         hist_list_bar[pl][br].Draw("COLZ")
      c.Write()

def draw_qdc_plane(hist_list, label):
    c  = ROOT.TCanvas(f"{label}",f"{label}",0,0,1000,1000)
    c.Divide(2,3)
    for i, plane in enumerate(hist_list):
      c.cd(i+1)
      hist_list[plane].Draw()
    c.Write()


def draw_qdc_vs_noh(h_qdc_hit, h_qdc_hit_norm):
    c  = ROOT.TCanvas("QDC vs. NOH","QDC vs. NOH",0,0,1000,1000)
    ROOT.gPad.SetLogz()
    h_qdc_hit.Draw('colz')
    c.Write()
    c  = ROOT.TCanvas("QDC/NOH vs. NOH","QDC/NOH vs. NOH",0,0,1000,1000)
    ROOT.gPad.SetLogz()
    h_qdc_hit_norm.Draw('colz')
    c.Write()

def draw_usshower(hist_usshower):
    c  = ROOT.TCanvas("US Shower","US Shower",0,0,1000,1000)
   #  ROOT.gPad.SetLogz()
    hist_usshower.Draw('colz')
    c.Write()


def draw_track_dist(File, hists):
   for hist in hists:
      c  = ROOT.TCanvas(hist.GetTitle(),hist.GetTitle(),0,0,1000,1000)
            # Set the marker color
      hist.SetMarkerSize(1.)

      # Set the error bar color
      hist.SetLineColor(ROOT.kRed)
      hist.SetLineWidth(2)
      hist.SetMarkerStyle(ROOT.kFullCircle)
      hist.Draw("PLC PMC")
      c.SetGrid(1,1)
      c.Write()

def draw_scifi_dist(File,hists, label):
   label = 'scifi_' + label
   c  = ROOT.TCanvas(f"{label}",f"{label}",0,0,1000,1000)
   c.Divide(4,2)
   for i, plane in enumerate(hists):
      tc = c.cd(i+1)
      tc.SetLogy(1)
      hists[plane].Draw()
   c.Write()


def fill_sipm(event, hist_list_sipm, sipm_cut_type, num_of_sipms, sqrt_on, offset, conv_mech, small_count, large_count):
      for aHit in event.Digi_MuFilterHits:
         detID = aHit.GetDetectorID()

         s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
         l  = (detID%10000)//1000  # plane number
         bar = (detID%1000)

         if s != 2 : continue
         qdc_value = ana.av_qdc(aHit, sipm_cut = sipm_cut_type, sat = False, cut = int(num_of_sipms), MPVs_sipm = {i:0 for i in range(16)}, sqrt = sqrt_on, offset = offset, conv_mech = conv_mech, small_count = small_count, large_count = large_count)
         if qdc_value == -999:
            continue
         if conv_mech == "dir":
           allChannels = {i:aHit.GetSignal(i) for i in range(16)}
         else:
           allChannels = map2Dict(aHit,'GetAllSignals')
         allTimes = aHit.GetAllTimes()
         for Si in allChannels:
            hist_list_sipm[s*10 + l][bar]
            if smallSiPMchannel(Si):
               qdc_1 = allChannels[Si]
            else:
               qdc_1 = allChannels[Si]
            # if qdc_1 == -1:
            #    continue
            if offset:
               qdc_1 += 7.
            hist_list_sipm[s*10 + l][bar][Si].Fill(qdc_1)


def track_reco_tool():
      optionTrack = "Scifi"
      for aTrack in Reco_MuonTracks: aTrack.Delete()
      Reco_MuonTracks.Clear()
      if optionTrack=='DS': rc = trackTask.ExecuteTask("DS")
      else                         : rc = trackTask.ExecuteTask("ScifiDS")
      #print("!!!!3")
      if not Reco_MuonTracks.GetEntries()==1:
         #print(f"{N}: {trackTask.fittedTracks.GetEntries()}") 
         return False, None
      theTrack = Reco_MuonTracks[0]
      if not theTrack.getFitStatus().isFitConverged():   # for H8 where only two planes / proj were avaiable
            return False, theTrack
      return True, theTrack
   # only take horizontal tracks
      state = theTrack.getFittedState(0)
      pos   = state.getPos()
      # print(f"Event: {N}", end = " ")
      # pos.Print()
      mom = state.getMom()
      slopeX= mom.X()/mom.Z()
      slopeY= mom.Y()/mom.Z()
      return True

from pathlib import Path

def MIP_study(Nev_st = 0, Nev_en = 1, list_of_events_key = False, title = "default", label_0 = "def", conv_mech = "dir", side_bar_key = False, mc = False, muon = True, small_count = False, large_count = False, fill_sipm_key = False, offset = 0., write_data = False):

 # veto system 2 layers with 7 bars and 8 sipm channels on both ends
 # US system 5 layers with 10 bars and 8 sipm channels on both ends
 # DS system horizontal(3) planes, 60 bars, readout on both sides, single channel
 #                         vertical(4) planes, 60 bar, readout on top, single channel
 

 wall_num = -999
 if list_of_events_key:
   wall_num = 3
   file_path = "/eos/user/u/ursovsnd/SWAN_projects/tests"
   file_list = list(map(str, Path().glob(f"./from_condor/{options.runNumber}_cal/*/*")))
   df_list = []
   for file in file_list:
      df_list.append(pd.read_csv(file))
   df = pd.concat(df_list)
   event_list = df.query(f"wall == {wall_num}")["event"].to_numpy()
 else:
   event_list = range(Nev_st, Nev_en)


 bin_min = -10
 bin_max = 250
 label = "[#sqrt{QDC_L #times QDC_R}]"
 label_sipm = "[QDC]"
 hist_list = {l: ROOT.TH1I("plane" + f"_{l}", "plane" + f"_{l}; {label};", 200, bin_min, bin_max) for l in range(20, 25)}
 hist_list_bar = {l: {bar: ROOT.TH1I("plane" + f"_{l}_{bar}", "plane" + f"_{l}_{bar}; {label};", 200, bin_min, bin_max) for bar in range(10)} for l in range(20, 25)}
 hist_list_sipm = {l: {bar: {sipm: ROOT.TH1I("plane" + f"_{l}_{bar}_{sipm}", "plane" + f"_{l}_{bar}_{sipm}; {label_sipm};", 200, bin_min, bin_max) for sipm in range(16)} for bar in range(10)} for l in range(20, 25)}
 l_us4_lausanne_check_hist = {bar: {sipm: ROOT.TH1I("PCB. " + f"{bar}_{sipm}", "PCB. " + f"_{bar}_{sipm}; {label_sipm};", 200, bin_min, bin_max) for sipm in filter(largeSiPMchannel, range(8))} for bar in range(10)}
 h_qdc_hit = ROOT.TH2I("qdc_vs_hit","QDC vs. Number of fired bars;Number of fired bars;QDC", 100,0,50.,100,0,200)
 h_qdc_hit_norm = ROOT.TH2I("qdc_vs_hit_norm","#sum QDC/Ntot vs. Number of fired bars;Number of fired bars; #sum QDC/Ntot", 100,0,50.,100,0,25)
 hist_usshower = ROOT.TH2I("US Shower", "US Shower; Plane; Bar", 5, -0.5, 4.5, 10, -0.5, 9.5)
 hist_usshower_fraction = ROOT.TH2I("US_Shower_fraction", "US Shower; Plane; Bar", 5, -0.5, 4.5, 10, -0.5, 9.5)
 l_vs_s_hist = {l: {bar: ROOT.TH2I(f"l_vs_s_{l}_{bar}", f"Large_vs_small. Plane_{l}_{bar}; Average signal of Lsipm/hit [QDC]; Average signal of Ssipm/hit [QDC]", 200, -10, 200, 200, -10, 50) for bar in range(10)} for l in range(20, 25)}
 l_vs_s_hist_left = {l: {bar: ROOT.TH2I(f"l_vs_s_{l}_{bar}", f"Large_vs_small. Left. Plane_{l}_{bar}; Average signal of Lsipm/hit [QDC]; Average signal of Ssipm/hit [QDC]", 200, -10, 200, 200, -10, 50) for bar in range(10)} for l in range(20, 25)}
 l_vs_s_hist_right = {l: {bar: ROOT.TH2I(f"l_vs_s_{l}_{bar}", f"Large_vs_small. Right. Plane_{l}_{bar}; Average signal of Lsipm/hit [QDC]; Average signal of Ssipm/hit [QDC]", 200, -10, 200, 200, -10, 50) for bar in range(10)} for l in range(20, 25)}
 l_vs_r_hist ={l: {bar: ROOT.TH2I(f"l_vs_r_{l}_{bar}", f"Left_vs_right. Plane_{l}_{bar}; Average signal of Left sipm/hit [QDC]; Average signal of Right sipm/hit [QDC]", 200, -10, 200, 200, -10, 200) for bar in range(10)} for l in range(20, 25)}
 l_vs_r_hist_large ={l: {bar: ROOT.TH2I(f"l_vs_r_{l}_{bar}_large", f"Left_vs_right. Large. Plane_{l}_{bar}; Average signal of Right sipm/hit [QDC]; Average signal of Left sipm/hit [QDC]", 200, -10, 200, 200, -10, 200) for bar in range(10)} for l in range(20, 25)}
 l_vs_r_hist_small ={l: {bar: ROOT.TH2I(f"l_vs_r_{l}_{bar}_small", f"Left_vs_right. Small. Plane_{l}_{bar}; Average signal of Left sipm/hit [QDC]; Average signal of Right sipm/hit [QDC]", 200, -10, 200, 200, -10, 200) for bar in range(10)} for l in range(20, 25)}
 h_xy_track =  ROOT.TH2I("track_xy","track_xy;x;y", 100,-100,100.,100,-100,100.)
 h_xy_track_res =  ROOT.TH2I("track_xy_res","track_xy_res;x;y", 100,-100,100.,100,-100,100.)
 hist_scifi_hit = {l: ROOT.TH1I(f"scifi_plane_{l}", f"scifi_plane_{'vertical' if l < 4 else 'horizontal'} {l+1 if l < 4 else l - 3}; Number of hits;", 100, 0, 750) for l in range(8)}
 hist_scifi_qdc = {l: ROOT.TH1I(f"scifi_plane_{l}", f"scifi_plane_{'vertical' if l < 4 else 'horizontal'} {l+1 if l < 4 else l - 3}; QDC;", 200, 0, 100) for l in range(8)}
 h_long = ROOT.TH1I("long_shower", f"long_shower; plane;{label}", 5, -0.5, 4.5,)
 h_small_sipm = ROOT.TH1I("small_signals", f"Small_sipms signals; Signal;", 200, -10, 30)
 event_freq_label = "Signal_small" if small_count else "Signal"
 event_freq = ROOT.TH1I("event_freq", f"Hit number. {event_freq_label} > 0; Event number [k]; Number of hits", Nev_en // 10000, 0., (Nev_en // 10000.)*10)
#  h_long_us_scifi = ROOT.TH1I("long_shower_us_scifi", f"long_shower; plane;{label}", 5, -0.5, 4.5,)
 l_vs_r_hist_large_all =  ROOT.TH2I("LRALL","LR;L;R", 100,0,200.,100,0,200.)
 N=-1
 if Nev_en < 0 : Nev_en = eventTree.GetEntries()
#  eventTree.GetEvent(0)
 us_scifi = {}
 us_qdc_list = []
 label = ""
 negative_signal = 0
 non_saturated_bars = [0, 1, 2, 3, 6, 7, 8, 9,]
 if muon:
     optionTrack = "Scifi"
     track_reco = True
     withReco = True
     time_cut = True
     sci_fi_cut = False
 else: 
     extrapolate_qdc = False
     optionTrack = "ScifiDS"
     track_reco = False
     withReco = True
     time_cut = False
     sci_fi_cut = False
 if mc:
   extrapolate_qdc = False
   track_reco = False
   time_cut = False
   offset = 0.
   sci_fi_cut = False
   dir = "any"


 sipm_cut = "large_11"
 sipm_cut_type, num_of_sipms = sipm_cut.split("_") 
 sqrt_on = False
#  sci_fi_cut = False
#  sci_fi_cut = False
 if sci_fi_cut:
    label += f"_{sci_fi_cut}_"
 if track_reco:
  label += "_with-reco_check-track-pos_slope_res-cut-strict"
 label += f"_{optionTrack}"
 label+= f"_{sipm_cut_type}-sipm-cut-{num_of_sipms}"
 if time_cut:
  label+= "time-cut"
 label+= f"_{conv_mech}"
 label+= f"_sqrt-{sqrt_on}"
 label+= f"_{side_bar_key}"
 
 if write_data:
    fout = open(str(options.runNumber), "w")
 if extrapolate_qdc:
   with open(f"/eos/user/u/ursovsnd/SWAN_projects/tests/large_small_cor_{options.runNumber}.json", "r") as read_content: 
      slope_info = pd.DataFrame(json.load(read_content))
      slope_info.columns = ["bar_id", "slopes", "slopes_err", "ch_small", "ch_other"]
#  print(slope_info)

 l_us4_integrated_lausanne_check = {i: {} for i in range(10)}
 for N in map(int, event_list):
#  for event in eventTree:
   #  N+=1
    eventTree.GetEvent(N)
   #  if N <= Nev_st - 1: continue 
    if N%10000 == 0: print('event ',N,' ',time.ctime())
   #  if N >= Nev_en + 1: break

    #if fill_sipm_key: fill_sipm(event, hist_list_sipm, sipm_cut_type, num_of_sipms, sqrt_on, offset, conv_mech, small_count, large_count)
#  if not ana.OneHitUS1(eventTree.Digi_MuFilterHits): continue  
    if sci_fi_cut:
      if not scifi_cut(eventTree.Digi_ScifiHits, "pion" if not muon else "muon", hist_scifi_hit, hist_scifi_qdc, us_scifi): pass
   #  if not ana.OneHitPerUS(eventTree.Digi_MuFilterHits): continue

    if muon:   
      if track_reco:
         track_des, theTrack = track_reco_tool()
         if not track_des:
            continue
      
    else:
      if track_reco:
         track_des, theTrack = track_reco_tool()
         if track_des:
            continue
      # if pos[0] > -37. or pos[0] < -40.:
      #    continue
      # # h_slope.Fill(slopeX, slopeY)
      # if abs(slopeX)>0.25: continue   # 4cm distance, 250mrad = 1cm
      # if abs(slopeY)>0.1: continue
      
      # h_xy_track_slope_abs.Fill(mom.x()/mom.Mag(),mom.y()/mom.Mag())

    number_of_hits = 0
    qdc_per_event = 0
   #  all_signal_num = {i: [0 for _ in range(10)] for i in range(20,25)}
   #  negative_signal_num = {i: [0 for _ in range(10)] for i in range(20,25)}
    all_signal_num = {}
    negative_signal_num = {}
    us_shower = {}
    l_us4_per_event_lausanne_check = {i: {} for i in range(10)}
    for hit_number, aHit in enumerate(eventTree.Digi_MuFilterHits):
      #   if not aHit.isValid(): continue
        detID = aHit.GetDetectorID()

        s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
        l  = (detID%10000)//1000  # plane number
        bar = (detID%1000)

        if s != 2 : continue
        if muon:
         if track_reco:
            
            resx, resy = ana.residual(theTrack, detID, MuFilter, h_xy_track)
            h_xy_track_res.Fill(resx, resy)
            #res = np.sqrt(resx**2 + resy**2)
            res = np.abs(resy)
            if res >= 2.95:
               continue
        else:
         if track_reco and track_des:   
            resx, resy = ana.residual(theTrack, detID, MuFilter, h_xy_track)
            h_xy_track_res.Fill(resx, resy)
            #res = np.sqrt(resx**2 + resy**2)
            res = np.abs(resy)
            if res < 2.95:
               continue           

      #   print("check before timing cut")

        if time_cut:
           if aHit.GetTime() >= 4. or aHit.GetTime() <= 0.:
              continue


        nSiPMs = aHit.GetnSiPMs()
        nSides  = aHit.GetnSides()

      #   print("check before side bar cut")
        if side_bar_key:
         if bar not in non_saturated_bars: continue
      #   print("check after side bar cut")
        if conv_mech == "dir":
          allChannels = {i: aHit.GetSignal(i) + offset if smallSiPMchannel(i) else aHit.GetSignal(i) for i in range(16)}
        else:
          allChannels = map2Dict(aHit,'GetAllSignals')
        allTimes = {i: aHit.GetTime(i) for i in range(16)}
      #   qdc_value = ana.av_qdc(aHit, sipm_cut = sipm_cut_type, sat = False, cut = int(num_of_sipms), MPVs_sipm = {i:0 for i in range(16)}, sqrt = sqrt_on, offset = offset, conv_mech = conv_mech, small_count = small_count, large_count = large_count)

      #   print(qdc_value)
      
      #   print(qdc_value)
      #   if qdc_value < 0:
      #      print("qdc value is negative!!!!")
      #   if qdc_value < 0. and conv_mech == "dir":
      #      continue
        # offset +7 QDC
      #   qdc_value += 7.
      #   print(f"Event: {N}", [aHit.GetSignal(i) for i in range(16)])

        #print(N, allChannels)

        small_list = [] #
        large_list = [] # 

        l_vs_r_list_right = [] #
        l_vs_r_list_left = [] #
        l_vs_r_hist_right_small_list = [] #
        l_vs_r_hist_left_small_list = [] #
        l_vs_r_hist_right_large_list = [] # 
        l_vs_r_hist_left_large_list = [] #
        l_us4_lausanne_check = {i: {} for i in range(10)}
      
        large_firing = 0
        for Si in allChannels:
            qdc_1 = allChannels[Si]
            if qdc_1 == -999:
               continue
            if not smallSiPMchannel(Si):
               large_firing += 1
      #   print("success")
        if large_firing < 11: continue

        if extrapolate_qdc:
            slope, slope_err = 1, 0
            def define_slope(bar, i):
                  try:
                     # print(s*100 + l*10 + bar)
                     slope, slope_err = slope_info.query(f"bar_id == {s*100 + l*10 + bar} & ch_other == {i}").apply(np.mean)["slopes"], slope_info.query(f"bar_id == {s*100 + l*10 + bar} & ch_other == {i}").apply(np.mean)["slopes_err"]
                     return slope, slope_err
                  except:
                     return None, None
            if bar not in non_saturated_bars:
                  small_list_extra = []
                  for Si in list(filter(smallSiPMchannel, allChannels)):
                     if allChannels[Si] != -999:
                        small_list_extra.append(allChannels[Si])
                  if len(small_list_extra) > 0:
                     for Si in allChannels:
                        if allChannels[Si] == -999 or allChannels[Si] > 100:
                           slope, slope_err = define_slope(bar, Si)

                           if slope > 0:
                                 # print(slope)
                                 # print("!!!!!!", bar, Si, slope, np.array(small_list_extra).mean())
                              if smallSiPMchannel(Si):
                                 allChannels[Si] = np.random.normal((np.array(small_list_extra).mean())/slope, np.abs(((np.array(small_list_extra).mean()))/slope)*slope_err/slope)
                              else:
                                 allChannels[Si] = np.random.normal((np.array(small_list_extra).mean()+7)/slope, np.abs(((np.array(small_list_extra).mean()+7))/slope)*slope_err/slope)                                 
                           else:

                              if smallSiPMchannel(Si):
                                 slope, slope_err = 1, 0.1
                                 allChannels[Si] = np.random.normal((np.array(small_list_extra).mean())/slope, np.abs(((np.array(small_list_extra).mean()))/slope)*slope_err*2/slope)
                              else:
                                 slope, slope_err = 0.01, 0.0025
                                 allChannels[Si] = np.random.normal((np.array(small_list_extra).mean()+7)/slope, np.abs(((np.array(small_list_extra).mean()+7))/slope)*slope_err*2/slope)                                                            
                        # print(allChannels)
                  else:
                     continue
      #   print(N, s*100 + l*10 + bar, hit_number, allChannels)
        if l not in all_signal_num:
           all_signal_num[l] = {}
           negative_signal_num[l] = {}
        if bar not in all_signal_num[l]:
           all_signal_num[l][bar] = 0
           negative_signal_num[l][bar] = 0
        for Si in allChannels:
            qdc_1 = allChannels[Si]
            if qdc_1 == -999:
               continue
            if smallSiPMchannel(Si):
               h_small_sipm.Fill(allChannels[Si])
            if small_count:
               if smallSiPMchannel(Si):
                  all_signal_num[l][bar] += 1 
                  if qdc_1 < 0:
                     negative_signal_num[l][bar] += 1
            else: 
               all_signal_num[l][bar] += 1 
               if qdc_1 < 0:
                     negative_signal_num[l][bar] += 1
            if smallSiPMchannel(Si):
               small_list.append(qdc_1)
               if Si >= 8:
                  l_vs_r_hist_right_small_list.append(qdc_1)
                  l_vs_r_list_right.append(qdc_1)
               else:
                  l_vs_r_hist_left_small_list.append(qdc_1)
                  l_vs_r_list_left.append(qdc_1)
            else:
               large_list.append(qdc_1)
               if Si >= 8:
                  l_vs_r_hist_right_large_list.append(qdc_1)
                  l_vs_r_list_right.append(qdc_1)
               else:
                  l_vs_r_hist_left_large_list.append(qdc_1)
                  l_vs_r_list_left.append(qdc_1)
                  if l == 3:
                     l_us4_lausanne_check[bar][Si] = qdc_1
                     if Si not in l_us4_per_event_lausanne_check[bar]:
                        l_us4_per_event_lausanne_check[bar][Si] = [qdc_1]
                        l_us4_integrated_lausanne_check[bar][Si] = [qdc_1]
                     else:
                        l_us4_per_event_lausanne_check[bar][Si].append(qdc_1)
                        l_us4_integrated_lausanne_check[bar][Si].append(qdc_1)
                     l_us4_lausanne_check_hist[bar][Si].Fill(qdc_1)
            hist_list_sipm[s*10 + l][bar][Si].Fill(qdc_1)
        l_vs_s_hist[s*10 + l][bar].Fill(np.array(large_list).mean(), np.array(small_list).mean())
        l_vs_r_hist[s*10 + l][bar].Fill(np.array(l_vs_r_list_right).mean(), np.array(l_vs_r_list_left).mean())
        l_vs_r_hist_large[s*10 + l][bar].Fill(np.array(l_vs_r_hist_right_large_list).mean(), np.array(l_vs_r_hist_left_large_list).mean())
        l_vs_r_hist_small[s*10 + l][bar].Fill(np.array(l_vs_r_hist_right_small_list).mean(), np.array(l_vs_r_hist_left_small_list).mean())
        l_vs_s_hist_left[s*10 + l][bar].Fill(np.array(l_vs_r_hist_left_large_list).mean(), np.array(l_vs_r_hist_left_small_list).mean())
        l_vs_s_hist_right[s*10 + l][bar].Fill(np.array(l_vs_r_hist_right_large_list).mean(), np.array(l_vs_r_hist_right_small_list).mean())
        l_vs_r_hist_large_all.Fill(np.array(l_vs_r_hist_left_large_list).mean(), np.array(l_vs_r_hist_right_large_list).mean())
      #   qdc_value_large = np.sqrt(sum(l_vs_r_hist_left_large_list)*sum(l_vs_r_hist_right_large_list))
        if write_data:
                  meta_data = [N, hit_number, detID]
                  string_output = ""
                  string_output += "\t".join(list(map(str, meta_data)))
                  string_output += "\t"
                  string_output += "\t".join(list(map(str, list(allChannels.values()))))
                  string_output += "\t"
                  string_output += "\t".join(list(map(str, list(allTimes.values()))))
                  string_output += "\n"
                  fout.write(string_output)



        qdc_value = np.sqrt(sum(l_vs_r_hist_left_small_list)*sum(l_vs_r_hist_right_small_list))


      #   if qdc_value == -999:
      #       continue

        hist_list[s*10 + l].Fill(qdc_value)
        hist_list_bar[s*10 + l][bar].Fill(qdc_value)



         

        if l not in us_shower:
            us_shower[l] = {}
      
        if bar not in us_shower[l]:
            us_shower[l][bar] = 0
            
        

        qdc_per_event += qdc_value
        number_of_hits += 1
        #print(N, qdc_value)
        us_shower[l][bar] += qdc_value
        
        
        
        
    if number_of_hits > 0:
      h_qdc_hit.Fill(number_of_hits, qdc_per_event)
      h_qdc_hit_norm.Fill(number_of_hits, qdc_per_event/(number_of_hits))
      us_qdc_list.append(qdc_per_event)

    event_freq.Fill((N // 10000)*10, number_of_hits)
    for plane in us_shower:
      sum_us = sum([us_shower[plane][i] for i in us_shower[plane]])
      h_long.Fill(plane, sum([us_shower[plane][i] for i in us_shower[plane]])) 
      if "us_" + str(plane + 1) not in us_scifi:
         us_scifi["us_" + str(plane + 1)] = []
      us_scifi["us_" + str(plane + 1)].append(sum_us) 
      for bar in us_shower[plane]:
         hist_usshower.Fill(plane, bar, us_shower[plane][bar])
         if all_signal_num[plane][bar]:
            hist_usshower_fraction.Fill(plane, bar, negative_signal_num[plane][bar]/all_signal_num[plane][bar])


#  for h in hist_scifi_qdc:
#    # if "scifi_" + h not in us_scifi:
#    #    us_scifi["scifi_" + h] = []
#    for bin_idx in range(1, hist_scifi_qdc[h].GetNbinsX() + 1):
#      us_scifi["scifi_" + h].append(hist_scifi_qdc[h].GetBinContent(bin_idx))

#  for 

 
 LIST = [[f"scifiX_{i}", f"scifiY_{i}"] for i in range(1, 5)]
 LIST_inter = []
 for x in LIST:
    for y in x:
       LIST_inter.append(y)
 us_scifi_sorted_keys = LIST_inter + [f"us_{i}" for i in range(1, 6)]
 us_scifi_sorted = {}
 for x in us_scifi_sorted_keys:
    us_scifi_sorted[x] = 0
 for x in us_scifi:
    us_scifi_sorted[x] = us_scifi[x]
 #print("!!!!!!", us_scifi, "!!!!!!!")
    
 fig, ax = plt.subplots(dpi = 200)
 plt.grid()
 y_data = [np.array(us_scifi_sorted[x]).mean() for x in us_scifi_sorted_keys]
 y_data_err = [np.array(us_scifi_sorted[x]).std() for x in us_scifi_sorted_keys]
 x_scifi = [1,2,3,4]
 x_us = [10,11,12,13,14]
 coef_scifi = np.polyfit(x_scifi,[(y_data[i+1] + y_data[i])/2 for i in range(0, len(y_data[:8]), 2)],1)
 coef_us = np.polyfit(x_us,y_data[8:],1)
 y_scifi = np.array([(y_data[i+1] + y_data[i])/2 for i in range(0, len(y_data[:8]), 2)])
 y_scifi_err = np.array([(y_data_err[i+1] + y_data_err[i])/2 for i in range(0, len(y_data_err[:8]), 2)])
 y_us = np.array(y_data[8:])
 y_us_err = np.array(y_data_err[8:])
#  print(x_us, y_data[8:])
#  print(coef_scifi, coef_us)

 poly1d_scifi = np.poly1d(coef_scifi)
 poly1d_us = np.poly1d(coef_us)
 #print(y_data_err)
#  a_fit,cov=curve_fit(linearFunc,x_scifi,y_scifi)
#  print("CHECK!!!!!!")
#  inter_s, slope_s  = a_fit[0], a_fit[1]
#  d_inter_s, d_slope_s = np.sqrt(cov[0][0]), np.sqrt(cov[1][1])
#  a_fit,cov=curve_fit(linearFunc,x_scifi,y_us)
#  inter_us, slope_us  = a_fit[0], a_fit[1]
#  d_inter_us, d_slope_us = np.sqrt(cov[0][0]), np.sqrt(cov[1][1])

 ax.errorbar(x = us_scifi_sorted_keys, y = y_data, yerr = y_data_err,  color = "red", fmt = "o")
#  ax.set_xticklabels(ax.get_xticks(), rotation=45)
 ax.plot(us_scifi_sorted_keys[1:9:2], poly1d_scifi(x_scifi))
 ax.plot(us_scifi_sorted_keys[8:], poly1d_us(x_us))
 props = dict(boxstyle='round', alpha=0.5)
 textstr = "$y_{scifi}$" + " = {:.1f}x + {:.1f}".format(coef_scifi[0], coef_scifi[1]) + \
    "\n" + "$y_{us}$" +  " = {:.1f}x + {:.1f}".format(coef_us[0], coef_us[1])
#  textstr = "Scifi slope: " + "$" + "{:.1f}".format(slope_s) + "\pm" + "{:.1f}".format(d_slope_s) + "$ \n" + "US slope: " + "$" + "{:.1f}".format(slope_us) + "\pm" + "{:.1f}".format(d_slope_us) + "$"
   
 # place a text box in upper left in axes coords
 ax.text(0.62, 0.95, textstr, transform=ax.transAxes, fontsize=12,
       verticalalignment='top', bbox=props)
 ax.set_xticklabels(us_scifi_sorted_keys, rotation=40)
 ax.set_xlabel("Plane")
 ax.set_ylabel("Average signal per event [a.u.]")
 ax.set_title(f"Longitudinal shower profile. {options.runNumber} {label_0} run")
 label_for_save = "_".join(label_0.split(" "))
 plt.savefig(f"us_scifi_{options.runNumber}_{label_for_save}.pdf")
 with open("/eos/user/u/ursovsnd/private/SND_Data/scifi_us_calibration/" + f"{options.runNumber}_side-bar_{side_bar_key}_small-sipm_{small_count}_large-sipm_{large_count}_{conv_mech}.json", "w") as f:
    data_to_save = {}
    for i, k in enumerate(us_scifi_sorted_keys):
       data_to_save[k] = y_data[i]
    json.dump(data_to_save, f)
 with open("/eos/user/u/ursovsnd/private/SND_Data/scifi_us_calibration/" + f"{options.runNumber}error_side-bar_{side_bar_key}_small-sipm_{small_count}_large-sipm_{large_count}_{conv_mech}.json", "w") as f:
    data_to_save = {}
    for i, k in enumerate(us_scifi_sorted_keys):
       data_to_save[k] = y_data_err[i]
    json.dump(data_to_save, f)
    
 
 #print(us_qdc_list)
 fig, ax = plt.subplots(dpi = 200)
 ax.hist(us_qdc_list, bins = 200)
 plt.savefig(f"hist_{options.runNumber}.pdf")
 print(f"{options.runNumber} run", f": {np.array(us_qdc_list).mean()}", np.array(us_qdc_list).std())
 
 
 File = ROOT.TFile.Open(f"{'muon' if muon else 'pion'}_{title}_{options.runNumber}{label}.root", "RECREATE")
 draw_qdc_sipm(hist_list_sipm, "US QDC distribution", wall_num)
 draw_pcb(l_us4_lausanne_check_hist, "PCB US4L")
 draw_qdc_bar(hist_list_bar, "US QDC distribution")
 draw_qdc_bar(l_vs_s_hist, "Large vs Small")
 draw_qdc_bar(l_vs_s_hist_right, "Large vs Small. Right")
 draw_qdc_bar(l_vs_s_hist_left, "Large vs Small. Large")
 draw_qdc_bar(l_vs_r_hist, "Left vs Right")
 draw_qdc_bar(l_vs_r_hist_small, "Left vs Right. Small")
 draw_qdc_bar(l_vs_r_hist_large, "Left vs Right. Large")
 draw_qdc_plane(hist_list, "US QDC distribution")
 draw_qdc_vs_noh(h_qdc_hit, h_qdc_hit_norm)
 draw_usshower(hist_usshower)
 draw_usshower(hist_usshower_fraction)
 draw_track_dist(File, [l_vs_r_hist_large_all])
 draw_track_dist(File, [h_xy_track, h_xy_track_res])
 draw_track_dist(File, [h_long])
 draw_track_dist(File, [h_small_sipm])
 draw_track_dist(File, [event_freq])
 draw_scifi_dist(File, hist_scifi_hit, "hit_distribution")
 draw_scifi_dist(File, hist_scifi_qdc, "qdc_distribution")
 File.Close()
#  File = ROOT.TFile.Open(f"{options.runNumber}_sipm_hists.root", "RECREATE") 
#  draw_scifi_dist(File, hist_scifi_qdc, "qdc_distribution")
 print(f"Data of {options.runNumber} run has been stored in {'muon' if muon else 'pion'}_{title}_{options.runNumber}{label}.root")
 print(f"Negative singals: {negative_signal}")

def run_energy(run):
   if run == 100631 or run == 100677:
      return 100
   if run == 100633 or run == 100673:
      return 140
   if run == 100648 or run == 100646:
      return 240
   if run == 100636:
      return 180
   if run == 100643 or run == 100645:
      return 300
import pickle

def signal_relation_digi_hists(signal_hit_scifi, event_info, energy):
    fig, ax = plt.subplots(figsize = (8,8), dpi = 300)
    ax.set_yscale("log")
   #  ax[1].set_yscale("log")
    ax.hist(signal_hit_scifi["signal_scifi_hit"], bins = 100, histtype = "step")
    ax.set_xlabel("Scifi hit signal [a.u.]")
   #  ax[1].set_xlabel("US signal sum/event [a.u.]")
   #  ax[1].hist(signal_hit_us["signal_us_hit"], bins = 100, histtype = "step")
    fig.savefig(f"signal_hists_signal_{energy}.pdf")      
    fig, ax = plt.subplots(figsize = (8,8), dpi = 300)
   #  ax[1].set_yscale("log")
    ax.hist(event_info["scifi_hit"], bins = 100, histtype = "step")
    ax.set_xlabel("Hits/event")

   #  ax[1].set_xlabel("US signal sum/event [a.u.]")
   #  ax[1].hist(signal_hit_us["signal_us_hit"], bins = 100, histtype = "step")
    fig.savefig(f"signal_hists_hit_{energy}.pdf")    



def def_tagging(Digi_ScifiHits, TY, Nsf_hits, SFqdcStation, Nsf_itofpet, signal_hit_scifi, N, file_info):
      for aHit in Digi_ScifiHits:        
         if (aHit.GetTime(0) - TY) < 0. or (aHit.GetTime(0) - TY) > 1.: continue
         Nsf_hits[aHit.GetStation()] += 1
         SFqdcStation[aHit.GetStation()] += aHit.GetSignal()
         X1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
         Y1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
         if aHit.isVertical():
               #if aHit.GetTofpetID(0) == 3:
                  #print("vertical: ", aHit.GetTofpetID(0), "station: ", aHit.GetStation(), "X,Y: ", X1, Y1)
               Nsf_itofpet[0][aHit.GetTofpetID(0)][aHit.GetStation()] += 1

         else:
               #if aHit.GetTofpetID(0) == 3:
                  #print("horizontal: ", aHit.GetTofpetID(0), "station: ", aHit.GetStation(), "X,Y: ", X1, Y1)
               Nsf_itofpet[1][aHit.GetTofpetID(0)][aHit.GetStation()] += 1
         signal_hit_scifi["signal_scifi_hit"].append(aHit.GetSignal())
         signal_hit_scifi["event_number"].append(N)
 
      Iwall = 0
      for t in range(4):
         if Iwall != 0: break
         Ishower=0
         NhitPerTofpetShowerId = 0
         for q in range(8):
            if Nsf_itofpet[0][q][t+1] > 35:
               NhitPerTofpetShowerId = Nsf_itofpet[0][q][t+1]
               Ishower = 1

         if Ishower == 1:
               for q in range(8):
                  if Nsf_itofpet[1][q][t+1] > 35:
                        Ishower=2
                        Iwall = t
      signal_sum = 0
      hit_sum = 0       
      for station in SFqdcStation:
         signal_sum += SFqdcStation[station]
         hit_sum += Nsf_hits[station]
         if f"scifi_{station}_hits" not in file_info:
            file_info[f"scifi_{station}_hits"] = [Nsf_hits[station]]    
            file_info[f"scifi_{station}_qdc"] = [SFqdcStation[station]]
         else:
            file_info[f"scifi_{station}_hits"].append(Nsf_hits[station])  
            file_info[f"scifi_{station}_qdc"].append(SFqdcStation[station])   
            
      file_info["energy"].append(run_energy(options.runNumber))
      file_info["event"].append(N)
      file_info["wall"].append(Iwall)
      file_info["scifi_signal"].append(signal_sum)
      file_info["scifi_hit"].append(hit_sum) 
        
def second_tagging(Digi_ScifiHits, TY, Nsf_hits, SFqdcStation, Nsf_itofpet, signal_hit_scifi, N, file_info, F = 0.10):

   for aHit in Digi_ScifiHits: 
      if (aHit.GetTime(0) - TY) < 0. or (aHit.GetTime(0) - TY) > 1.: continue
      Nsf_hits[aHit.GetStation()] += 1
      SFqdcStation[aHit.GetStation()] += aHit.GetSignal()
   Iwall = 0
   for st in [1,2,3]:
      if Iwall != 0: break
      #if Nsf_hits[st]/(Nsf_hits[st] + Nsf_hits[st+1]) < F:
      if Nsf_hits[st+1] > (1/F - 1)*Nsf_hits[st] and Nsf_hits[st+1] > 15:
            Iwall = st

   signal_sum = 0
   hit_sum = 0       
   for station in SFqdcStation:
      signal_sum += SFqdcStation[station]
      hit_sum += Nsf_hits[station]
      if f"scifi_{station}_hits" not in file_info:
         file_info[f"scifi_{station}_hits"] = [Nsf_hits[station]]    
         file_info[f"scifi_{station}_qdc"] = [SFqdcStation[station]]
      else:
         file_info[f"scifi_{station}_hits"].append(Nsf_hits[station])  
         file_info[f"scifi_{station}_qdc"].append(SFqdcStation[station])   
         
   file_info["energy"].append(run_energy(options.runNumber))
   file_info["event"].append(N)
   file_info["wall"].append(Iwall)
   file_info["scifi_signal"].append(signal_sum)
   file_info["scifi_hit"].append(hit_sum) 
   

def gen_events_scifi(Nev_st, Nev_en, save = True, condor = False):
 print(Nev_st, Nev_en)
 file_info = {"energy": [], "event": [], "wall": [], "scifi_signal": [], "scifi_hit": [], "us_signal": [], "us_hit": [], "sat_sipms": [], "notsat_sipms": [], "ds_vert_signal": [], "ds_vert_hit": [], "ds_hor_signal": [], "ds_hor_hit": []}
 signal_hit_scifi = {"event_number": [], "signal_scifi_hit": []}
 with open(f'/eos/user/u/ursovsnd/SWAN_projects/tests/scifi_events_bar_level_{options.runNumber}_first_alg_35_1-1.pkl', 'wb') as f:  # open a text file
   good = 0
   for N in range(Nev_st, Nev_en + 1):
   #  for event in eventTree:
      #  N+=1
      eventTree.GetEvent(N)
      #  if N <= Nev_st - 1: continue 
      if N%10000 == 0: print('event ',N,' ',time.ctime())   
      
      N1Y = 0
      N1X = 0
      Nsf_hits = {i: 0 for i in range(1,5)}
      SFqdcStation = {i: 0 for i in range(1,5)}
      Nsf_itofpet = {0: [{t:0 for t in range(1,5)} for _ in range(8)], 1: [{t:0 for t in range(1,5)} for _ in range(8)]}
      Nus_hits = {i: 0 for i in range(5)}
      Nus_qdc = {i: 0 for i in range(5)}
      Nds_hits = {i: 0 for i in range(2)}
      Nds_qdc = {i: 0 for i in range(2)}
      Nusbars_hits = {i: {j: 0 for j in range(10)} for i in range(5)}
      Nusbars_qdc = {i: {j: 0 for j in range(10)} for i in range(5)}
      # signal_sum = 0
      # hit_number_1 = 0
      for aHit in eventTree.Digi_ScifiHits:
         # if aHit.GetSignal():
         #    print(aHit.GetSignal())
         # signal_sum += aHit.GetSignal()
         # hit_number_1 += 1
         if aHit.GetStation()==1:
            if aHit.isVertical():
               N1Y += 1
               X1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
               TY = aHit.GetTime(0)
            else:
               N1X += 1
               Y1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
               TX = aHit.GetTime(0)
      # print("CHECK!!!", hit_number_1, signal_sum)             
      if N1X == 1 and N1Y == 1:
         good += 1
      else:
         continue
      # print("CHECK_1!!!")            
      if np.abs(TX - TY) > 1:
         continue
      # print("CHECK_2!!!")            
      def_tagging(eventTree.Digi_ScifiHits, TY, Nsf_hits, SFqdcStation, Nsf_itofpet, signal_hit_scifi, N, file_info)
      #second_tagging(eventTree.Digi_ScifiHits, TY, Nsf_hits, SFqdcStation, Nsf_itofpet, signal_hit_scifi, N, file_info)

      # if (Iwall == 0): continue
      ######### 
      num_of_sipms_sat = 0  
      num_of_sipms_not_sat = 0        
      sat_thr = 130
      for aHit in eventTree.Digi_MuFilterHits:
         detID = aHit.GetDetectorID()
         s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
         l  = (detID%10000)//1000  # plane number
         bar = (detID%1000)
         if aHit.GetSystem() == 3:
            if aHit.isVertical():
               if aHit.isVertical() > -100:
                  Nds_qdc[0] += aHit.GetSignal()
                  Nds_hits[0] += 1
            else:
               if aHit.isVertical() > -100:
                  for ch in range(2):
                     Nds_qdc[1] += aHit.GetSignal(ch)
                  Nds_hits[1] += 1                              
            # print("!!!!!! DS:", aHit.GetSystem(), s, aHit.isVertical())
            # if aHit.GetSignal() < 0:
            #    print("NEGATIVE: ", aHit.GetSignal())
            # print([[sig.first, sig.second] for sig in aHit.GetAllSignals()])
            # print("indir: ", [aHit.GetSignal(ks) for ks in range(len(aHit.GetAllSignals()))])
         else:
            hit_check = 0
            for ch in range(16):
               if aHit.GetSignal(ch)>-100 and not aHit.isMasked(ch) and not aHit.isShort(ch) and (aHit.GetTime(ch)-TY) < 3.:
                  if aHit.GetSignal(ch) > sat_thr:
                     num_of_sipms_sat += 1 
                  else:
                     num_of_sipms_not_sat += 1 
                  hit_check = 1
                  Nus_qdc[l] += aHit.GetSignal(ch)
                  Nusbars_qdc[l][bar] += aHit.GetSignal(ch)
            if hit_check:
               Nus_hits[l] += 1
               Nusbars_hits[l][bar] += 1       

      signal_us_sum = 0
      hit_us_sum = 0       
      for l in Nus_qdc:
         signal_us_sum += Nus_qdc[l]
         hit_us_sum += Nus_hits[l] 
         if f"us_{l}_hits" not in file_info:
            file_info[f"us_{l}_hits"] = [Nus_hits[l]]
            file_info[f"us_{l}_qdc"] = [Nus_qdc[l]]
         else:
            file_info[f"us_{l}_hits"].append(Nus_hits[l])
            file_info[f"us_{l}_qdc"].append(Nus_qdc[l])
         for bar in Nusbars_qdc[l]:
            if f"us_{l}_{bar}_hits" not in file_info:
               file_info[f"us_{l}_{bar}_hits"] = [Nusbars_hits[l][bar]]
               file_info[f"us_{l}_{bar}_qdc"] = [Nusbars_qdc[l][bar]]
            else:
               file_info[f"us_{l}_{bar}_hits"].append(Nusbars_hits[l][bar])
               file_info[f"us_{l}_{bar}_qdc"].append(Nusbars_qdc[l][bar])                 
   
      #########
        
      
      # Iwall_signal[Iwall].append(signal_sum)
      # print(N)
      file_info["us_signal"].append(signal_us_sum)
      file_info["us_hit"].append(hit_us_sum)
      file_info["ds_vert_signal"].append(Nds_qdc[0])
      file_info["ds_vert_hit"].append(Nds_hits[0])
      file_info["ds_hor_signal"].append(Nds_qdc[1])
      file_info["ds_hor_hit"].append(Nds_hits[1])
      file_info["sat_sipms"].append(num_of_sipms_sat)
      file_info["notsat_sipms"].append(num_of_sipms_not_sat)
   if save:
      pickle.dump(file_info, f) # serialize the list
   if condor:
      pd.DataFrame(file_info).to_csv("output.csv")
   print(good, "number of events")
   signal_relation_digi_hists(signal_hit_scifi, file_info, 300)
   signal_hit_scifi_out = np.array(signal_hit_scifi["signal_scifi_hit"])
   
   

# Open the ROOT file
def read_good_events(run_number):
    with uproot.open("/afs/cern.ch/work/f/fmei/public/TB_showertags.root") as file:

        # Access the 'ShowerTags' tree
        tree = file["ShowerTags"]

        # Get the branches from the tree
        branch_names = tree.keys()

        # Create an empty dictionary to store branch data
        data_dict = {}

        # Access each branch and retrieve data
        for branch_name in branch_names:
            branch_data = tree[branch_name].array()
            data_dict[branch_name] = branch_data

        # Create a pandas DataFrame from the dictionary
        df = pd.DataFrame(data_dict).query(f"run_number == {run_number}")
        # print(df)
        return df
     
import statistics as st


def draw_cors(hists):
   c  = ROOT.TCanvas(f"CORS",f"CORS",0,0,1600,1600)
   c.Divide(2,2)
   for i, plane in enumerate(hists):
      tc = c.cd(i+1)
      hists[plane].Draw("COLZ")
   c.Write()



def gen_events_scifi_filippo(Nev_st, Nev_en, save = True, condor = False): 
 hist_list_3walls = {f"Wall {i}": ROOT.TH2I("Scifi_US","Scifi_US;US;SciFi", 200,0,20000.,200,0,8000) for i in [1,2,3]}
 hist_time_scifinocut = ROOT.TH1I("SciFiTimeNocut", "SciFi Time;Time [clk cycles]", 100, 0, 5)
 hist_time_usnocut = ROOT.TH2I("USQDCTimeNocut", "US Signal-Time;Time - tref [clk cycles]; Signal [QDC]", 100, 0, 15, 100, 0, 180)
 hist_time_scifi2dnocut = ROOT.TH2I("SciFiQDCTimeNocut", "SciFi Signal-Time;Time [clk cycles]; Signal [QDC]", 100, 0, 15, 100, -2, 10)
 hist_time_scifi = ROOT.TH1I("SciFiTime", "SciFi Time;Time [clk cycles]", 100, 0, 5)
 hist_time_us = ROOT.TH2I("USQDCTime", "US Signal-Time;Time - tref [clk cycles]; Signal [QDC]", 100, 0, 15, 100, 0, 180)
 hist_time_scifi2d = ROOT.TH2I("SciFiQDCTime", "SciFi Signal-Time;Time[clk cycles]; Signal [QDC]", 100, 0, 15, 100, -2, 10)
 hist_list_3walls["All"] = ROOT.TH2I("Scifi_US","Scifi_US;US;SciFi", 200,0,20000.,200,0,8000)
 file_info = {"energy": [], "event": [], "wall": [], "scifi_signal": [], "scifi_hit": [], "us_signal": [], "us_hit": [], "sat_sipms": [], "notsat_sipms": [], "ds_vert_signal": [], "ds_vert_hit": [], "ds_hor_signal": [], "ds_hor_hit": []}
 signal_hit_scifi = {"event_number": [], "signal_scifi_hit": []}
 
 


 with open(f'/eos/user/u/ursovsnd/SWAN_projects/tests/scifi_events_bar_level_{options.runNumber}_filippo_ready_events.pkl', 'wb') as f:  # open a text file 
   df_events = read_good_events(options.runNumber)
   sss = 0
   if condor:
      df_loop = df_events.query(f"event_number > {Nev_st} & event_number < {Nev_en}")["event_number"]
   else:
      df_loop = df_events
   #.query(f"event_number > {Nev_st} & event_number < {Nev_en}")
   ke = 0 
   print(df_loop)
   for N in df_loop["event_number"]:
      # print(df_loop.query(f"event_number == {N}")["wall"].to_list())
      wall_num = df_loop.query(f"event_number == {N}")["wall"].to_list()[0]
      if wall_num == 0: continue
      # print(wall_num.to_list()[0])
   #  for event in eventTree:
      #  N+=1
      if not condor:
         if sss > Nev_en:
            break
         sss += 1    
         if sss%5000 == 0: print('event ',sss,' ',time.ctime())       
      eventTree.GetEvent(N)
      #  if N <= Nev_st - 1: continue  
      
      N1Y = 0
      N1X = 0
      Nsf_hits = {i: 0 for i in range(1,5)}
      SFqdcStation = {i: 0 for i in range(1,5)}
      Nsf_itofpet = {0: [{t:0 for t in range(1,5)} for _ in range(8)], 1: [{t:0 for t in range(1,5)} for _ in range(8)]}
      Nus_hits = {i: 0 for i in range(5)}
      Nus_qdc = {i: 0 for i in range(5)}
      Nds_hits = {i: 0 for i in range(2)}
      Nds_qdc = {i: 0 for i in range(2)}
      Nusbars_hits = {i: {j: 0 for j in range(10)} for i in range(5)}
      Nusbars_qdc = {i: {j: 0 for j in range(10)} for i in range(5)}
      good = 0
      time_dist = []
      for aHit in eventTree.Digi_ScifiHits: 
         time_dist.append(aHit.GetTime(0))
         if not ke and sss == 10:
            hist_time_scifinocut.Fill(aHit.GetTime(0))
         if aHit.GetStation()==1:
            if aHit.isVertical():
               N1Y += 1
               X1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
               TY = aHit.GetTime(0)
               #print(X1)
            else:
               N1X += 1
               Y1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
               TX = aHit.GetTime(0)
         else:
            pass
      if sss == 10:
         ke = 1     
      # if N1X and N1Y:
      #    good += 1
      # else:
      #    continue
      # # print("CHECK_1!!!")            
      # if np.abs(TX - TY) > 1:
      #    continue
      print(time_dist)
      MOD = st.mode(np.array(time_dist))
      signal_sum = 0
      hit_sum = 0
      for aHit in eventTree.Digi_ScifiHits:        
         # if (aHit.GetTime(0) - TY) < 0. or (aHit.GetTime(0) - TY) > 1.: continue
         # if aHit.GetStation()==1: continue
         hist_time_scifi2dnocut.Fill(aHit.GetSignal(), aHit.GetTime(0))
         if np.abs(aHit.GetTime(0) - MOD) > 0.5: continue
         hist_time_scifi2d.Fill(aHit.GetSignal(), aHit.GetTime(0))
         hist_time_scifi.Fill(aHit.GetTime(0))
         Nsf_hits[aHit.GetStation()] += 1
         SFqdcStation[aHit.GetStation()] += aHit.GetSignal()
         signal_sum += aHit.GetSignal()
         hit_sum += 1
         X1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
         Y1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
   

      for station in SFqdcStation:
         if f"scifi_{station}_hits" not in file_info:
            file_info[f"scifi_{station}_hits"] = [Nsf_hits[station]]    
            file_info[f"scifi_{station}_qdc"] = [SFqdcStation[station]]
         else:
            file_info[f"scifi_{station}_hits"].append(Nsf_hits[station])  
            file_info[f"scifi_{station}_qdc"].append(SFqdcStation[station])   
            
      file_info["energy"].append(run_energy(options.runNumber))
      file_info["event"].append(N)
      file_info["wall"].append(wall_num)
      file_info["scifi_signal"].append(signal_sum)
      file_info["scifi_hit"].append(hit_sum) 
      
      # if (Iwall == 0): continue
      ######### 
      num_of_sipms_sat = 0  
      num_of_sipms_not_sat = 0        
      sat_thr = 130
      signal_us_sum = 0
      hit_us_sum = 0     
      for aHit in eventTree.Digi_MuFilterHits:
         detID = aHit.GetDetectorID()
         s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
         l  = (detID%10000)//1000  # plane number
         bar = (detID%1000)
         if aHit.GetSystem() == 3:
            if aHit.isVertical():
               if aHit.isVertical() > -100:
                  Nds_qdc[0] += aHit.GetSignal()
                  Nds_hits[0] += 1
            else:
               if aHit.isVertical() > -100:
                  for ch in range(2):
                     Nds_qdc[1] += aHit.GetSignal(ch)
                  Nds_hits[1] += 1                              
            # print("!!!!!! DS:", aHit.GetSystem(), s, aHit.isVertical())
            # if aHit.GetSignal() < 0:
            #    print("NEGATIVE: ", aHit.GetSignal())
            # print([[sig.first, sig.second] for sig in aHit.GetAllSignals()])
            # print("indir: ", [aHit.GetSignal(ks) for ks in range(len(aHit.GetAllSignals()))])
         else:
            hit_check = 0
            for ch in range(16):
               if aHit.GetSignal(ch) < -100 or smallSiPMchannel(ch) or aHit.isMasked(ch): continue
               hist_time_usnocut.Fill(aHit.GetTime(ch) - MOD, aHit.GetSignal(ch))
               if np.abs(aHit.GetTime(ch) - MOD) > 3: continue         
               hist_time_us.Fill(aHit.GetTime(ch) - MOD, aHit.GetSignal(ch))
               if aHit.GetSignal(ch) > sat_thr:
                  num_of_sipms_sat += 1 
               else:
                  num_of_sipms_not_sat += 1 
               hit_check = 1
               Nus_qdc[l] += aHit.GetSignal(ch)
               Nusbars_qdc[l][bar] += aHit.GetSignal(ch)
               signal_us_sum += aHit.GetSignal(ch)
                  
            if hit_check:
               Nus_hits[l] += 1
               Nusbars_hits[l][bar] += 1       
               hit_us_sum += 1
      for l in Nus_qdc:
         if f"us_{l}_hits" not in file_info:
            file_info[f"us_{l}_hits"] = [Nus_hits[l]]
            file_info[f"us_{l}_qdc"] = [Nus_qdc[l]]
         else:
            file_info[f"us_{l}_hits"].append(Nus_hits[l])
            file_info[f"us_{l}_qdc"].append(Nus_qdc[l])
         for bar in Nusbars_qdc[l]:
            if f"us_{l}_{bar}_hits" not in file_info:
               file_info[f"us_{l}_{bar}_hits"] = [Nusbars_hits[l][bar]]
               file_info[f"us_{l}_{bar}_qdc"] = [Nusbars_qdc[l][bar]]
            else:
               file_info[f"us_{l}_{bar}_hits"].append(Nusbars_hits[l][bar])
               file_info[f"us_{l}_{bar}_qdc"].append(Nusbars_qdc[l][bar])                 
   
      #########
      hist_list_3walls[f"Wall {wall_num}"].Fill(signal_us_sum, signal_sum) 
      hist_list_3walls[f"All"].Fill(signal_us_sum, signal_sum)       
      # Iwall_signal[Iwall].append(signal_sum)
      # print(N)
      file_info["us_signal"].append(signal_us_sum)
      file_info["us_hit"].append(hit_us_sum)
      file_info["ds_vert_signal"].append(Nds_qdc[0])
      file_info["ds_vert_hit"].append(Nds_hits[0])
      file_info["ds_hor_signal"].append(Nds_qdc[1])
      file_info["ds_hor_hit"].append(Nds_hits[1])
      file_info["sat_sipms"].append(num_of_sipms_sat)
      file_info["notsat_sipms"].append(num_of_sipms_not_sat)
   
   
   File = ROOT.TFile.Open(f"COR_{options.runNumber}_new_filippo_w_small.root", "RECREATE")
   draw_cors(hist_list_3walls)
   
   
   def write_dict_to_file(dict_obj, File):
      for key in dict_obj.keys():
         File.WriteObject(dict_obj[key], dict_obj[key].GetTitle())
         
      
   hist_dict = {"us_time": hist_time_us, "scifi_time": hist_time_scifi, "us_timenocut": hist_time_usnocut, "scifi_timenocut": hist_time_scifinocut, "scifi_time2d": hist_time_scifi2d, "scifi_time2dnocut": hist_time_scifi2dnocut}
   for key in hist_dict.keys():
      File.WriteObject(hist_dict[key], hist_dict[key].GetTitle())
   
   if save:
      pickle.dump(file_info, f) # serialize the list
   if condor:
      pd.DataFrame(file_info).to_csv("output.csv")
   # signal_relation_digi_hists(signal_hit_scifi, file_info, 300)
   #print(signal_hit_scifi_out[signal_hit_scifi_out<0].size/signal_hit_scifi_out.size)
   print(pd.DataFrame(file_info))






def gen_events_scifi_mc(Nev_st, Nev_en, energy = 100, save = True, condor = False): 
 hist_list_3walls = {f"Wall {i}": ROOT.TH2I("Scifi_US","Scifi_US;US;SciFi", 200,0,20000.,200,0,50000) for i in [1,2,3]}
 hist_time_scifinocut = ROOT.TH1I("SciFiTimeNocut", "SciFi Time;Time [clk cycles]", 100, 0, 10)
 hist_signal_scifinocut = ROOT.TH1I("SciFiSignalpereventNocut", "SciFi Signal per event;SciFi Signal [QDC]", 200, 0, 50000)
 hist_time_usnocut = ROOT.TH2I("USQDCTimeNocut", "US Signal-Time;Time - tref [clk cycles]; Signal [QDC]", 100, 0, 15, 100, 0, 180)
 hist_time_scifi2dnocut = ROOT.TH2I("SciFiQDCTimeNocut", "SciFi Signal-Time;Time [clk cycles]; Signal [QDC]", 100, 0, 15, 100, -2, 10)
 hist_signal_scifi = ROOT.TH1I("SciFiSignalperevent", "SciFi Signal per event;SciFi Signal [QDC]", 200, 0, 50000) 
 hist_time_scifi = ROOT.TH1I("SciFiTime", "SciFi Time;Time [clk cycles]", 100, 0, 10)
 hist_time_us = ROOT.TH2I("USQDCTime", "US Signal-Time;Time - tref [clk cycles]; Signal [QDC]", 100, 0, 15, 100, 0, 180)
 hist_time_scifi2d = ROOT.TH2I("SciFiQDCTime", "SciFi Signal-Time;Time[clk cycles]; Signal [QDC]", 100, 0, 15, 100, -2, 10)
 hist_list_3walls["All"] = ROOT.TH2I("Scifi_US","Scifi_US;US;SciFi", 200,0,20000.,200,0,50000)
 file_info = {"energy": [], "event": [], "wall": [], "scifi_signal": [], "scifi_hit": [], "us_signal": [], "us_hit": [], "sat_sipms": [], "notsat_sipms": [], "ds_vert_signal": [], "ds_vert_hit": [], "ds_hor_signal": [], "ds_hor_hit": []}
 signal_hit_scifi = {"event_number": [], "signal_scifi_hit": []}
 


 with open(f'/eos/user/u/ursovsnd/SWAN_projects/tests/scifi_events_bar_level_{options.runNumber}_filippo_ready_events.pkl', 'wb') as f:  # open a text file 
   df_events = pd.read_csv(f"/eos/user/u/ursovsnd/SWAN_projects/tests/shower_tagging/merged_{300}.csv")
   sss = 0
   #.query(f"event_number > {Nev_st} & event_number < {Nev_en}")
   ke = 0
   for N in df_events["event_number"]:
      wall_num = df_events.loc[df_events["event_number"] == N]["wall"].to_list()[0]
      if wall_num == 0: continue
      # print(wall_num.to_list()[0])
   #  for event in eventTree:
      #  N+=1
      if not condor:
         if sss > Nev_en or N > 10000:
            break
         sss += 1    
         if sss%5000 == 0: print('event ',sss,' ',time.ctime())       
      eventTree.GetEvent(N)
      #  if N <= Nev_st - 1: continue  
      
      N1Y = 0
      N1X = 0
      Nsf_hits = {i: 0 for i in range(1,5)}
      SFqdcStation = {i: 0 for i in range(1,5)}
      Nsf_itofpet = {0: [{t:0 for t in range(1,5)} for _ in range(8)], 1: [{t:0 for t in range(1,5)} for _ in range(8)]}
      Nus_hits = {i: 0 for i in range(5)}
      Nus_qdc = {i: 0 for i in range(5)}
      Nds_hits = {i: 0 for i in range(2)}
      Nds_qdc = {i: 0 for i in range(2)}
      Nusbars_hits = {i: {j: 0 for j in range(10)} for i in range(5)}
      Nusbars_qdc = {i: {j: 0 for j in range(10)} for i in range(5)}
      good = 0
      time_dist = []
      for aHit in eventTree.Digi_ScifiHits: 
            time_dist.append(aHit.GetTime(0))
            if not ke and sss == 10:
               hist_time_scifinocut.Fill(aHit.GetTime(0))
            if aHit.GetStation()==1:
               if aHit.isVertical():
                  N1Y += 1
                  X1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
                  TY = aHit.GetTime(0)
                  #print(X1)
               else:
                  N1X += 1
                  Y1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
                  TX = aHit.GetTime(0)
            else:
               pass
      if sss == 10:
         ke = 1     
      # if N1X and N1Y:
      #    good += 1
      # else:
      #    continue
      # # print("CHECK_1!!!")            
      # if np.abs(TX - TY) > 1:
      #    continue
      #print(time_dist)
      MOD = st.mode(np.array(time_dist))
      signal_sum = 0
      hit_sum = 0
      for aHit in eventTree.Digi_ScifiHits:        
         # if (aHit.GetTime(0) - TY) < 0. or (aHit.GetTime(0) - TY) > 1.: continue
         # if aHit.GetStation()==1: continue
         hist_time_scifi2dnocut.Fill(aHit.GetSignal(), aHit.GetTime(0))
         #if np.abs(aHit.GetTime(0) - MOD) > 0.5: continue
         hist_time_scifi2d.Fill(aHit.GetSignal(), aHit.GetTime(0))
         hist_time_scifi.Fill(aHit.GetTime(0))
         Nsf_hits[aHit.GetStation()] += 1
         SFqdcStation[aHit.GetStation()] += aHit.GetSignal()
         signal_sum += aHit.GetSignal()
         hit_sum += 1
         X1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
         Y1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
   

      for station in SFqdcStation:
         if f"scifi_{station}_hits" not in file_info:
            file_info[f"scifi_{station}_hits"] = [Nsf_hits[station]]    
            file_info[f"scifi_{station}_qdc"] = [SFqdcStation[station]]
         else:
            file_info[f"scifi_{station}_hits"].append(Nsf_hits[station])  
            file_info[f"scifi_{station}_qdc"].append(SFqdcStation[station])   
            
      file_info["energy"].append(run_energy(options.runNumber))
      file_info["event"].append(N)
      file_info["wall"].append(wall_num)
      file_info["scifi_signal"].append(signal_sum)
      file_info["scifi_hit"].append(hit_sum) 
      hist_signal_scifi.Fill(signal_sum)
      # if (Iwall == 0): continue
      ######### 
      num_of_sipms_sat = 0  
      num_of_sipms_not_sat = 0        
      sat_thr = 130
      signal_us_sum = 0
      hit_us_sum = 0     
      for aHit in eventTree.Digi_MuFilterHits:
         detID = aHit.GetDetectorID()
         s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
         l  = (detID%10000)//1000  # plane number
         bar = (detID%1000)
         if aHit.GetSystem() == 3:
            if aHit.isVertical():
               if aHit.isVertical() > -100:
                  Nds_qdc[0] += aHit.GetSignal()
                  Nds_hits[0] += 1
            else:
               if aHit.isVertical() > -100:
                  for ch in range(2):
                     Nds_qdc[1] += aHit.GetSignal(ch)
                  Nds_hits[1] += 1                              
            # print("!!!!!! DS:", aHit.GetSystem(), s, aHit.isVertical())
            # if aHit.GetSignal() < 0:
            #    print("NEGATIVE: ", aHit.GetSignal())
            # print([[sig.first, sig.second] for sig in aHit.GetAllSignals()])
            # print("indir: ", [aHit.GetSignal(ks) for ks in range(len(aHit.GetAllSignals()))])
         else:
            hit_check = 0
            for ch in range(16):
               if aHit.GetSignal(ch) < -100 or smallSiPMchannel(ch) or aHit.isMasked(ch): continue
               hist_time_usnocut.Fill(aHit.GetTime(ch) - MOD, aHit.GetSignal(ch))
               #if np.abs(aHit.GetTime(ch) - MOD) > 3: continue         
               hist_time_us.Fill(aHit.GetTime(ch) - MOD, aHit.GetSignal(ch))
               if aHit.GetSignal(ch) > sat_thr:
                  num_of_sipms_sat += 1 
               else:
                  num_of_sipms_not_sat += 1 
               hit_check = 1
               Nus_qdc[l] += aHit.GetSignal(ch)
               Nusbars_qdc[l][bar] += aHit.GetSignal(ch)
               signal_us_sum += aHit.GetSignal(ch)
                  
            if hit_check:
               Nus_hits[l] += 1
               Nusbars_hits[l][bar] += 1       
               hit_us_sum += 1
      for l in Nus_qdc:
         if f"us_{l}_hits" not in file_info:
            file_info[f"us_{l}_hits"] = [Nus_hits[l]]
            file_info[f"us_{l}_qdc"] = [Nus_qdc[l]]
         else:
            file_info[f"us_{l}_hits"].append(Nus_hits[l])
            file_info[f"us_{l}_qdc"].append(Nus_qdc[l])
         for bar in Nusbars_qdc[l]:
            if f"us_{l}_{bar}_hits" not in file_info:
               file_info[f"us_{l}_{bar}_hits"] = [Nusbars_hits[l][bar]]
               file_info[f"us_{l}_{bar}_qdc"] = [Nusbars_qdc[l][bar]]
            else:
               file_info[f"us_{l}_{bar}_hits"].append(Nusbars_hits[l][bar])
               file_info[f"us_{l}_{bar}_qdc"].append(Nusbars_qdc[l][bar])                 
   
      #########
      hist_list_3walls[f"Wall {wall_num}"].Fill(signal_us_sum, signal_sum) 
      hist_list_3walls[f"All"].Fill(signal_us_sum, signal_sum)       
      # Iwall_signal[Iwall].append(signal_sum)
      # print(N)
      file_info["us_signal"].append(signal_us_sum)
      file_info["us_hit"].append(hit_us_sum)
      file_info["ds_vert_signal"].append(Nds_qdc[0])
      file_info["ds_vert_hit"].append(Nds_hits[0])
      file_info["ds_hor_signal"].append(Nds_qdc[1])
      file_info["ds_hor_hit"].append(Nds_hits[1])
      file_info["sat_sipms"].append(num_of_sipms_sat)
      file_info["notsat_sipms"].append(num_of_sipms_not_sat)
   
   
   File = ROOT.TFile.Open(f"COR_{options.runNumber}_new_filippo_w_small.root", "RECREATE")
   draw_cors(hist_list_3walls)
   
   
   def write_dict_to_file(dict_obj, File):
      for key in dict_obj.keys():
         File.WriteObject(dict_obj[key], dict_obj[key].GetTitle())
         
      
   hist_dict = {"us_time": hist_time_us, "scifi_time": hist_time_scifi, 
                "us_timenocut": hist_time_usnocut, "scifi_timenocut": hist_time_scifinocut, "scifi_time2d": hist_time_scifi2d, "scifi_time2dnocut": hist_time_scifi2dnocut,
                "scifi-signal-per-event": hist_signal_scifi}
   for key in hist_dict.keys():
      File.WriteObject(hist_dict[key], hist_dict[key].GetTitle())
   
   if save:
      pickle.dump(file_info, f) # serialize the list
   if condor:
      pd.DataFrame(file_info).to_csv("output.csv")
   # signal_relation_digi_hists(signal_hit_scifi, file_info, 300)
   #print(signal_hit_scifi_out[signal_hit_scifi_out<0].size/signal_hit_scifi_out.size)
   print(pd.DataFrame(file_info))





def gen_events_scifi_all(Nev_st, Nev_en, energy = 100, is_mc = False, save = True, condor = False, tagged_events = True): 
 hist_list_3walls = {f"Wall {i}": ROOT.TH2I("Scifi_US","Scifi_US;US;SciFi", 200,0,20000.,200,0,50000) for i in [1,2,3]}
 hist_time_scifinocut = ROOT.TH1I("SciFiTimeNocut", "SciFi Time;Time [clk cycles]", 100, 0, 10)
 hist_signal_scifinocut = ROOT.TH1I("SciFiSignalpereventNocut", "SciFi Signal per event;SciFi Signal [QDC]", 200, 0, 50000)
 hist_time_usnocut = ROOT.TH2I("USQDCTimeNocut", "US Signal-Time;Time - tref [clk cycles]; Signal [QDC]", 100, 0, 15, 100, 0, 180)
 hist_time_scifi2dnocut = ROOT.TH2I("SciFiQDCTimeNocut", "SciFi Signal-Time;Time [clk cycles]; Signal [QDC]", 100, 0, 15, 100, -2, 10)
 hist_signal_scifi = ROOT.TH1I("SciFiSignalperevent", "SciFi Signal per event;SciFi Signal [QDC]", 200, 0, 50000) 
 hist_time_scifi = ROOT.TH1I("SciFiTime", "SciFi Time;Time [clk cycles]", 100, 0, 10)
 hist_time_us = ROOT.TH2I("USQDCTime", "US Signal-Time;Time - tref [clk cycles]; Signal [QDC]", 100, 0, 15, 100, 0, 180)
 hist_time_scifi2d = ROOT.TH2I("SciFiQDCTime", "SciFi Signal-Time;Time[clk cycles]; Signal [QDC]", 100, 0, 15, 100, -2, 10)
 hist_list_3walls["All"] = ROOT.TH2I("Scifi_US","Scifi_US;US;SciFi", 200,0,20000.,200,0,50000)
 file_info = {"energy": [], "event": [], "wall": [], "scifi_signal": [], "scifi_hit": [], "us_signal": [], "us_hit": [], "sat_sipms": [], "notsat_sipms": [], "ds_vert_signal": [], "ds_vert_hit": [], "ds_hor_signal": [], "ds_hor_hit": []}
 signal_hit_scifi = {"event_number": [], "signal_scifi_hit": []}
 
 if options.runNumber != 0:
    is_mc = False
 else:
    is_mc = True

 with open(f'/eos/user/u/ursovsnd/SWAN_projects/tests/scifi_events_bar_level_{options.runNumber}_filippo_ready_events.pkl', 'wb') as f:  # open a text file 
   
   if tagged_events:
      if is_mc:
         df_events = pd.read_csv(f"/eos/user/u/ursovsnd/SWAN_projects/tests/shower_tagging/merged_{300}.csv")
      else:
         df_events = read_good_events(options.runNumber)
   else:
      df_events = pd.DataFrame({"event_number": np.arange(Nev_en)})
      # df_events = pd.DataFrame({"event_number": np.arangerange(100000)})
   sss = 0
   #.query(f"event_number > {Nev_st} & event_number < {Nev_en}")
   
   
   if condor:
      df_loop = df_events.query(f"event_number > {Nev_st} & event_number < {Nev_en}")["event_number"]
   else:
      df_loop = df_events
   #.query(f"event_number > {Nev_st} & event_number < {Nev_en}")
   ke = 0 
   print(df_loop)
   for N in df_loop["event_number"]:
      # print(df_loop.query(f"event_number == {N}")["wall"].to_list())
      if tagged_events:
         wall_num = df_loop.query(f"event_number == {N}")["wall"].to_list()[0]
      else:
         wall_num = 1
      #wall_num = 1
      if wall_num == 0: continue
      # print(wall_num.to_list()[0])
   #  for event in eventTree:
      #  N+=1
      if not condor:
         if sss > Nev_en:
            break
         if is_mc and N > 10000:
            break
         sss += 1    
         if sss%5000 == 0: print('event ',sss,' ',time.ctime())
      eventTree.GetEvent(N)
      #  if N <= Nev_st - 1: continue  
      
      N1Y = 0
      N1X = 0
      Nsf_hits = {i: 0 for i in range(1,5)}
      SFqdcStation = {i: 0 for i in range(1,5)}
      Nsf_itofpet = {0: [{t:0 for t in range(1,5)} for _ in range(8)], 1: [{t:0 for t in range(1,5)} for _ in range(8)]}
      Nus_hits = {i: 0 for i in range(5)}
      Nus_qdc = {i: 0 for i in range(5)}
      Nds_hits = {i: 0 for i in range(2)}
      Nds_qdc = {i: 0 for i in range(2)}
      Nusbars_hits = {i: {j: 0 for j in range(10)} for i in range(5)}
      Nusbars_qdc = {i: {j: 0 for j in range(10)} for i in range(5)}
      good = 0
      time_dist = []
      for aHit in eventTree.Digi_ScifiHits: 
            time_dist.append(aHit.GetTime(0))
            if not ke and sss == 10:
               hist_time_scifinocut.Fill(aHit.GetTime(0))
            if aHit.GetStation()==1:
               if aHit.isVertical():
                  N1Y += 1
                  X1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
                  TY = aHit.GetTime(0)
                  #print(X1)
               else:
                  N1X += 1
                  Y1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
                  TX = aHit.GetTime(0)
            else:
               pass
      if sss == 10:
         ke = 1     
      # if N1X and N1Y:
      #    good += 1
      # else:
      #    continue
      # # print("CHECK_1!!!")            
      # if np.abs(TX - TY) > 1:
      #    continue
      #print(time_dist)
      MOD = st.mode(np.array(time_dist))
      signal_sum = 0
      hit_sum = 0
      for aHit in eventTree.Digi_ScifiHits:        
         # if (aHit.GetTime(0) - TY) < 0. or (aHit.GetTime(0) - TY) > 1.: continue
         # if aHit.GetStation()==1: continue
         hist_time_scifi2dnocut.Fill(aHit.GetTime(0), aHit.GetSignal())
         if np.abs(aHit.GetTime(0) - MOD) > 0.5 and not is_mc: continue
         hist_time_scifi2d.Fill(aHit.GetTime(0), aHit.GetSignal())
         hist_time_scifi.Fill(aHit.GetTime(0))
         Nsf_hits[aHit.GetStation()] += 1
         SFqdcStation[aHit.GetStation()] += aHit.GetSignal()
         signal_sum += aHit.GetSignal()
         hit_sum += 1
         X1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
         Y1 = (64*aHit.GetTofpetID(0)+63-aHit.Getchannel(0))*0.025
   

      for station in SFqdcStation:
         if f"scifi_{station}_hits" not in file_info:
            file_info[f"scifi_{station}_hits"] = [Nsf_hits[station]]    
            file_info[f"scifi_{station}_qdc"] = [SFqdcStation[station]]
         else:
            file_info[f"scifi_{station}_hits"].append(Nsf_hits[station])  
            file_info[f"scifi_{station}_qdc"].append(SFqdcStation[station])   
         
      file_info["energy"].append(energy)
      file_info["event"].append(N)
      file_info["wall"].append(wall_num)
      file_info["scifi_signal"].append(signal_sum)
      file_info["scifi_hit"].append(hit_sum) 
      hist_signal_scifi.Fill(signal_sum)
      # if (Iwall == 0): continue
      ######### 
      num_of_sipms_sat = 0  
      num_of_sipms_not_sat = 0        
      sat_thr = 130
      signal_us_sum = 0
      hit_us_sum = 0     
      for aHit in eventTree.Digi_MuFilterHits:
         detID = aHit.GetDetectorID()
         s = detID//10000 # det number (s = 1 -- veto, s = 2 -- US, s = 3 -- DS)
         l  = (detID%10000)//1000  # plane number
         bar = (detID%1000)
         if aHit.GetSystem() == 3:
            if aHit.isVertical():
               if aHit.isVertical() > -100:
                  Nds_qdc[0] += aHit.GetSignal()
                  Nds_hits[0] += 1
            else:
               if aHit.isVertical() > -100:
                  for ch in range(2):
                     Nds_qdc[1] += aHit.GetSignal(ch)
                  Nds_hits[1] += 1                              
            # print("!!!!!! DS:", aHit.GetSystem(), s, aHit.isVertical())
            # if aHit.GetSignal() < 0:
            #    print("NEGATIVE: ", aHit.GetSignal())
            # print([[sig.first, sig.second] for sig in aHit.GetAllSignals()])
            # print("indir: ", [aHit.GetSignal(ks) for ks in range(len(aHit.GetAllSignals()))])
         else:
            hit_check = 0
            for ch in range(16):
               if aHit.GetSignal(ch) < -100 or smallSiPMchannel(ch) or aHit.isMasked(ch): continue
               hist_time_usnocut.Fill(aHit.GetTime(ch) - MOD, aHit.GetSignal(ch))
               if np.abs(aHit.GetTime(ch) - MOD) > 3 and not is_mc: continue         
               hist_time_us.Fill(aHit.GetTime(ch) - MOD, aHit.GetSignal(ch))
               if aHit.GetSignal(ch) > sat_thr:
                  num_of_sipms_sat += 1 
               else:
                  num_of_sipms_not_sat += 1 
               hit_check = 1
               Nus_qdc[l] += aHit.GetSignal(ch)
               Nusbars_qdc[l][bar] += aHit.GetSignal(ch)
               signal_us_sum += aHit.GetSignal(ch)
                  
            if hit_check:
               Nus_hits[l] += 1
               Nusbars_hits[l][bar] += 1       
               hit_us_sum += 1
      for l in Nus_qdc:
         if f"us_{l}_hits" not in file_info:
            file_info[f"us_{l}_hits"] = [Nus_hits[l]]
            file_info[f"us_{l}_qdc"] = [Nus_qdc[l]]
         else:
            file_info[f"us_{l}_hits"].append(Nus_hits[l])
            file_info[f"us_{l}_qdc"].append(Nus_qdc[l])
         for bar in Nusbars_qdc[l]:
            if f"us_{l}_{bar}_hits" not in file_info:
               file_info[f"us_{l}_{bar}_hits"] = [Nusbars_hits[l][bar]]
               file_info[f"us_{l}_{bar}_qdc"] = [Nusbars_qdc[l][bar]]
            else:
               file_info[f"us_{l}_{bar}_hits"].append(Nusbars_hits[l][bar])
               file_info[f"us_{l}_{bar}_qdc"].append(Nusbars_qdc[l][bar])                 
   
      #########
      hist_list_3walls[f"Wall {wall_num}"].Fill(signal_us_sum, signal_sum) 
      hist_list_3walls[f"All"].Fill(signal_us_sum, signal_sum)       
      # Iwall_signal[Iwall].append(signal_sum)
      # print(N)
      file_info["us_signal"].append(signal_us_sum)
      file_info["us_hit"].append(hit_us_sum)
      file_info["ds_vert_signal"].append(Nds_qdc[0])
      file_info["ds_vert_hit"].append(Nds_hits[0])
      file_info["ds_hor_signal"].append(Nds_qdc[1])
      file_info["ds_hor_hit"].append(Nds_hits[1])
      file_info["sat_sipms"].append(num_of_sipms_sat)
      file_info["notsat_sipms"].append(num_of_sipms_not_sat)
   
   
   File = ROOT.TFile.Open(f"COR_{options.runNumber}_new_filippo_w_small_check_1.root", "RECREATE")
   #draw_cors(hist_list_3walls)
   for key in hist_list_3walls.keys():
      File.WriteObject(hist_list_3walls[key], hist_list_3walls[key].GetTitle())   
   
   def write_dict_to_file(dict_obj, File):
      for key in dict_obj.keys():
         File.WriteObject(dict_obj[key], dict_obj[key].GetTitle())
         
      
   hist_dict = {"us_time": hist_time_us, "scifi_time": hist_time_scifi, 
                "us_timenocut": hist_time_usnocut, "scifi_timenocut": hist_time_scifinocut, "scifi_time2d": hist_time_scifi2d, "scifi_time2dnocut": hist_time_scifi2dnocut,
                "scifi-signal-per-event": hist_signal_scifi}
   for key in hist_dict.keys():
      File.WriteObject(hist_dict[key], hist_dict[key].GetTitle())
   
   if save:
      pickle.dump(file_info, f) # serialize the list
   if condor:
      pd.DataFrame(file_info).to_csv("output.csv")
   # signal_relation_digi_hists(signal_hit_scifi, file_info, 300)
   #print(signal_hit_scifi_out[signal_hit_scifi_out<0].size/signal_hit_scifi_out.size)
   print(pd.DataFrame(file_info))

    

# title = "muon_august_2023"
# label = "muon 100 GeV"
title = "august"
label = "pion 300 GeV"
# title = "pion_august_2023_after_all_sipm_cut_onehit_us_offset-7-QDC_off"
# title = "muon_H8"
title = "pion_august_2023_after_all_sipm_cut_onehit_scifi_antitag"
# MIP_study(Nev_st = options.Nstart,
#    Nev_en = options.Nstart + options.nEvents, 
#    list_of_events_key = False,
#    title = title,
#    label_0 = label,
#    conv_mech = "dir",
#    offset = 0.,
#    side_bar_key = False,
#    mc = True,
#    muon = False,
#    small_count = True,
#    large_count = True,
#    fill_sipm_key = False,
#    write_data = False)

#gen_events_scifi(options.Nstart, options.Nstart + options.nEvents, save = True, condor = False)


#gen_events_scifi_filippo(options.Nstart, options.Nstart + options.nEvents, save = False, condor = False)
#file_info = gen_events_scifi_mc(options.Nstart, options.Nstart + options.nEvents, energy = 300, save = False, condor = False)
file_info = gen_events_scifi_all(options.Nstart, options.Nstart + options.nEvents, energy = 300, is_mc = True, save = False, condor = False, tagged_events = True)

# Mufi_hitMaps(100000)
# smallVsLargeSiPMs(1000000)
# Scifi_hitMaps(100)
# Scifi_slopes(100000)
