// -*- C++ -*-
//
// Package:    Test/QCDEventShape
// Class:      QCDEventShape
// 
/**\class QCDEventShape QCDEventShape.cc Test/QCDEventShape/plugins/QCDEventShape.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Soumyadip Barman
//         Created:  Tue, 02 Feb 2021 10:24:21 GMT
//
//


// system include files

#define DIJETAVE 

////for data
//#define JETENERGY
//#define TRIGGER

////for Madgraph
//#define LHAPDF
//#define JETRESO
//#define TRACKSYS
//#define TRIGGER

////for Pythia8 & Herwig7
#define JETRESO
#define TRIGGER

////For Flat
//#define FLAT

//#define LUMIWEIGHT

//#define PREFIRE


//#define MERGE_BIN
//#define PUNOMI
//#define PUUP
//#define PUDOWN

////For GenParticle only
//#define GENPART


#include <memory>
#include <map>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TFormula.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include <cmath>
#include "TMath.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "TRandom.h"
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"
#include "TUnfoldBinningXML.h"
#include "TUnfold.h"
#include "TUnfoldSys.h"

#include "TH2F.h"
#include "TProfile.h"
#include <fstream>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <time.h>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtTrigReport.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtTrigReportEntry.h"
#include "CondFormats/DataRecord/interface/L1GtStableParametersRcd.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "FWCore/Framework/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"

#include "Test/QCDEventShape/plugins/EventShape_vector.h" 


#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <CondFormats/DataRecord/interface/JetResolutionRcd.h>
#include <CondFormats/DataRecord/interface/JetResolutionScaleFactorRcd.h>
#include "FWCore/Utilities/interface/typelookup.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "TProfile.h"
using namespace edm;
using namespace reco;
using namespace std;
using namespace CLHEP;
using namespace trigger;
using namespace math;

static const int nvar=32;
//static const int njcvar=3; // 3 definitions of jet charge observables
static const int nhist=10;
static const int typen=2; 
static const int nHLTmx=9;

const int nkappa=10;
double kappa[nkappa]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

double recojet1_pt; 
double recojet1g_pt, recojet1qg_pt, recojet1ag_pt;
double recojet1b_pt, recojet1qb_pt, recojet1ab_pt;
double recojet1c_pt, recojet1qc_pt, recojet1ac_pt;
double recojet1s_pt, recojet1qs_pt, recojet1as_pt;
double recojet1u_pt, recojet1qu_pt, recojet1au_pt;
double recojet1d_pt, recojet1qd_pt, recojet1ad_pt; 

double recojet2_pt; 
double recojet2g_pt, recojet2qg_pt, recojet2ag_pt; 
double recojet2b_pt, recojet2qb_pt, recojet2ab_pt;
double recojet2c_pt, recojet2qc_pt, recojet2ac_pt;
double recojet2s_pt, recojet2qs_pt, recojet2as_pt;
double recojet2u_pt, recojet2qu_pt, recojet2au_pt;
double recojet2d_pt, recojet2qd_pt, recojet2ad_pt;

double genrecojet1_pt, genrecojet1g_pt, genrecojet1b_pt, genrecojet1c_pt, genrecojet1s_pt, genrecojet1u_pt, genrecojet1d_pt; 
double genrecojet2_pt, genrecojet2g_pt, genrecojet2b_pt, genrecojet2c_pt, genrecojet2s_pt, genrecojet2u_pt, genrecojet2d_pt; 
//const char* jcvarname[njcvar]={"Default", "Longitudinal", "Transverse"}; // jet charge

//-------------------------------------------
const char* varname[nvar]={"y3anti", "y3ceanti", "y3cranti", "thrustc", "thrustce", "thrustcr",
                           "minorc", "minorce", "minorcr", "tmass", "tmasse", "tmassr",
                           "hmass", "hmasse", "hmassr", "y3c", "y3ce", "y3cr",
                           "broadt", "broadte", "broadtr", "broadw", "broadwe", "broadwr",
                           "ttmass", "ttmasse", "ttmassr", "htmass", "htmasse", "htmassr",
                           "sphericity", "cparameter"};

const char* vartitle[nvar]={"Anti-Y_{23,C} ", "Anti-Y_{23,E} ", "Anti-Y_{23,R} ",
                            "#tau_{_{#perp} _{   ,C}} ", "#tau_{_{#perp} _{   ,E}} ", "#tau_{_{#perp} _{   ,R}} ",
                            "T_{ m,C} ", "T_{ m,E} ", "T_{ m,R} ",
                            "#rho_{Tot,C} ", "#rho_{Tot,E} ", "#rho_{Tot,R} ",
                            "#rho_{H,C} ", "#rho_{H,E} ", "#rho_{H,R} ",
                            "Y_{23,C} ", "Y_{23,E} ", "Y_{23,R} ",
                            "B_{ T,C} ", "B_{ T,E} ", "B_{ T,R} ",
                            "B_{ W,C} ", "B_{ W,E} ", "B_{ W,R} ",
                            "#rho^{T}_{Tot,C} ", "#rho^{T}_{Tot,E} ", "#rho^{T}_{Tot,R} ",
                            "#rho^{T}_{H,C} ", "#rho^{T}_{H,E} ", "#rho^{T}_{H,R} ",
                            "S_{_{#perp} _{   ,C}}", "C-parameter_{C}"};

//--------------------------For fixed Binning
int nbinsx[nvar]={120, 120, 120, 120, 120, 120,  //Number of Bins 
                   120, 120, 120, 120, 120, 120, 
                   120, 120, 120, 120, 120, 120, 
                   120, 120, 120, 120, 120, 120, 
                   120, 120, 120, 120, 120, 120, 
                   120, 120};

double endxt[nvar]={8.0, 8.0, 6.0, 8.0, 6.0,  5.0, //Lower Edges
                    4.0, 4.0,  2.0, 7.0, 7.0, 4.0,
                    7.0, 7.0, 4.0, 8.0, 8.0, 6.0,
                    5.0, 5.0, 4.0, 5.0, 5.0, 4.0,
                    10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
                    0.0, 0.0};
double startx[nvar]={2., 0.5, 1., 1.0, 0.0, 2.,   //Lower Edges
                     0., 0., 1., 0.2, -1., -1.,
                     0., 0., 1., 1.0, 2., 0.5,
                     0.2, -0.5, -1., 0., -1., -1.,
                     1.0, 0., -0.5, 0., 0., -0.5,
                    -1., -1.};
double endx[nvar]={10.0, 10.0, 10.0, 8.0, 12.0, 12.0,
                   19.0, 19.0, 10.0, 7.0, 10.0, 10.0,
                   13.0, 13.0,  10.0, 8.0, 10.0, 6.0,
                   6.0, 5.0, 5.0, 8.0, 8.0, 8.0,
                   7.0, 8.0, 8.0, 12.0, 12.0, 8.0,
                   0.0, 0.0};
/*
//-------------------------No. of Bins For Reco Level for Fixed Bin
//For Jet 
const int rnmxbins=26;
int rnbinsx0[nvar]={0,0,0,26,0,0,
                   0,0,0,20,0,0,
                   0,0,0,12,0,0,
                   20,0,0,0,0,0,
                   22,0,0,0,0,0,0,0};

//For charge Particles
int rnbinsx1[nvar]={0,0,0,20,0,0,
                   0,0,0,16,0,0,
                   0,0,0,8,0,0,
                   22,0,0,0,0,0,
                   14,0,0,0,0,0,0,0};
*/
//-----------------------------------For Reco Level
/*
int jcbins[60]={40,40,40,30,20,20,20,20,20,20,
		40,40,40,30,20,20,20,20,20,20,
		20,20,20,20,20,20,20,20,20,20,
		20,20,20,20,20,20,20,20,20,20,
		20,20,20,20,20,20,20,20,20,20,
		20,20,20,20,20,20,20,20,20,20};

double jcminran[60]={-8.0,-5.0,-4.0,-3.0,-2.0,-2.0,-2.0,-2.0,-1.0,-1.0,
		     -8.0,-5.0,-4.0,-3.0,-2.0,-2.0,-2.0,-2.0,-1.0,-1.0,
		     -1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,
		     -1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,
		     -1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,
		     -1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,};

double jcmaxran[60]={8.0,5.0,4.0,3.0,2.0,2.0,2.0,2.0,1.0,1.0,
		     8.0,5.0,4.0,3.0,2.0,2.0,2.0,2.0,1.0,1.0,
		     1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
		     1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
		     1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
	             1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,};
*/
//----------------------------------------------

/*
int recojcd1bins[nkappa]={40,40,40,40,20,20,20,20,20,20};
double recojcd1minran[nkappa]={-8.0,-5.0,-4.0,-3.0,-2.0,-2.0,-2.0,-2.0,-1.0,-1.0};
double recojcd1maxran[nkappa]={8.0,5.0,4.0,3.0,2.0,2.0,2.0,2.0,1.0,1.0};

int recojcd23bins[nkappa]={20,20,20,20,20,20,20,20,20,20};
double recojcd23minran[nkappa]={-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
double recojcd23maxran[nkappa]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};

int genjcd1bins[nkappa]={20,20,20,20,10,10,10,10,10,10};
double genjcd1minran[nkappa]={-8.0,-5.0,-4.0,-3.0,-2.0,-2.0,-2.0,-2.0,-1.0,-1.0};
double genjcd1maxran[nkappa]={8.0,5.0,4.0,3.0,2.0,2.0,2.0,2.0,1.0,1.0};

int genjcd23bins[nkappa]={10,10,10,10,10,10,10,10,10,10};
double genjcd23minran[nkappa]={-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
double genjcd23maxran[nkappa]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
*/

//Reco level
int recojcd1bins[nkappa]={32,32,32,20,20,20,20,20,20,20};
double recojcd1minran[nkappa]={-6.0,-4.0,-3.0,-2.0,-2.0,-1.0,-1.0,-1.0,-1.0,-1.0};
double recojcd1maxran[nkappa]={6.0,4.0,3.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0};

int recojcd23bins[nkappa]={20,20,20,20,20,20,20,20,20,20};
double recojcd23minran[nkappa]={-0.4,-0.6,-0.6,-0.8,-0.8,-1.0,-1.0,-1.0,-1.0,-1.0};
double recojcd23maxran[nkappa]={0.4,0.6,0.6,0.8,0.8,1.0,1.0,1.0,1.0,1.0};

//Gen level
int genjcd1bins[nkappa]={16,16,16,10,10,10,10,10,10,10};
double genjcd1minran[nkappa]={-6.0,-4.0,-3.0,-2.0,-2.0,-1.0,-1.0,-1.0,-1.0,-1.0};
double genjcd1maxran[nkappa]={6.0,4.0,3.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0};

int genjcd23bins[nkappa]={10,10,10,10,10,10,10,10,10,10};
double genjcd23minran[nkappa]={-0.4,-0.6,-0.6,-0.8,-0.8,-1.0,-1.0,-1.0,-1.0,-1.0};
double genjcd23maxran[nkappa]={0.4,0.6,0.6,0.8,0.8,1.0,1.0,1.0,1.0,1.0};
//-------------------------------------------
const int rnmxbins=32;  //Maximum Bins in bellow array
int rnbinsx0[nvar]={0,0,0,32,0,0,0,0,0,
                  32,0,0,0,0,0,14,0,
                  0,24,0,0,0,0,0,26,
                  0,0,0,0,0,0,0};

//For charge Particles
int rnbinsx1[nvar]= {0,0,0,18,0,0,0,0,
                  0,18,0,0,0,0,0,8,
                  0,0,16,0,0,0,0,0,
                  14,0,0,0,0,0,0,0};

//-----------------------------------For Gen Level
const int nmxbins=16;  //Maximum Bins in bellow array
int nbinsx0[nvar]={0,0,0,16,0,0,0,0,
                  0,16,0,0,0,0,0,7,
                  0,0,12,0,0,0,0,0,
                  13,0,0,0,0,0,0,0};

 int  nbinsx1[nvar]={0,0,0,9,0,0,0,0,
                  0,9,0,0,0,0,0,4,
                  0,0,8,0,0,0,0,0,
                  7,0,0,0,0,0,0,0};

//-------------------------------------------
//Reco level  /Jet
double rbinrngs0[nvar][rnmxbins+1] ={{},{},{},
                                      {-6.71, -6.4, -6.11, -5.83, -5.56, -5.3, -5.05, -4.8, -4.56, -4.33, -4.11, -3.89, -3.68, -3.48, -3.28, -3.09, -2.91, -2.74, -2.58, -2.43, -2.29, -2.15, -2.02, -1.9, -1.79, -1.69, -1.6, -1.51, -1.43, -1.36, -1.3, -1.24, -1.19},//32
                                      {},{},{},{},{},
                                      {-5.75, -5.45, -5.18, -4.92, -4.68, -4.45, -4.22, -4, -3.79, -3.58, -3.37, -3.17, -2.97, -2.77, -2.57, -2.38, -2.19, -2.01, -1.84, -1.67, -1.51, -1.36, -1.22, -1.09, -0.97, -0.86, -0.76, -0.67, -0.59, -0.52, -0.45, -0.39, -0.34}, //32
                                      {},{},{},{},{},
                                      {-6.71, -6.18, -5.67, -5.18, -4.73, -4.31, -3.92, -3.56, -3.22, -2.9, -2.59, -2.3, -2.01, -1.73, -1.44}, //14
                                      {},{},
                                      {-3.58, -3.34, -3.1, -2.87, -2.65, -2.44, -2.24, -2.05, -1.88, -1.72, -1.57, -1.43, -1.3, -1.18, -1.06, -0.95, -0.85, -0.76, -0.67, -0.59, -0.51, -0.43, -0.36, -0.29, -0.23}, //24
                                      {},{},{},{},{},
                                      {-5.89, -5.5, -5.17, -4.88, -4.61, -4.36, -4.12, -3.89, -3.67, -3.45, -3.24, -3.03, -2.83, -2.64, -2.45, -2.27, -2.1, -1.94, -1.79, -1.66, -1.54, -1.43, -1.34, -1.26, -1.19, -1.13, -1.08},//26
                                     {},{},{},{},{},{},{}};

double rbinrngs1[nvar][rnmxbins+1] = {{},{},{},
                                        {-5.66, -5.36, -5.06, -4.76, -4.46, -4.16, -3.86, -3.56, -3.27, -2.98, -2.7, -2.44, -2.19, -1.96, -1.75, -1.56, -1.4, -1.26, -1.14}, //18
                                        {},{},{},{},{},
                                        {-5.75, -5.3, -4.87, -4.45, -4.04, -3.65, -3.27, -2.9, -2.54, -2.2, -1.87, -1.56, -1.28, -1.02, -0.79, -0.59, -0.41, -0.26, -0.13 }, //18
                                        {},{},{},{},{},
                                        {-6.71, -6.04, -5.36, -4.69, -4.04, -3.41, -2.79, -2.15, -1.58 }, //8
                                        {},{},
                                        {-3.77, -3.49, -3.21, -2.93, -2.66, -2.39, -2.13, -1.87, -1.62, -1.39, -1.17, -0.97, -0.79, -0.62, -0.47, -0.34, -0.23 }, //16
                                        {},{},{},{},{},
                                        {-5.89, -5.4, -4.95, -4.52, -4.11, -3.72, -3.35, -2.99, -2.65, -2.33, -2.03, -1.75, -1.49, -1.25, -1.04},
                                        {},{},{},{},{},{},{}};//14

//Gen Level
double binrngs0[nvar][nmxbins+1] = {{},{},{},
                                      {-6.71,-6.11,-5.56,-5.05,-4.56,-4.11,-3.68,-3.28,-2.91,-2.58,-2.29,-2.02,-1.79,-1.6,-1.43,-1.3,-1.19}, //16
                                      {},{},{},{},{},
                                      {-5.75,-5.18,-4.68,-4.22,-3.79,-3.37,-2.97,-2.57,-2.19,-1.84,-1.51,-1.22,-0.97,-0.76,-0.59,-0.45,-0.34}, //16
                                      {},{},{},{},{},
                                      {-6.71,-5.67,-4.73,-3.92,-3.22,-2.59,-2.01,-1.44},
                                      {},{},
                                      {-3.58,-3.1,-2.65,-2.24,-1.88,-1.57,-1.3,-1.06,-0.85,-0.67,-0.51,-0.36,-0.23},//12
                                      {},{},{},{},{},
                                      {-5.89,-5.17,-4.61,-4.12,-3.67,-3.24,-2.83,-2.45,-2.1,-1.79,-1.54,-1.34,-1.19,-1.08},
                                      {},{},{},{},{},{},{}};

double binrngs1[nvar][nmxbins+1]={{},{},{},
                                    {-5.66,-5.06,-4.46,-3.86,-3.27,-2.7,-2.19,-1.75,-1.4,-1.14},
                                    {},{},{},{},{},
                                    {-5.75,-4.87,-4.04,-3.27,-2.54,-1.87,-1.28,-0.79,-0.41,-0.13},
                                    {},{},{},{},{},
                                    {-6.71,-5.36,-4.04,-2.79,-1.58},
                                    {},{},
                                    {-3.77,-3.21,-2.66,-2.13,-1.62,-1.17,-0.79,-0.47,-0.23},
                                    {},{},{},{},{},
                                    {-5.89,-4.95,-4.11,-3.35,-2.65,-2.03,-1.49,-1.04},
                                    {},{},{},{},{},{},{}};

//----------------------------------------HT2 Binning For 2D unfold 
//double recohtbins[nHLTmx+1] = {83, 109, 172, 241, 309, 377, 462, 570, 3000.0};
//double recohtbins[nHLTmx+1] = {89, 118, 186, 259, 333, 403, 488, 528, 584, 688, 3000.0};
double recohtbins[nHLTmx+1] = {76, 101, 164, 232, 295, 373, 449, 507, 612, 3000.0};

//int recohtnbins[nkappa]={10,10,10,10,10,10,10,10,10,10};
//double recohtbinsmin[nkappa] = {92, 119, 185, 251, 319, 388, 467, 518, 579, 669};
//double recohtbinsmax[nkappa] = {119, 185, 251, 319, 388, 467, 518, 579, 669, 3000};

/*
#ifdef MERGE_BIN
double recohtbins[nHLTmx+1] = {83, 109, 176, 247, 318, 387, 477, 573, 3000.0}; //For 2018
#else
double recohtbins[nHLTmx+1] = {66, 95, 155, 227, 286, 350, 444, 557, 3000.0}; // For 2016 Trigger
#endif
*/

//---------------------------------------------Lumi weight
/*
int iera = 2;// 0 for Run 2016, 1 for Run 2017 , 2 for Run 2018
double lumi[3] = {36330, 41480, 59830};
double total_lumi = lumi[0]+lumi[1]+lumi[2];
double lumiwtt = lumi[iera]/total_lumi;
*/
//---------------------------------------------
const int nusedvar=5;
double usedvars[nusedvar]={3, 9, 15, 18, 24};

int isItUsed(int ival) {
	for (int ij=0; ij<nusedvar; ij++) {
		if (ival==usedvars[ij]) {return 1;}
	}
	return 0;
} 
//-------------------------------------------
const int npileupmx=99; //49;
double rat_pileup[nHLTmx][npileupmx]={{0}};
//clock_t t1,t2;
//UL PU
double mcpileup[npileupmx] ={1.00402e-05, 5.76499e-05, 7.37891e-05, 0.000110933, 0.000158858, 0.000368637, 0.000893114, 0.00189701, 0.0035888, 0.00636053, 0.0104174, 0.0158123, 0.0223786, 0.0299187, 0.0380276, 0.0454314, 0.0511181, 0.0547435, 0.0567906, 0.0577145, 0.0578177, 0.0571252, 0.0555457, 0.0531344, 0.0501519, 0.0466816, 0.0429245, 0.0389567, 0.0348507, 0.0307357, 0.0267712, 0.022972, 0.0193389, 0.0159603, 0.0129311, 0.0102889, 0.00798783, 0.00606652, 0.00447821, 0.0032159, 0.00224504, 0.00151447, 0.000981184, 0.00060967, 0.000362193, 0.000211573, 0.000119152, 6.49134e-05, 3.57796e-05, 1.99044e-05, 1.13639e-05, 6.49624e-06, 3.96626e-06, 2.3791e-06, 1.50997e-06, 1.09817e-06, 7.31299e-07, 6.10399e-07, 3.74846e-07, 2.65177e-07, 2.01924e-07, 1.39348e-07, 8.326e-08, 6.04932e-08, 6.52537e-08, 5.90575e-08, 2.29162e-08, 1.97295e-08, 1.77311e-08, 3.57548e-09, 1.3504e-09, 8.50071e-09, 5.02792e-09, 4.93737e-10, 8.1392e-10, 5.62779e-09, 5.15141e-10, 8.21677e-10, 0, 1.49167e-09, 8.43518e-09, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

//HLT Path PileUP
//20APril Ultra Legacy true for HLT_DiPFJet with 8 trigger paths
double datpileup[nHLTmx][npileupmx] ={{5.2228e-06, 4.78698e-05, 0.00031893, 0.00015203, 0.000324416, 0.00151584, 0.00516684, 0.00772539, 0.00854466, 0.0114081, 0.0168755, 0.0239617, 0.0320735, 0.0408648, 0.0499775, 0.0580075, 0.0630094, 0.0644992, 0.0632582, 0.0603569, 0.0568499, 0.0534058, 0.0500768, 0.0465926, 0.042772, 0.0386465, 0.0343659, 0.0301007, 0.0260019, 0.0221828, 0.0187114, 0.0156125, 0.0128799, 0.0104923, 0.00842414, 0.00665175, 0.00515313, 0.00390718, 0.0028921, 0.0020846, 0.0014596, 0.000990586, 0.000650469, 0.000412764, 0.000252968, 0.000149756, 8.57179e-05, 4.75242e-05, 2.55932e-05, 1.34397e-05, 6.91719e-06, 3.51216e-06, 1.77333e-06, 8.98744e-07, 4.61913e-07, 2.43173e-07, 1.32199e-07, 7.45422e-08, 4.35795e-08, 2.62928e-08, 1.6251e-08, 1.02079e-08, 6.46951e-09, 4.11319e-09, 2.61233e-09, 1.65257e-09, 1.03933e-09, 6.49073e-10, 4.02229e-10, 2.4724e-10, 1.50709e-10, 9.10956e-11, 5.45988e-11, 3.24493e-11, 1.91242e-11, 1.11773e-11, 6.47884e-12, 3.72464e-12, 2.12385e-12, 1.20126e-12, 6.73976e-13, 3.75116e-13, 2.07118e-13, 1.13453e-13, 6.16558e-14, 3.32436e-14, 1.7784e-14, 9.43963e-15, 4.9716e-15, 2.59815e-15, 1.34733e-15, 6.93321e-16, 3.5404e-16, 1.79404e-16, 9.02143e-17, 4.50171e-17, 2.22908e-17, 1.09521e-17, 5.33933e-18},
{4.89492e-06, 4.76567e-05, 0.000263877, 0.000188908, 0.000418299, 0.00171571, 0.00593946, 0.0112904, 0.0168787, 0.0240734, 0.0314472, 0.0370343, 0.0411121, 0.0454532, 0.0513282, 0.0576159, 0.0618456, 0.0629171, 0.0613305, 0.0581588, 0.0543907, 0.05055, 0.0466628, 0.0425628, 0.0381995, 0.0337011, 0.0292685, 0.0250796, 0.0212552, 0.0178546, 0.0148856, 0.0123194, 0.0101101, 0.0082103, 0.00658029, 0.00519021, 0.00401753, 0.00304362, 0.00225084, 0.0016208, 0.00113373, 0.000768708, 0.000504346, 0.0003198, 0.000195861, 0.000115874, 6.62793e-05, 3.67184e-05, 1.97554e-05, 1.03619e-05, 5.32531e-06, 2.69906e-06, 1.35991e-06, 6.87579e-07, 3.52505e-07, 1.85139e-07, 1.00455e-07, 5.65722e-08, 3.3059e-08, 1.99515e-08, 1.23426e-08, 7.76261e-09, 4.92681e-09, 3.13702e-09, 1.99524e-09, 1.26392e-09, 7.95901e-10, 4.97621e-10, 3.08696e-10, 1.89927e-10, 1.15872e-10, 7.00926e-11, 4.20403e-11, 2.50018e-11, 1.47439e-11, 8.62222e-12, 5.00056e-12, 2.87634e-12, 1.64102e-12, 9.28672e-13, 5.21328e-13, 2.90323e-13, 1.60395e-13, 8.7914e-14, 4.78078e-14, 2.57945e-14, 1.38088e-14, 7.33499e-15, 3.86607e-15, 2.02199e-15, 1.04938e-15, 5.40438e-16, 2.76197e-16, 1.40073e-16, 7.04936e-17, 3.52045e-17, 1.74457e-17, 8.57807e-18, 4.185e-18},
{1.33191e-05, 8.41562e-05, 0.000336949, 0.000247009, 0.000854888, 0.00689466, 0.0245049, 0.0281708, 0.0163041, 0.0147222, 0.0200221, 0.0265182, 0.0330448, 0.0396851, 0.0467685, 0.0533749, 0.0577262, 0.0592521, 0.0585386, 0.0564231, 0.0537013, 0.0508362, 0.0478055, 0.0443781, 0.0404634, 0.0361914, 0.0317973, 0.0275112, 0.0235081, 0.0198888, 0.0166845, 0.0138762, 0.0114237, 0.00928747, 0.00743786, 0.00585427, 0.00451961, 0.00341603, 0.00252284, 0.00181673, 0.00127275, 0.000865498, 0.000570098, 0.000363126, 0.00022341, 0.00013271, 7.61446e-05, 4.22573e-05, 2.27382e-05, 1.19075e-05, 6.10035e-06, 3.07864e-06, 1.54393e-06, 7.77484e-07, 3.97689e-07, 2.08932e-07, 1.13713e-07, 6.43742e-08, 3.78546e-08, 2.29865e-08, 1.42947e-08, 9.02658e-09, 5.74521e-09, 3.66472e-09, 2.33321e-09, 1.47862e-09, 9.31075e-10, 5.8194e-10, 3.60809e-10, 2.21846e-10, 1.35255e-10, 8.17682e-11, 4.90207e-11, 2.91468e-11, 1.71902e-11, 1.00581e-11, 5.83928e-12, 3.36418e-12, 1.92369e-12, 1.09188e-12, 6.15244e-13, 3.44179e-13, 1.91168e-13, 1.05427e-13, 5.77308e-14, 3.1389e-14, 1.69454e-14, 9.08272e-15, 4.83332e-15, 2.55336e-15, 1.339e-15, 6.96971e-16, 3.60058e-16, 1.84589e-16, 9.39012e-17, 4.73931e-17, 2.37293e-17, 1.17847e-17, 5.80466e-18},
{4.16221e-06, 4.20045e-05, 0.000224165, 0.000160465, 0.000332433, 0.00117499, 0.00383611, 0.00717808, 0.0107503, 0.0159182, 0.0226394, 0.0297432, 0.0367484, 0.0439936, 0.0518253, 0.0591955, 0.0641381, 0.0659295, 0.065037, 0.0623062, 0.0586976, 0.0548553, 0.0508768, 0.0465898, 0.041915, 0.0369678, 0.0319664, 0.0271353, 0.02266, 0.0186647, 0.0152061, 0.0122793, 0.00983689, 0.00781234, 0.00613877, 0.00475928, 0.00362906, 0.00271327, 0.00198322, 0.0014134, 0.000979767, 0.000659182, 0.00042967, 0.00027098, 0.000165231, 9.74029e-05, 5.55513e-05, 3.07016e-05, 1.6486e-05, 8.63383e-06, 4.43255e-06, 2.2456e-06, 1.13185e-06, 5.73032e-07, 2.94462e-07, 1.55134e-07, 8.44596e-08, 4.77086e-08, 2.79384e-08, 1.68768e-08, 1.04382e-08, 6.55754e-09, 4.15476e-09, 2.63991e-09, 1.67527e-09, 1.05878e-09, 6.65204e-10, 4.14986e-10, 2.5689e-10, 1.57737e-10, 9.60534e-11, 5.8004e-11, 3.47353e-11, 2.06288e-11, 1.21506e-11, 7.09863e-12, 4.11381e-12, 2.36506e-12, 1.34898e-12, 7.63414e-13, 4.28685e-13, 2.38871e-13, 1.32086e-13, 7.24827e-14, 3.94739e-14, 2.13351e-14, 1.14445e-14, 6.09285e-15, 3.21936e-15, 1.68827e-15, 8.78695e-16, 4.5389e-16, 2.32685e-16, 1.18381e-16, 5.97683e-17, 2.99443e-17, 1.48862e-17, 7.34246e-18, 3.59306e-18},
{4.19508e-06, 4.54036e-05, 0.000294051, 0.000155697, 0.000292895, 0.000836147, 0.00260713, 0.00564786, 0.00970735, 0.0149791, 0.0218178, 0.0289569, 0.0354815, 0.04158, 0.0480195, 0.0543231, 0.0588081, 0.0606926, 0.060364, 0.0585843, 0.0561192, 0.0534228, 0.0504891, 0.0471161, 0.0432365, 0.0389841, 0.0345794, 0.0302318, 0.0261045, 0.0223035, 0.0188796, 0.0158385, 0.0131593, 0.0108113, 0.00876475, 0.0069947, 0.00548052, 0.00420426, 0.00314862, 0.00229557, 0.00162496, 0.00111422, 0.000738753, 0.000473053, 0.000292403, 0.000174502, 0.000100637, 5.61783e-05, 3.04302e-05, 1.60479e-05, 8.27602e-06, 4.19719e-06, 2.10821e-06, 1.05797e-06, 5.35915e-07, 2.77058e-07, 1.47681e-07, 8.17357e-08, 4.70595e-08, 2.80849e-08, 1.7244e-08, 1.07962e-08, 6.83529e-09, 4.34712e-09, 2.76377e-09, 1.75078e-09, 1.10276e-09, 6.8974e-10, 4.28074e-10, 2.63508e-10, 1.60848e-10, 9.73512e-11, 5.84192e-11, 3.47587e-11, 2.05059e-11, 1.19955e-11, 6.95829e-12, 4.00266e-12, 2.28336e-12, 1.29181e-12, 7.24831e-13, 4.03368e-13, 2.22641e-13, 1.21889e-13, 6.61894e-14, 3.56528e-14, 1.90499e-14, 1.00973e-14, 5.30938e-15, 2.76968e-15, 1.43344e-15, 7.3606e-16, 3.75013e-16, 1.89581e-16, 9.50975e-17, 4.73344e-17, 2.33787e-17, 1.14574e-17, 5.57155e-18},
{4.356e-06, 4.44041e-05, 0.000283907, 0.000154777, 0.000305286, 0.00104224, 0.00347718, 0.00693546, 0.01064, 0.0150729, 0.02068, 0.0268104, 0.033288, 0.0403348, 0.0479969, 0.0550832, 0.0597774, 0.061566, 0.0610849, 0.0591776, 0.0566076, 0.0538174, 0.0507929, 0.0473282, 0.0433585, 0.039024, 0.0345511, 0.0301509, 0.0259852, 0.0221578, 0.0187173, 0.0156678, 0.0129867, 0.0106418, 0.00860174, 0.00684094, 0.00533838, 0.004076, 0.00303621, 0.00220028, 0.00154714, 0.00105316, 0.0006928, 0.000439929, 0.000269546, 0.000159403, 9.10872e-05, 5.03923e-05, 2.70691e-05, 1.41746e-05, 7.27316e-06, 3.68106e-06, 1.8526e-06, 9.36045e-07, 4.79841e-07, 2.52168e-07, 1.37e-07, 7.72906e-08, 4.52558e-08, 2.73636e-08, 1.69541e-08, 1.06754e-08, 6.78113e-09, 4.32017e-09, 2.74889e-09, 1.74192e-09, 1.09726e-09, 6.86269e-10, 4.25881e-10, 2.62131e-10, 1.5999e-10, 9.68217e-11, 5.80953e-11, 3.45622e-11, 2.03875e-11, 1.19246e-11, 6.916e-12, 3.97757e-12, 2.26853e-12, 1.28308e-12, 7.19701e-13, 4.00366e-13, 2.20892e-13, 1.20874e-13, 6.56036e-14, 3.53166e-14, 1.88582e-14, 9.98875e-15, 5.24842e-15, 2.73572e-15, 1.4147e-15, 7.25812e-16, 3.69465e-16, 1.86607e-16, 9.35203e-17, 4.65067e-17, 2.29489e-17, 1.12367e-17, 5.45945e-18},
{4.24547e-06, 3.80708e-05, 0.000159027, 0.000147613, 0.000277688, 0.00088616, 0.00283265, 0.00544837, 0.00846445, 0.0128465, 0.0189894, 0.0259604, 0.0330516, 0.040244, 0.0478646, 0.0551157, 0.0602123, 0.0624065, 0.0621714, 0.0603534, 0.0577819, 0.054945, 0.0518418, 0.0482671, 0.0441596, 0.0396701, 0.0350384, 0.0304891, 0.0261943, 0.022264, 0.0187483, 0.0156489, 0.0129382, 0.0105782, 0.00853271, 0.00677207, 0.00527306, 0.0040163, 0.00298352, 0.00215541, 0.00151038, 0.00102426, 0.000671052, 0.000424272, 0.000258769, 0.00015231, 8.66201e-05, 4.76974e-05, 2.55096e-05, 1.33072e-05, 6.80862e-06, 3.44087e-06, 1.73235e-06, 8.77557e-07, 4.52065e-07, 2.392e-07, 1.30988e-07, 7.44832e-08, 4.39109e-08, 2.66898e-08, 1.65961e-08, 1.04732e-08, 6.66081e-09, 4.24601e-09, 2.70228e-09, 1.7124e-09, 1.07856e-09, 6.74493e-10, 4.18516e-10, 2.57565e-10, 1.57185e-10, 9.51131e-11, 5.70639e-11, 3.39447e-11, 2.00207e-11, 1.17083e-11, 6.78935e-12, 3.90388e-12, 2.22592e-12, 1.25857e-12, 7.05687e-13, 3.92394e-13, 2.16381e-13, 1.18334e-13, 6.41814e-14, 3.45246e-14, 1.84198e-14, 9.74752e-15, 5.11656e-15, 2.66414e-15, 1.37612e-15, 7.0518e-16, 3.58519e-16, 1.80848e-16, 9.05161e-17, 4.49537e-17, 2.21534e-17, 1.0833e-17, 5.25658e-18},
{2.6023e-06, 2.0996e-05, 6.51172e-05, 9.41871e-05, 0.000148034, 0.000291241, 0.000704223, 0.00141856, 0.00252469, 0.00439193, 0.00751877, 0.0120581, 0.017926, 0.0248897, 0.0326516, 0.0404338, 0.0470294, 0.0518213, 0.0548712, 0.0565006, 0.0571406, 0.0570952, 0.0563517, 0.05476, 0.0523059, 0.0491505, 0.0455113, 0.0415816, 0.0375179, 0.0334442, 0.0294538, 0.0256111, 0.0219607, 0.018538, 0.0153747, 0.0125008, 0.00994118, 0.00771318, 0.00582371, 0.00426785, 0.00302825, 0.00207592, 0.00137259, 0.000874405, 0.000536469, 0.000317084, 0.000180748, 9.95605e-05, 5.31458e-05, 2.7603e-05, 1.4023e-05, 7.0154e-06, 3.48514e-06, 1.7364e-06, 8.77271e-07, 4.54445e-07, 2.43627e-07, 1.35898e-07, 7.88916e-08, 4.74368e-08, 2.93109e-08, 1.84476e-08, 1.17316e-08, 7.49018e-09, 4.77875e-09, 3.0369e-09, 1.91841e-09, 1.20307e-09, 7.48437e-10, 4.61708e-10, 2.8239e-10, 1.71234e-10, 1.02945e-10, 6.13675e-11, 3.62783e-11, 2.1271e-11, 1.23716e-11, 7.13865e-12, 4.08711e-12, 2.32207e-12, 1.30928e-12, 7.32698e-13, 4.06981e-13, 2.24387e-13, 1.22802e-13, 6.67108e-14, 3.59722e-14, 1.92533e-14, 1.02282e-14, 5.39289e-15, 2.82199e-15, 1.46546e-15, 7.55173e-16, 3.86138e-16, 1.95897e-16, 9.85978e-17, 4.92291e-17, 2.43805e-17, 1.19758e-17},
{2.6023e-06, 2.0996e-05, 6.51172e-05, 9.41871e-05, 0.000148034, 0.000291241, 0.000704223, 0.00141856, 0.00252469, 0.00439193, 0.00751877, 0.0120581, 0.017926, 0.0248897, 0.0326516, 0.0404338, 0.0470294, 0.0518213, 0.0548712, 0.0565006, 0.0571406, 0.0570952, 0.0563517, 0.05476, 0.0523059, 0.0491505, 0.0455113, 0.0415816, 0.0375179, 0.0334442, 0.0294538, 0.0256111, 0.0219607, 0.018538, 0.0153747, 0.0125008, 0.00994118, 0.00771318, 0.00582371, 0.00426785, 0.00302825, 0.00207592, 0.00137259, 0.000874405, 0.000536469, 0.000317084, 0.000180748, 9.95605e-05, 5.31458e-05, 2.7603e-05, 1.4023e-05, 7.0154e-06, 3.48514e-06, 1.7364e-06, 8.77271e-07, 4.54445e-07, 2.43627e-07, 1.35898e-07, 7.88916e-08, 4.74368e-08, 2.93109e-08, 1.84476e-08, 1.17316e-08, 7.49018e-09, 4.77875e-09, 3.0369e-09, 1.91841e-09, 1.20307e-09, 7.48437e-10, 4.61708e-10, 2.8239e-10, 1.71234e-10, 1.02945e-10, 6.13675e-11, 3.62783e-11, 2.1271e-11, 1.23716e-11, 7.13865e-12, 4.08711e-12, 2.32207e-12, 1.30928e-12, 7.32698e-13, 4.06981e-13, 2.24387e-13, 1.22802e-13, 6.67108e-14, 3.59722e-14, 1.92533e-14, 1.02282e-14, 5.39289e-15, 2.82199e-15, 1.46546e-15, 7.55173e-16, 3.86138e-16, 1.95897e-16, 9.85978e-17, 4.92291e-17, 2.43805e-17, 1.19758e-17}};

//-------------------------------------------
static const int nsrc = 27;   // Change form 26 as for 2015 data .  See JEC for 2017 94X
const char* srcnames[nsrc] = {"AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation", "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF","RelativePtBB", "RelativePtEC1", "RelativePtEC2","RelativePtHF","RelativeBal", "RelativeSample", "RelativeFSR", "RelativeStatFSR", "RelativeStatEC", "RelativeStatHF", "PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF"};
//-------------------------------------------
double intlumi[nHLTmx]={1., 1, 1, 1, 1, 1,1,1,1};
double lumiwt[nHLTmx]={1., 1, 1, 1, 1, 1,1,1,1};
//unsigned int l1trg[4], hlttr[8], tetrg[2];
unsigned int mypow_2[32];

//std::ofstream myfile;
//myfile.open("txt.log");

//const bool m_trigeff = true;
const int njetptmn=nHLTmx; // 8; //10
const int njetptbin=120;

#ifdef DIJETAVE
//const char* jethlt_name[nHLTmx]={"HLT_DiPFJetAve60_v","HLT_DiPFJetAve80_v", "HLT_DiPFJetAve140_v", "HLT_DiPFJetAve200_v", "HLT_DiPFJetAve260_v", "HLT_DiPFJetAve320_v", "HLT_DiPFJetAve400_v", "HLT_DiPFJetAve500_v"};
const char* jethlt_name[nHLTmx]={"HLT_PFJet60_v","HLT_PFJet80_v","HLT_PFJet140_v","HLT_PFJet200_v","HLT_PFJet260_v","HLT_PFJet320_v","HLT_PFJet400_v","HLT_PFJet450_v","HLT_PFJet500_v"};

//double leadingPtThreshold[njetptmn+1] ={90.0, 120.0, 180.0, 250.0, 320.0, 400.0, 480.0, 600.0, 2000.0};
//double leadingPtThreshold[njetptmn+1] ={83, 109, 172, 241, 309, 377, 462, 570, 3000.0}; //Fit Value dijet trigger
//double leadingPtThreshold[njetptmn+1] ={89, 118, 186, 259, 333, 403, 488, 528, 584, 688, 3000.0}; //trigger turn on for 2017 UL JetHT sample
double leadingPtThreshold[njetptmn+1] = {76, 101, 164, 232, 295, 373, 449, 507, 612, 3000.0};

/*
#ifdef MERGE_BIN
double leadingPtThreshold[njetptmn+1] ={83, 109, 176, 247, 318, 387, 477, 573, 3000.0}; //Fit Value dijet trigger 2018
#else
double leadingPtThreshold[njetptmn+1] ={66, 95, 155, 227, 286, 350, 444, 557, 3000.0}; //Fit Value dijet trigger 2016
#endif
//double leadingPtThreshold[njetptmn+1] ={83, 109, 172, 241, 309, 377, 462, 570, 3000.0}; //Fit Value dijet trigger 2017
*/

//double compres[njetptmn] = {1630, 5320, 62.1, 38.9, 27.0, 4.33, 1.23, 1.0};
//double compres[njetptmn] = {1630, 5320, 62.1, 38.9, 27.0, 4.33, 1.23, 1.0};

//const char* jethlt_lowest={"HLT_DiPFJetAve40_v"};
const char* jethlt_lowest={"HLT_PFJet40_v"};

//#else

#endif

#ifdef DIJETAVE
//double jethlt_thr[nHLTmx]={60,80,140,200,260,320,400,500};
double jethlt_thr[nHLTmx]={60,80,140,200,260,320,400,450,500};
//#else

#endif
double prescl[nHLTmx];

#ifdef TRACKSYS
const int ntype=3;
//const int ijet=2; 
#else
const int ntype=2;
//const int ijet=2;
#endif

const int njetetamn=1; // GMA 4;
#ifdef  LHAPDF
const int nnnmx=101;
double pdfwt[nnnmx];
  TH1F* h_genevtvarpdf[ntype][njetptmn][njetetamn][nvar][nnnmx];
  TH1* h_genevtvarpdf_2D[ntype][njetetamn][nvar][nnnmx];
#endif

#ifdef  JETENERGY//const int nsrc = 26;
//const int nsrc = 26;
const int njecmx=2*nsrc+1;
  TH1F* h_recoevtvarjec[ntype][njetptmn][njetetamn][nvar][njecmx];
  TH1* h_recoevtvarjec_2D[ntype][njetetamn][nvar][njecmx]; //For 2D
#elif defined(JETRESO)
const int njecmx = 3;
  TH1F* h_recoevtvarres[ntype][njetptmn][njetetamn][nvar][njecmx];
  TH1* h_recoevtvarres_2D[ntype][njetetamn][nvar][njecmx]; //For 2D
#else
const int njecmx=1;
#endif

//#ifdef  JETRESO
//  const int nGenReso = 3;
//  TH1F* h_genevtvarres[ntype][njetptmn][njetetamn][nvar][nGenReso];
//#else
const int nGenReso=1;
//const int nGenReso=njecmx;
//#endif

//int trgbit[nHLTmx]={10,11,12,13,14,16};
//double trgpas[nHLTmx+1]={0,0,0,0,0,0,0,0,0};

//const int njetetamn=3;
//double etarange[njetetamn] ={2.4}; //{3.0, 2.4, 1.8, 1.3};
double etarange[njetetamn] ={2.5};
double resetarange[njetetamn+4] ={0, 0.5, 1.0, 1.5}; //, 2.0, 2.5, 3.0, 3.5};
double par0[njetetamn+4]={1.02, 1.02, 1.022, 1.017, 0.98}; //, 0.9327};
double par1[njetetamn+4]={7.3e-6, -7.3e-6, -5.66e-6, -9.9e-6, 1.41e-4}; //, 4.6e-4};
double par2[njetetamn+4]={-8.2e-9, -8.2e-9, -3.58e-9, -4.18e-9, -6.104e-8}; //, -4.041e-7};
double particlept[4]={0.0, 0.25, 0.50, 1.00};

//const int ntype=4;
//const char* typname[ntype]={"Jets", "All particle", "All particle: P_{T}>0.25", "All particle: P_{T}>0.50", "All particle: P_{T}>1"};

#ifdef TRACKSYS
const char* typname[ntype]={"Jets", "Charged Particles"};
#else
const char* typname[ntype]={"Jets", "Charged Particles"};
#endif
static const int njetmx =30;

int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -2;
  for (int ix=0; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return -3;
}

double Phi_0_2pi(double x) {
  while (x >= 2*M_PI) x -= 2*M_PI;
  while (x <     0.)  x += 2*M_PI;
  return x;
}

double Phi_mpi_pi(double x) {
  while (x >= M_PI) x -= 2*M_PI;
  while (x < -M_PI) x += 2*M_PI;
  return x;
}

double dPhi(double phi1,double phi2){
  phi1=Phi_0_2pi(phi1);
  phi2=Phi_0_2pi(phi2);
  return Phi_mpi_pi(phi1-phi2);
}

 int sbitx(unsigned ival, int ibit) {
 unsigned den = mypow_2[ibit]; // unsigned(pow(2., double(ibit)));
 int isel = unsigned(ival/den)%2;
 //  int isel = unsigned(ival/den);
 //cout <<"iv "<< ival<<" "<<ibit<<" "<<den<<" "<<ival/den<<" "<<unsigned(ival/den)<<" "<<isel<<endl;

 return isel;
}

double respfun(double a, double b, double c, double x){
  double func=a+b*x+c*x*x;
  return func;
}

double JetCharge1(int charge, double candspt, double jpt, double k);

double candsmom(int charge, double candspt, double k);

double dotproduct(double candspx, double candspy, double candspz, double jpx, double jpy, double jpz, double jpt, double k);

double crossproduct(double candspx, double candspy, double candspz, double jpx, double jpy, double jpz, double jpt, double k);

struct triggervar{
  HepLorentzVector trg4v;
  bool		  both;
  bool            level1;
  bool            highl;
  int             ihlt;
  int             prescl;
};

//
// class declaration
//

class QCDEventShape : public edm::EDAnalyzer {
   public:
      explicit QCDEventShape(const edm::ParameterSet&);
      ~QCDEventShape();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 
  // int sbitx(unsigned ival, int ibit);

  bool isHistFill;
  bool isTrigger;
  bool isRECO[ntype][njetetamn];
  bool isMC;
  //bool isRECOJC;
  //bool isGENJC;
  //  bool isParticle; //Do we want particle level informations, other than jets ?
  //  bool isGenParticle; //Do we want Simulated particle level informations, other than jets ?
  bool isReconstruct; // otherwise Only generator level informations  
  //  bool isPartQCD; //For tracker variables, recosntruct QCD EVT variables
  bool isJetQCD;  //For Jet variables, recosntruct QCD EVT variables
  bool isGenJET; // Genjet information or note (for herwig/alpgen, donot store this ?)
  //  double trackPtThreshold; //Threshold of track Pt to store it in root file, -ve implies don't store

  //  double etarange; //Eta range of all jets
  double ptthreshold; //Pt threshold of JEC jets
  double leadingPtthreshold; //Pt threshold of JEC leading jet
  bool   isOtherAlgo; // store Kt4 and ak7 variables or not
  double weight=1; //weight for histogramme fit
  double weight2=1;
 
  std::string m_resolutions_file;
  std::string scalefile;

  std::string theHLTTag;
  //unsigned int mypow_2[32];
  int nevt;

  std::string theRootFileName;
  //TFile* //theFile;
  //TTree* //T1;

  //ifstream myfile ("example.txt");
  //std::ofstream myfile;
  //myfile.open("txt.log");
  TDirectoryFile *TUnfoldBinng2D =new TDirectoryFile("analyzeBasicPat2D","2D Unfolding Historgams");
  //TDirectoryFile *binning2D =new TDirectoryFile("Binning2D","Binning");
  
  //ESV
  //TH1F* h_recoevtvar[10][ntype][njetptmn][njetetamn][nvar];
  TH1F* h_recoevtvar[ntype][njetptmn][njetetamn][nvar];
  TH1F* h_recoevtfake[ntype][njetptmn][njetetamn][nvar];  //For fake
  TH1F* h_genevtmiss[ntype][njetptmn][njetetamn][nvar];   //For miss

  TH1F* h_genevtvar[ntype][njetptmn][njetetamn][nvar];    //For Gen
  TH1F* h_genevtvar2[ntype][njetptmn][njetetamn][nvar];
  
  TH2F* h_2devtvar[ntype][njetptmn][njetetamn][nvar];     //For RM
  TH2F* h_2ht;

  TH1F* vec_anglex[nhist];

  //Jet Charge 
  //Reco
  TH1F* h_recojc_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_D3J2[nkappa][njetptmn][njetetamn];

  //Reco Gluon jet
  TH1F* h_recojc_gjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_gjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_gjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_gjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_gjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_gjt_D3J2[nkappa][njetptmn][njetetamn];
  //
  TH1F* h_recojc_qgjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qgjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_qgjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qgjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_qgjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qgjt_D3J2[nkappa][njetptmn][njetetamn];
  //
  TH1F* h_recojc_agjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_agjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_agjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_agjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_agjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_agjt_D3J2[nkappa][njetptmn][njetetamn];
  //Reco b jet
  TH1F* h_recojc_bjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_bjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_bjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_bjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_bjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_bjt_D3J2[nkappa][njetptmn][njetetamn];
  //
  TH1F* h_recojc_qbjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qbjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_qbjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qbjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_qbjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qbjt_D3J2[nkappa][njetptmn][njetetamn];
  //
  TH1F* h_recojc_abjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_abjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_abjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_abjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_abjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_abjt_D3J2[nkappa][njetptmn][njetetamn];

  //Reco c jet
  TH1F* h_recojc_cjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_cjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_cjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_cjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_cjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_cjt_D3J2[nkappa][njetptmn][njetetamn];
  //
  TH1F* h_recojc_qcjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qcjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_qcjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qcjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_qcjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qcjt_D3J2[nkappa][njetptmn][njetetamn];
  //
  TH1F* h_recojc_acjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_acjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_acjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_acjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_acjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_acjt_D3J2[nkappa][njetptmn][njetetamn];

  //Reco s jet
  TH1F* h_recojc_sjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_sjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_sjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_sjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_sjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_sjt_D3J2[nkappa][njetptmn][njetetamn];
  //
  TH1F* h_recojc_qsjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qsjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_qsjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qsjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_qsjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qsjt_D3J2[nkappa][njetptmn][njetetamn];
  //
  TH1F* h_recojc_asjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_asjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_asjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_asjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_asjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_asjt_D3J2[nkappa][njetptmn][njetetamn];
  //Reco u jet
  TH1F* h_recojc_ujt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_ujt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_ujt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_ujt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_ujt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_ujt_D3J2[nkappa][njetptmn][njetetamn];
  //
  TH1F* h_recojc_qujt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qujt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_qujt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qujt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_qujt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qujt_D3J2[nkappa][njetptmn][njetetamn];
  //
  TH1F* h_recojc_aujt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_aujt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_aujt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_aujt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_aujt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_aujt_D3J2[nkappa][njetptmn][njetetamn];
  //Reco d jet
  TH1F* h_recojc_djt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_djt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_djt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_djt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_djt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_djt_D3J2[nkappa][njetptmn][njetetamn];
  //
  TH1F* h_recojc_qdjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qdjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_qdjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qdjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_qdjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_qdjt_D3J2[nkappa][njetptmn][njetetamn];
  //
  TH1F* h_recojc_adjt_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_adjt_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_adjt_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_adjt_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recojc_adjt_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recojc_adjt_D3J2[nkappa][njetptmn][njetetamn];

  //Gen
  TH1F* h_genjc_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_genjc_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_genjc_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_genjc_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_genjc_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_genjc_D3J2[nkappa][njetptmn][njetetamn];

  //RM
  TH2F* h_RM_D1J1[nkappa][njetptmn][njetetamn];
  TH2F* h_RM_D1J2[nkappa][njetptmn][njetetamn];

  TH2F* h_RM_D2J1[nkappa][njetptmn][njetetamn];
  TH2F* h_RM_D2J2[nkappa][njetptmn][njetetamn];

  TH2F* h_RM_D3J1[nkappa][njetptmn][njetetamn];
  TH2F* h_RM_D3J2[nkappa][njetptmn][njetetamn];

  //Recofake
  TH1F* h_recofake_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recofake_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recofake_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recofake_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_recofake_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_recofake_D3J2[nkappa][njetptmn][njetetamn];

  //Genmiss
  TH1F* h_genmiss_D1J1[nkappa][njetptmn][njetetamn];
  TH1F* h_genmiss_D1J2[nkappa][njetptmn][njetetamn];

  TH1F* h_genmiss_D2J1[nkappa][njetptmn][njetetamn];
  TH1F* h_genmiss_D2J2[nkappa][njetptmn][njetetamn];

  TH1F* h_genmiss_D3J1[nkappa][njetptmn][njetetamn];
  TH1F* h_genmiss_D3J2[nkappa][njetptmn][njetetamn];
  
  //Profile histogram
  TProfile* hprof;
  TProfile* hchpt; // jetpt vs. charge

  //static const int njetmx =30;
  //int npfjets; 
  int nchg;
  int pfjetmul[njetmx];
  float pfjetpx[njetmx], pfjetpy[njetmx], pfjetpz[njetmx], pfjeten[njetmx],  pfjetenuc[njetmx], neuemf[njetmx], neuhad[njetmx];
  float pfjetenscl[njetmx], pfjetensmr[njetmx];
  float jetpt, jeteta, jetphi; 
  int nallpf, ncharged;
  float thphi[nhist], thrust[nhist], anglex[nhist];
  float jtthan;
  int irunhlt, l1pres[nHLTmx],  hltpres[nHLTmx], compres[nHLTmx]; 
  static const int nprimx=150;
  int nprim, ntkpm[nprimx];
  //float  primdx[nprimx], primdy[nprimx], primdz[nprimx], 
  float primpr[nprimx];
  int irun, ilumi, ibrnc;
  unsigned int ievt;
  float inslumi;
  int nsicls, ntottrk;
//#ifdef FLAT 
  bool isFlat=1;
//#else 
  //bool isFlat=0;
//#endif

   float defweight=1.0, weighttrg=1., qlow=-10., qhigh=100000.;
//-------------------------------------------TunfoldBinning

//-------------------------------------------2D Bining using TUnfoldBinning ESV
  TUnfoldBinning *binsRec2D[ntype][njetetamn][nvar];
  TUnfoldBinning *RecoBinning2D[ntype][njetetamn][nvar];
  TUnfoldBinning *binsGen2D[ntype][njetetamn][nvar];
  TUnfoldBinning *GenBinning2D[ntype][njetetamn][nvar];
//-------------------------------------------2D Bining using TUnfoldBinning Jet Charge
  //Reco
  TUnfoldBinning *binsRec2D_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_D3J2[nkappa][njetetamn];

  // g jets
  TUnfoldBinning *binsRec2D_gjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_gjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_gjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_gjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_gjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_gjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *binsRec2D_qgjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qgjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_qgjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qgjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_qgjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qgjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *binsRec2D_agjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_agjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_agjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_agjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_agjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_agjt_D3J2[nkappa][njetetamn];

  // b jets
  TUnfoldBinning *binsRec2D_bjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_bjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_bjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_bjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_bjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_bjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *binsRec2D_qbjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qbjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_qbjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qbjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_qbjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qbjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *binsRec2D_abjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_abjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_abjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_abjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_abjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_abjt_D3J2[nkappa][njetetamn];

  // c jets
  TUnfoldBinning *binsRec2D_cjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_cjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_cjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_cjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_cjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_cjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *binsRec2D_qcjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qcjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_qcjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qcjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_qcjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qcjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *binsRec2D_acjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_acjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_acjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_acjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_acjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_acjt_D3J2[nkappa][njetetamn];

  // s jets
  TUnfoldBinning *binsRec2D_sjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_sjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_sjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_sjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_sjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_sjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *binsRec2D_qsjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qsjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_qsjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qsjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_qsjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qsjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *binsRec2D_asjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_asjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_asjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_asjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_asjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_asjt_D3J2[nkappa][njetetamn];

  // u jets
  TUnfoldBinning *binsRec2D_ujt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_ujt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_ujt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_ujt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_ujt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_ujt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *binsRec2D_qujt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qujt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_qujt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qujt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_qujt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qujt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *binsRec2D_aujt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_aujt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_aujt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_aujt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_aujt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_aujt_D3J2[nkappa][njetetamn];

  // d jets
  TUnfoldBinning *binsRec2D_djt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_djt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_djt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_djt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_djt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_djt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *binsRec2D_qdjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qdjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_qdjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qdjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_qdjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_qdjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *binsRec2D_adjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_adjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_adjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_adjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsRec2D_adjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsRec2D_adjt_D3J2[nkappa][njetetamn];

  //Gen
  TUnfoldBinning *binsGen2D_D1J1[nkappa][njetetamn];
  TUnfoldBinning *binsGen2D_D1J2[nkappa][njetetamn];

  TUnfoldBinning *binsGen2D_D2J1[nkappa][njetetamn];
  TUnfoldBinning *binsGen2D_D2J2[nkappa][njetetamn];

  TUnfoldBinning *binsGen2D_D3J1[nkappa][njetetamn];
  TUnfoldBinning *binsGen2D_D3J2[nkappa][njetetamn]; 

  //Reco
  TUnfoldBinning *RecoBinning2D_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_D3J2[nkappa][njetetamn];

  // g jets
  TUnfoldBinning *RecoBinning2D_gjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_gjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_gjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_gjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_gjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_gjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *RecoBinning2D_qgjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qgjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_qgjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qgjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_qgjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qgjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *RecoBinning2D_agjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_agjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_agjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_agjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_agjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_agjt_D3J2[nkappa][njetetamn];


  // b jets
  TUnfoldBinning *RecoBinning2D_bjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_bjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_bjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_bjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_bjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_bjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *RecoBinning2D_qbjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qbjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_qbjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qbjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_qbjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qbjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *RecoBinning2D_abjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_abjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_abjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_abjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_abjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_abjt_D3J2[nkappa][njetetamn];

  // c jets
  TUnfoldBinning *RecoBinning2D_cjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_cjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_cjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_cjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_cjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_cjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *RecoBinning2D_qcjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qcjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_qcjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qcjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_qcjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qcjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *RecoBinning2D_acjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_acjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_acjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_acjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_acjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_acjt_D3J2[nkappa][njetetamn];

  // s jets
  TUnfoldBinning *RecoBinning2D_sjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_sjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_sjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_sjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_sjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_sjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *RecoBinning2D_qsjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qsjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_qsjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qsjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_qsjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qsjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *RecoBinning2D_asjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_asjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_asjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_asjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_asjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_asjt_D3J2[nkappa][njetetamn];

  // u jets
  TUnfoldBinning *RecoBinning2D_ujt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_ujt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_ujt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_ujt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_ujt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_ujt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *RecoBinning2D_qujt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qujt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_qujt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qujt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_qujt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qujt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *RecoBinning2D_aujt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_aujt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_aujt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_aujt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_aujt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_aujt_D3J2[nkappa][njetetamn];

  // d jets
  TUnfoldBinning *RecoBinning2D_djt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_djt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_djt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_djt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_djt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_djt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *RecoBinning2D_qdjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qdjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_qdjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qdjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_qdjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_qdjt_D3J2[nkappa][njetetamn];
  //
  TUnfoldBinning *RecoBinning2D_adjt_D1J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_adjt_D1J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_adjt_D2J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_adjt_D2J2[nkappa][njetetamn];

  TUnfoldBinning *RecoBinning2D_adjt_D3J1[nkappa][njetetamn];
  TUnfoldBinning *RecoBinning2D_adjt_D3J2[nkappa][njetetamn];

  //Gen
  TUnfoldBinning *GenBinning2D_D1J1[nkappa][njetetamn];
  TUnfoldBinning *GenBinning2D_D1J2[nkappa][njetetamn];

  TUnfoldBinning *GenBinning2D_D2J1[nkappa][njetetamn];
  TUnfoldBinning *GenBinning2D_D2J2[nkappa][njetetamn];

  TUnfoldBinning *GenBinning2D_D3J1[nkappa][njetetamn];
  TUnfoldBinning *GenBinning2D_D3J2[nkappa][njetetamn];

//-------------------------------------------ESV 2D TUnfoldBinning Histograms
  TH1* h_recovar_2D[ntype][njetetamn][nvar]; //Reco
  TH1* h_recofake_2D[ntype][njetetamn][nvar];//For fake
  //TH1* h_recofakeOutE_2D[type][njetetamn][nvar];//For fake
  //TH1* h_recofakeOutHT_2D[type][njetetamn][nvar];//For fake
   
  TH1* h_genvar_2D[ntype][njetetamn][nvar]; //Gen
  TH1* h_genmiss_2D[ntype][njetetamn][nvar]; //For Miss
  //TH1* h_genmissOutE_2D[type][njetetamn][nvar]; //For Miss
  //TH1* h_genmissOutHT_2D[type][njetetamn][nvar]; //For Miss

  TH2* RM_2D[ntype][njetetamn][nvar];
//-------------------------------------------Jet Charge 2D TUnfoldBinning Histograms
  //Reco
  TH1* h_recovar_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_2D_D3J2[nkappa][njetetamn];
  
  //Gluon Jet
  TH1* h_recovar_gjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_gjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_gjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_gjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_gjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_gjt_2D_D3J2[nkappa][njetetamn];

  //
  TH1* h_recovar_qgjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_qgjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_qgjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_qgjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_qgjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_qgjt_2D_D3J2[nkappa][njetetamn];

  //
  TH1* h_recovar_agjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_agjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_agjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_agjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_agjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_agjt_2D_D3J2[nkappa][njetetamn];

  // b jets
  TH1* h_recovar_bjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_bjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_bjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_bjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_bjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_bjt_2D_D3J2[nkappa][njetetamn];
  
  //
  TH1* h_recovar_qbjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_qbjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_qbjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_qbjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_qbjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_qbjt_2D_D3J2[nkappa][njetetamn];
  //
  TH1* h_recovar_abjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_abjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_abjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_abjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_abjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_abjt_2D_D3J2[nkappa][njetetamn];

  // c jets
  TH1* h_recovar_cjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_cjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_cjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_cjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_cjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_cjt_2D_D3J2[nkappa][njetetamn];
  //
  TH1* h_recovar_qcjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_qcjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_qcjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_qcjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_qcjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_qcjt_2D_D3J2[nkappa][njetetamn];
  //
  TH1* h_recovar_acjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_acjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_acjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_acjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_acjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_acjt_2D_D3J2[nkappa][njetetamn];

  // s jets
  TH1* h_recovar_sjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_sjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_sjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_sjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_sjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_sjt_2D_D3J2[nkappa][njetetamn];
  //
  TH1* h_recovar_qsjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_qsjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_qsjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_qsjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_qsjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_qsjt_2D_D3J2[nkappa][njetetamn];
  //
  TH1* h_recovar_asjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_asjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_asjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_asjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_asjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_asjt_2D_D3J2[nkappa][njetetamn];

  // d jets
  TH1* h_recovar_djt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_djt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_djt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_djt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_djt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_djt_2D_D3J2[nkappa][njetetamn];
  //
  TH1* h_recovar_qdjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_qdjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_qdjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_qdjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_qdjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_qdjt_2D_D3J2[nkappa][njetetamn];
  //
  TH1* h_recovar_adjt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_adjt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_adjt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_adjt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_adjt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_adjt_2D_D3J2[nkappa][njetetamn];

  // u jets
  TH1* h_recovar_ujt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_ujt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_ujt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_ujt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_ujt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_ujt_2D_D3J2[nkappa][njetetamn];
  //
  TH1* h_recovar_qujt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_qujt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_qujt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_qujt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_qujt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_qujt_2D_D3J2[nkappa][njetetamn];
  //
  TH1* h_recovar_aujt_2D_D1J1[nkappa][njetetamn];
  TH1* h_recovar_aujt_2D_D1J2[nkappa][njetetamn];

  TH1* h_recovar_aujt_2D_D2J1[nkappa][njetetamn];
  TH1* h_recovar_aujt_2D_D2J2[nkappa][njetetamn];

  TH1* h_recovar_aujt_2D_D3J1[nkappa][njetetamn];
  TH1* h_recovar_aujt_2D_D3J2[nkappa][njetetamn];

  //Gen
  TH1* h_genvar_2D_D1J1[nkappa][njetetamn];
  TH1* h_genvar_2D_D1J2[nkappa][njetetamn];

  TH1* h_genvar_2D_D2J1[nkappa][njetetamn];
  TH1* h_genvar_2D_D2J2[nkappa][njetetamn];

  TH1* h_genvar_2D_D3J1[nkappa][njetetamn];
  TH1* h_genvar_2D_D3J2[nkappa][njetetamn];

  //Recofake
  TH1* h_recofake_2D_D1J1[nkappa][njetetamn];
  TH1* h_recofake_2D_D1J2[nkappa][njetetamn];

  TH1* h_recofake_2D_D2J1[nkappa][njetetamn];
  TH1* h_recofake_2D_D2J2[nkappa][njetetamn];

  TH1* h_recofake_2D_D3J1[nkappa][njetetamn];
  TH1* h_recofake_2D_D3J2[nkappa][njetetamn];
  
  //Genmiss
  TH1* h_genmiss_2D_D1J1[nkappa][njetetamn];
  TH1* h_genmiss_2D_D1J2[nkappa][njetetamn];

  TH1* h_genmiss_2D_D2J1[nkappa][njetetamn];
  TH1* h_genmiss_2D_D2J2[nkappa][njetetamn];

  TH1* h_genmiss_2D_D3J1[nkappa][njetetamn];
  TH1* h_genmiss_2D_D3J2[nkappa][njetetamn];

  //RM
  TH2* RM_2D_D1J1[nkappa][njetetamn];
  TH2* RM_2D_D1J2[nkappa][njetetamn];

  TH2* RM_2D_D2J1[nkappa][njetetamn];
  TH2* RM_2D_D2J2[nkappa][njetetamn];

  TH2* RM_2D_D3J1[nkappa][njetetamn];
  TH2* RM_2D_D3J2[nkappa][njetetamn];
//-------------------------------------------
  //TH1F* recojt_hist;
  //TH1F* recojt_pt[njetetamn][nHLTmx];
  TH1F* recojt_pt[njetetamn];
  TH1F* recojt_eta;
  TH1F* recojt_phi;

  TH1F* recojtg_pt[njetetamn];
  TH1F* recojtg_eta;
  TH1F* recojtg_phi;

  TH1F* recojtqg_pt[njetetamn];
  TH1F* recojtqg_eta;
  TH1F* recojtqg_phi;

  TH1F* recojtag_pt[njetetamn];
  TH1F* recojtag_eta;
  TH1F* recojtag_phi;

  TH1F* recojtb_pt[njetetamn];
  TH1F* recojtb_eta;
  TH1F* recojtb_phi;

  TH1F* recojtqb_pt[njetetamn];
  TH1F* recojtqb_eta;
  TH1F* recojtqb_phi;

  TH1F* recojtab_pt[njetetamn];
  TH1F* recojtab_eta;
  TH1F* recojtab_phi;

  TH1F* recojtc_pt[njetetamn];
  TH1F* recojtc_eta;
  TH1F* recojtc_phi;

  TH1F* recojtqc_pt[njetetamn];
  TH1F* recojtqc_eta;
  TH1F* recojtqc_phi;

  TH1F* recojtac_pt[njetetamn];
  TH1F* recojtac_eta;
  TH1F* recojtac_phi;

  TH1F* recojts_pt[njetetamn];
  TH1F* recojts_eta;
  TH1F* recojts_phi;

  TH1F* recojtqs_pt[njetetamn];
  TH1F* recojtqs_eta;
  TH1F* recojtqs_phi;

  TH1F* recojtas_pt[njetetamn];
  TH1F* recojtas_eta;
  TH1F* recojtas_phi;

  TH1F* recojtu_pt[njetetamn];
  TH1F* recojtu_eta;
  TH1F* recojtu_phi;

  TH1F* recojtqu_pt[njetetamn];
  TH1F* recojtqu_eta;
  TH1F* recojtqu_phi;

  TH1F* recojtau_pt[njetetamn];
  TH1F* recojtau_eta;
  TH1F* recojtau_phi;

  TH1F* recojtd_pt[njetetamn];
  TH1F* recojtd_eta;
  TH1F* recojtd_phi;

  TH1F* recojtqd_pt[njetetamn];
  TH1F* recojtqd_eta;
  TH1F* recojtqd_phi;

  TH1F* recojtad_pt[njetetamn];
  TH1F* recojtad_eta;
  TH1F* recojtad_phi;

  TH1F* recojtallave_pt[njetetamn];
  TH1F* recojtallavewt1_pt[njetetamn];

  TH1F* recojtave_pt[njetetamn][nHLTmx];
  TH1F* recojtavewt1_pt[njetetamn][nHLTmx];
  TH1F* recojt1_pt[njetetamn];
  TH1F* recojt1_eta;
  TH1F* recojt1_phi;
 
  TH1F* recojt1g_pt[njetetamn];
  TH1F* recojt1g_eta;
  TH1F* recojt1g_phi;

  TH1F* recojt1qg_pt[njetetamn];
  TH1F* recojt1qg_eta;
  TH1F* recojt1qg_phi;

  TH1F* recojt1ag_pt[njetetamn];
  TH1F* recojt1ag_eta;
  TH1F* recojt1ag_phi;

  TH1F* recojt1b_pt[njetetamn];
  TH1F* recojt1b_eta;
  TH1F* recojt1b_phi;

  TH1F* recojt1qb_pt[njetetamn];
  TH1F* recojt1qb_eta;
  TH1F* recojt1qb_phi;

  TH1F* recojt1ab_pt[njetetamn];
  TH1F* recojt1ab_eta;
  TH1F* recojt1ab_phi;

  TH1F* recojt1c_pt[njetetamn];
  TH1F* recojt1c_eta;
  TH1F* recojt1c_phi;

  TH1F* recojt1qc_pt[njetetamn];
  TH1F* recojt1qc_eta;
  TH1F* recojt1qc_phi;

  TH1F* recojt1ac_pt[njetetamn];
  TH1F* recojt1ac_eta;
  TH1F* recojt1ac_phi;

  TH1F* recojt1s_pt[njetetamn];
  TH1F* recojt1s_eta;
  TH1F* recojt1s_phi;

  TH1F* recojt1qs_pt[njetetamn];
  TH1F* recojt1qs_eta;
  TH1F* recojt1qs_phi;
  
  TH1F* recojt1as_pt[njetetamn];
  TH1F* recojt1as_eta;
  TH1F* recojt1as_phi;

  TH1F* recojt1u_pt[njetetamn];
  TH1F* recojt1u_eta;
  TH1F* recojt1u_phi;

  TH1F* recojt1qu_pt[njetetamn];
  TH1F* recojt1qu_eta;
  TH1F* recojt1qu_phi;
 
  TH1F* recojt1au_pt[njetetamn];
  TH1F* recojt1au_eta;
  TH1F* recojt1au_phi;

  TH1F* recojt1d_pt[njetetamn];
  TH1F* recojt1d_eta;
  TH1F* recojt1d_phi;
  
  TH1F* recojt1qd_pt[njetetamn];
  TH1F* recojt1qd_eta;
  TH1F* recojt1qd_phi;

  TH1F* recojt1ad_pt[njetetamn];
  TH1F* recojt1ad_eta;
  TH1F* recojt1ad_phi;


  TH1F* recojt2_pt[njetetamn];
  TH1F* recojt2_eta;
  TH1F* recojt2_phi;
  
  TH1F* recojt2g_pt[njetetamn];
  TH1F* recojt2g_eta;
  TH1F* recojt2g_phi;

  TH1F* recojt2qg_pt[njetetamn];
  TH1F* recojt2qg_eta;
  TH1F* recojt2qg_phi;

  TH1F* recojt2ag_pt[njetetamn];
  TH1F* recojt2ag_eta;
  TH1F* recojt2ag_phi;

  TH1F* recojt2b_pt[njetetamn];
  TH1F* recojt2b_eta;
  TH1F* recojt2b_phi;

  TH1F* recojt2qb_pt[njetetamn];
  TH1F* recojt2qb_eta;
  TH1F* recojt2qb_phi;

  TH1F* recojt2ab_pt[njetetamn];
  TH1F* recojt2ab_eta;
  TH1F* recojt2ab_phi;

  TH1F* recojt2c_pt[njetetamn];
  TH1F* recojt2c_eta;
  TH1F* recojt2c_phi;

  TH1F* recojt2qc_pt[njetetamn];
  TH1F* recojt2qc_eta;
  TH1F* recojt2qc_phi;

  TH1F* recojt2ac_pt[njetetamn];
  TH1F* recojt2ac_eta;
  TH1F* recojt2ac_phi;

  TH1F* recojt2s_pt[njetetamn];
  TH1F* recojt2s_eta;
  TH1F* recojt2s_phi;

  TH1F* recojt2qs_pt[njetetamn];
  TH1F* recojt2qs_eta;
  TH1F* recojt2qs_phi;

  TH1F* recojt2as_pt[njetetamn];
  TH1F* recojt2as_eta;
  TH1F* recojt2as_phi;

  TH1F* recojt2u_pt[njetetamn];
  TH1F* recojt2u_eta;
  TH1F* recojt2u_phi;

  TH1F* recojt2qu_pt[njetetamn];
  TH1F* recojt2qu_eta;
  TH1F* recojt2qu_phi;

  TH1F* recojt2au_pt[njetetamn];
  TH1F* recojt2au_eta;
  TH1F* recojt2au_phi;

  TH1F* recojt2d_pt[njetetamn];
  TH1F* recojt2d_eta;
  TH1F* recojt2d_phi;

  TH1F* recojt2qd_pt[njetetamn];
  TH1F* recojt2qd_eta;
  TH1F* recojt2qd_phi;

  TH1F* recojt2ad_pt[njetetamn];
  TH1F* recojt2ad_eta;
  TH1F* recojt2ad_phi;

  TH1F* recojt3_pt[njetetamn];
  TH1F* recojt3_eta;
  TH1F* recojt3_phi;

  TH1F* recoht2_pt[njetetamn];

  TH1F* hjetdpt[njetetamn];
  TH1F* hjetdphi[njetetamn];
  TH1F* hjetptbypl[njetetamn];
  TH1F* hjetpt2bypt1[njetetamn];
  TH1F* hjetpt3bypt2[njetetamn];
  //TH1F* recochg_hist;
  TH1F* recochg_pt;
  TH1F* recochg_eta;
  TH1F* recochg_phi;

  TH1F* recochg1_pt;
  TH1F* recochg1_eta;
  TH1F* recochg1_phi;

  TH1F* recochg2_pt;
  TH1F* recochg2_eta;
  TH1F* recochg2_phi;

  TH1F* recochg3_pt;
  TH1F* recochg3_eta;
  TH1F* recochg3_phi;
//-------------------------------------------
  //TH1F* genjt_hist;
  TH1F* genjt_pt[njetetamn];
  TH1F* genjt_eta;
  TH1F* genjt_phi;
  TH1F* genjtallave_pt[njetetamn];

  TH1F* genjt1_pt[njetetamn];
  TH1F* genjt1_eta;
  TH1F* genjt1_phi;

  TH1F* genjt2_pt[njetetamn];
  TH1F* genjt2_eta;
  TH1F* genjt2_phi;

  TH1F* genjt3_pt[njetetamn];
  TH1F* genjt3_eta;
  TH1F* genjt3_phi;

  TH1F* genjetdpt[njetetamn];
  TH1F* genjetdphi[njetetamn];
  TH1F* genjetptbypl[njetetamn];
  TH1F* genjetpt2bypt1[njetetamn];
  TH1F* genjetpt3bypt2[njetetamn];

  //TH1F* genchg_hist;
  TH1F* genchg_pt;
  TH1F* genchg_eta;
  TH1F* genchg_phi;

  TH1F* genchg1_pt;
  TH1F* genchg1_eta;
  TH1F* genchg1_phi;

  TH1F* genchg2_pt;
  TH1F* genchg2_eta;
  TH1F* genchg2_phi;

  TH1F* genchg3_pt;
  TH1F* genchg3_eta;
  TH1F* genchg3_phi;
/*
  TH1F* genneu_hist;
  TH1F* genneu_pt;
  TH1F* genneu_eta;
  TH1F* genneu_phi;

  TH1F* genjt_oth_pt[njetetamn];
  TH1F* genjt_oth_eta;
  TH1F* genjt_oth_phi;

  TH1F* genchg_oth_hist;
  TH1F* genchg_oth_pt;
  TH1F* genchg_oth_eta;
  TH1F* genchg_oth_phi;

  TH1F* genneu_oth_hist;
  TH1F* genneu_oth_pt;
  TH1F* genneu_oth_eta;
  TH1F* genneu_oth_phi;
*/
  
  //Response hist
  TH2F* resp_jet[njetetamn+2];
  TH1F* resp_jet1[njetetamn+2];

  TH1F* prim_hist[nHLTmx+1];
  TH1F* prim_sel[nHLTmx+1];

  TH1F* prim_hist_rewt[nHLTmx+1];
  TH1F* prim_sel_rewt[nHLTmx+1];

  TH2F* prim_correl;

  TH1F* prim_alltrk[2];
  TH1F* prim_seltrk[2];
  TH1F* prim_goodtrk[2];
  TH1F* prim_dx[2];
  TH1F* prim_dy[2];
  TH2F* prim_dxy[2];
  TH1F* prim_dz[2];  
  TH1F* prim_prob[2];

  TH1F* h_jetpt[nHLTmx][njetetamn];
  TH1F* h_jeteta[nHLTmx];
  TH1F* h_jetphi[nHLTmx][njetetamn];
  TH1F* h_njets[njetetamn];
  TH1F* h_nchg[njetetamn];

  TH1F* gen_njets[njetetamn];

  TH1F* trgjet_angle[nHLTmx][2];
  TH2F* trgjet_2dangle[nHLTmx][2];
  TH1F* trgjet_pt[nHLTmx][2];
  TH1F* trgjet_eta[nHLTmx][2];
  TH1F* trgjet_phi[nHLTmx][2];
  TH1F* prbjet_pt[nHLTmx][2];
  TH1F* prbjet_eta[nHLTmx][2];
  TH1F* prbjet_phi[nHLTmx][2];

  //Dijet trigger efficiency
  TH1F* hlt_dijettag[nHLTmx][njetetamn];
  TH1F* hlt_dijetprob[nHLTmx][njetetamn];

  //Trigger Normal case
  TH1F* counthist; 

//-------------------------------------------Member data
  edm::EDGetTokenT<GenEventInfoProduct> generator1_;
  edm::EDGetTokenT<pat::JetCollection> jetSrcToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > genSrcToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> PFSrcToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  edm::EDGetTokenT<reco::GenJetCollection> genjetToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileup_;
  edm::EDGetTokenT<reco::PFJetCollection> ak5PFjetToken_;
  edm::EDGetTokenT<reco::GenJetCollection> ak5GenJetToken_;
  const edm::EDGetTokenT<std::vector<double> > pdfCTEQWeightsInputToken_;
  const edm::EDGetTokenT<std::vector<double> > pdfMMTHWeightsInputToken_;
  const edm::EDGetTokenT<std::vector<double> > pdfNNPDFWeightsInputToken_;
  const edm::EDGetTokenT<LHERunInfoProduct> LHERunInfoToken_;
  const edm::EDGetTokenT<LHEEventProduct> lheEventProductToken_;
  edm::EDGetTokenT<double> m_rho_token;

  edm::EDGetTokenT< double > prefweight_token;
  edm::EDGetTokenT< double > prefweightup_token;
  edm::EDGetTokenT< double > prefweightdown_token;


  float qscale;
  float wtfact; //MC : eventinfo->weight(); data : hltpres[ihltfill]*l1pres[ihltfill];
  int procid, npilup1, npilup2; //1:-5 to -1, 2:0 to 3

  int idall;
  float xfrac1, xfrac2, xpdf1, xpdf2;  

  //HLTConfigProvider hltConfig_;
  HLTPrescaleProvider hltPrescaleProvider_;
  int nreco, naa, nbb, ncc;

std::vector<JetCorrectionUncertainty*> vsrc; // (nsrc);
reweight::PoissonMeanShifter PShiftUp_;
reweight::PoissonMeanShifter PShiftDown_;
edm::LumiReWeighting *LumiWeights_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
  QCDEventShape::QCDEventShape(const edm::ParameterSet& iConfig):
  generator1_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("evtinfo"))),
  jetSrcToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
  genSrcToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genSrc"))),
  PFSrcToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfSrc"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metSrc"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  beamSpot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bsSrc"))),
  genjetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genjetSrc"))),
  pileup_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSrc"))),
  ak5PFjetToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("ak5pfJetSrc"))),
  ak5GenJetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("ak5genJetSrc"))),
  pdfCTEQWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFCTEQWeightsInputTag"))),
  pdfMMTHWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFMMTHWeightsInputTag"))),
  pdfNNPDFWeightsInputToken_(consumes<std::vector<double> >(iConfig.getParameter<edm::InputTag>("PDFNNPDFWeightsInputTag"))),
  LHERunInfoToken_(consumes<LHERunInfoProduct, edm::InRun >(iConfig.getParameter<edm::InputTag>("LHERunInfoProductInputTag"))),
  lheEventProductToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEEventProductInputTag"))),
  prefweight_token(consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"))),
  prefweightup_token(consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"))),
  prefweightdown_token(consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"))),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  m_rho_token = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
  //m_resolutions_file = iConfig.getParameter<edm::FileEEInPath>("resolutionsFile").fullPath();
  //scalefile = iConfig.getParameter<edm::FileInPath>("scaleFactorsFile").fullPath();
  isHistFill = iConfig.getUntrackedParameter<bool>("HistFill", true);
  //isHistFill2 = pset.getUntrackedParameter<bool>("HistFill2", false);                                            
  isTrigger = iConfig.getUntrackedParameter<bool>("Trigger", true);
  //isRECO = iConfig.getUntrackedParameter<bool>("RECO", false);
  isMC = iConfig.getUntrackedParameter<bool>("MonteCarlo", false);
  isReconstruct = iConfig.getUntrackedParameter<bool>("Reconstruct", true);
  isJetQCD = iConfig.getUntrackedParameter<bool>("JetQCD", false);
  isGenJET = iConfig.getUntrackedParameter<bool>("GenJET", false);
  //etarange = iConfig.getUntrackedParameter<double>("EtaRange", 5.0);
  ptthreshold = iConfig.getUntrackedParameter<double>("PtThreshold", 10.0);
  //leadingPtthreshold = iConfig.getUntrackedParameter<double>("LeadingPtThreshold", 40.0);
  isOtherAlgo = iConfig.getUntrackedParameter<bool>("OtherAlgo", false);
  weight2 = iConfig.getUntrackedParameter<double>("HistWeight", 1.0);
  weight = weight2;
  theHLTTag = iConfig.getUntrackedParameter<string>("HLTTag", "HLT");
  theRootFileName = iConfig.getUntrackedParameter<string>("RootFileName");
  

  //theFile = new TFile(theRootFileName.c_str(), "RECREATE");
  //theFile->cd();
  //T1 = new TTree("T1", "QCDEvt");

  //T1->Branch("irun", &irun, "irun/I");  
  //T1->Branch("ilumi", &ilumi, "ilumi/I");  
  //T1->Branch("ievt", &ievt, "ievt/i");
  //T1->Branch("ibrnc", &ibrnc, "ibrnc/I");  
  //T1->Branch("nsicls", &nsicls, "nsicls/I");  //to pf neutral paritcle (excluding HF)
  //T1->Branch("ntottrk", &ntottrk, "ntottrk/I");  //total pfcharged particle (HF)

  //T1->Branch("jetpt", &jetpt, "jetpt/F");
  //T1->Branch("jeteta", &jeteta, "jeteta/F");
  //T1->Branch("jetphi", &jetphi, "jetphi/F");

  //T1->Branch("nallpf", &nallpf, "nallpf/I");
  //T1->Branch("ncharged", &ncharged, "ncharged/I");
  //T1->Branch("jtthan",&jtthan,"jtthan/F");
  //T1->Branch("thphi",thphi,"thphi[10]/F");

  //T1->Branch("thrust",thrust,"thrust[10]/F");
  //T1->Branch("anglex",anglex,"anglex[10]/F");

  //cout << "Testing 1 ==== " <<njecmx<< endl;
                                                                                                                                                                        
//-------------------------------------------Define TUnfoldBinning ESV
  char RecoBinName[100], GenBinName[100], Axisname[100]; 
   for (int ityp=0; ityp<ntype; ityp++) {
      for (int iet=0; iet<njetetamn; iet++) {
        for (int ij=0; ij<nvar; ij++) {
           if (isItUsed(ij)) {
//-------------------------------------------2D Binning Reco
         sprintf(RecoBinName, "Detector2d_typ_%i_eta%i_%i", ityp, iet, ij);
         binsRec2D[ityp][iet][ij] = new TUnfoldBinning(RecoBinName);

         sprintf(RecoBinName, "Recobin2d_typ_%i_eta%i_%i",ityp, iet, ij);
         RecoBinning2D[ityp][iet][ij]= binsRec2D[ityp][iet][ij]->AddBinning(RecoBinName);

         sprintf(Axisname, "var_%i", ij);
         RecoBinning2D[ityp][iet][ij]->AddAxis(Axisname,(ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij],(ityp==0) ? rbinrngs0[ij]: rbinrngs1[ij],false,false);
         sprintf(Axisname, "ht");
         RecoBinning2D[ityp][iet][ij]->AddAxis(Axisname, nHLTmx, recohtbins, false,false);
//-------------------------------------------2D Binning Gen
         sprintf(GenBinName, "Generator2d_typ_%i_eta%i_%i",ityp, iet, ij);
         binsGen2D[ityp][iet][ij] = new TUnfoldBinning(GenBinName);

         sprintf(GenBinName, "Genbin2d_typ_%i_eta%i_%i", ityp, iet, ij);
         GenBinning2D[ityp][iet][ij]= binsGen2D[ityp][iet][ij]->AddBinning(GenBinName);

         sprintf(Axisname, "var_%i", ij);
         GenBinning2D[ityp][iet][ij]->AddAxis(Axisname,(ityp==0) ? nbinsx0[ij] : nbinsx1[ij],(ityp==0) ? binrngs0[ij]: binrngs1[ij],false,false);
         sprintf(Axisname, "ht");
         GenBinning2D[ityp][iet][ij]->AddAxis(Axisname, nHLTmx, recohtbins,false,false);
	   }
	 }
      }
    }
//-------------------------------------------
/*
	for (int iet=0; iet<njetetamn; iet++) {
		for (int id=1; id<2; id++){
			for (int ij=1; ij<2; ij++){
				for (int ik=1; ik<12; ik++){
					sprintf(RecoBinName, "Detector2d_D%i_j%i_k%i_eta%i", id, ij, ik, iet);
					binsRectest2D[id][ij][ik][iet] = new TUnfoldBinning(RecoBinName);
					cout <<"RecoBinName"<<RecoBinName<<endl;

					sprintf(RecoBinName, "Recobin2d_D%i_j%i_k%i_eta%i",id, ij, ik, iet);
					RecoBinningtest2D[id][ij][ik][iet]= binsRectest2D[id][ij][ik][iet]->AddBinning(RecoBinName);
					cout<<"RecoBinName"<<RecoBinName<<endl;				
			
					sprintf(Axisname, "var_%i_%i_%i", id, ij, ik);
					cout<<"Axisname"<<Axisname<<endl;
					RecoBinningtest2D[id][ij][ik][iet]->AddAxis(Axisname, jcbins[ik], jcminran[ik], jcmaxran[ik] ,false,false);
					sprintf(Axisname, "ht");
		        		RecoBinningtest2D[id][ij][ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false,false);		
					}
				}
			}
		}

char name1[200], title1[200];
	for (int iet=0; iet<njetetamn; iet++) {
		for (int id=1; id<2; id++){
			for (int ij=1; ij<2; ij++){
				for (int ik=1; ik<12; ik++){
                			if (isReconstruct) {
                        			sprintf(name1, "dd_reco_D%i_j%i_k%i_eta%i", id, ij, ik, iet);
						cout<<"name1"<<name1<<endl;
                        			sprintf(title1, "2D Reco Q_{D%i,}_{j%i}^{k=%i} %g", id, ij, ik, etarange[iet]);
						cout<<"title1"<<title1<<endl;
                        			h_recovartest_2D[id][ij][ik][iet] = binsRectest2D[id][ij][ik][iet]->CreateHistogram(name1,false,0,title1); //false : global bin ID
                        			h_recovartest_2D[id][ij][ik][iet]->Sumw2();
						}
					}
				}
			}
		}
*/
//-------------------------------------------
  char name[200];
  char title[200];

  for (int ityp=0; ityp<ntype; ityp++) {
      for (int iet=0; iet<njetetamn; iet++) {
        for (int ij=0; ij<nvar; ij++) {
             if (isItUsed(ij)) {
	     if (isReconstruct) { 
              sprintf(name, "dd_reco_typ_%i_eta%i_%i", ityp, iet, ij);
	      sprintf(title, "2D Reco %s %g %s", typname[ityp], etarange[iet], vartitle[ij]);
              h_recovar_2D[ityp][iet][ij] = binsRec2D[ityp][iet][ij]->CreateHistogram(name,false,0,title); //false : global bin ID
              h_recovar_2D[ityp][iet][ij]->Sumw2();              

              sprintf(name, "dd_fake_reco_typ_%i_eta%i_%i", ityp, iet, ij);
              sprintf(title, "2D Fake Reco %s %g %s", typname[ityp], etarange[iet], vartitle[ij]);
              h_recofake_2D[ityp][iet][ij] = binsRec2D[ityp][iet][ij]->CreateHistogram(name,false,0,title);//false : global bin ID
              h_recofake_2D[ityp][iet][ij]->Sumw2();             
	     }
#ifdef  LHAPDF
            for (int ix=1; ix<nnnmx; ix++) {
             sprintf(name, "dd_genpdf_typ_%i_eta%i_%i_%i", ityp, iet, ij, ix);
             sprintf(title, "2D Genpdf %s %g %s %i", typname[ityp], etarange[iet], vartitle[ij], ix);
             h_genevtvarpdf_2D[ityp][iet][ij][ix] = binsGen2D[ityp][iet][ij]->CreateHistogram(name,false,0,title);
             h_genevtvarpdf_2D[ityp][iet][ij][ix]->Sumw2();
            }
#endif
#ifdef  JETENERGY
            for (int ix=1; ix<njecmx; ix++) {
              sprintf(name, "dd_recojec_typ_%i_eta%i_%i_%i", ityp, iet, ij, ix);
              sprintf(title, "2D Recojec %s %g %s %i", typname[ityp], etarange[iet], vartitle[ij], ix);
              h_recoevtvarjec_2D[ityp][iet][ij][ix] = binsRec2D[ityp][iet][ij]->CreateHistogram(name,false,0,title);
              h_recoevtvarjec_2D[ityp][iet][ij][ix]->Sumw2();
            }
#elif defined(JETRESO)
            for (int ix=1; ix<njecmx; ix++ ) {
              sprintf(name, "dd_recoreso_typ_%i_eta%i_%i_%i", ityp, iet, ij, ix);
              sprintf(title, "2D Recoreso %s %g %s %i", typname[ityp], etarange[iet], vartitle[ij], ix);
              h_recoevtvarres_2D[ityp][iet][ij][ix] = binsRec2D[ityp][iet][ij]->CreateHistogram(name,false,0,title);
              h_recoevtvarres_2D[ityp][iet][ij][ix]->Sumw2();
            }
#endif
	      sprintf(name, "dd_gen_typ_%i_eta%i_%i", ityp, iet, ij);
              sprintf(title, "2D Gen %s %g %s", typname[ityp] , etarange[iet], vartitle[ij]);
              h_genvar_2D[ityp][iet][ij] = binsGen2D[ityp][iet][ij]->CreateHistogram(name,false,0,title);//false : global bin ID
              h_genvar_2D[ityp][iet][ij]->Sumw2();              

              sprintf(title, "Miss Gen %s %g %s", typname[ityp], etarange[iet], vartitle[ij]);
	      h_genmiss_2D[ityp][iet][ij] = binsGen2D[ityp][iet][ij]->CreateHistogram(name,false,0,title); //false : global bin ID
              h_genmiss_2D[ityp][iet][ij]->Sumw2();              

	      if (isReconstruct) {
              sprintf(name, "dd_corr_typ_%i_eta%i_%i", ityp , iet, ij);
              sprintf(title, "Gen_Reco %s %g %s", typname[ityp], etarange[iet], vartitle[ij]);
	      RM_2D[ityp][iet][ij] = TUnfoldBinning::CreateHistogramOfMigrations(binsRec2D[ityp][iet][ij], binsGen2D[ityp][iet][ij], name ,0,0, title);
              RM_2D[ityp][iet][ij]->Sumw2();              
             
	      } //if (isReconstruct)
	    }//if (isItUsed(ij))
          }
        }
      }
//-------------------------------------------Define TUnfoldBinning Jet Charge
for (int iet=0; iet<njetetamn; iet++) {                         //Default 1
	for (int ik=0; ik<nkappa; ik++){
		//Reco
                sprintf(RecoBinName, "Detector2d_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);	
		sprintf(RecoBinName, "Recobin2d_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_D1J1[ik][iet]= binsRec2D_D1J1[ik][iet]->AddBinning(RecoBinName);
		sprintf(Axisname, "d0_j0_k%i", ik);
                RecoBinning2D_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
		sprintf(Axisname, "pt");
                RecoBinning2D_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		//RecoBinning2D_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		sprintf(RecoBinName, "Detector2d_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_D1J2[ik][iet]= binsRec2D_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "d0_j1_k%i", ik);
                RecoBinning2D_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
                RecoBinning2D_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		//RecoBinning2D_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);
		
		sprintf(RecoBinName, "Detector2d_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_D2J1[ik][iet]= binsRec2D_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "d1_j0_k%i", ik);
                RecoBinning2D_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
                RecoBinning2D_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_D2J2[ik][iet]= binsRec2D_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "d1_j1_k%i", ik);
                RecoBinning2D_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
                RecoBinning2D_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		//RecoBinning2D_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		sprintf(RecoBinName, "Detector2d_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_D3J1[ik][iet]= binsRec2D_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "d2_j0_k%i", ik);
                RecoBinning2D_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
                RecoBinning2D_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_D3J2[ik][iet]= binsRec2D_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "d2_j1_k%i", ik);
                RecoBinning2D_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
                RecoBinning2D_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		// gluon jets
		sprintf(RecoBinName, "Detector2d_gjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_gjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_gjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_gjt_D1J1[ik][iet]= binsRec2D_gjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "gjt_d0_j0_k%i", ik);
                RecoBinning2D_gjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_gjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_gjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		sprintf(RecoBinName, "Detector2d_gjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_gjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_gjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_gjt_D1J2[ik][iet]= binsRec2D_gjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "gjt_d0_j1_k%i", ik);
                RecoBinning2D_gjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_gjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_gjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		sprintf(RecoBinName, "Detector2d_gjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_gjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_gjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_gjt_D2J1[ik][iet]= binsRec2D_gjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "gjt_d1_j0_k%i", ik);
                RecoBinning2D_gjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_gjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_gjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);
          
		sprintf(RecoBinName, "Detector2d_gjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_gjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_gjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_gjt_D2J2[ik][iet]= binsRec2D_gjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "gjt_d1_j1_k%i", ik);
                RecoBinning2D_gjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_gjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_gjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);
                
		sprintf(RecoBinName, "Detector2d_gjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_gjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_gjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_gjt_D3J1[ik][iet]= binsRec2D_gjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "gjt_d2_j0_k%i", ik);
                RecoBinning2D_gjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_gjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_gjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		sprintf(RecoBinName, "Detector2d_gjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_gjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_gjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_gjt_D3J2[ik][iet]= binsRec2D_gjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "gjt_d2_j1_k%i", ik);
                RecoBinning2D_gjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_gjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_gjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

      		//
      		sprintf(RecoBinName, "Detector2d_qgjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_qgjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qgjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qgjt_D1J1[ik][iet]= binsRec2D_qgjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qgjt_d0_j0_k%i", ik);
                RecoBinning2D_qgjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qgjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qgjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qgjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_qgjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qgjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qgjt_D1J2[ik][iet]= binsRec2D_qgjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qgjt_d0_j1_k%i", ik);
                RecoBinning2D_qgjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qgjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qgjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qgjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_qgjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qgjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qgjt_D2J1[ik][iet]= binsRec2D_qgjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qgjt_d1_j0_k%i", ik);
                RecoBinning2D_qgjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qgjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qgjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qgjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_qgjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qgjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qgjt_D2J2[ik][iet]= binsRec2D_qgjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qgjt_d1_j1_k%i", ik);
                RecoBinning2D_qgjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qgjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qgjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qgjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_qgjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qgjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qgjt_D3J1[ik][iet]= binsRec2D_qgjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qgjt_d2_j0_k%i", ik);
                RecoBinning2D_qgjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qgjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qgjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qgjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_qgjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qgjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qgjt_D3J2[ik][iet]= binsRec2D_qgjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qgjt_d2_j1_k%i", ik);
                RecoBinning2D_qgjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qgjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qgjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		//
		sprintf(RecoBinName, "Detector2d_agjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_agjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_agjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_agjt_D1J1[ik][iet]= binsRec2D_agjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "agjt_d0_j0_k%i", ik);
                RecoBinning2D_agjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_agjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_agjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_agjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_agjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_agjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_agjt_D1J2[ik][iet]= binsRec2D_agjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "agjt_d0_j1_k%i", ik);
                RecoBinning2D_agjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_agjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_agjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_agjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_agjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_agjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_agjt_D2J1[ik][iet]= binsRec2D_agjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "agjt_d1_j0_k%i", ik);
                RecoBinning2D_agjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_agjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_agjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_agjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_agjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_agjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_agjt_D2J2[ik][iet]= binsRec2D_agjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "agjt_d1_j1_k%i", ik);
                RecoBinning2D_agjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_agjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_agjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_agjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_agjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_agjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_agjt_D3J1[ik][iet]= binsRec2D_agjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "agjt_d2_j0_k%i", ik);
                RecoBinning2D_agjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_agjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_agjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_agjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_agjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_agjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_agjt_D3J2[ik][iet]= binsRec2D_agjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "agjt_d2_j1_k%i", ik);
                RecoBinning2D_agjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_agjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_agjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		
		// b jets
		sprintf(RecoBinName, "Detector2d_bjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_bjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_bjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_bjt_D1J1[ik][iet]= binsRec2D_bjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "bjt_d0_j0_k%i", ik);
                RecoBinning2D_bjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_bjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_bjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_bjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_bjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_bjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_bjt_D1J2[ik][iet]= binsRec2D_bjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "bjt_d0_j1_k%i", ik);
                RecoBinning2D_bjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_bjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_bjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_bjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_bjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_bjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_bjt_D2J1[ik][iet]= binsRec2D_bjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "bjt_d1_j0_k%i", ik);
                RecoBinning2D_bjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_bjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_bjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_bjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_bjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_bjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_bjt_D2J2[ik][iet]= binsRec2D_bjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "bjt_d1_j1_k%i", ik);
                RecoBinning2D_bjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_bjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_bjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_bjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_bjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_bjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_bjt_D3J1[ik][iet]= binsRec2D_bjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "bjt_d2_j0_k%i", ik);
                RecoBinning2D_bjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_bjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_bjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_bjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_bjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_bjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_bjt_D3J2[ik][iet]= binsRec2D_bjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "bjt_d2_j1_k%i", ik);
                RecoBinning2D_bjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_bjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_bjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		//
		sprintf(RecoBinName, "Detector2d_qbjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_qbjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qbjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qbjt_D1J1[ik][iet]= binsRec2D_qbjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qbjt_d0_j0_k%i", ik);
                RecoBinning2D_qbjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qbjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qbjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qbjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_qbjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qbjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qbjt_D1J2[ik][iet]= binsRec2D_qbjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qbjt_d0_j1_k%i", ik);
                RecoBinning2D_qbjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qbjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qbjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qbjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_qbjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qbjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qbjt_D2J1[ik][iet]= binsRec2D_qbjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qbjt_d1_j0_k%i", ik);
                RecoBinning2D_qbjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qbjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qbjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qbjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_qbjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qbjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qbjt_D2J2[ik][iet]= binsRec2D_qbjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qbjt_d1_j1_k%i", ik);
                RecoBinning2D_qbjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qbjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qbjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qbjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_qbjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qbjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qbjt_D3J1[ik][iet]= binsRec2D_qbjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qbjt_d2_j0_k%i", ik);
                RecoBinning2D_qbjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qbjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qbjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qbjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_qbjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qbjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qbjt_D3J2[ik][iet]= binsRec2D_qbjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "bjt_d2_j1_k%i", ik);
                RecoBinning2D_qbjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qbjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qbjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		//
		sprintf(RecoBinName, "Detector2d_abjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_abjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_abjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_abjt_D1J1[ik][iet]= binsRec2D_abjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "abjt_d0_j0_k%i", ik);
                RecoBinning2D_abjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_abjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_abjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_abjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_abjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_abjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_abjt_D1J2[ik][iet]= binsRec2D_abjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "abjt_d0_j1_k%i", ik);
                RecoBinning2D_abjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_abjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_abjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_abjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_abjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_abjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_abjt_D2J1[ik][iet]= binsRec2D_abjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "abjt_d1_j0_k%i", ik);
                RecoBinning2D_abjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_abjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_abjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_abjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_abjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_abjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_abjt_D2J2[ik][iet]= binsRec2D_abjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "abjt_d1_j1_k%i", ik);
                RecoBinning2D_abjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_abjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_abjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_abjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_abjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_abjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_abjt_D3J1[ik][iet]= binsRec2D_abjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "abjt_d2_j0_k%i", ik);
                RecoBinning2D_abjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_abjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_abjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_abjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_abjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_abjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_abjt_D3J2[ik][iet]= binsRec2D_abjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "abjt_d2_j1_k%i", ik);
                RecoBinning2D_abjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_abjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_abjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		
		// c jets
		sprintf(RecoBinName, "Detector2d_cjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_cjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_cjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_cjt_D1J1[ik][iet]= binsRec2D_cjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "cjt_d0_j0_k%i", ik);
                RecoBinning2D_cjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_cjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_cjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_cjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_cjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_cjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_cjt_D1J2[ik][iet]= binsRec2D_cjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "cjt_d0_j1_k%i", ik);
                RecoBinning2D_cjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_cjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_cjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_cjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_cjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_cjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_cjt_D2J1[ik][iet]= binsRec2D_cjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "cjt_d1_j0_k%i", ik);
                RecoBinning2D_cjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_cjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_cjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_cjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_cjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_cjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_cjt_D2J2[ik][iet]= binsRec2D_cjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "cjt_d1_j1_k%i", ik);
                RecoBinning2D_cjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_cjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_cjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_cjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_cjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_cjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_cjt_D3J1[ik][iet]= binsRec2D_cjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "cjt_d2_j0_k%i", ik);
                RecoBinning2D_cjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_cjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_cjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_cjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_cjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_cjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_cjt_D3J2[ik][iet]= binsRec2D_cjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "cjt_d2_j1_k%i", ik);
                RecoBinning2D_cjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_cjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_cjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		//
		sprintf(RecoBinName, "Detector2d_qcjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_qcjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qcjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qcjt_D1J1[ik][iet]= binsRec2D_qcjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qcjt_d0_j0_k%i", ik);
                RecoBinning2D_qcjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qcjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qcjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qcjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_qcjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qcjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qcjt_D1J2[ik][iet]= binsRec2D_qcjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qcjt_d0_j1_k%i", ik);
                RecoBinning2D_qcjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qcjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qcjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qcjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_qcjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qcjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qcjt_D2J1[ik][iet]= binsRec2D_qcjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qcjt_d1_j0_k%i", ik);
                RecoBinning2D_qcjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qcjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qcjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qcjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_qcjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qcjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qcjt_D2J2[ik][iet]= binsRec2D_qcjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qcjt_d1_j1_k%i", ik);
                RecoBinning2D_qcjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qcjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qcjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qcjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_qcjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qcjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qcjt_D3J1[ik][iet]= binsRec2D_qcjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qcjt_d2_j0_k%i", ik);
                RecoBinning2D_qcjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qcjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qcjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qcjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_qcjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qcjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qcjt_D3J2[ik][iet]= binsRec2D_qcjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qcjt_d2_j1_k%i", ik);
                RecoBinning2D_qcjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qcjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qcjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		//
		sprintf(RecoBinName, "Detector2d_acjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_acjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_acjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_acjt_D1J1[ik][iet]= binsRec2D_acjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "acjt_d0_j0_k%i", ik);
                RecoBinning2D_acjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_acjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_acjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_acjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_acjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_acjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_acjt_D1J2[ik][iet]= binsRec2D_acjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "acjt_d0_j1_k%i", ik);
                RecoBinning2D_acjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_acjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_acjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_acjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_acjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_acjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_acjt_D2J1[ik][iet]= binsRec2D_acjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "acjt_d1_j0_k%i", ik);
                RecoBinning2D_acjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_acjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_acjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_acjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_acjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_acjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_acjt_D2J2[ik][iet]= binsRec2D_acjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "acjt_d1_j1_k%i", ik);
                RecoBinning2D_acjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_acjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_acjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_acjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_acjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_acjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_acjt_D3J1[ik][iet]= binsRec2D_acjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "acjt_d2_j0_k%i", ik);
                RecoBinning2D_acjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_acjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_acjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_acjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_acjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_acjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_acjt_D3J2[ik][iet]= binsRec2D_acjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "acjt_d2_j1_k%i", ik);
                RecoBinning2D_acjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_acjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_acjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);


		// s jets
		sprintf(RecoBinName, "Detector2d_sjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_sjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_sjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_sjt_D1J1[ik][iet]= binsRec2D_sjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "sjt_d0_j0_k%i", ik);
                RecoBinning2D_sjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_sjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_sjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_sjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_sjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_sjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_sjt_D1J2[ik][iet]= binsRec2D_sjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "sjt_d0_j1_k%i", ik);
                RecoBinning2D_sjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_sjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_sjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_sjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_sjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_sjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_sjt_D2J1[ik][iet]= binsRec2D_sjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "sjt_d1_j0_k%i", ik);
                RecoBinning2D_sjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_sjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_sjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_sjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_sjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_sjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_sjt_D2J2[ik][iet]= binsRec2D_sjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "sjt_d1_j1_k%i", ik);
                RecoBinning2D_sjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_sjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_sjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_sjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_sjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_sjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_sjt_D3J1[ik][iet]= binsRec2D_sjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "sjt_d2_j0_k%i", ik);
                RecoBinning2D_sjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_sjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_sjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_sjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_sjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_sjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_sjt_D3J2[ik][iet]= binsRec2D_sjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "sjt_d2_j1_k%i", ik);
                RecoBinning2D_sjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_sjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_sjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		//
		sprintf(RecoBinName, "Detector2d_qsjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_qsjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qsjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qsjt_D1J1[ik][iet]= binsRec2D_qsjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qsjt_d0_j0_k%i", ik);
                RecoBinning2D_qsjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qsjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qsjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qsjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_qsjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qsjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qsjt_D1J2[ik][iet]= binsRec2D_qsjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qsjt_d0_j1_k%i", ik);
                RecoBinning2D_qsjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qsjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qsjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qsjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_qsjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qsjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qsjt_D2J1[ik][iet]= binsRec2D_qsjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qsjt_d1_j0_k%i", ik);
                RecoBinning2D_qsjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qsjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qsjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qsjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_qsjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qsjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qsjt_D2J2[ik][iet]= binsRec2D_qsjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qsjt_d1_j1_k%i", ik);
                RecoBinning2D_qsjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qsjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qsjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qsjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_qsjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qsjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qsjt_D3J1[ik][iet]= binsRec2D_qsjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qsjt_d2_j0_k%i", ik);
                RecoBinning2D_qsjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qsjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qsjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qsjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_qsjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qsjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qsjt_D3J2[ik][iet]= binsRec2D_qsjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qsjt_d2_j1_k%i", ik);
                RecoBinning2D_qsjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qsjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qsjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		//
		sprintf(RecoBinName, "Detector2d_asjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_asjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_asjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_asjt_D1J1[ik][iet]= binsRec2D_asjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "asjt_d0_j0_k%i", ik);
                RecoBinning2D_asjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_asjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_asjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_asjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_asjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_asjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_asjt_D1J2[ik][iet]= binsRec2D_asjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "asjt_d0_j1_k%i", ik);
                RecoBinning2D_asjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_asjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_asjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_asjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_asjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_asjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_asjt_D2J1[ik][iet]= binsRec2D_asjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "asjt_d1_j0_k%i", ik);
                RecoBinning2D_asjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_asjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_asjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_asjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_asjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_asjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_asjt_D2J2[ik][iet]= binsRec2D_asjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "asjt_d1_j1_k%i", ik);
                RecoBinning2D_asjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_asjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_asjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_asjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_asjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_asjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_asjt_D3J1[ik][iet]= binsRec2D_asjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "asjt_d2_j0_k%i", ik);
                RecoBinning2D_asjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_asjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_asjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_asjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_asjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_asjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_asjt_D3J2[ik][iet]= binsRec2D_asjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "asjt_d2_j1_k%i", ik);
                RecoBinning2D_asjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_asjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_asjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		// d jets
		sprintf(RecoBinName, "Detector2d_djt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_djt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_djt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_djt_D1J1[ik][iet]= binsRec2D_djt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "djt_d0_j0_k%i", ik);
                RecoBinning2D_djt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_djt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_djt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_djt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_djt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_djt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_djt_D1J2[ik][iet]= binsRec2D_djt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "djt_d0_j1_k%i", ik);
                RecoBinning2D_djt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_djt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_djt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_djt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_djt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_djt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_djt_D2J1[ik][iet]= binsRec2D_djt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "djt_d1_j0_k%i", ik);
                RecoBinning2D_djt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_djt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_djt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_djt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_djt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_djt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_djt_D2J2[ik][iet]= binsRec2D_djt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "djt_d1_j1_k%i", ik);
                RecoBinning2D_djt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_djt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_djt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_djt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_djt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_djt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_djt_D3J1[ik][iet]= binsRec2D_djt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "djt_d2_j0_k%i", ik);
                RecoBinning2D_djt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_djt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_djt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_djt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_djt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_djt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_djt_D3J2[ik][iet]= binsRec2D_djt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "djt_d2_j1_k%i", ik);
                RecoBinning2D_djt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_djt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_djt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		//
		sprintf(RecoBinName, "Detector2d_qdjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_qdjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qdjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qdjt_D1J1[ik][iet]= binsRec2D_qdjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qdjt_d0_j0_k%i", ik);
                RecoBinning2D_qdjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qdjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qdjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qdjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_qdjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qdjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qdjt_D1J2[ik][iet]= binsRec2D_qdjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qdjt_d0_j1_k%i", ik);
                RecoBinning2D_qdjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qdjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qdjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qdjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_qdjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qdjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qdjt_D2J1[ik][iet]= binsRec2D_qdjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qdjt_d1_j0_k%i", ik);
                RecoBinning2D_qdjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qdjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qdjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qdjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_qdjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qdjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qdjt_D2J2[ik][iet]= binsRec2D_qdjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qdjt_d1_j1_k%i", ik);
                RecoBinning2D_qdjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qdjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qdjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qdjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_qdjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qdjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qdjt_D3J1[ik][iet]= binsRec2D_qdjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qdjt_d2_j0_k%i", ik);
                RecoBinning2D_qdjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qdjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qdjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qdjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_qdjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qdjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qdjt_D3J2[ik][iet]= binsRec2D_qdjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qdjt_d2_j1_k%i", ik);
                RecoBinning2D_qdjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qdjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qdjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		//
		sprintf(RecoBinName, "Detector2d_adjt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_adjt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_adjt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_adjt_D1J1[ik][iet]= binsRec2D_adjt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "adjt_d0_j0_k%i", ik);
                RecoBinning2D_adjt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_adjt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_adjt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_adjt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_adjt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_adjt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_adjt_D1J2[ik][iet]= binsRec2D_adjt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "adjt_d0_j1_k%i", ik);
                RecoBinning2D_adjt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_adjt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_adjt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_adjt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_adjt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_adjt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_adjt_D2J1[ik][iet]= binsRec2D_adjt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "adjt_d1_j0_k%i", ik);
                RecoBinning2D_adjt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_adjt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_adjt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_adjt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_adjt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_adjt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_adjt_D2J2[ik][iet]= binsRec2D_adjt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "adjt_d1_j1_k%i", ik);
                RecoBinning2D_adjt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_adjt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_adjt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_adjt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_adjt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_adjt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_adjt_D3J1[ik][iet]= binsRec2D_adjt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "adjt_d2_j0_k%i", ik);
                RecoBinning2D_adjt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_adjt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_adjt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_adjt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_adjt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_adjt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_adjt_D3J2[ik][iet]= binsRec2D_adjt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "adjt_d2_j1_k%i", ik);
                RecoBinning2D_adjt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_adjt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_adjt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		
		// u jets
		sprintf(RecoBinName, "Detector2d_ujt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_ujt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_ujt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_ujt_D1J1[ik][iet]= binsRec2D_ujt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "ujt_d0_j0_k%i", ik);
                RecoBinning2D_ujt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_ujt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_ujt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_ujt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_ujt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_ujt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_ujt_D1J2[ik][iet]= binsRec2D_ujt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "ujt_d0_j1_k%i", ik);
                RecoBinning2D_ujt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_ujt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_ujt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_ujt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_ujt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_ujt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_ujt_D2J1[ik][iet]= binsRec2D_ujt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "ujt_d1_j0_k%i", ik);
                RecoBinning2D_ujt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_ujt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_ujt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_ujt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_ujt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_ujt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_ujt_D2J2[ik][iet]= binsRec2D_ujt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "ujt_d1_j1_k%i", ik);
                RecoBinning2D_ujt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_ujt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_ujt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_ujt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_ujt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_ujt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_ujt_D3J1[ik][iet]= binsRec2D_ujt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "ujt_d2_j0_k%i", ik);
                RecoBinning2D_ujt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_ujt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_ujt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_ujt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_ujt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_ujt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_ujt_D3J2[ik][iet]= binsRec2D_ujt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "ujt_d2_j1_k%i", ik);
                RecoBinning2D_ujt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_ujt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_ujt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		//
		sprintf(RecoBinName, "Detector2d_qujt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_qujt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qujt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qujt_D1J1[ik][iet]= binsRec2D_qujt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qujt_d0_j0_k%i", ik);
                RecoBinning2D_qujt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qujt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qujt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qujt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_qujt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qujt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qujt_D1J2[ik][iet]= binsRec2D_qujt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qujt_d0_j1_k%i", ik);
                RecoBinning2D_qujt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qujt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qujt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qujt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_qujt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qujt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qujt_D2J1[ik][iet]= binsRec2D_qujt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qujt_d1_j0_k%i", ik);
                RecoBinning2D_qujt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qujt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qujt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qujt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_qujt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qujt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qujt_D2J2[ik][iet]= binsRec2D_qujt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qujt_d1_j1_k%i", ik);
                RecoBinning2D_qujt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qujt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qujt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qujt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_qujt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qujt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_qujt_D3J1[ik][iet]= binsRec2D_qujt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qujt_d2_j0_k%i", ik);
                RecoBinning2D_qujt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qujt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qujt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_qujt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_qujt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_qujt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_qujt_D3J2[ik][iet]= binsRec2D_qujt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "qujt_d2_j1_k%i", ik);
                RecoBinning2D_qujt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_qujt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_qujt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		//
		sprintf(RecoBinName, "Detector2d_aujt_d0_j0_k%i_eta%i", ik, iet);
                binsRec2D_aujt_D1J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_aujt_d0_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_aujt_D1J1[ik][iet]= binsRec2D_aujt_D1J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "aujt_d0_j0_k%i", ik);
                RecoBinning2D_aujt_D1J1[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_aujt_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_aujt_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_aujt_d0_j1_k%i_eta%i", ik, iet);
                binsRec2D_aujt_D1J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_aujt_d0_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_aujt_D1J2[ik][iet]= binsRec2D_aujt_D1J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "aujt_d0_j1_k%i", ik);
                RecoBinning2D_aujt_D1J2[ik][iet]->AddAxis(Axisname, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_aujt_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_aujt_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_aujt_d1_j0_k%i_eta%i", ik, iet);
                binsRec2D_aujt_D2J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_aujt_d1_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_aujt_D2J1[ik][iet]= binsRec2D_aujt_D2J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "aujt_d1_j0_k%i", ik);
                RecoBinning2D_aujt_D2J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_aujt_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_aujt_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_aujt_d1_j1_k%i_eta%i", ik, iet);
                binsRec2D_aujt_D2J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_aujt_d1_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_aujt_D2J2[ik][iet]= binsRec2D_aujt_D2J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "aujt_d1_j1_k%i", ik);
                RecoBinning2D_aujt_D2J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_aujt_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_aujt_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_aujt_d2_j0_k%i_eta%i", ik, iet);
                binsRec2D_aujt_D3J1[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_aujt_d2_j0_k%i_eta%i",ik, iet);
                RecoBinning2D_aujt_D3J1[ik][iet]= binsRec2D_aujt_D3J1[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "aujt_d2_j0_k%i", ik);
                RecoBinning2D_aujt_D3J1[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_aujt_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_aujt_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(RecoBinName, "Detector2d_aujt_d2_j1_k%i_eta%i", ik, iet);
                binsRec2D_aujt_D3J2[ik][iet] = new TUnfoldBinning(RecoBinName);
                sprintf(RecoBinName, "Recobin2d_aujt_d2_j1_k%i_eta%i",ik, iet);
                RecoBinning2D_aujt_D3J2[ik][iet]= binsRec2D_aujt_D3J2[ik][iet]->AddBinning(RecoBinName);
                sprintf(Axisname, "aujt_d2_j1_k%i", ik);
                RecoBinning2D_aujt_D3J2[ik][iet]->AddAxis(Axisname, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
		RecoBinning2D_aujt_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //RecoBinning2D_aujt_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);


		//Gen
		sprintf(GenBinName, "Generator2d_d0_j0_k%i_eta%i", ik, iet);
                binsGen2D_D1J1[ik][iet] = new TUnfoldBinning(GenBinName);
                sprintf(GenBinName, "Genbin2d_d0_j0_k%i_eta%i",ik, iet);
                GenBinning2D_D1J1[ik][iet]= binsGen2D_D1J1[ik][iet]->AddBinning(GenBinName);
                sprintf(Axisname, "d0_j0_k%i", ik);
                GenBinning2D_D1J1[ik][iet]->AddAxis(Axisname, genjcd1bins[ik], genjcd1minran[ik], genjcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
                GenBinning2D_D1J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		//GenBinning2D_D1J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(GenBinName, "Generator2d_d0_j1_k%i_eta%i", ik, iet);
                binsGen2D_D1J2[ik][iet] = new TUnfoldBinning(GenBinName);
                sprintf(GenBinName, "Genbin2d_d0_j1_k%i_eta%i",ik, iet);
                GenBinning2D_D1J2[ik][iet]= binsGen2D_D1J2[ik][iet]->AddBinning(GenBinName);
                sprintf(Axisname, "d0_j1_k%i", ik);
                GenBinning2D_D1J2[ik][iet]->AddAxis(Axisname, genjcd1bins[ik], genjcd1minran[ik], genjcd1maxran[ik], false, false);
                sprintf(Axisname, "pt");
                GenBinning2D_D1J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
                //GenBinning2D_D1J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		sprintf(GenBinName, "Generator2d_d1_j0_k%i_eta%i", ik, iet);
                binsGen2D_D2J1[ik][iet] = new TUnfoldBinning(GenBinName);
                sprintf(GenBinName, "Genbin2d_d1_j0_k%i_eta%i",ik, iet);
                GenBinning2D_D2J1[ik][iet]= binsGen2D_D2J1[ik][iet]->AddBinning(GenBinName);
                sprintf(Axisname, "d1_j0_k%i", ik);
                GenBinning2D_D2J1[ik][iet]->AddAxis(Axisname, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
                GenBinning2D_D2J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		//GenBinning2D_D2J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

                sprintf(GenBinName, "Generator2d_d1_j1_k%i_eta%i", ik, iet);
                binsGen2D_D2J2[ik][iet] = new TUnfoldBinning(GenBinName);
                sprintf(GenBinName, "Genbin2d_d1_j1_k%i_eta%i",ik, iet);
                GenBinning2D_D2J2[ik][iet]= binsGen2D_D2J2[ik][iet]->AddBinning(GenBinName);
                sprintf(Axisname, "d1_j1_k%i", ik);
                GenBinning2D_D2J2[ik][iet]->AddAxis(Axisname, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
                GenBinning2D_D2J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		//GenBinning2D_D2J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);

		sprintf(GenBinName, "Generator2d_d2_j0_k%i_eta%i", ik, iet);
                binsGen2D_D3J1[ik][iet] = new TUnfoldBinning(GenBinName);
                sprintf(GenBinName, "Genbin2d_d2_j0_k%i_eta%i",ik, iet);
                GenBinning2D_D3J1[ik][iet]= binsGen2D_D3J1[ik][iet]->AddBinning(GenBinName);
                sprintf(Axisname, "d2_j0_k%i", ik);
                GenBinning2D_D3J1[ik][iet]->AddAxis(Axisname, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
                GenBinning2D_D3J1[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		//GenBinning2D_D3J1[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);
        
                sprintf(GenBinName, "Generator2d_d2_j1_k%i_eta%i", ik, iet);
                binsGen2D_D3J2[ik][iet] = new TUnfoldBinning(GenBinName);
                sprintf(GenBinName, "Genbin2d_d2_j1_k%i_eta%i",ik, iet);
                GenBinning2D_D3J2[ik][iet]= binsGen2D_D3J2[ik][iet]->AddBinning(GenBinName);
                sprintf(Axisname, "d2_j1_k%i", ik);
                GenBinning2D_D3J2[ik][iet]->AddAxis(Axisname, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik], false, false);
                sprintf(Axisname, "pt");
                GenBinning2D_D3J2[ik][iet]->AddAxis(Axisname, nHLTmx, recohtbins, false, false);
		//GenBinning2D_D3J2[ik][iet]->AddAxis(Axisname, recohtnbins[ik], recohtbinsmin[ik],recohtbinsmax[ik], false, false);
		}
}
//-------------------------------------------
for (int iet=0; iet<njetetamn; iet++) {
	for (int ik=0; ik<10; ik++){
                if (isReconstruct) {
			//Reco
			sprintf(name, "dd_reco_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_2D_D1J1[ik][iet] = binsRec2D_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_2D_D1J1[ik][iet]->Sumw2();

			sprintf(name, "dd_reco_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_2D_D1J2[ik][iet] = binsRec2D_D1J2[ik][iet]->CreateHistogram(name,false,0,title); 
                        h_recovar_2D_D1J2[ik][iet]->Sumw2();

			sprintf(name, "dd_reco_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_2D_D2J1[ik][iet] = binsRec2D_D2J1[ik][iet]->CreateHistogram(name,false,0,title); 
                        h_recovar_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_2D_D2J2[ik][iet] = binsRec2D_D2J2[ik][iet]->CreateHistogram(name,false,0,title); 
                        h_recovar_2D_D2J2[ik][iet]->Sumw2();

			sprintf(name, "dd_reco_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_2D_D3J1[ik][iet] = binsRec2D_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_2D_D3J2[ik][iet] = binsRec2D_D3J2[ik][iet]->CreateHistogram(name,false,0,title); 
                        h_recovar_2D_D3J2[ik][iet]->Sumw2();

                        //Gluon Jet
                        sprintf(name, "dd_reco_gjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco gjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_gjt_2D_D1J1[ik][iet] = binsRec2D_gjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_gjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_gjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco gJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_gjt_2D_D1J2[ik][iet] = binsRec2D_gjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_gjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_gjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco gJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_gjt_2D_D2J1[ik][iet] = binsRec2D_gjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_gjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_gjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco gJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_gjt_2D_D2J2[ik][iet] = binsRec2D_gjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_gjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_gjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco gJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_gjt_2D_D3J1[ik][iet] = binsRec2D_gjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_gjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_gjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco gJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_gjt_2D_D3J2[ik][iet] = binsRec2D_gjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_gjt_2D_D3J2[ik][iet]->Sumw2();

			//
			sprintf(name, "dd_reco_qgjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qgjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qgjt_2D_D1J1[ik][iet] = binsRec2D_qgjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_qgjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qgjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qgJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qgjt_2D_D1J2[ik][iet] = binsRec2D_qgjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qgjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qgjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qgJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qgjt_2D_D2J1[ik][iet] = binsRec2D_qgjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qgjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qgjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qgJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qgjt_2D_D2J2[ik][iet] = binsRec2D_qgjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qgjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qgjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qgJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qgjt_2D_D3J1[ik][iet] = binsRec2D_qgjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qgjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qgjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qgJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qgjt_2D_D3J2[ik][iet] = binsRec2D_qgjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qgjt_2D_D3J2[ik][iet]->Sumw2();

			//
			sprintf(name, "dd_reco_agjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco agjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_agjt_2D_D1J1[ik][iet] = binsRec2D_agjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_agjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_agjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco agJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_agjt_2D_D1J2[ik][iet] = binsRec2D_agjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_agjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_agjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco agJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_agjt_2D_D2J1[ik][iet] = binsRec2D_agjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_agjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_agjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco agJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_agjt_2D_D2J2[ik][iet] = binsRec2D_agjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_agjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_agjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco agJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_agjt_2D_D3J1[ik][iet] = binsRec2D_agjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_agjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_agjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco agJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_agjt_2D_D3J2[ik][iet] = binsRec2D_agjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_agjt_2D_D3J2[ik][iet]->Sumw2();
	
                        // b jets
                        sprintf(name, "dd_reco_bjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco bjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_bjt_2D_D1J1[ik][iet] = binsRec2D_bjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_bjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_bjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco bJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_bjt_2D_D1J2[ik][iet] = binsRec2D_bjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_bjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_bjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco bJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_bjt_2D_D2J1[ik][iet] = binsRec2D_bjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_bjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_bjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco bJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_bjt_2D_D2J2[ik][iet] = binsRec2D_bjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_bjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_bjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco bJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_bjt_2D_D3J1[ik][iet] = binsRec2D_bjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_bjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_bjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco bJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_bjt_2D_D3J2[ik][iet] = binsRec2D_bjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_bjt_2D_D3J2[ik][iet]->Sumw2();

			//
			sprintf(name, "dd_reco_qbjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qbjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qbjt_2D_D1J1[ik][iet] = binsRec2D_qbjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_qbjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qbjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qbJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qbjt_2D_D1J2[ik][iet] = binsRec2D_qbjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qbjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qbjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qbJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qbjt_2D_D2J1[ik][iet] = binsRec2D_qbjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qbjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qbjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qbJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qbjt_2D_D2J2[ik][iet] = binsRec2D_qbjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qbjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qbjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qbJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qbjt_2D_D3J1[ik][iet] = binsRec2D_qbjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qbjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qbjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qbJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qbjt_2D_D3J2[ik][iet] = binsRec2D_qbjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qbjt_2D_D3J2[ik][iet]->Sumw2();

			//
			sprintf(name, "dd_reco_abjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco abjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_abjt_2D_D1J1[ik][iet] = binsRec2D_abjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_abjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_abjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco abJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_abjt_2D_D1J2[ik][iet] = binsRec2D_abjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_abjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_abjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco abJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_abjt_2D_D2J1[ik][iet] = binsRec2D_abjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_abjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_abjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco abJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_abjt_2D_D2J2[ik][iet] = binsRec2D_abjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_abjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_abjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco abJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_abjt_2D_D3J1[ik][iet] = binsRec2D_abjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_abjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_abjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco abJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_abjt_2D_D3J2[ik][iet] = binsRec2D_abjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_abjt_2D_D3J2[ik][iet]->Sumw2();

			// c jets
			sprintf(name, "dd_reco_cjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco cjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_cjt_2D_D1J1[ik][iet] = binsRec2D_cjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_cjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_cjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco cJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_cjt_2D_D1J2[ik][iet] = binsRec2D_cjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_cjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_cjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco cJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_cjt_2D_D2J1[ik][iet] = binsRec2D_cjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_cjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_cjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco cJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_cjt_2D_D2J2[ik][iet] = binsRec2D_cjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_cjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_cjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco cJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_cjt_2D_D3J1[ik][iet] = binsRec2D_cjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_cjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_cjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco cJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_cjt_2D_D3J2[ik][iet] = binsRec2D_cjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_cjt_2D_D3J2[ik][iet]->Sumw2();

			//
			sprintf(name, "dd_reco_qcjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qcjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qcjt_2D_D1J1[ik][iet] = binsRec2D_qcjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_qcjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qcjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qcJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qcjt_2D_D1J2[ik][iet] = binsRec2D_qcjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qcjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qcjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qcJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qcjt_2D_D2J1[ik][iet] = binsRec2D_qcjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qcjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qcjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qcJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qcjt_2D_D2J2[ik][iet] = binsRec2D_qcjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qcjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qcjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qcJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qcjt_2D_D3J1[ik][iet] = binsRec2D_qcjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qcjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qcjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qcJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qcjt_2D_D3J2[ik][iet] = binsRec2D_qcjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qcjt_2D_D3J2[ik][iet]->Sumw2();

			//
			sprintf(name, "dd_reco_acjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco acjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_acjt_2D_D1J1[ik][iet] = binsRec2D_acjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_acjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_acjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco acJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_acjt_2D_D1J2[ik][iet] = binsRec2D_acjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_acjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_acjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco acJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_acjt_2D_D2J1[ik][iet] = binsRec2D_acjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_acjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_acjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco acJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_acjt_2D_D2J2[ik][iet] = binsRec2D_acjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_acjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_acjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco acJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_acjt_2D_D3J1[ik][iet] = binsRec2D_acjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_acjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_acjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco acJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_acjt_2D_D3J2[ik][iet] = binsRec2D_acjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_acjt_2D_D3J2[ik][iet]->Sumw2();

			// s jets
			sprintf(name, "dd_reco_sjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco sjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_sjt_2D_D1J1[ik][iet] = binsRec2D_sjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_sjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_sjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco sJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_sjt_2D_D1J2[ik][iet] = binsRec2D_sjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_sjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_sjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco sJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_sjt_2D_D2J1[ik][iet] = binsRec2D_sjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_sjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_sjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco sJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_sjt_2D_D2J2[ik][iet] = binsRec2D_sjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_sjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_sjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco sJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_sjt_2D_D3J1[ik][iet] = binsRec2D_sjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_sjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_sjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco sJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_sjt_2D_D3J2[ik][iet] = binsRec2D_sjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_sjt_2D_D3J2[ik][iet]->Sumw2();

			//
			sprintf(name, "dd_reco_qsjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qsjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qsjt_2D_D1J1[ik][iet] = binsRec2D_qsjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_qsjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qsjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qsJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qsjt_2D_D1J2[ik][iet] = binsRec2D_qsjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qsjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qsjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qsJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qsjt_2D_D2J1[ik][iet] = binsRec2D_qsjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qsjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qsjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qsJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qsjt_2D_D2J2[ik][iet] = binsRec2D_qsjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qsjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qsjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qsJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qsjt_2D_D3J1[ik][iet] = binsRec2D_qsjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qsjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qsjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qsJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qsjt_2D_D3J2[ik][iet] = binsRec2D_qsjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qsjt_2D_D3J2[ik][iet]->Sumw2();

			//
			sprintf(name, "dd_reco_asjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco asjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_asjt_2D_D1J1[ik][iet] = binsRec2D_asjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_asjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_asjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco asJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_asjt_2D_D1J2[ik][iet] = binsRec2D_asjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_asjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_asjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco asJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_asjt_2D_D2J1[ik][iet] = binsRec2D_asjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_asjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_asjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco asJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_asjt_2D_D2J2[ik][iet] = binsRec2D_asjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_asjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_asjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco asJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_asjt_2D_D3J1[ik][iet] = binsRec2D_asjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_asjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_asjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco asJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_asjt_2D_D3J2[ik][iet] = binsRec2D_asjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_asjt_2D_D3J2[ik][iet]->Sumw2();

			// d jets
			sprintf(name, "dd_reco_djt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco djet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_djt_2D_D1J1[ik][iet] = binsRec2D_djt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_djt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_djt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco dJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_djt_2D_D1J2[ik][iet] = binsRec2D_djt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_djt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_djt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco dJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_djt_2D_D2J1[ik][iet] = binsRec2D_djt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_djt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_djt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco dJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_djt_2D_D2J2[ik][iet] = binsRec2D_djt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_djt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_djt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco dJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_djt_2D_D3J1[ik][iet] = binsRec2D_djt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_djt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_djt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco dJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_djt_2D_D3J2[ik][iet] = binsRec2D_djt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_djt_2D_D3J2[ik][iet]->Sumw2();

			//
			sprintf(name, "dd_reco_qdjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qdjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qdjt_2D_D1J1[ik][iet] = binsRec2D_qdjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_qdjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qdjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qdJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qdjt_2D_D1J2[ik][iet] = binsRec2D_qdjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qdjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qdjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qdJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qdjt_2D_D2J1[ik][iet] = binsRec2D_qdjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qdjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qdjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qdJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qdjt_2D_D2J2[ik][iet] = binsRec2D_qdjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qdjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qdjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qdJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qdjt_2D_D3J1[ik][iet] = binsRec2D_qdjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qdjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qdjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qdJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qdjt_2D_D3J2[ik][iet] = binsRec2D_qdjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qdjt_2D_D3J2[ik][iet]->Sumw2();

			//
			sprintf(name, "dd_reco_adjt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco adjet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_adjt_2D_D1J1[ik][iet] = binsRec2D_adjt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_adjt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_adjt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco adJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_adjt_2D_D1J2[ik][iet] = binsRec2D_adjt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_adjt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_adjt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco adJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_adjt_2D_D2J1[ik][iet] = binsRec2D_adjt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_adjt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_adjt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco adJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_adjt_2D_D2J2[ik][iet] = binsRec2D_adjt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_adjt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_adjt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco adJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_adjt_2D_D3J1[ik][iet] = binsRec2D_adjt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_adjt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_adjt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco adJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_adjt_2D_D3J2[ik][iet] = binsRec2D_adjt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_adjt_2D_D3J2[ik][iet]->Sumw2();

			// u jets
			sprintf(name, "dd_reco_ujt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco ujet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_ujt_2D_D1J1[ik][iet] = binsRec2D_ujt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_ujt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_ujt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco uJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_ujt_2D_D1J2[ik][iet] = binsRec2D_ujt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_ujt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_ujt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco uJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_ujt_2D_D2J1[ik][iet] = binsRec2D_ujt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_ujt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_ujt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco uJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_ujt_2D_D2J2[ik][iet] = binsRec2D_ujt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_ujt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_ujt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco uJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_ujt_2D_D3J1[ik][iet] = binsRec2D_ujt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_ujt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_ujt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco uJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_ujt_2D_D3J2[ik][iet] = binsRec2D_ujt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_ujt_2D_D3J2[ik][iet]->Sumw2();

			//
			sprintf(name, "dd_reco_qujt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco qujet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qujt_2D_D1J1[ik][iet] = binsRec2D_qujt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_qujt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qujt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco quJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qujt_2D_D1J2[ik][iet] = binsRec2D_qujt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qujt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qujt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco quJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qujt_2D_D2J1[ik][iet] = binsRec2D_qujt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qujt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qujt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco quJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qujt_2D_D2J2[ik][iet] = binsRec2D_qujt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qujt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qujt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco quJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qujt_2D_D3J1[ik][iet] = binsRec2D_qujt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qujt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_qujt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco quJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_qujt_2D_D3J2[ik][iet] = binsRec2D_qujt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_qujt_2D_D3J2[ik][iet]->Sumw2();

			//
			sprintf(name, "dd_reco_aujt_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco aujet Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_aujt_2D_D1J1[ik][iet] = binsRec2D_aujt_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recovar_aujt_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_aujt_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco auJet Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_aujt_2D_D1J2[ik][iet] = binsRec2D_aujt_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_aujt_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_aujt_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco auJet Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_aujt_2D_D2J1[ik][iet] = binsRec2D_aujt_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_aujt_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_aujt_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco auJet Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_aujt_2D_D2J2[ik][iet] = binsRec2D_aujt_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_aujt_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_aujt_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco auJet Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_aujt_2D_D3J1[ik][iet] = binsRec2D_aujt_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_aujt_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_reco_aujt_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco auJet Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recovar_aujt_2D_D3J2[ik][iet] = binsRec2D_aujt_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recovar_aujt_2D_D3J2[ik][iet]->Sumw2();

			//recofake
	     		sprintf(name, "dd_recofake_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco Fake Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recofake_2D_D1J1[ik][iet] = binsRec2D_D1J1[ik][iet]->CreateHistogram(name,false,0,title); //false : global bin ID
                        h_recofake_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_recofake_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco Fake Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recofake_2D_D1J2[ik][iet] = binsRec2D_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recofake_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_recofake_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco Fake Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recofake_2D_D2J1[ik][iet] = binsRec2D_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recofake_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_recofake_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco Fake Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recofake_2D_D2J2[ik][iet] = binsRec2D_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recofake_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_recofake_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco Fake Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recofake_2D_D3J1[ik][iet] = binsRec2D_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recofake_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_recofake_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Reco Fake Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_recofake_2D_D3J2[ik][iet] = binsRec2D_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_recofake_2D_D3J2[ik][iet]->Sumw2();
			}

			//Gen
			sprintf(name, "dd_gen_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Gen Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_genvar_2D_D1J1[ik][iet] = binsGen2D_D1J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_genvar_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_gen_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Gen Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_genvar_2D_D1J2[ik][iet] = binsGen2D_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_genvar_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_gen_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Gen Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_genvar_2D_D2J1[ik][iet] = binsGen2D_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_genvar_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_gen_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Gen Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_genvar_2D_D2J2[ik][iet] = binsGen2D_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_genvar_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_gen_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Gen Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_genvar_2D_D3J1[ik][iet] = binsGen2D_D3J1[ik][iet]->CreateHistogram(name,false,0,title); 
                        h_genvar_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_gen_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Gen Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_genvar_2D_D3J2[ik][iet] = binsGen2D_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_genvar_2D_D3J2[ik][iet]->Sumw2();

			//genmiss
			sprintf(name, "dd_genmiss_d0_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Gen Miss Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
                        h_genmiss_2D_D1J1[ik][iet] = binsGen2D_D1J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_genmiss_2D_D1J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_genmiss_d0_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Gen Miss Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                        h_genmiss_2D_D1J2[ik][iet] = binsGen2D_D1J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_genmiss_2D_D1J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_genmiss_d1_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Gen Miss Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_genmiss_2D_D2J1[ik][iet] = binsGen2D_D2J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_genmiss_2D_D2J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_genmiss_d1_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Gen Miss Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                        h_genmiss_2D_D2J2[ik][iet] = binsGen2D_D2J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_genmiss_2D_D2J2[ik][iet]->Sumw2();

                        sprintf(name, "dd_genmiss_d2_j0_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Gen Miss Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_genmiss_2D_D3J1[ik][iet] = binsGen2D_D3J1[ik][iet]->CreateHistogram(name,false,0,title);
                        h_genmiss_2D_D3J1[ik][iet]->Sumw2();

                        sprintf(name, "dd_genmiss_d2_j1_k%i_eta%i", ik, iet);
                        sprintf(title, "2D Gen Miss Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                        h_genmiss_2D_D3J2[ik][iet] = binsGen2D_D3J2[ik][iet]->CreateHistogram(name,false,0,title);
                        h_genmiss_2D_D3J2[ik][iet]->Sumw2();

			//RM
			if (isReconstruct) {
              			sprintf(name, "dd_corr_d0_j0_k%i_eta%i", ik , iet);
              			sprintf(title, "Gen_Reco Q_{1}^{%g} %g", kappa[ik], etarange[iet]);
              			RM_2D_D1J1[ik][iet] = TUnfoldBinning::CreateHistogramOfMigrations(binsRec2D_D1J1[ik][iet], binsGen2D_D1J1[ik][iet], name ,0,0, title);
              			RM_2D_D1J1[ik][iet]->Sumw2();

				sprintf(name, "dd_corr_d0_j1_k%i_eta%i", ik , iet);
                                sprintf(title, "Gen_Reco Q_{2}^{%g} %g", kappa[ik], etarange[iet]);
                                RM_2D_D1J2[ik][iet] = TUnfoldBinning::CreateHistogramOfMigrations(binsRec2D_D1J2[ik][iet], binsGen2D_D1J2[ik][iet], name ,0,0, title);
                                RM_2D_D1J2[ik][iet]->Sumw2();

				sprintf(name, "dd_corr_d1_j0_k%i_eta%i", ik , iet);
                                sprintf(title, "Gen_Reco Q_{1L}^{%g} %g", kappa[ik], etarange[iet]);
                                RM_2D_D2J1[ik][iet] = TUnfoldBinning::CreateHistogramOfMigrations(binsRec2D_D2J1[ik][iet], binsGen2D_D2J1[ik][iet], name ,0,0, title);
                                RM_2D_D2J1[ik][iet]->Sumw2();

                                sprintf(name, "dd_corr_d1_j1_k%i_eta%i", ik , iet);
                                sprintf(title, "Gen_Reco Q_{2L}^{%g} %g", kappa[ik], etarange[iet]);
                                RM_2D_D2J2[ik][iet] = TUnfoldBinning::CreateHistogramOfMigrations(binsRec2D_D2J2[ik][iet], binsGen2D_D2J2[ik][iet], name ,0,0, title);
                                RM_2D_D2J2[ik][iet]->Sumw2();

				sprintf(name, "dd_corr_d2_j0_k%i_eta%i", ik , iet);
                                sprintf(title, "Gen_Reco Q_{1T}^{%g} %g", kappa[ik], etarange[iet]);
                                RM_2D_D3J1[ik][iet] = TUnfoldBinning::CreateHistogramOfMigrations(binsRec2D_D3J1[ik][iet], binsGen2D_D3J1[ik][iet], name ,0,0, title);
                                RM_2D_D3J1[ik][iet]->Sumw2();

                                sprintf(name, "dd_corr_d2_j1_k%i_eta%i", ik , iet);
                                sprintf(title, "Gen_Reco Q_{2T}^{%g} %g", kappa[ik], etarange[iet]);
                                RM_2D_D3J2[ik][iet] = TUnfoldBinning::CreateHistogramOfMigrations(binsRec2D_D3J2[ik][iet], binsGen2D_D3J2[ik][iet], name ,0,0, title);
                                RM_2D_D3J2[ik][iet]->Sumw2();

		              }
		}
	}
//-------------------------------------------Jet Charge Histogram
for (int ipt=0; ipt<njetptmn; ipt++) {
	for (int iet=0; iet<njetetamn; iet++) {
		//for (int ij=1; ij<11; ij++){
			for (int ik=0; ik<10; ik++){
        		if (isReconstruct) {
				//Reco
				sprintf(name, "reco_jc_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
        	        	sprintf(title, "Reco Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
				h_recojc_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);	
				h_recojc_D1J1[ik][ipt][iet]->Sumw2();

				sprintf(name, "reco_jc_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_D1J2[ik][ipt][iet]->Sumw2();

				sprintf(name, "reco_jc_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_D2J2[ik][ipt][iet]->Sumw2();

				sprintf(name, "reco_jc_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_D3J2[ik][ipt][iet]->Sumw2();

				//Gluon Jets
				sprintf(name, "reco_jc_gjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco gJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_gjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_gjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_gjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco gJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_gjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_gjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_gjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco gJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_gjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_gjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_gjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco gJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_gjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_gjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_gjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco gJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_gjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_gjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_gjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco gJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_gjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_gjt_D3J2[ik][ipt][iet]->Sumw2();

				//
				sprintf(name, "reco_jc_qgjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qgJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qgjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_qgjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qgjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qgJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qgjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_qgjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qgjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qgJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qgjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qgjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qgjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qgJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qgjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qgjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qgjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qgJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qgjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qgjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qgjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qgJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qgjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qgjt_D3J2[ik][ipt][iet]->Sumw2();

				//
				sprintf(name, "reco_jc_agjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco agJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_agjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_agjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_agjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco agJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_agjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_agjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_agjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco agJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_agjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_agjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_agjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco agJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_agjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_agjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_agjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco agJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_agjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_agjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_agjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco agJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_agjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_agjt_D3J2[ik][ipt][iet]->Sumw2();


				//b Jets
				sprintf(name, "reco_jc_bjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco bJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_bjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_bjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_bjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco bJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_bjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_bjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_bjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco bJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_bjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_bjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_bjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco bJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_bjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_bjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_bjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco bJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_bjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_bjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_bjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco bJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_bjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_bjt_D3J2[ik][ipt][iet]->Sumw2();

				//
				sprintf(name, "reco_jc_qbjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qbJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qbjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_qbjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qbjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qbJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qbjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_qbjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qbjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qbJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qbjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qbjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qbjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qbJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qbjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qbjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qbjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qbJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qbjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qbjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qbjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qbJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qbjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qbjt_D3J2[ik][ipt][iet]->Sumw2();

				//
				sprintf(name, "reco_jc_abjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco abJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_abjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_abjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_abjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco abJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_abjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_abjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_abjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco abJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_abjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_abjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_abjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco abJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_abjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_abjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_abjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco abJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_abjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_abjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_abjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco abJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_abjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_abjt_D3J2[ik][ipt][iet]->Sumw2();

				//c Jets
				sprintf(name, "reco_jc_cjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco cJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_cjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_cjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_cjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco cJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_cjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_cjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_cjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco cJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_cjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_cjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_cjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco cJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_cjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_cjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_cjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco cJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_cjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_cjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_cjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco cJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_cjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_cjt_D3J2[ik][ipt][iet]->Sumw2();

				//
				sprintf(name, "reco_jc_qcjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qcJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qcjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_qcjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qcjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qcJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qcjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_qcjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qcjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qcJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qcjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qcjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qcjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qcJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qcjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qcjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qcjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qcJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qcjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qcjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qcjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qcJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qcjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qcjt_D3J2[ik][ipt][iet]->Sumw2();

				//
				sprintf(name, "reco_jc_acjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco acJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_acjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_acjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_acjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco acJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_acjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_acjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_acjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco acJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_acjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_acjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_acjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco acJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_acjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_acjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_acjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco acJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_acjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_acjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_acjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco acJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_acjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_acjt_D3J2[ik][ipt][iet]->Sumw2();


				//s Jets
				sprintf(name, "reco_jc_sjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco sJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_sjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_sjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_sjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco sJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_sjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_sjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_sjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco sJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_sjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_sjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_sjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco sJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_sjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_sjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_sjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco sJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_sjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_sjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_sjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco sJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_sjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_sjt_D3J2[ik][ipt][iet]->Sumw2();

				//
				sprintf(name, "reco_jc_qsjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qsJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qsjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_qsjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qsjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qsJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qsjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_qsjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qsjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qsJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qsjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qsjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qsjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qsJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qsjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qsjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qsjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qsJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qsjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qsjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qsjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qsJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qsjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qsjt_D3J2[ik][ipt][iet]->Sumw2();

				//
				sprintf(name, "reco_jc_asjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco asJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_asjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_asjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_asjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco asJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_asjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_asjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_asjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco asJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_asjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_asjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_asjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco asJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_asjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_asjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_asjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco asJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_asjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_asjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_asjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco asJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_asjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_asjt_D3J2[ik][ipt][iet]->Sumw2();


				//u Jet
				sprintf(name, "reco_jc_ujt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco uJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_ujt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_ujt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_ujt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco uJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_ujt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_ujt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_ujt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco uJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_ujt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_ujt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_ujt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco uJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_ujt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_ujt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_ujt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco uJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_ujt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_ujt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_ujt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco uJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_ujt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_ujt_D3J2[ik][ipt][iet]->Sumw2();

				//
				sprintf(name, "reco_jc_qujt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco quJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qujt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_qujt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qujt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco quJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qujt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_qujt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qujt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco quJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qujt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qujt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qujt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco quJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qujt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qujt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qujt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco quJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qujt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qujt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qujt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco quJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qujt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qujt_D3J2[ik][ipt][iet]->Sumw2();

				//
				sprintf(name, "reco_jc_aujt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco auJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_aujt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_aujt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_aujt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco auJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_aujt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_aujt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_aujt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco auJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_aujt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_aujt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_aujt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco auJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_aujt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_aujt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_aujt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco auJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_aujt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_aujt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_aujt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco auJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_aujt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_aujt_D3J2[ik][ipt][iet]->Sumw2();

				//d Jet
				sprintf(name, "reco_jc_djt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco dJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_djt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_djt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_djt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco dJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_djt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_djt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_djt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco dJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_djt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_djt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_djt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco dJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_djt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_djt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_djt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco dJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_djt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_djt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_djt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco dJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_djt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_djt_D3J2[ik][ipt][iet]->Sumw2();

				//
				sprintf(name, "reco_jc_qdjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qdJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qdjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_qdjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qdjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qdJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qdjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_qdjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qdjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qdJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qdjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qdjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qdjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qdJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qdjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qdjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qdjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qdJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qdjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qdjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_qdjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco qdJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_qdjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_qdjt_D3J2[ik][ipt][iet]->Sumw2();

				//
				sprintf(name, "reco_jc_adjt_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco adJet Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_adjt_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_adjt_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_adjt_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco adJet Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_adjt_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recojc_adjt_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_adjt_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco adJet Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_adjt_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_adjt_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_adjt_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco adJet Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_adjt_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_adjt_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_adjt_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco adJet Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_adjt_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_adjt_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "reco_jc_adjt_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco adJet Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_recojc_adjt_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recojc_adjt_D3J2[ik][ipt][iet]->Sumw2();

				}
				
				//Gen
				sprintf(name, "gen_jc_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Gen Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_genjc_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, genjcd1bins[ik], genjcd1minran[ik], genjcd1maxran[ik]);
                                h_genjc_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "gen_jc_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Gen Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_genjc_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, genjcd1bins[ik], genjcd1minran[ik], genjcd1maxran[ik]);
                                h_genjc_D1J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "gen_jc_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Gen Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_genjc_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik]);
                                h_genjc_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "gen_jc_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Gen Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_genjc_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik]);
                                h_genjc_D2J2[ik][ipt][iet]->Sumw2();

                                sprintf(name, "gen_jc_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Gen Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_genjc_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik]);
                                h_genjc_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "gen_jc_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Gen Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet]);
                                h_genjc_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik]);
                                h_genjc_D3J2[ik][ipt][iet]->Sumw2();

			}	
		}
	}
//}
//-------------------------------------------ESV Root Histogram
  for (int ityp=0; ityp<ntype; ityp++) {
    for (int ipt=0; ipt<njetptmn; ipt++) {
      for (int iet=0; iet<njetetamn; iet++) {
	for (int ij=0; ij<nvar; ij++) {
	  if (isItUsed(ij)) { 
	    if (isReconstruct) { 
	      sprintf(name, "reco_typ_%i_pt%i_eta%i_%i", ityp, ipt, iet, ij);
	      sprintf(title, "Reco %s %i %g %s", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij]);
	      h_recoevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij], (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij] );  
	      //h_recoevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij] : nbinsx1[ij], (ityp==0) ? binrngs0[ij] : binrngs1[ij] );  
              //h_recoevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt] );
              //h_recoevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, nbinsx[ij], -endx[ij], -startx[ij]);
	      h_recoevtvar[ityp][ipt][iet][ij]->Sumw2();
	     //For Fake
	     sprintf(name, "fake_reco_typ_%i_pt%i_eta%i_%i", ityp, ipt, iet, ij);
             sprintf(title, "Fake Reco %s %i %g %s", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij]);
             h_recoevtfake[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij], (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij]);  //Reco : only Var
	     // h_recoevtfake[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt] );//Reco each HT,var
             h_recoevtfake[ityp][ipt][iet][ij]->Sumw2();
	    }
#ifdef  LHAPDF
	    for (int ix=1; ix<nnnmx; ix++) {
	     sprintf(name, "genpdf_typ_%i_pt%i_eta%i_%i_%i", ityp, ipt, iet, ij, ix);
	     sprintf(title, "Genpdf %s %i %g %s %i", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij], ix);
	     h_genevtvarpdf[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title,  (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij] , (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij]);
	     // h_genevtvarpdf[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title,  (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt] , (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]);//Reco each HT, var
	     h_genevtvarpdf[ityp][ipt][iet][ij][ix]->Sumw2();
	    }
#endif
	    
#ifdef  JETENERGY
	    for (int ix=1; ix<njecmx; ix++) {
	      sprintf(name, "recojec_typ_%i_pt%i_eta%i_%i_%i", ityp, ipt, iet, ij, ix);
	      sprintf(title, "Recojec %s %i %g %s %i", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij], ix);
	      h_recoevtvarjec[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title,  (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij] , (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij]);
	      //h_recoevtvarjec[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title,  (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt] , (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]);
	      h_recoevtvarjec[ityp][ipt][iet][ij][ix]->Sumw2();
	    }
#elif defined(JETRESO)
	    for (int ix=1; ix<njecmx; ix++ ) {
	      sprintf(name, "recoreso_typ_%i_pt%i_eta%i_%i_%i", ityp, ipt, iet, ij, ix);
	      sprintf(title, "Recoreso %s %i %g %s %i", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij], ix);
              h_recoevtvarres[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title, (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij] , (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij]);
	      sprintf(title, "Recoreso %s %i %g %s %i", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij], ix);
              h_recoevtvarres[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title, (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij] , (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij]);
              //h_recoevtvarres[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij] : nbinsx1[ij] , (ityp==0) ? binrngs0[ij] : binrngs1[ij]);
              //h_recoevtvarres[ityp][ipt][iet][ij][ix] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt] , (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]);
	      h_recoevtvarres[ityp][ipt][iet][ij][ix]->Sumw2();
	    }
#endif	    
	    sprintf(name, "gen_typ_%i_pt%i_eta%i_%i", ityp, ipt, iet, ij);
	    sprintf(title, "Gen %s %i %g %s", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij]);
            //h_genevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt] );
            h_genevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij] : nbinsx1[ij], (ityp==0) ? binrngs0[ij] : binrngs1[ij] );
            //h_genevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, nbinsxgen[ij], -endx[ij], -startx[ij]);
            //h_genevtvar[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, nbinsx[ij], -endx[ij], -startx[ij]); //For Equal Binning
            h_genevtvar[ityp][ipt][iet][ij]->Sumw2();
            //for miss
            sprintf(name, "miss_gen_typ_%i_pt%i_eta%i_%i", ityp, ipt, iet, ij);
            sprintf(title, "Miss Gen %s %i %g %s", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij]);
            h_genevtmiss[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij] : nbinsx1[ij], (ityp==0) ? binrngs0[ij] : binrngs1[ij] );
            //h_genevtmiss[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]);
            h_genevtmiss[ityp][ipt][iet][ij]->Sumw2();
	    	    
	    sprintf(name, "gen2_typ_%i_pt%i_eta%i_%i", ityp, ipt, iet, ij);
	    sprintf(title, "Gen2 %s %i %g %s", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij]);
            h_genevtvar2[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij] : nbinsx1[ij] , (ityp==0) ? binrngs0[ij] : binrngs1[ij]);	   
            //h_genevtvar2[ityp][ipt][iet][ij] = fs->make<TH1F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt] , (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]);	   
	    h_genevtvar2[ityp][ipt][iet][ij]->Sumw2();

	    if (isReconstruct) { 
	      sprintf(name, "corr_typ_%i_pt%i_eta%i_%i", ityp , ipt, iet, ij);
	      sprintf(title, "Gen_Reco %s %i %g %s", typname[ityp], int(leadingPtThreshold[ipt]), etarange[iet], vartitle[ij]);
              h_2devtvar[ityp][ipt][iet][ij] = fs->make<TH2F>(name, title, (ityp==0) ? rnbinsx0[ij] : rnbinsx1[ij], (ityp==0) ? rbinrngs0[ij] : rbinrngs1[ij], (ityp==0) ? nbinsx0[ij] : nbinsx1[ij], (ityp==0) ? binrngs0[ij] : binrngs1[ij]);
             //h_2devtvar[ityp][ipt][iet][ij] = fs->make<TH2F>(name, title, (ityp==0) ? nbinsx0[ij] : nbinsx1[ij], (ityp==0) ? binrngs0[ij] : binrngs1[ij], (ityp==0) ? nbinsx0[ij] : nbinsx1[ij], (ityp==0) ? binrngs0[ij] : binrngs1[ij]);
             //h_2devtvar[ityp][ipt][iet][ij] = fs->make<TH2F>(name, title, (ityp==0) ? rnbinsx0[ij][ipt] : rnbinsx1[ij][ipt], (ityp==0) ? rbinrngs0[ij][ipt] : rbinrngs1[ij][ipt], (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]);
             //h_2devtvar[ityp][ipt][iet][ij] = fs->make<TH2F>(name, title, (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt], (ityp==0) ? nbinsx0[ij][ipt] : nbinsx1[ij][ipt], (ityp==0) ? binrngs0[ij][ipt] : binrngs1[ij][ipt]); //Reco =Gen bin Each HT, Var
             //h_2devtvar[ityp][ipt][iet][ij] = fs->make<TH2F>(name, title, nbinsxgen[ij], -endx[ij], -startx[ij], nbinsx[ij], -endx[ij], -startx[ij]);
	     //h_2devtvar[ityp][ipt][iet][ij] = fs->make<TH2F>(name, title, nbinsx[ij], -endx[ij], -startx[ij], nbinsx[ij], -endx[ij], -startx[ij]);
	     h_2devtvar[ityp][ipt][iet][ij]->Sumw2();
	     }
	   }
	     //cout <<"ijx "<< ityp<<" "<< ipt<<" "<< iet<<" "<<ij<<endl;	  
	}
      }
    }
  }
//-------------------------------------------
 sprintf(name, "corr_jet");
 sprintf(title, "Gen_Reco_HT2");
 h_2ht=fs->make<TH2F>(name, title, 9, leadingPtThreshold, 9, leadingPtThreshold);
//-------------------------------------------
//RM 
for(int ipt=0; ipt<njetptmn; ipt++){
        for(int iet=0; iet<njetetamn; iet++){
		for (int ik=0; ik<10; ik++){
                	if (isReconstruct){
				sprintf(name, "RM_jc_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
	                        sprintf(title, "RM Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
        	                h_RM_D1J1[ik][ipt][iet] = fs->make<TH2F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], genjcd1bins[ik], genjcd1minran[ik], genjcd1maxran[ik]);
                	        h_RM_D1J1[ik][ipt][iet]->Sumw2();

				sprintf(name, "RM_jc_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "RM Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_RM_D1J2[ik][ipt][iet] = fs->make<TH2F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik], genjcd1bins[ik], genjcd1minran[ik], genjcd1maxran[ik]);
                                h_RM_D1J2[ik][ipt][iet]->Sumw2();

				sprintf(name, "RM_jc_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "RM Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_RM_D2J1[ik][ipt][iet] = fs->make<TH2F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik]);
                                h_RM_D2J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "RM_jc_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "RM Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_RM_D2J2[ik][ipt][iet] = fs->make<TH2F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik]);
                                h_RM_D2J2[ik][ipt][iet]->Sumw2();

				sprintf(name, "RM_jc_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "RM Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_RM_D3J1[ik][ipt][iet] = fs->make<TH2F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik]);
                                h_RM_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "RM_jc_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "RM Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_RM_D3J2[ik][ipt][iet] = fs->make<TH2F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik], genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik]);
                                h_RM_D3J2[ik][ipt][iet]->Sumw2();

			}
		}
	}
}
//-------------------------------------------
for(int ipt=0; ipt<njetptmn; ipt++){
	for(int iet=0; iet<njetetamn; iet++){
		for (int ik=0; ik<10; ik++){
        		if (isReconstruct){
				//recofake
				sprintf(name, "recofake_jc_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
	                        sprintf(title, "Reco Fake Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
        	                h_recofake_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                	        h_recofake_D1J1[ik][ipt][iet]->Sumw2();

				sprintf(name, "recofake_jc_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco Fake Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_recofake_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd1bins[ik], recojcd1minran[ik], recojcd1maxran[ik]);
                                h_recofake_D1J2[ik][ipt][iet]->Sumw2();

				sprintf(name, "recofake_jc_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco Fake Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_recofake_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recofake_D2J1[ik][ipt][iet]->Sumw2();

				sprintf(name, "recofake_jc_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco Fake Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_recofake_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recofake_D2J2[ik][ipt][iet]->Sumw2();

				sprintf(name, "recofake_jc_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco Fake Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_recofake_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recofake_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "recofake_jc_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Reco Fake Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_recofake_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, recojcd23bins[ik], recojcd23minran[ik], recojcd23maxran[ik]);
                                h_recofake_D3J2[ik][ipt][iet]->Sumw2();
			}
				//genmiss
				sprintf(name, "genmiss_jc_d0_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Gen Miss Q_{1}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_genmiss_D1J1[ik][ipt][iet] = fs->make<TH1F>(name, title, genjcd1bins[ik], genjcd1minran[ik], genjcd1maxran[ik]);
                                h_genmiss_D1J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "genmiss_jc_d0_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Gen Miss Q_{2}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_genmiss_D1J2[ik][ipt][iet] = fs->make<TH1F>(name, title, genjcd1bins[ik], genjcd1minran[ik], genjcd1maxran[ik]);
                                h_genmiss_D1J2[ik][ipt][iet]->Sumw2();

				sprintf(name, "genmiss_jc_d1_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Gen Miss Q_{1L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_genmiss_D2J1[ik][ipt][iet] = fs->make<TH1F>(name, title, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik]);
                                h_genmiss_D2J1[ik][ipt][iet]->Sumw2();
                                
				sprintf(name, "genmiss_jc_d1_j1_k%i_pt%i_eta%i", ik, ipt, iet);
				sprintf(title, "Gen Miss Q_{2L}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_genmiss_D2J2[ik][ipt][iet] = fs->make<TH1F>(name, title, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik]);
                                h_genmiss_D2J2[ik][ipt][iet]->Sumw2();

				sprintf(name, "genmiss_jc_d2_j0_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Gen Miss Q_{1T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_genmiss_D3J1[ik][ipt][iet] = fs->make<TH1F>(name, title, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik]);
                                h_genmiss_D3J1[ik][ipt][iet]->Sumw2();

                                sprintf(name, "genmiss_jc_d2_j1_k%i_pt%i_eta%i", ik, ipt, iet);
                                sprintf(title, "Gen Miss Q_{2T}^{%g} %i %g", kappa[ik], int(leadingPtThreshold[ipt]), etarange[iet] );
                                h_genmiss_D3J2[ik][ipt][iet] = fs->make<TH1F>(name, title, genjcd23bins[ik], genjcd23minran[ik], genjcd23maxran[ik]);
                                h_genmiss_D3J2[ik][ipt][iet]->Sumw2();				
				}
			}
		}
//-------------------------------------------
#ifndef GENPART                     
  //recojt_hist = fs->make<TH1F>("recojt_hist","# of recojets",20,-0.5, 19.5);
  //recojt_hist->Sumw2();
  //recojt_pt = fs->make<TH1F>("recojt_pt","Et_{recojets}",100,20., 2020.);
  //recojt_pt->Sumw2();
  recojt_eta = fs->make<TH1F>("recojt_eta","#eta_{recojets}",100,-2.5, 2.5);
  recojt_eta->Sumw2();
  recojt_phi = fs->make<TH1F>("recojt_phi","#phi_{recojets}",100,-M_PI, M_PI);
  recojt_phi->Sumw2();

  recojtg_eta = fs->make<TH1F>("recojtg_eta","#eta_{recojetsg}",100,-2.5, 2.5);
  recojtg_eta->Sumw2();
  recojtg_phi = fs->make<TH1F>("recojtg_phi","#phi_{recojetsg}",100,-M_PI, M_PI);
  recojtg_phi->Sumw2();

  recojtqg_eta = fs->make<TH1F>("recojtqg_eta","#eta_{recojetsqg}",100,-2.5, 2.5);
  recojtqg_eta->Sumw2();
  recojtqg_phi = fs->make<TH1F>("recojtqg_phi","#phi_{recojetsqg}",100,-M_PI, M_PI);
  recojtqg_phi->Sumw2();

  recojtag_eta = fs->make<TH1F>("recojtag_eta","#eta_{recojetsag}",100,-2.5, 2.5);
  recojtag_eta->Sumw2();
  recojtag_phi = fs->make<TH1F>("recojtag_phi","#phi_{recojetsag}",100,-M_PI, M_PI);
  recojtag_phi->Sumw2();

  recojtb_eta = fs->make<TH1F>("recojtb_eta","#eta_{recojetsb}",100,-2.5, 2.5);
  recojtb_eta->Sumw2();
  recojtb_phi = fs->make<TH1F>("recojtb_phi","#phi_{recojetsb}",100,-M_PI, M_PI);
  recojtb_phi->Sumw2();

  recojtqb_eta = fs->make<TH1F>("recojtqb_eta","#eta_{recojetsqb}",100,-2.5, 2.5);
  recojtqb_eta->Sumw2();
  recojtqb_phi = fs->make<TH1F>("recojtqb_phi","#phi_{recojetsqb}",100,-M_PI, M_PI);
  recojtqb_phi->Sumw2();

  recojtab_eta = fs->make<TH1F>("recojtab_eta","#eta_{recojetsab}",100,-2.5, 2.5);
  recojtab_eta->Sumw2();
  recojtab_phi = fs->make<TH1F>("recojtab_phi","#phi_{recojetsab}",100,-M_PI, M_PI);
  recojtab_phi->Sumw2();

  recojtc_eta = fs->make<TH1F>("recojtc_eta","#eta_{recojetsc}",100,-2.5, 2.5);
  recojtc_eta->Sumw2();
  recojtc_phi = fs->make<TH1F>("recojtc_phi","#phi_{recojetsc}",100,-M_PI, M_PI);
  recojtc_phi->Sumw2();

  recojtqc_eta = fs->make<TH1F>("recojtqc_eta","#eta_{recojetsqc}",100,-2.5, 2.5);
  recojtqc_eta->Sumw2();
  recojtqc_phi = fs->make<TH1F>("recojtqc_phi","#phi_{recojetsqc}",100,-M_PI, M_PI);
  recojtqc_phi->Sumw2();

  recojtac_eta = fs->make<TH1F>("recojtac_eta","#eta_{recojetsac}",100,-2.5, 2.5);
  recojtac_eta->Sumw2();
  recojtac_phi = fs->make<TH1F>("recojtac_phi","#phi_{recojetsac}",100,-M_PI, M_PI);
  recojtac_phi->Sumw2();

  recojts_eta = fs->make<TH1F>("recojts_eta","#eta_{recojetss}",100,-2.5, 2.5);
  recojts_eta->Sumw2();
  recojts_phi = fs->make<TH1F>("recojts_phi","#phi_{recojetss}",100,-M_PI, M_PI);
  recojts_phi->Sumw2();

  recojtqs_eta = fs->make<TH1F>("recojtqs_eta","#eta_{recojetsqs}",100,-2.5, 2.5);
  recojtqs_eta->Sumw2();
  recojtqs_phi = fs->make<TH1F>("recojtqs_phi","#phi_{recojetsqs}",100,-M_PI, M_PI);
  recojtqs_phi->Sumw2();

  recojtas_eta = fs->make<TH1F>("recojtas_eta","#eta_{recojetsas}",100,-2.5, 2.5);
  recojtas_eta->Sumw2();
  recojtas_phi = fs->make<TH1F>("recojtas_phi","#phi_{recojetsas}",100,-M_PI, M_PI);
  recojtas_phi->Sumw2();

  recojtu_eta = fs->make<TH1F>("recojtu_eta","#eta_{recojetsu}",100,-2.5, 2.5);
  recojtu_eta->Sumw2();
  recojtu_phi = fs->make<TH1F>("recojtu_phi","#phi_{recojetsu}",100,-M_PI, M_PI);
  recojtu_phi->Sumw2();

  recojtqu_eta = fs->make<TH1F>("recojtqu_eta","#eta_{recojetsqu}",100,-2.5, 2.5);
  recojtqu_eta->Sumw2();
  recojtqu_phi = fs->make<TH1F>("recojtqu_phi","#phi_{recojetsqu}",100,-M_PI, M_PI);
  recojtqu_phi->Sumw2();

  recojtau_eta = fs->make<TH1F>("recojtau_eta","#eta_{recojetsau}",100,-2.5, 2.5);
  recojtau_eta->Sumw2();
  recojtau_phi = fs->make<TH1F>("recojtau_phi","#phi_{recojetsau}",100,-M_PI, M_PI);
  recojtau_phi->Sumw2();

  recojtd_eta = fs->make<TH1F>("recojtd_eta","#eta_{recojetsd}",100,-2.5, 2.5);
  recojtd_eta->Sumw2();
  recojtd_phi = fs->make<TH1F>("recojtd_phi","#phi_{recojetsd}",100,-M_PI, M_PI);
  recojtd_phi->Sumw2();

  recojtqd_eta = fs->make<TH1F>("recojtqd_eta","#eta_{recojetsqd}",100,-2.5, 2.5);
  recojtqd_eta->Sumw2();
  recojtqd_phi = fs->make<TH1F>("recojtqd_phi","#phi_{recojetsqd}",100,-M_PI, M_PI);
  recojtqd_phi->Sumw2();

  recojtad_eta = fs->make<TH1F>("recojtad_eta","#eta_{recojetsad}",100,-2.5, 2.5);
  recojtad_eta->Sumw2();
  recojtad_phi = fs->make<TH1F>("recojtad_phi","#phi_{recojetsad}",100,-M_PI, M_PI);
  recojtad_phi->Sumw2();

  //recojt1_pt = fs->make<TH1F>("recojet1_pt","Et_{recojets}",100,20., 2020.);
  //recojt1_pt->Sumw2();
  recojt1_eta = fs->make<TH1F>("recojet1_eta","#eta_{recojets}",100,-2.5, 2.5);
  recojt1_eta->Sumw2();
  recojt1_phi = fs->make<TH1F>("recojet1_phi","#phi_{recojets}",100,-M_PI, M_PI);
  recojt1_phi->Sumw2();

  recojt1g_eta = fs->make<TH1F>("recojet1g_eta","#eta_{recojetsg}",100,-2.5, 2.5);
  recojt1g_eta->Sumw2();
  recojt1g_phi = fs->make<TH1F>("recojet1g_phi","#phi_{recojetsg}",100,-M_PI, M_PI);
  recojt1g_phi->Sumw2();

  recojt1qg_eta = fs->make<TH1F>("recojet1qg_eta","#eta_{recojetsqg}",100,-2.5, 2.5);
  recojt1qg_eta->Sumw2();
  recojt1qg_phi = fs->make<TH1F>("recojet1qg_phi","#phi_{recojetsqg}",100,-M_PI, M_PI);
  recojt1qg_phi->Sumw2();

  recojt1ag_eta = fs->make<TH1F>("recojet1ag_eta","#eta_{recojetsag}",100,-2.5, 2.5);
  recojt1ag_eta->Sumw2();
  recojt1ag_phi = fs->make<TH1F>("recojet1ag_phi","#phi_{recojetsag}",100,-M_PI, M_PI);
  recojt1ag_phi->Sumw2();

  recojt1b_eta = fs->make<TH1F>("recojet1b_eta","#eta_{recojetsb}",100,-2.5, 2.5);
  recojt1b_eta->Sumw2();
  recojt1b_phi = fs->make<TH1F>("recojet1b_phi","#phi_{recojetsb}",100,-M_PI, M_PI);
  recojt1b_phi->Sumw2();

  recojt1qb_eta = fs->make<TH1F>("recojet1qb_eta","#eta_{recojetsqb}",100,-2.5, 2.5);
  recojt1qb_eta->Sumw2();
  recojt1qb_phi = fs->make<TH1F>("recojet1qb_phi","#phi_{recojetsqb}",100,-M_PI, M_PI);
  recojt1qb_phi->Sumw2();

  recojt1ab_eta = fs->make<TH1F>("recojet1ab_eta","#eta_{recojetsab}",100,-2.5, 2.5);
  recojt1ab_eta->Sumw2();
  recojt1ab_phi = fs->make<TH1F>("recojet1ab_phi","#phi_{recojetsab}",100,-M_PI, M_PI);
  recojt1ab_phi->Sumw2();

  recojt1c_eta = fs->make<TH1F>("recojet1c_eta","#eta_{recojetsc}",100,-2.5, 2.5);
  recojt1c_eta->Sumw2();
  recojt1c_phi = fs->make<TH1F>("recojet1c_phi","#phi_{recojetsc}",100,-M_PI, M_PI);
  recojt1c_phi->Sumw2();

  recojt1qc_eta = fs->make<TH1F>("recojet1qc_eta","#eta_{recojetsqc}",100,-2.5, 2.5);
  recojt1qc_eta->Sumw2();
  recojt1qc_phi = fs->make<TH1F>("recojet1qc_phi","#phi_{recojetsqc}",100,-M_PI, M_PI);
  recojt1qc_phi->Sumw2();

  recojt1ac_eta = fs->make<TH1F>("recojet1ac_eta","#eta_{recojetsac}",100,-2.5, 2.5);
  recojt1ac_eta->Sumw2();
  recojt1ac_phi = fs->make<TH1F>("recojet1ac_phi","#phi_{recojetsac}",100,-M_PI, M_PI);
  recojt1ac_phi->Sumw2();

  recojt1s_eta = fs->make<TH1F>("recojet1s_eta","#eta_{recojetss}",100,-2.5, 2.5);
  recojt1s_eta->Sumw2();
  recojt1s_phi = fs->make<TH1F>("recojet1s_phi","#phi_{recojetss}",100,-M_PI, M_PI);
  recojt1s_phi->Sumw2();

  recojt1qs_eta = fs->make<TH1F>("recojet1qs_eta","#eta_{recojetsqs}",100,-2.5, 2.5);
  recojt1qs_eta->Sumw2();
  recojt1qs_phi = fs->make<TH1F>("recojet1qs_phi","#phi_{recojetsqs}",100,-M_PI, M_PI);
  recojt1qs_phi->Sumw2();

  recojt1as_eta = fs->make<TH1F>("recojet1as_eta","#eta_{recojetsas}",100,-2.5, 2.5);
  recojt1as_eta->Sumw2();
  recojt1as_phi = fs->make<TH1F>("recojet1as_phi","#phi_{recojetsas}",100,-M_PI, M_PI);
  recojt1as_phi->Sumw2();

  recojt1u_eta = fs->make<TH1F>("recojet1u_eta","#eta_{recojetsu}",100,-2.5, 2.5);
  recojt1u_eta->Sumw2();
  recojt1u_phi = fs->make<TH1F>("recojet1u_phi","#phi_{recojetsu}",100,-M_PI, M_PI);
  recojt1u_phi->Sumw2();

  recojt1qu_eta = fs->make<TH1F>("recojet1qu_eta","#eta_{recojetsqu}",100,-2.5, 2.5);
  recojt1qu_eta->Sumw2();
  recojt1qu_phi = fs->make<TH1F>("recojet1qu_phi","#phi_{recojetsqu}",100,-M_PI, M_PI);
  recojt1qu_phi->Sumw2();

  recojt1au_eta = fs->make<TH1F>("recojet1au_eta","#eta_{recojetsau}",100,-2.5, 2.5);
  recojt1au_eta->Sumw2();
  recojt1au_phi = fs->make<TH1F>("recojet1au_phi","#phi_{recojetsau}",100,-M_PI, M_PI);
  recojt1au_phi->Sumw2();

  recojt1d_eta = fs->make<TH1F>("recojet1d_eta","#eta_{recojetsd}",100,-2.5, 2.5);
  recojt1d_eta->Sumw2();
  recojt1d_phi = fs->make<TH1F>("recojet1d_phi","#phi_{recojetsd}",100,-M_PI, M_PI);
  recojt1d_phi->Sumw2();

  recojt1qd_eta = fs->make<TH1F>("recojet1qd_eta","#eta_{recojetsqd}",100,-2.5, 2.5);
  recojt1qd_eta->Sumw2();
  recojt1qd_phi = fs->make<TH1F>("recojet1qd_phi","#phi_{recojetsqd}",100,-M_PI, M_PI);
  recojt1qd_phi->Sumw2();

  recojt1ad_eta = fs->make<TH1F>("recojet1ad_eta","#eta_{recojetsad}",100,-2.5, 2.5);
  recojt1ad_eta->Sumw2();
  recojt1ad_phi = fs->make<TH1F>("recojet1ad_phi","#phi_{recojetsad}",100,-M_PI, M_PI);
  recojt1ad_phi->Sumw2();

  //recojt2_pt = fs->make<TH1F>("recojet2_pt","Et_{recojets}",100,20., 2020.);
  //recojt2_pt->Sumw2();
  recojt2_eta = fs->make<TH1F>("recojet2_eta","#eta_{recojets}",100,-2.5, 2.5);
  recojt2_eta->Sumw2();
  recojt2_phi = fs->make<TH1F>("recojet2_phi","#phi_{recojets}",100,-M_PI, M_PI);
  recojt2_phi->Sumw2();
  
  recojt2g_eta = fs->make<TH1F>("recojet2g_eta","#eta_{recojetsg}",100,-2.5, 2.5);
  recojt2g_eta->Sumw2();
  recojt2g_phi = fs->make<TH1F>("recojet2g_phi","#phi_{recojetsg}",100,-M_PI, M_PI);
  recojt2g_phi->Sumw2();

  recojt2qg_eta = fs->make<TH1F>("recojet2qg_eta","#eta_{recojetsqg}",100,-2.5, 2.5);
  recojt2qg_eta->Sumw2();
  recojt2qg_phi = fs->make<TH1F>("recojet2qg_phi","#phi_{recojetsqg}",100,-M_PI, M_PI);
  recojt2qg_phi->Sumw2();

  recojt2ag_eta = fs->make<TH1F>("recojet2ag_eta","#eta_{recojetsag}",100,-2.5, 2.5);
  recojt2ag_eta->Sumw2();
  recojt2ag_phi = fs->make<TH1F>("recojet2ag_phi","#phi_{recojetsag}",100,-M_PI, M_PI);
  recojt2ag_phi->Sumw2();

  recojt2b_eta = fs->make<TH1F>("recojet2b_eta","#eta_{recojetsb}",100,-2.5, 2.5);
  recojt2b_eta->Sumw2();
  recojt2b_phi = fs->make<TH1F>("recojet2b_phi","#phi_{recojetsb}",100,-M_PI, M_PI);
  recojt2b_phi->Sumw2();

  recojt2qb_eta = fs->make<TH1F>("recojet2qb_eta","#eta_{recojetsqb}",100,-2.5, 2.5);
  recojt2qb_eta->Sumw2();
  recojt2qb_phi = fs->make<TH1F>("recojet2qb_phi","#phi_{recojetsqb}",100,-M_PI, M_PI);
  recojt2qb_phi->Sumw2();

  recojt2ab_eta = fs->make<TH1F>("recojet2ab_eta","#eta_{recojetsab}",100,-2.5, 2.5);
  recojt2ab_eta->Sumw2();
  recojt2ab_phi = fs->make<TH1F>("recojet2ab_phi","#phi_{recojetsab}",100,-M_PI, M_PI);
  recojt2ab_phi->Sumw2();

  recojt2c_eta = fs->make<TH1F>("recojet2c_eta","#eta_{recojetsc}",100,-2.5, 2.5);
  recojt2c_eta->Sumw2();
  recojt2c_phi = fs->make<TH1F>("recojet2c_phi","#phi_{recojetsc}",100,-M_PI, M_PI);
  recojt2c_phi->Sumw2();

  recojt2qc_eta = fs->make<TH1F>("recojet2qc_eta","#eta_{recojetsqc}",100,-2.5, 2.5);
  recojt2qc_eta->Sumw2();
  recojt2qc_phi = fs->make<TH1F>("recojet2qc_phi","#phi_{recojetsqc}",100,-M_PI, M_PI);
  recojt2qc_phi->Sumw2();

  recojt2ac_eta = fs->make<TH1F>("recojet2ac_eta","#eta_{recojetsac}",100,-2.5, 2.5);
  recojt2ac_eta->Sumw2();
  recojt2ac_phi = fs->make<TH1F>("recojet2ac_phi","#phi_{recojetsac}",100,-M_PI, M_PI);
  recojt2ac_phi->Sumw2();

  recojt2s_eta = fs->make<TH1F>("recojet2s_eta","#eta_{recojetss}",100,-2.5, 2.5);
  recojt2s_eta->Sumw2();
  recojt2s_phi = fs->make<TH1F>("recojet2s_phi","#phi_{recojetss}",100,-M_PI, M_PI);
  recojt2s_phi->Sumw2();

  recojt2qs_eta = fs->make<TH1F>("recojet2qs_eta","#eta_{recojetsqs}",100,-2.5, 2.5);
  recojt2qs_eta->Sumw2();
  recojt2qs_phi = fs->make<TH1F>("recojet2qs_phi","#phi_{recojetsqs}",100,-M_PI, M_PI);
  recojt2qs_phi->Sumw2();

  recojt2as_eta = fs->make<TH1F>("recojet2as_eta","#eta_{recojetsas}",100,-2.5, 2.5);
  recojt2as_eta->Sumw2();
  recojt2as_phi = fs->make<TH1F>("recojet2as_phi","#phi_{recojetsas}",100,-M_PI, M_PI);
  recojt2as_phi->Sumw2();

  recojt2u_eta = fs->make<TH1F>("recojet2u_eta","#eta_{recojetsu}",100,-2.5, 2.5);
  recojt2u_eta->Sumw2();
  recojt2u_phi = fs->make<TH1F>("recojet2u_phi","#phi_{recojetsu}",100,-M_PI, M_PI);
  recojt2u_phi->Sumw2();

  recojt2qu_eta = fs->make<TH1F>("recojet2qu_eta","#eta_{recojetsqu}",100,-2.5, 2.5);
  recojt2qu_eta->Sumw2();
  recojt2qu_phi = fs->make<TH1F>("recojet2qu_phi","#phi_{recojetsqu}",100,-M_PI, M_PI);
  recojt2qu_phi->Sumw2();

  recojt2au_eta = fs->make<TH1F>("recojet2au_eta","#eta_{recojetsau}",100,-2.5, 2.5);
  recojt2au_eta->Sumw2();
  recojt2au_phi = fs->make<TH1F>("recojet2au_phi","#phi_{recojetsau}",100,-M_PI, M_PI);
  recojt2au_phi->Sumw2();

  recojt2d_eta = fs->make<TH1F>("recojet2d_eta","#eta_{recojetsd}",100,-2.5, 2.5);
  recojt2d_eta->Sumw2();
  recojt2d_phi = fs->make<TH1F>("recojet2d_phi","#phi_{recojetsd}",100,-M_PI, M_PI);
  recojt2d_phi->Sumw2();

  recojt2qd_eta = fs->make<TH1F>("recojet2qd_eta","#eta_{recojetsqd}",100,-2.5, 2.5);
  recojt2qd_eta->Sumw2();
  recojt2qd_phi = fs->make<TH1F>("recojet2qd_phi","#phi_{recojetsqd}",100,-M_PI, M_PI);
  recojt2qd_phi->Sumw2();

  recojt2ad_eta = fs->make<TH1F>("recojet2ad_eta","#eta_{recojetsad}",100,-2.5, 2.5);
  recojt2ad_eta->Sumw2();
  recojt2ad_phi = fs->make<TH1F>("recojet2ad_phi","#phi_{recojetsad}",100,-M_PI, M_PI);
  recojt2ad_phi->Sumw2();

  //recojt3_pt = fs->make<TH1F>("recojet2_pt","Et_{recojets}",100,20., 2020.);
  //recojt3_pt->Sumw2();
  recojt3_eta = fs->make<TH1F>("recojet3_eta","#eta_{recojets}",100,-2.5, 2.5);
  recojt3_eta->Sumw2();
  recojt3_phi = fs->make<TH1F>("recojet3_phi","#phi_{recojets}",100,-M_PI, M_PI);
  recojt3_phi->Sumw2();

  hprof = fs->make<TProfile>("hprof","hprof",40,60,2000);
  hprof->Sumw2();

  hchpt = fs->make<TProfile>("hchpt","hchpt",40,60,2000);
  hchpt->Sumw2();

  for(int jk=0; jk<njetetamn; jk++){
    sprintf(name, "recojetallave_pt_%i",jk);
    sprintf(title, "Et_{recojetsallave}_%g", etarange[jk]);
    recojtallave_pt[jk] = fs->make<TH1F>(name,title,400, 20., 2020.);
    recojtallave_pt[jk]->Sumw2();

    sprintf(name, "recojetallavewt1_pt_%i",jk);
    sprintf(title, "Et_{recojetsallavewt1}_%g", etarange[jk]);
    recojtallavewt1_pt[jk] = fs->make<TH1F>(name,title,400, 20., 2020.);
    recojtallavewt1_pt[jk]->Sumw2();

    sprintf(name, "recojt_pt_%i",jk);
    sprintf(title, "Et_{recojets}_%g", etarange[jk]);
    recojt_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt_pt[jk]->Sumw2();

    sprintf(name, "recojtg_pt_%i",jk);
    sprintf(title, "Et_{recojetsg}_%g", etarange[jk]);
    recojtg_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtg_pt[jk]->Sumw2();

    sprintf(name, "recojtqg_pt_%i",jk);
    sprintf(title, "Et_{recojetsqg}_%g", etarange[jk]);
    recojtqg_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtqg_pt[jk]->Sumw2();

    sprintf(name, "recojtag_pt_%i",jk);
    sprintf(title, "Et_{recojetsag}_%g", etarange[jk]);
    recojtag_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtag_pt[jk]->Sumw2();

    sprintf(name, "recojtb_pt_%i",jk);
    sprintf(title, "Et_{recojetsb}_%g", etarange[jk]);
    recojtb_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtb_pt[jk]->Sumw2();

    sprintf(name, "recojtqb_pt_%i",jk);
    sprintf(title, "Et_{recojetsqb}_%g", etarange[jk]);
    recojtqb_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtqb_pt[jk]->Sumw2();

    sprintf(name, "recojtab_pt_%i",jk);
    sprintf(title, "Et_{recojetsab}_%g", etarange[jk]);
    recojtab_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtab_pt[jk]->Sumw2();

    sprintf(name, "recojtc_pt_%i",jk);
    sprintf(title, "Et_{recojetsc}_%g", etarange[jk]);
    recojtc_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtc_pt[jk]->Sumw2();

    sprintf(name, "recojtqc_pt_%i",jk);
    sprintf(title, "Et_{recojetsqc}_%g", etarange[jk]);
    recojtqc_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtqc_pt[jk]->Sumw2();

    sprintf(name, "recojtac_pt_%i",jk);
    sprintf(title, "Et_{recojetsac}_%g", etarange[jk]);
    recojtac_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtac_pt[jk]->Sumw2();

    sprintf(name, "recojts_pt_%i",jk);
    sprintf(title, "Et_{recojetss}_%g", etarange[jk]);
    recojts_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojts_pt[jk]->Sumw2();

    sprintf(name, "recojtqs_pt_%i",jk);
    sprintf(title, "Et_{recojetsqs}_%g", etarange[jk]);
    recojtqs_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtqs_pt[jk]->Sumw2();

    sprintf(name, "recojtas_pt_%i",jk);
    sprintf(title, "Et_{recojetsas}_%g", etarange[jk]);
    recojtas_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtas_pt[jk]->Sumw2();

    sprintf(name, "recojtu_pt_%i",jk);
    sprintf(title, "Et_{recojetsu}_%g", etarange[jk]);
    recojtu_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtu_pt[jk]->Sumw2();

    sprintf(name, "recojtqu_pt_%i",jk);
    sprintf(title, "Et_{recojetsqu}_%g", etarange[jk]);
    recojtqu_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtqu_pt[jk]->Sumw2();

    sprintf(name, "recojtau_pt_%i",jk);
    sprintf(title, "Et_{recojetsau}_%g", etarange[jk]);
    recojtau_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtau_pt[jk]->Sumw2();

    sprintf(name, "recojtd_pt_%i",jk);
    sprintf(title, "Et_{recojetsd}_%g", etarange[jk]);
    recojtd_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtd_pt[jk]->Sumw2();

    sprintf(name, "recojtqd_pt_%i",jk);
    sprintf(title, "Et_{recojetsqd}_%g", etarange[jk]);
    recojtqd_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtqd_pt[jk]->Sumw2();

    sprintf(name, "recojtad_pt_%i",jk);
    sprintf(title, "Et_{recojetsad}_%g", etarange[jk]);
    recojtad_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojtad_pt[jk]->Sumw2();

    sprintf(name, "recojet1_pt_%i",jk);
    sprintf(title, "Et_{recojets1}_%g", etarange[jk]);
    recojt1_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1_pt[jk]->Sumw2();

    sprintf(name, "recojet1g_pt_%i",jk);
    sprintf(title, "Et_{recojets1g}_%g", etarange[jk]);
    recojt1g_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1g_pt[jk]->Sumw2();

    sprintf(name, "recojet1qg_pt_%i",jk);
    sprintf(title, "Et_{recojets1qg}_%g", etarange[jk]);
    recojt1qg_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1qg_pt[jk]->Sumw2();

    sprintf(name, "recojet1ag_pt_%i",jk);
    sprintf(title, "Et_{recojets1ag}_%g", etarange[jk]);
    recojt1ag_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1ag_pt[jk]->Sumw2();

    sprintf(name, "recojet1b_pt_%i",jk);
    sprintf(title, "Et_{recojets1b}_%g", etarange[jk]);
    recojt1b_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1b_pt[jk]->Sumw2();

    sprintf(name, "recojet1qb_pt_%i",jk);
    sprintf(title, "Et_{recojets1qb}_%g", etarange[jk]);
    recojt1qb_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1qb_pt[jk]->Sumw2();

    sprintf(name, "recojet1ab_pt_%i",jk);
    sprintf(title, "Et_{recojets1ab}_%g", etarange[jk]);
    recojt1ab_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1ab_pt[jk]->Sumw2();

    sprintf(name, "recojet1c_pt_%i",jk);
    sprintf(title, "Et_{recojets1c}_%g", etarange[jk]);
    recojt1c_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1c_pt[jk]->Sumw2();

    sprintf(name, "recojet1qc_pt_%i",jk);
    sprintf(title, "Et_{recojets1qc}_%g", etarange[jk]);
    recojt1qc_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1qc_pt[jk]->Sumw2();

    sprintf(name, "recojet1ac_pt_%i",jk);
    sprintf(title, "Et_{recojets1ac}_%g", etarange[jk]);
    recojt1ac_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1ac_pt[jk]->Sumw2();

    sprintf(name, "recojet1s_pt_%i",jk);
    sprintf(title, "Et_{recojets1s}_%g", etarange[jk]);
    recojt1s_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1s_pt[jk]->Sumw2();

    sprintf(name, "recojet1qs_pt_%i",jk);
    sprintf(title, "Et_{recojets1qs}_%g", etarange[jk]);
    recojt1qs_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1qs_pt[jk]->Sumw2();

    sprintf(name, "recojet1as_pt_%i",jk);
    sprintf(title, "Et_{recojets1as}_%g", etarange[jk]);
    recojt1as_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1as_pt[jk]->Sumw2();

    sprintf(name, "recojet1u_pt_%i",jk);
    sprintf(title, "Et_{recojets1u}_%g", etarange[jk]);
    recojt1u_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1u_pt[jk]->Sumw2();

    sprintf(name, "recojet1qu_pt_%i",jk);
    sprintf(title, "Et_{recojets1qu}_%g", etarange[jk]);
    recojt1qu_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1qu_pt[jk]->Sumw2();

    sprintf(name, "recojet1au_pt_%i",jk);
    sprintf(title, "Et_{recojets1au}_%g", etarange[jk]);
    recojt1au_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1au_pt[jk]->Sumw2();

    sprintf(name, "recojet1d_pt_%i",jk);
    sprintf(title, "Et_{recojets1d}_%g", etarange[jk]);
    recojt1d_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1d_pt[jk]->Sumw2();

    sprintf(name, "recojet1qd_pt_%i",jk);
    sprintf(title, "Et_{recojets1qd}_%g", etarange[jk]);
    recojt1qd_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1qd_pt[jk]->Sumw2();

    sprintf(name, "recojet1ad_pt_%i",jk);
    sprintf(title, "Et_{recojets1ad}_%g", etarange[jk]);
    recojt1ad_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt1ad_pt[jk]->Sumw2();



    sprintf(name, "recojet2_pt_%i",jk);
    sprintf(title, "Et_{recojets2}_%g", etarange[jk]);
    recojt2_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2_pt[jk]->Sumw2();

    sprintf(name, "recojet2g_pt_%i",jk);
    sprintf(title, "Et_{recojets2g}_%g", etarange[jk]);
    recojt2g_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2g_pt[jk]->Sumw2();

    sprintf(name, "recojet2qg_pt_%i",jk);
    sprintf(title, "Et_{recojets2qg}_%g", etarange[jk]);
    recojt2qg_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2qg_pt[jk]->Sumw2();

    sprintf(name, "recojet2ag_pt_%i",jk);
    sprintf(title, "Et_{recojets2ag}_%g", etarange[jk]);
    recojt2ag_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2ag_pt[jk]->Sumw2();

    sprintf(name, "recojet2b_pt_%i",jk);
    sprintf(title, "Et_{recojets2b}_%g", etarange[jk]);
    recojt2b_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2b_pt[jk]->Sumw2();

    sprintf(name, "recojet2qb_pt_%i",jk);
    sprintf(title, "Et_{recojets2qb}_%g", etarange[jk]);
    recojt2qb_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2qb_pt[jk]->Sumw2();

    sprintf(name, "recojet2ab_pt_%i",jk);
    sprintf(title, "Et_{recojets2ab}_%g", etarange[jk]);
    recojt2ab_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2ab_pt[jk]->Sumw2();

    sprintf(name, "recojet2c_pt_%i",jk);
    sprintf(title, "Et_{recojets2c}_%g", etarange[jk]);
    recojt2c_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2c_pt[jk]->Sumw2();

    sprintf(name, "recojet2qc_pt_%i",jk);
    sprintf(title, "Et_{recojets2qc}_%g", etarange[jk]);
    recojt2qc_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2qc_pt[jk]->Sumw2();

    sprintf(name, "recojet2ac_pt_%i",jk);
    sprintf(title, "Et_{recojets2ac}_%g", etarange[jk]);
    recojt2ac_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2ac_pt[jk]->Sumw2();

    sprintf(name, "recojet2s_pt_%i",jk);
    sprintf(title, "Et_{recojets2s}_%g", etarange[jk]);
    recojt2s_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2s_pt[jk]->Sumw2();

    sprintf(name, "recojet2qs_pt_%i",jk);
    sprintf(title, "Et_{recojets2qs}_%g", etarange[jk]);
    recojt2qs_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2qs_pt[jk]->Sumw2();

    sprintf(name, "recojet2as_pt_%i",jk);
    sprintf(title, "Et_{recojets2as}_%g", etarange[jk]);
    recojt2as_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2as_pt[jk]->Sumw2();

    sprintf(name, "recojet2u_pt_%i",jk);
    sprintf(title, "Et_{recojets2u}_%g", etarange[jk]);
    recojt2u_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2u_pt[jk]->Sumw2();

    sprintf(name, "recojet2qu_pt_%i",jk);
    sprintf(title, "Et_{recojets2qu}_%g", etarange[jk]);
    recojt2qu_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2qu_pt[jk]->Sumw2();

    sprintf(name, "recojet2au_pt_%i",jk);
    sprintf(title, "Et_{recojets2au}_%g", etarange[jk]);
    recojt2au_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2au_pt[jk]->Sumw2();

    sprintf(name, "recojet2d_pt_%i",jk);
    sprintf(title, "Et_{recojets2d}_%g", etarange[jk]);
    recojt2d_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2d_pt[jk]->Sumw2();

    sprintf(name, "recojet2qd_pt_%i",jk);
    sprintf(title, "Et_{recojets2qd}_%g", etarange[jk]);
    recojt2qd_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2qd_pt[jk]->Sumw2();

    sprintf(name, "recojet2ad_pt_%i",jk);
    sprintf(title, "Et_{recojets2ad}_%g", etarange[jk]);
    recojt2ad_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt2ad_pt[jk]->Sumw2();


    sprintf(name, "recojet3_pt_%i",jk);
    sprintf(title, "Et_{recojets3}_%g", etarange[jk]);
    recojt3_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt3_pt[jk]->Sumw2();
    sprintf(title, "Et_{recojets3}_%g", etarange[jk]);
    recojt3_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    recojt3_pt[jk]->Sumw2();

    for (int kl=0; kl<nHLTmx; kl++) { 
      //sprintf(name, "recojt_pt_%i_%i",jk, kl);
      //sprintf(title, "Et_{recojets}_%g_%i", etarange[jk], kl);
      //recojt_pt[jk][kl] = fs->make<TH1F>(name,title, 400, 20., 2020.);
      //recojt_pt[jk][kl]->Sumw2();

      sprintf(name, "recojetave_pt_%i_%i",jk, kl);
      sprintf(title, "Et_{recojetsave}_%g_%i", etarange[jk], kl);
      recojtave_pt[jk][kl] = fs->make<TH1F>(name,title, 400, 20., 2020.);
      recojtave_pt[jk][kl]->Sumw2();

      sprintf(name, "recojetavewt1_pt_%i_%i",jk, kl);
      sprintf(title, "Et_{recojetsavewt1}_%g_%i", etarange[jk], kl);
      recojtavewt1_pt[jk][kl] = fs->make<TH1F>(name,title, 400, 20., 2020.);
      recojtavewt1_pt[jk][kl]->Sumw2();
    }

    sprintf(name, "recojetHT2_%i",jk);
    sprintf(title, "recojetsHT2_%g", etarange[jk]);

    recoht2_pt[jk] = fs->make<TH1F>(name, title, 400,20., 1500.);
    recoht2_pt[jk]->Sumw2();


    sprintf(name, "hjetdpt_%i",jk);
    sprintf(title, "dpt_{recojets12}_%g", etarange[jk]);

    hjetdpt[jk] = fs->make<TH1F>(name, title, 100,20., 500.);
    hjetdpt[jk]->Sumw2();

    sprintf(name, "hjetpt2bypt1_%i",jk);
    sprintf(title, "hjetpt2bypt1 reco jet_%g", etarange[jk]);

    hjetpt2bypt1[jk] = fs->make<TH1F>(name, title, 60, 0., 1.0);
    hjetpt2bypt1[jk]->Sumw2();

    sprintf(name, "hjetpt3bypt2_%i",jk);
    sprintf(title, "hjetpt3bypt2 reco jet_%g", etarange[jk]);
    hjetpt3bypt2[jk] = fs->make<TH1F>(name, title, 60, 0., 1.0);
    hjetpt3bypt2[jk]->Sumw2();

    sprintf(name, "hjetdphi_%i",jk);
    sprintf(title, "#phi_{recojets}_%g", etarange[jk]);
    hjetdphi[jk] = fs->make<TH1F>(name,title,100,-M_PI, M_PI);
    hjetdphi[jk]->Sumw2();
    sprintf(name, "hjetptbypl_%i",jk);
    sprintf(title, "1st recojet Pt*sin/1st Recojet_%g", etarange[jk]);
    hjetptbypl[jk] = fs->make<TH1F>(name, title, 60, 0., 1.0);
    hjetptbypl[jk]->Sumw2();

    //hjetpt2bypt1 = fs->make<TH1F>("hjetpt2bypt1", "hjetpt2bypt1 reco jet", 60, 0., 1.0);
    //hjetpt2bypt1->Sumw2();
    //hjetpt3bypt2 = fs->make<TH1F>("hjetpt2bypt1", "hjetpt2bypt1 reco jet", 60, 0., 1.0);
    //hjetpt3bypt2->Sumw2();
  }

  recochg_pt = fs->make<TH1F>("recochg_pt","Et_{recocharge_alljet}",100, 1., 101.);
  recochg_pt->Sumw2();
  recochg_eta = fs->make<TH1F>("recochg_eta","#eta_{recocharge_alljet}",100,-3., 3.);
  recochg_eta->Sumw2();
  recochg_phi = fs->make<TH1F>("recochg_phi","#phi_{recocharge_alljet}",100,-M_PI, M_PI);
  recochg_phi->Sumw2();

  recochg1_pt = fs->make<TH1F>("recochg1_pt","Et_{recocharge_jet1}",100, 1., 101.);
  recochg1_pt->Sumw2();
  recochg1_eta = fs->make<TH1F>("recochg1_eta","#eta_{recocharge_jet1}",100,-3., 3.);
  recochg1_eta->Sumw2();
  recochg1_phi = fs->make<TH1F>("recochg1_phi","#phi_{recocharge_jet1}",100,-M_PI, M_PI);
  recochg1_phi->Sumw2();

  recochg2_pt = fs->make<TH1F>("recochg2_pt","Et_{recocharge_jet2}",100, 1., 101.);
  recochg2_pt->Sumw2();
  recochg2_eta = fs->make<TH1F>("recochg2_eta","#eta_{recocharge_jet2}",100,-3., 3.);
  recochg2_eta->Sumw2();
  recochg2_phi = fs->make<TH1F>("recochg2_phi","#phi_{recocharge_jet2}",100,-M_PI, M_PI);
  recochg2_phi->Sumw2();

  recochg3_pt = fs->make<TH1F>("recochg3_pt","Et_{recocharge_jet3}",100, 1., 101.);
  recochg3_pt->Sumw2();
  recochg3_eta = fs->make<TH1F>("recochg3_eta","#eta_{recocharge_jet3}",100,-3., 3.);
  recochg3_eta->Sumw2();
  recochg3_phi = fs->make<TH1F>("recochg3_phi","#phi_{recocharge_jet3}",100,-M_PI, M_PI);
  recochg3_phi->Sumw2();

#endif
//-------------------------------------------
  for (int ij=0; ij<nhist; ij++) {
    sprintf(name, "anglex_%i", ij);
    vec_anglex[ij] = fs->make<TH1F>(name, name, 240, 0.7, 1.0);
  }
  //genjt_hist = fs->make<TH1F>("genjt_hist","# of genjets",20,-0.5, 19.5);
  //genjt_hist->Sumw2();
  for(int jk=0; jk<njetetamn; jk++){
    sprintf(name, "genjetallave_pt_%i",jk);
    sprintf(title, "Et_{genjetsallave}_%g", etarange[jk]);
    genjtallave_pt[jk] = fs->make<TH1F>(name,title,400, 20., 2020.);
    genjtallave_pt[jk]->Sumw2();

    sprintf(name, "genjt_pt_%i",jk);
    sprintf(title, "Et_{genjets}_%g", etarange[jk]);
    genjt_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    genjt_pt[jk]->Sumw2();

    sprintf(name, "genjet1_pt_%i",jk);
    sprintf(title, "Et_{genjets1}_%g", etarange[jk]);
    genjt1_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    genjt1_pt[jk]->Sumw2();

    sprintf(name, "genjet2_pt_%i",jk);
    sprintf(title, "Et_{genjets2}_%g", etarange[jk]);
    genjt2_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    genjt2_pt[jk]->Sumw2();

    sprintf(name, "genjet3_pt_%i",jk);
    sprintf(title, "Et_{genjets3}_%g", etarange[jk]);
    genjt3_pt[jk] = fs->make<TH1F>(name,title, 400, 20., 2020.);
    genjt3_pt[jk]->Sumw2();

    /*sprintf(name, "genjt_oth_pt_%i",jk);
    sprintf(title, "#Et_{genjets_oth}_%g", etarange[jk]);

    genjt_oth_pt[jk] = fs->make<TH1F>(name,title,100, 20., 2020.);
    genjt_oth_pt[jk]->Sumw2();
   */
    sprintf(name, "genjetdpt_%i",jk);
    sprintf(title, "dpt_{genjets12}_%g", etarange[jk]);

    genjetdpt[jk] = fs->make<TH1F>(name, title, 100,20., 500.);
    genjetdpt[jk]->Sumw2();

    sprintf(name, "genjetpt2bypt1_%i",jk);
    sprintf(title, "jetpt2bypt1 gen jet_%g", etarange[jk]);

    genjetpt2bypt1[jk] = fs->make<TH1F>(name, title, 60, 0., 1.0);
    genjetpt2bypt1[jk]->Sumw2();

    sprintf(name, "genjetpt3bypt2_%i",jk);
    sprintf(title, "hjetpt3bypt2 gen jet_%g", etarange[jk]);
    genjetpt3bypt2[jk] = fs->make<TH1F>(name, title, 60, 0., 1.0);
    genjetpt3bypt2[jk]->Sumw2();

    sprintf(name, "genjetdphi_%i",jk);
    sprintf(title, "#phi_{genjets}_%g", etarange[jk]);
    genjetdphi[jk] = fs->make<TH1F>(name,title,100,-M_PI, M_PI);
    genjetdphi[jk]->Sumw2();
    sprintf(name, "genjetptbypl_%i",jk);
    sprintf(title, "1st genjet Pt*sin/1st genjet_%g", etarange[jk]);
    genjetptbypl[jk] = fs->make<TH1F>(name, title, 60, 0., 1.0);
    genjetptbypl[jk]->Sumw2();
  }

  for(int jk=0; jk<njetetamn+2; jk++){
    sprintf(name, "response_jet_distribution_%i",jk);
    sprintf(title, "response_2d_ratio_eta_%g_%g", resetarange[jk], resetarange[jk+1]);
    resp_jet[jk] = fs->make<TH2F>(name,title,400, 20., 2020., 50 , 0.0, 5.0);
    resp_jet[jk]->Sumw2();
    sprintf(name, "response_jet_distribution_resolution%i",jk);
    sprintf(title, "response_1d_resulation_eta_%g_%g", resetarange[jk], resetarange[jk+1]);
    resp_jet1[jk] = fs->make<TH1F>(name,title, 100 , 0.0, 1);
    resp_jet1[jk]->Sumw2();
  }

  //genjt_pt = fs->make<TH1F>("genjt_pt","Et_{genjets}",100,20., 2020.);
  //genjt_pt->Sumw2();
  genjt_eta = fs->make<TH1F>("genjt_eta","#eta_{genjets}",100,-2.5, 2.5);
  genjt_eta->Sumw2();
  genjt_phi = fs->make<TH1F>("genjt_phi","#phi_{genjets}",100,-M_PI, M_PI);
  genjt_phi->Sumw2();

  genjt1_eta = fs->make<TH1F>("genjet1_eta","#eta_{genjets}",100,-2.5, 2.5);
  genjt1_eta->Sumw2();
  genjt1_phi = fs->make<TH1F>("genjet1_phi","#phi_{genjets}",100,-M_PI, M_PI);
  genjt1_phi->Sumw2();

  genjt2_eta = fs->make<TH1F>("genjet2_eta","#eta_{genjets}",100,-2.5, 2.5);
  genjt2_eta->Sumw2();
  genjt2_phi = fs->make<TH1F>("genjet2_phi","#phi_{genjets}",100,-M_PI, M_PI);
  genjt2_phi->Sumw2();

  genjt3_eta = fs->make<TH1F>("genjet3_eta","#eta_{genjets}",100,-2.5, 2.5);
  genjt3_eta->Sumw2();
  genjt3_phi = fs->make<TH1F>("genjet3_phi","#phi_{genjets}",100,-M_PI, M_PI);
  genjt3_phi->Sumw2();
  //genjt_oth_pt = fs->make<TH1F>("genjt_oth_pt","Et_{genjets_oth}",100, 20., 2020.);
  //genjt_oth_pt->Sumw2();
  /*  
  genjt_oth_eta = fs->make<TH1F>("genjt_oth_eta","#eta_{genjets_oth}",100,-5., 5.);
  genjt_oth_eta->Sumw2();
  genjt_oth_phi = fs->make<TH1F>("genjt_oth_phi","#phi_{genjets_oth}",100,-M_PI, M_PI);
  genjt_oth_phi->Sumw2();
  */
  //genchg_hist = fs->make<TH1F>("genchg_hist","# of genchargeds",120,-0.5, 239.5);
  //genchg_hist->Sumw2();
  genchg_pt = fs->make<TH1F>("genchg_pt","Et_{gencharge_alljet}",100, 1., 101.);
  genchg_pt->Sumw2();
  genchg_eta = fs->make<TH1F>("genchg_eta","#eta_{gencharge_alljet)",100,-3., 3.);
  genchg_eta->Sumw2();
  genchg_phi = fs->make<TH1F>("genchg_phi","#phi_{gencharge_alljet}",100,-M_PI, M_PI);
  genchg_phi->Sumw2();

  genchg1_pt = fs->make<TH1F>("genchg1_pt","Et_{gencharge_jet1}",100, 1., 101.);
  genchg1_pt->Sumw2();
  genchg1_eta = fs->make<TH1F>("genchg1_eta","#eta_{gencharge_jet1}",100,-3., 3.);
  genchg1_eta->Sumw2();
  genchg1_phi = fs->make<TH1F>("genchg1_phi","#phi_{gencharge_jet1}",100,-M_PI, M_PI);
  genchg1_phi->Sumw2(); 

  genchg2_pt = fs->make<TH1F>("genchg2_pt","Et_{gencharge_jet2}",100, 1., 101.);
  genchg2_pt->Sumw2();
  genchg2_eta = fs->make<TH1F>("genchg2_eta","#eta_{gencharge_jet2}",100,-3., 3.);
  genchg2_eta->Sumw2();
  genchg2_phi = fs->make<TH1F>("genchg2_phi","#phi_{gencharge_jet2}",100,-M_PI, M_PI);
  genchg2_phi->Sumw2();

  genchg3_pt = fs->make<TH1F>("genchg3_pt","Et_{gencharge_jet3}",100, 1., 101.);
  genchg3_pt->Sumw2();
  genchg3_eta = fs->make<TH1F>("genchg3_eta","#eta_{gencharge_jet3}",100,-3., 3.);
  genchg3_eta->Sumw2();
  genchg3_phi = fs->make<TH1F>("genchg3_phi","#phi_{gencharge_jet3}",100,-M_PI, M_PI);
  genchg3_phi->Sumw2();

  //genchg_oth_hist = fs->make<TH1F>("genchg_oth_hist","# of genchargeds (others)",120,-0.5, 239.5);
  //genchg_oth_hist->Sumw2();
  /* 
  genchg_oth_pt = fs->make<TH1F>("genchg_oth_pt","Et_{genchargeds_oth}",100,1., 101.);
  genchg_oth_pt->Sumw2();
  genchg_oth_eta = fs->make<TH1F>("genchg_oth_eta","#eta_{genchargeds_oth}",100,-5., 5.);
  genchg_oth_eta->Sumw2();
  genchg_oth_phi = fs->make<TH1F>("genchg_oth_phi","#phi_{genchargeds_oth}",100,-M_PI, M_PI);
  genchg_oth_phi->Sumw2();
  genneu_hist = fs->make<TH1F>("genneu_hist","# of genneutrals",120,-0.5, 239.5);
  genneu_hist->Sumw2();
  genneu_pt = fs->make<TH1F>("genneu_pt","Et_{genneutrals}",100,1., 101.);
  genneu_pt->Sumw2();
  genneu_eta = fs->make<TH1F>("genneu_eta","#eta_{genneutrals}",100,-3., 3.);
  genneu_eta->Sumw2();
  genneu_phi = fs->make<TH1F>("genneu_phi","#phi_{genneutrals}",100,-M_PI, M_PI);
  genneu_phi->Sumw2();

  genneu_oth_hist = fs->make<TH1F>("genneu_oth_hist","# of genneutrals (others)",120,-0.5, 239.5);
  genneu_oth_hist->Sumw2();
  genneu_oth_pt = fs->make<TH1F>("genneu_oth_pt","Et_{genneutrals_oth}",100, 1., 101.);
  genneu_oth_pt->Sumw2();
  genneu_oth_eta = fs->make<TH1F>("genneu_oth_eta","#eta_{genneutrals_oth}",100,-5., 5.);
  genneu_oth_eta->Sumw2();
  genneu_oth_phi = fs->make<TH1F>("genneu_oth_phi","#phi_{genneutrals_oth}",100,-M_PI, M_PI);
  genneu_oth_phi->Sumw2();
  */
  for (int ij=0; ij<nHLTmx; ij++) { 
    sprintf(name, "nprimall_%i", ij);
    sprintf(title, "# of primary vtx (%s)", (ij==0) ? "ALL" : jethlt_name[ij-1]);
    prim_hist[ij] = fs->make<TH1F>(name, title, 60, -0.5, 59.5);
    prim_hist[ij]->Sumw2();

    sprintf(name, "nprimsel_%i", ij);
    sprintf(title, "Selected # of primary vtx (%s)", (ij==0) ? "ALL" : jethlt_name[ij-1]);
    prim_sel[ij] = fs->make<TH1F>(name, title, 60, -0.5, 59.5);
    prim_sel[ij]->Sumw2();

    sprintf(name, "nprimall_rewt_%i", ij);
    sprintf(title, "# of rewighted primary vtx (%s)", (ij==0) ? "ALL" : jethlt_name[ij-1]);
    prim_hist_rewt[ij] = fs->make<TH1F>(name, title, 60, -0.5, 59.5);
    prim_hist_rewt[ij]->Sumw2();

    sprintf(name, "nprimsel_rewt_%i", ij);
    sprintf(title, "Selected # of reweighted primary vtx (%s)", (ij==0) ? "ALL" : jethlt_name[ij-1]);
    prim_sel_rewt[ij] = fs->make<TH1F>(name, title, 60, -0.5, 59.5);
    prim_sel_rewt[ij]->Sumw2();
  }

  prim_correl = fs->make<TH2F>("correl", "Correlation of all and Selected # of primary vtx", 60, -0.5, 59.5, 60, -0.5, 59.5);
  const char* namex[2]={"Selected", "Rejected"};
  for (int ij=0; ij<2; ij++) {
    sprintf(name, "primalltrk_%i", ij);
    sprintf(title, "All tracks in primary vtx (%s)", namex[ij]);
    prim_alltrk[ij] = fs->make<TH1F>(name, title, 240, -0.5, 239.5);

    sprintf(name, "primgoodtrk_%i", ij);
    sprintf(title, "Good tracks in primary vtx (%s)", namex[ij]);
    prim_goodtrk[ij] = fs->make<TH1F>(name, title, 240, -0.5, 239.5);

    sprintf(name, "primseltrk_%i", ij);
    sprintf(title, "Selected tracks in primary vtx (%s)", namex[ij]);
    prim_seltrk[ij] = fs->make<TH1F>(name, title, 240, -0.5, 239.5);

    sprintf(name, "primdx_%i", ij);
    sprintf(title, "#Delta x of prim wrt beam spot (%s)", namex[ij]);
    prim_dx[ij] = fs->make<TH1F>(name, title, 120, -2.4, 2.4);

    sprintf(name, "primdy_%i", ij);
    sprintf(title, "#Delta y of prim wrt beam spot (%s)", namex[ij]);
    prim_dy[ij] = fs->make<TH1F>(name, title, 120, -2.4, 2.4);

    sprintf(name, "primdxy_%i", ij);
    sprintf(title, "#Delta y vs #Delta x of prim (%s)", namex[ij]);
    prim_dxy[ij] = fs->make<TH2F>(name, title, 60, -0.15, 0.15, 60, -0.15, 0.15);

    sprintf(name, "primdz_%i", ij);
    sprintf(title, "#Delta z of prim wrt beam spo (%s)", namex[ij]);
    prim_dz[ij] = fs->make<TH1F>(name, title, 120, -30.0, 30.0); 

    sprintf(name, "primprob_%i", ij);
    sprintf(title, "log10(vertex fit prob) (%s)", namex[ij]);
    prim_prob[ij] = fs->make<TH1F>(name, title, 120, -20.0, 0.0);   
  }

  for(int ij=0; ij<njetetamn; ij++){
    sprintf(name, "njets_%i",ij);
    sprintf(title, "No of Jets_eta range_%gs", etarange[ij]);
    h_njets[ij] = fs->make<TH1F>(name, title, 9, 1, 10);
    h_njets[ij]->Sumw2();
  }

  for(int ij=0; ij<njetetamn; ij++){
    sprintf(name, "ncharges_%i",ij);
    sprintf(title, "No of charge particles_eta range_%gs", etarange[ij]);
    h_nchg[ij] = fs->make<TH1F>(name, title, 139, 1, 140);
    h_nchg[ij]->Sumw2();
  }


  for(int ij=0; ij<njetetamn; ij++){
    sprintf(name, "gennjets_%i",ij);
    sprintf(title, "No of GenJets_eta range_%gs", etarange[ij]);
    gen_njets[ij] = fs->make<TH1F>(name, title, 9, 1, 10);
    gen_njets[ij]->Sumw2();
  }
//-------------------------------------------
#ifdef TRIGGER
  const char* trigvar[2]={"L1", "HLT"};
  for(int ij=0; ij<nHLTmx; ij++){
    for(int jk=0; jk<2; jk++){
      sprintf(name, "trgjet_pt_%i_%i", ij, jk);
      sprintf(title, "trgjet_pt_%s_%s", jethlt_name[ij], trigvar[jk]);
      trgjet_pt[ij][jk] = fs->make<TH1F>(name, title, njetptbin, 20,1500);
      trgjet_pt[ij][jk]->Sumw2();

      sprintf(name, "trgjet_eta_%i_%i", ij, jk);
      sprintf(title, "trgjet_eta_%s_%s", jethlt_name[ij], trigvar[jk]);
      trgjet_eta[ij][jk] = fs->make<TH1F>(name, title, njetptbin, -5., 5.);
      trgjet_eta[ij][jk]->Sumw2();

      sprintf(name, "trgjet_phi_%i_%i", ij, jk);
      sprintf(title, "trgjet_phi_%s_%s", jethlt_name[ij], trigvar[jk]);
      trgjet_phi[ij][jk] = fs->make<TH1F>(name, title, 180,-M_PI, M_PI);
      trgjet_phi[ij][jk]->Sumw2();

      sprintf(name, "prbjet_pt_%i_%i", ij, jk);
      sprintf(title, "prbjet_pt_%s_%s", jethlt_name[ij], trigvar[jk]);
      prbjet_pt[ij][jk] = fs->make<TH1F>(name, title, njetptbin, 20,1500);
      prbjet_pt[ij][jk]->Sumw2();

      sprintf(name, "prbjet_eta_%i_%i", ij, jk);
      sprintf(title, "prbjet_eta_%s_%s", jethlt_name[ij], trigvar[jk]);
      prbjet_eta[ij][jk] = fs->make<TH1F>(name, title, 100,-5., 5.);
      prbjet_eta[ij][jk]->Sumw2();

      sprintf(name, "prbjet_phi_%i_%i", ij, jk);
      sprintf(title, "prbjet_phi_%s_%s", jethlt_name[ij], trigvar[jk]);
      prbjet_phi[ij][jk] = fs->make<TH1F>(name, title, 180,-M_PI, M_PI);
      prbjet_phi[ij][jk]->Sumw2();
    }
  } 
#endif
//Trigger special
//-------------------------------------------
	if (isReconstruct) { 
		for(int ij=0; ij<nHLTmx; ij++){
			for(int jk=0; jk<njetetamn; jk++){
				sprintf(name, "jetpt_%i_%i",jk,ij);
				sprintf(title, "jetpt_%s_%g", jethlt_name[ij], etarange[jk]);
				h_jetpt[ij][jk] = fs->make<TH1F>(name, title, 300, 50, 1550);
				h_jetpt[ij][jk]->Sumw2();
				
				sprintf(name, "jetphi_%i_%i",jk, ij);
				sprintf(title, "jetphi_%s_%g", jethlt_name[ij],etarange[jk]);
				h_jetphi[ij][jk] = fs->make<TH1F>(name, title, 180,-M_PI, M_PI);
				h_jetphi[ij][jk]->Sumw2();				
			}
		}
	}
#ifdef TRIGGER
  for(int ij=0; ij<nHLTmx; ij++){

    sprintf(name, "jeteta_%i", ij);
    sprintf(title, "jetphi_%s", jethlt_name[ij]);//, jetvar[ij]);
    h_jeteta[ij] = fs->make<TH1F>(name, title, 100, -5, 5);
    h_jeteta[ij]->Sumw2();

    for (int jk=0; jk<2; jk++){ 
      sprintf(name, "angle1d_%s_%i", jethlt_name[ij], jk);
      sprintf(title, "Angle%s_%i", jethlt_name[ij], jk);
      trgjet_angle[ij][jk] = fs->make<TH1F>(name, title, 90 , 0.1, 2.5);

      sprintf(name, "angle2d_%s_%i", jethlt_name[ij], jk);
      sprintf(title, "Angle_2d_hist%s_%i", jethlt_name[ij], jk);
      trgjet_2dangle[ij][jk] = fs->make<TH2F>(name, title, njetptbin, 20, 1500, 30 , 0.1, 2.5);
    }
  }

  for (int ij=0; ij<nHLTmx; ij++) {
    for (int jk=0; jk<njetetamn; jk++) {
      sprintf(name, "hlt_dijettag_%i_%i", ij, jk);
      sprintf(title, "dijet tagged P_T : (%s) |i#eta|<%g", jethlt_name[ij], etarange[jk]);
      hlt_dijettag[ij][jk] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
      hlt_dijettag[ij][jk]->Sumw2();

      sprintf(name, "hlt_dijetprob_%i_%i", ij, jk);
      sprintf(title, "dijet probed P_T : (%s) |i#eta|<%g", jethlt_name[ij], etarange[jk]);
      hlt_dijetprob[ij][jk] = fs->make<TH1F>(name, title, 60, 0.4*leadingPtThreshold[ij], 2.5*leadingPtThreshold[ij]);
      hlt_dijetprob[ij][jk]->Sumw2();
    }
  }
#endif
  counthist = fs->make<TH1F>("count","No of events",2,0,2); 

  for (int ix=0; ix<32; ix++) { mypow_2[ix] = pow(2,ix);}
  nevt = 0;
  //irun_old=-1;
  //trig_init=0;
  nreco=naa= nbb= ncc=0;
}

QCDEventShape::~QCDEventShape()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//member functions

// ------------ method called for each event  ------------
void QCDEventShape::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //t1=clock();
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  //gRandom->SetSeed(19919925);
  //float rn=gRandom->Uniform();
  //cout << " Random Number ini = " << rn << endl;
  //if (rn >0.90) return;
  //cout << " Random Number = " << rn << endl;
  //cout << "Time = " << t1 << "; " << t2 << endl;
  nevt++;
  //int ievt = iEvent.id().event();
  counthist->Fill(1); 
  //if (nevt%100==1)   std::cout<<"QCDEventShape::analyze "<< nevt<<" IRUN= "<<iEvent.id().run()<<" ievt= "<< iEvent.id().event()<<" "<<ievt<<endl;
  if (nevt%100==1)   std::cout<<"Jet Charge Analysis Run No. =  "<< nevt<<endl;
  //std::cout << "ok1"<<endl;

  //" ilumi" <<
  //iEvent.luminosityBlock() << " ibunch " << iEvent.bunchCrossing() <<std::endl;
  //cout << "NEvent = " <<  nevt << endl;
  //if(iEvent.luminosityBlock()==9881 || iEvent.luminosityBlock()==23185 || iEvent.luminosityBlock()==25334 || iEvent.luminosityBlock()== 26584 ||iEvent.luminosityBlock()== 35674 || iEvent.luminosityBlock()==32764 || iEvent.luminosityBlock()== 35675 || iEvent.luminosityBlock()==53681) return ;
  //if(iEvent.luminosityBlock()==2 || iEvent.luminosityBlock()==7175 || iEvent.luminosityBlock()==41151 || iEvent.luminosityBlock()==7389697 || iEvent.luminosityBlock()==60334 || iEvent.luminosityBlock()==51317 || iEvent.luminosityBlock()==53654 || iEvent.luminosityBlock()==10333 || iEvent.luminosityBlock()==54778 || iEvent.luminosityBlock()==10082 || iEvent.luminosityBlock()==54322 || iEvent.luminosityBlock()==64667 || iEvent.luminosityBlock()==65977 || iEvent.luminosityBlock()==55534 || iEvent.luminosityBlock()==55781 || iEvent.luminosityBlock()==55782 || iEvent.luminosityBlock()==55783 || iEvent.luminosityBlock()==61360 || iEvent.luminosityBlock()==61370 ||iEvent.luminosityBlock()==68258 || iEvent.luminosityBlock()==62147 || iEvent.luminosityBlock()==67194 || iEvent.luminosityBlock()==43070 || iEvent.luminosityBlock()==49429 || iEvent.luminosityBlock()==15102 || iEvent.luminosityBlock()==23306 || iEvent.luminosityBlock()==14242|| iEvent.luminosityBlock()==19080 || iEvent.luminosityBlock()==9312025) return;
  //npfjets = 0;
  //nchg =0;
  //if(iEvent.luminosityBlock()<4401) return; 
  //if(nevt<3442) return;
  //if(nevt!=3080) return;
  //cout << "Write test 1 = ok " << endl;
//-------------------------------------------
  vector<double> recovar;
  vector<double> recovar1;
  std::vector<HepLorentzVector> recomom[njecmx][ntype][njetetamn];
  std::vector<HepLorentzVector> tmpjt4v; 
  //std::vector<HepLorentzVector> tmpcand4v;                
  std::vector<HepLorentzVector> tmpgen4v;                                                                                             
  
  std::vector<HepLorentzVector> genmom[nGenReso][ntype][njetetamn];
  vector<double> genvar;
//-------------------------------------------
  //vector<double> ijet1candsmom[nkappa];
  //vector<double> ijet1candsmom;

  double ijet1candsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
///
  double ijet1gcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1g_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1g_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1g_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1g_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qgcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qg_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1qg_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qg_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1qg_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1agcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1ag_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1ag_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1ag_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1ag_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
///
  double ijet1bcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1b_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1b_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1b_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1b_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qbcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qb_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1qb_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qb_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1qb_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1abcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1ab_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1ab_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1ab_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1ab_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
///
  double ijet1ccandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1c_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1c_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1c_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1c_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qccandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qc_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1qc_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qc_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1qc_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1accandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1ac_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1ac_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1ac_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1ac_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};


///
  double ijet1scandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1s_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1s_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1s_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1s_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qscandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qs_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1qs_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qs_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1qs_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1ascandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1as_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1as_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1as_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1as_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
///
  double ijet1ucandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1u_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1u_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1u_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1u_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qucandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qu_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1qu_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qu_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1qu_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1aucandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1au_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1au_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1au_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1au_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

///
  double ijet1dcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1d_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1d_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1d_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1d_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qdcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qd_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1qd_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1qd_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1qd_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1adcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1ad_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1ad_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet1ad_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet1ad_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
/////////////////

  double ijet2candsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
///
  double ijet2gcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2g_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2g_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2g_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2g_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qgcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qg_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2qg_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qg_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2qg_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2agcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2ag_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2ag_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2ag_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2ag_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

///
  double ijet2bcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2b_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2b_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2b_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2b_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qbcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qb_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2qb_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qb_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2qb_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  
  double ijet2abcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2ab_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2ab_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2ab_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2ab_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
///
  double ijet2ccandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2c_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2c_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2c_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2c_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qccandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qc_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2qc_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qc_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2qc_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2accandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2ac_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2ac_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2ac_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2ac_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

///
  double ijet2scandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2s_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2s_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2s_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2s_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qscandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qs_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2qs_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qs_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2qs_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2ascandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2as_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2as_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2as_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2as_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

///
  double ijet2ucandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2u_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2u_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2u_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2u_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qucandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qu_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2qu_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qu_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2qu_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2aucandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2au_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2au_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2au_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2au_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

///
  double ijet2dcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2d_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2d_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2d_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2d_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qdcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qd_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2qd_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2qd_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2qd_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2adcandsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2ad_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2ad_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double ijet2ad_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double ijet2ad_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

/*
  double ijet1candsmom[nkappa];

  double ijet1_long_num[nkappa];
  double ijet1_long_den[nkappa];

  double ijet1_tran_num[nkappa];
  double ijet1_tran_den[nkappa];

  double ijet2candsmom[nkappa];

  double ijet2_long_num[nkappa];
  double ijet2_long_den[nkappa];

  double ijet2_tran_num[nkappa];
  double ijet2_tran_den[nkappa];
*/

  double igenjet1candsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double igenjet1_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double igenjet1_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double igenjet1_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double igenjet1_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double igenjet2candsmom[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double igenjet2_long_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double igenjet2_long_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  double igenjet2_tran_num[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double igenjet2_tran_den[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//-------------------------------------------
  wtfact=1.0;
  //double px=0;
  //double py=0;
  //double ptxy=0;

  //int ncount=0;
  unsigned ncount=0;
  //double recterm=0;
  //int ithird=-1;
  int irecoht=-1;
	//#ifdef JETENERGY
	int irecohtjec[njecmx];
	for (int ij=0; ij<njecmx; ij++) { irecohtjec[ij]=-1;}
	//#endif	
  double aveleadingptjec[njecmx] ={0};//14Sep20
  double leadingptjec[njecmx] ={0};   // jet charge 

  int igenht=-1;
	//#ifdef  JETRESO
	int igenhtres[nGenReso];
	for (int ij=0; ij<nGenReso; ij++) { igenhtres[ij]=-1;}
	//#endif
        double avegenptres[nGenReso]={0};//14Sep20
	double leadgenptres[nGenReso]={0};// jet charge

#ifdef TRIGGER
  const char* variab1;
#endif
#ifndef DIJETAVE
  const char* variab2; 
#endif
  if (isMC) {
#ifdef LHAPDF
    edm::Handle<LHEEventProduct> EvtHandle ;
    iEvent.getByToken( lheEventProductToken_ , EvtHandle ) ;
		
		for ( unsigned int weightIndex = 0; weightIndex < EvtHandle->weights().size(); ++weightIndex ) {
			//cout<< EvtHandle->weights()[weightIndex].wgt <<endl;
                        //systematicWeightIDs->push_back( atoi(EvtHandle->weights()[weightIndex].id.c_str()) );
			if (weightIndex>=9 && weightIndex<=109) {
				pdfwt[weightIndex-9] = EvtHandle->weights()[weightIndex].wgt/EvtHandle->originalXWGTUP(); 
				//std::cout << weightIndex << " " << EvtHandle->weights()[weightIndex].id << " " << EvtHandle->weights()[weightIndex].wgt <<" "<<pdfwt[weightIndex-9]<< std::endl;
				}
    			}
#endif
    //cout<<"AAAAAAAAAAAAAA"<<endl;
    edm::Handle<GenEventInfoProduct> eventinfo;
    iEvent.getByToken(generator1_, eventinfo);
    if (eventinfo.isValid()) { 
      qscale = eventinfo->qScale(); 
      wtfact = eventinfo->weight();
      //weight = weight2*wtfact;
      procid = eventinfo->signalProcessID();
      //cout << " qscale = " <<setw(14)<< qscale << " ; wtfact = " << wtfact << " ; procid = " << procid  << endl;

      if (eventinfo->hasPDF()) {
	const gen::PdfInfo* xpdf = eventinfo->pdf();
	
	int id1 = xpdf->id.first;
	int id2 = xpdf->id.second;
	
	idall = 100*(id1+50)+ (id2+50); 
	
	qscale = xpdf->scalePDF;
	
	xfrac1 = xpdf->x.first;
	xfrac2 = xpdf->x.second;
	
	xpdf1 = xfrac1*xpdf->xPDF.first;
	xpdf2 = xfrac2*xpdf->xPDF.second; 
      }
    }
  }
 
/*
edm::Handle< double > theprefweight;
iEvent.getByToken(prefweight_token, theprefweight ) ;
double _prefiringweight =(*theprefweight);

edm::Handle< double > theprefweightup;
iEvent.getByToken(prefweightup_token, theprefweightup ) ;
double _prefiringweightup =(*theprefweightup);

edm::Handle< double > theprefweightdown;
iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
double _prefiringweightdown =(*theprefweightdown);
*/
//cout << "Prefire  wt : " << _prefiringweight << endl;
//cout << "Prefire up wt : " << _prefiringweightup << endl;
//cout << "Prefire down wt : " << _prefiringweightdown << endl;
 
#ifdef TRIGGER
  edm::Handle<edm::TriggerResults> trigRes;
  iEvent.getByToken(triggerBits_, trigRes);
  
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
//-------------------------------------------Trigger
  const edm::TriggerNames &names = iEvent.triggerNames(*trigRes);
  //int ihltfill = -1;
#endif
  tmpjt4v.clear();
  //tmpcand4v.clear();
  tmpgen4v.clear();

  double aveleadingpt =0;  //ESVs (dijet)
  double leadingpt = 0;    //Jet charge (single jet)
  bool isInEtaRange[njetetamn]={0}; //GMA{0,0,0,0};
  
  recojet1_pt = 0.0; 
  recojet1g_pt = 0.0, recojet1qg_pt = 0.0, recojet1ag_pt = 0.0;
  recojet1b_pt = 0.0, recojet1qb_pt = 0.0, recojet1ab_pt = 0.0;
  recojet1c_pt = 0.0, recojet1qc_pt = 0.0, recojet1ac_pt = 0.0;
  recojet1s_pt = 0.0, recojet1qs_pt = 0.0, recojet1as_pt = 0.0;
  recojet1u_pt = 0.0, recojet1qu_pt = 0.0, recojet1au_pt = 0.0;
  recojet1d_pt = 0.0, recojet1qd_pt = 0.0, recojet1ad_pt = 0.0;
  
  recojet2_pt = 0.0; 
  recojet2g_pt = 0.0, recojet2qg_pt = 0.0, recojet2ag_pt = 0.0;
  recojet2b_pt = 0.0, recojet2qb_pt = 0.0, recojet2ab_pt = 0.0;
  recojet2c_pt = 0.0, recojet2qc_pt = 0.0, recojet2ac_pt = 0.0;
  recojet2s_pt = 0.0, recojet2qs_pt = 0.0, recojet2as_pt = 0.0;
  recojet2u_pt = 0.0, recojet2qu_pt = 0.0, recojet2au_pt = 0.0;
  recojet2d_pt = 0.0, recojet2qd_pt = 0.0, recojet2ad_pt = 0.0;
 
  genrecojet1_pt = 0.0, genrecojet1g_pt = 0.0, genrecojet1b_pt = 0.0, genrecojet1c_pt = 0.0, genrecojet1s_pt = 0.0, genrecojet1u_pt = 0.0, genrecojet1d_pt = 0.0;
  genrecojet2_pt = 0.0, genrecojet2g_pt = 0.0, genrecojet2b_pt = 0.0, genrecojet2c_pt = 0.0, genrecojet2s_pt = 0.0, genrecojet2u_pt = 0.0, genrecojet2d_pt = 0.0;
#ifndef GENPART
  edm::Handle<pat::JetCollection> ak4PFJets;
  if (isReconstruct) { 
    iEvent.getByToken(jetSrcToken_, ak4PFJets);
  }
  //cout<<"1 aveleadingpt"<<endl;
  if (isReconstruct && ((!ak4PFJets.isValid()) || ak4PFJets->size() <2)) return; //GMA, do we use this
  
  if (ak4PFJets.isValid() && ak4PFJets->size()>=2) {
#ifdef DIJETAVE
    //aveleadingpt = 0.5*((*ak4PFJets)[0].pt() + (*ak4PFJets)[1].pt());
    //cout<<"1 aveleadingpt"<<aveleadingpt<<endl;

    for (int iet=0; iet<njetetamn; iet++) {
      isInEtaRange[iet] = true;
      }
    
    for (int ij=0; ij<2; ij++) { 
      for (int iet=0; iet<njetetamn; iet++) {
	if (abs((*ak4PFJets)[ij].eta())>etarange[iet]) { isInEtaRange[iet] = false;}
      }
      
      //Jet ID ================= 2017 & 2018 UL jetID recomendation 
      double NHF = (*ak4PFJets)[ij].neutralHadronEnergyFraction();
      double NEMF = (*ak4PFJets)[ij].neutralEmEnergyFraction();
      double CHF = (*ak4PFJets)[ij].chargedHadronEnergyFraction();
      //double MUF = (*ak4PFJets)[ij].muonEnergyFraction();
      //double CEMF = (*ak4PFJets)[ij].chargedEmEnergyFraction();
      int NumConst = (*ak4PFJets)[ij].chargedMultiplicity()+(*ak4PFJets)[ij].neutralMultiplicity();
      //int NumNeutralParticles =(*ak4PFJets)[ij].neutralMultiplicity();
      int CHM = (*ak4PFJets)[ij].chargedMultiplicity();

      bool TightJetID =false;
      if(abs((*ak4PFJets)[ij].eta())<=2.7){
      if(NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 && abs((*ak4PFJets)[ij].eta())<=2.4 )  TightJetID =true;
      if(NHF<0.90 && NEMF<0.99 && abs((*ak4PFJets)[ij].eta())>2.4 )  TightJetID =true;} 
      else {TightJetID =false;}
      if (abs((*ak4PFJets)[ij].eta())>2.7) {TightJetID = false;}  //2.5 or 2.6
      if ((*ak4PFJets)[ij].pt()<30.0) {TightJetID = false;}

      if (TightJetID) { aveleadingpt +=(*ak4PFJets)[ij].pt();
			leadingpt = (*ak4PFJets)[0].pt();
      //std::cout<<"ok2"<<endl;
      //std::cout << "Leading jet pt :"<<leadingpt<<endl;
      //std::cout << "Avg. Leading jet pt :"<<aveleadingpt<<endl;
			} else {aveleadingpt -=100000;
			        leadingpt -=100000;}
    }
    //std::cout << "Leading jet pt at end of the loop :"<<leadingpt<<endl;
    //std::cout << "Avg. Leading jet pt at end of the loop :"<<aveleadingpt<<endl;
    aveleadingpt /=2.0;
    //leadingpt /=2.0;
#else

#endif
  }//if (ak4PFJets.isValid() && ak4PFJets->size()>=2) {
#endif
  //cout<<"leadingpt:"<<leadingpt<<endl;
  //cout<<"aveleadingpt:"<<aveleadingpt<<endl;
  //if (isReconstruct && isMC && aveleadingpt>3*qscale) return;
  if (isReconstruct && isMC && leadingpt>3*qscale) return;
  //irecoht = getbinid(aveleadingpt, nHLTmx, leadingPtThreshold);
  irecoht = getbinid(leadingpt, nHLTmx, leadingPtThreshold);

#ifdef TRIGGER
  bool trgpas[nHLTmx]={0,0,0,0,0,0,0,0,0};  //8 or 10

  //  if (isMC && ak4PFJets.isValid() &&  ak4PFJets->size()>=2) {
  //    aveleadingpt = (*ak4PFJets)[0].pt();
  //    irecoht = getbinid(aveleadingpt, nHLTmx, leadingPtThreshold);
  //  }

  //cout<<"ave"<<aveleadingpt<<" ; Jet1 Pt= " <<   (*ak4PFJets)[0].pt() <<" ; Jet2 Pt= " <<  (*ak4PFJets)[1].pt()<<endl; 
  //std::pair<std::vector<std::pair<std::string,int> >,int> prescaleValuesInDetail(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& trigger) const;
/* 
  //Calcualte Trigger Efficiency for dijet events Manas Sir
 bool trg_fired1=false;
  for (int jk=0; jk<nHLTmx-1; jk++) {
      bool trg_fired=false;
    for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
      std::string name = names.triggerName(ij);
      variab1 = name.c_str();
      //bool trg_fired=false;
      for (int iet=0; iet<njetetamn; iet++) {
	if(aveleadingpt*2>jethlt_thr[jk] && (strstr(variab1,jethlt_name[jk]) && (strlen(variab1)-strlen(jethlt_name[jk])<5)) ){ 
	//if(aveleadingpt*2>jethlt_thr[jk] && aveleadingpt*2<jethlt_thr[jk] && (strstr(variab1,jethlt_name[jk]) && (strlen(variab1)-strlen(jethlt_name[jk])<5)) ){ 
	  //for (int iet=0; iet<njetetamn; iet++) {
	  if (isInEtaRange[iet]) {
	    if (trigRes->accept(ij)) {
                hlt_dijettag[jk][iet]->Fill(aveleadingpt);
	        trg_fired=true;
               } 
	    if(trg_fired) trg_fired1=true;
            else trg_fired1=false;
          if(iet==0)  cout << jk << " ; ij" <<ij <<" ; Prb name "  << variab1 << endl;
	 // }
	//}
	
	if(trg_fired1 && (strstr(variab1,jethlt_name[jk+1]) && (strlen(variab1)-strlen(jethlt_name[jk+1])<5))){
	  if (trigRes->accept(ij) && isInEtaRange[iet]) {hlt_dijetprob[jk][iet]->Fill(aveleadingpt);}
         if(iet==0) cout << jk <<" Tag name "  << variab1 << endl;
	}
       }
      } //for (int iet=0; iet<njetetamn; iet++) {
    }
 }
}
*/
//Calcualte Trigger Efficiency for dijet events Test
/*      
 bool trg_fired1=false; 
  for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
    std::string name = names.triggerName(ij);
    variab1 = name.c_str();
//   std::cout << "Trigger " << names.triggerName(ij) << 
  //              ", prescale " << triggerPrescales->getPrescaleForIndex(ij) <<
    //            ": " << (trigRes->accept(ij) ? "PASS" : "fail (or not run)") 
      //          << std::endl;
    bool trg_fired=false;
    for (int iet=0; iet<njetetamn-2; iet++) {
      if (isInEtaRange[iet]) {
	if((strstr(variab1,jethlt_name[2]) && (strlen(variab1)-strlen(jethlt_name[2])<5))){
	  if (trigRes->accept(ij)) {
	  //if (isL3) {
	    hlt_dijettag[2][iet]->Fill(aveleadingpt);
	    trg_fired=true;
	  if(iet==0) cout <<" Tag name "  << variab1 << endl;
	  }
	  if(trg_fired) trg_fired1=true;
	  else trg_fired1=false;
	}
	if(trg_fired1 && (strstr(variab1,jethlt_name[3]) && (strlen(variab1)-strlen(jethlt_name[3])<5))){
	  //if (trigRes->accept(ij) && isInEtaRange[iet]) {hlt_dijetprob[2][iet]->Fill(aveleadingpt);}
	   if (trigRes->accept(ij)) {
          hlt_dijetprob[2][iet]->Fill(aveleadingpt);
	  if(iet==0) cout <<"Prob name "  << variab1 << endl;
         }
	}   
      }    
    }
  }
*/

  //Calcualte Trigger Efficiency for dijet events
  bool trg_prev=false;

  //if (!isMC) {
  for (int jk=-1; jk<nHLTmx; jk++) {
    for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
      std::string name = names.triggerName(ij);
      variab1 = name.c_str(); 
      if ((jk<0 && strstr(variab1,jethlt_lowest) && strlen(variab1)-strlen(jethlt_lowest)<5) || 
	  (jk>=0 && strstr(variab1,jethlt_name[jk]) && strlen(variab1)-strlen(jethlt_name[jk])<5)) {
	
	 //const std::pair<std::vector<std::pair<std::string,int> >,int> prescalesInDetail(hltConfig_.prescaleValuesInDetail(iEvent,iSetup, variab1));
	 const std::pair<std::vector<std::pair<std::string,int> >,int> prescalesInDetail(hltPrescaleProvider_.prescaleValuesInDetail(iEvent,iSetup,variab1));
	 if (jk>=0) { 
          //cout<<variab1<<endl;
	  //==============================================================================
	  //double tmpp1= prescalesInDetail.first[0].second;
	  //double tmpp2 = prescalesInDetail.first[1].second;
	  //l1pres[jk] =min(tmpp1, tmpp2);
	  //==============================================================================
	  l1pres[jk] = prescalesInDetail.first[0].second;
	  
	 //if (jk>=3 && l1pres[jk]>1) { l1pres[jk]=1.0;}
	 if(l1pres[jk]<=0){l1pres[jk]=1.0;}
         hltpres[jk] = prescalesInDetail.second;	  

	 //compres[jk] = (l1pres[jk])*(triggerPrescales->getPrescaleForIndex(ij)); 
	 //compres[jk] = triggerPrescales->getPrescaleForIndex(ij);
	 compres[jk] = (l1pres[jk])*(hltpres[jk]);
         //if (jk==9){compres[jk]=1.0;} // fix for 2017E sample pt spike
	         	//cout<<"Print Trigger : "<<variab1<< ":"<< compres[jk]<<endl;} 
	 //cout << "Prescale" << "Path " << variab1 <<" " <<   compres[jk] << endl;
	 //cout<<"Run NO= "<< iEvent.id().run()<<" ; Event No = "<< iEvent.id().event()<< " ; ilumi = " << iEvent.luminosityBlock() << 
	 //	" ; ibunch = " << iEvent.bunchCrossing()<<" ; L1 Pres0 = " << l1pres[jk] <<" "<<
	 //            " ; HLT Path= "<<name <<" ; HLT Pres = " <<hltpres[jk]<<" ; compres ="<<compres[jk] <<"; irecoht = "<< irecoht <<"; Pt=" <<aveleadingpt<<endl;
	 if (trigRes->accept(ij)) {trgpas[jk] = true;
		//cout << "Prescale" << "Path " << variab1 <<" " <<   compres[jk] << endl;} // ihltfill = jk;}
	 }
	 //if (trg_prev && compres[jk]>0.99) {
	 if (trg_prev){
	    for (int iet=0; iet<njetetamn; iet++) {
	      if (isInEtaRange[iet]) { 
		//hlt_dijettag[jk][iet]->Fill(aveleadingpt,compres[jk]);
                hlt_dijettag[jk][iet]->Fill(leadingpt,compres[jk]);
		//if (trigRes->accept(ij)) {hlt_dijetprob[jk][iet]->Fill(aveleadingpt, compres[jk]);} //{, (isMC) ? 1.0 : compres[jk]);}
                if (trigRes->accept(ij)) {hlt_dijetprob[jk][iet]->Fill(leadingpt, compres[jk]);}
	      }
	    }
	  }
	  /*
 	    for (int iet=0; iet<njetetamn; iet++) {
	    if (isInEtaRange[iet]) { 
	    if(trg_prev) hlt_dijettag[jk][iet]->Fill(aveleadingpt);
	    if (trg_prev && trigRes->accept(ij)) {hlt_dijetprob[jk][iet]->Fill(aveleadingpt);} 
	    }    
	    }*/
	  //if (trg_prev) cout << "Accept =" << " name = " <<name <<endl;
	  trg_prev = trigRes->accept(ij);
	  //trg_prev = trg_prev|trigRes->accept(ij);
	  //if (!trg_prev) { trg_prev = trigRes->accept(ij);}
	  break;
	} else {
	  trg_prev = trigRes->accept(ij);
	  break;
	}
      }
    }
  }
#endif
  //cout<<"ihltfill "<<ihltfill<<endl;
  //cout<<"3 aveleadingpt "<<aveleadingpt<< " ; "<<ihltfill<<" "<<irecoht<<endl;
  
  //if ((irecoht <0 || irecoht >=nHLTmx) || ((!isMC) && (!trgpas[irecoht]))) return; //GMA remopve this condition
  //cout <<"irecoht = "<<irecoht<<endl;
  //if (irecoht==-3) return;
#ifdef TRIGGER
  if (irecoht>=0 && ((!isMC) && (!trgpas[irecoht]))) return;
  if (irecoht==-2 && ((!isMC) && (!trgpas[0]))) return;
#endif
  
  if (!isMC) {
    if (irecoht>=0) {
      wtfact = compres[irecoht];
    } else if (irecoht==-2) {
      wtfact = compres[0];
    } else {
      return ;
    }
  }
  //for (int ij=0; ij<nHLTmx; ij++) {lumiwt[ij]=intlumi[nHLTmx-1]/intlumi[ij];}// cout<<"nt "<<datpileup[ij][0]<<endl;}
  if (isMC) {
#ifndef GENPART
    //Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    //iEvent.getByLabel("addPileupInfo", PupInfo);
    edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(pileup_, PupInfo);
    int npu = -1;
    //int tnpv  = -1;
    PShiftDown_ = reweight::PoissonMeanShifter(-0.5);
    PShiftUp_ = reweight::PoissonMeanShifter(0.5);
    if (PupInfo.isValid()) {
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
	if (PVI->getBunchCrossing()==0) {
          //npu = PVI->getPU_NumInteractions();
	  npu = PVI->getTrueNumInteractions();
          //tnpv  = PVI->getTrueNumInteractions();
	  break;
	}
      }
    }
    //double MyWeight = LumiWeights_->weight(npu);
    
    //cout << "Main weight = " <<MyWeight << endl;
    //double TotalWeight_plus = MyWeight*PShiftUp_.ShiftWeight( npu );
    //double TotalWeight_plus = PShiftUp_.ShiftWeight( npu );
    //double TotalWeight_minus = PShiftDown_.ShiftWeight( npu ); 

    //cout << "Plus " << wtfact*TotalWeight_plus << " Mi = " << endl;
    //cout << "wt= " <<  wtfact << " : weightmi" <<wtfact*TotalWeight_minus << " Mi = " << endl;
    //wtfact=wtfact*TotalWeight_plus; 
    //wtfact=wtfact-TotalWeight_minus; 
    //cout << "npu Number of interactions : " << npu << endl; 
    //cout << "tnpv Number of true interactions : " << tnpv << endl; 
    if (npu<0) return; //GMA  
    if (isFlat) {
      weight =weight2*wtfact; // for flat MC sample
    } else {
      weight =weight2;
    }
#endif
    defweight = weight;
//-------------------------------------------
#ifndef GENPART
    int tmprecht = (irecoht>=0) ? irecoht : 0; //GMA
    
    if (npu<npileupmx) {
          weight *=rat_pileup[tmprecht][npu]; //GMA
    } else {
            weight *=rat_pileup[tmprecht][npileupmx-1]; //GMA
    }
#endif
    
    weighttrg = weight;
    //cout <<"weight  "<<weight<<" "<< weight2<<endl;
    //sar 3D PU reweighting 111028
  } else {
    weight = weight2;
    defweight = weight2;
    weighttrg = weight*wtfact; // *lumiwt[irecoht];
    //weighttrg = weight*lumiwt[3];
    //cout <<"TEST2  weighttrg "<< weighttrg<<" ; weight "<<weight<<" ; "<< wtfact<<endl;
  }

/*
#ifdef PREFIRE
double tmpwt = weighttrg;
weighttrg = tmpwt*_prefiringweight;
//weighttrg = tmpwt*_prefiringweightup;
//weighttrg = tmpwt*_prefiringweightdown;
#endif
*/

//---------------Lumiweight------------------
/*
//cout << "before lumi " << weighttrg ;
#ifdef LUMIWEIGHT
double tmpwt = weighttrg;
weighttrg = tmpwt*lumiwtt;
#endif
//cout << "   After lumi " << weighttrg << endl;
*/

//-------------------------------------------
#ifndef GENPART
  if(!isMC){
    reco::TrackBase::Point beamPoint(0,0, 0);
    //math::XYZPoint beamPoint(0,0, 0); 
    
    edm::Handle<reco::BeamSpot> beamSpotH;
    iEvent.getByToken(beamSpot_,beamSpotH);
    if (beamSpotH.isValid()){
      beamPoint = beamSpotH->position();
    }
    
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(vtxToken_, primaryVertices);
    
    int tmpvert=0;
    nprim=0;
    if (primaryVertices.isValid()) {
      tmpvert = primaryVertices->size();
      //cout<<"temp"<<tmpvert<<endl;
      for (reco::VertexCollection::const_iterator vert=primaryVertices->begin(); vert<primaryVertices->end(); vert++) {
	int isel = (vert->isValid() && !vert->isFake()) ? 1 : 0;
	int ngoodtrk = 0;
	int nseltrk = 0;
	double prob = ChiSquaredProbability(vert->chi2(),vert->ndof());
	for (reco::Vertex::trackRef_iterator reftrk =vert->tracks_begin(); reftrk<vert->tracks_end(); reftrk++) {
	  if ((*reftrk)->quality(TrackBase::highPurity) && vert->trackWeight(*reftrk)>0) {
	    ngoodtrk++; 
	    if ((*reftrk)->normalizedChi2()<100000 && 
		abs((*reftrk)->dxy()) < 10000 && 
		(*reftrk)->pt() >0.50) {nseltrk++; } 
	  }
	}
	prim_alltrk[isel]->Fill(vert->tracksSize());
	prim_goodtrk[isel]->Fill(ngoodtrk);
	prim_seltrk[isel]->Fill(nseltrk);
	prim_dx[isel]->Fill(vert->position().x() - beamPoint.x());
	prim_dy[isel]->Fill(vert->position().y() - beamPoint.y());
	prim_dxy[isel]->Fill(vert->position().x() - beamPoint.x(), vert->position().y() - beamPoint.y());
	prim_dz[isel]->Fill(vert->position().z() - beamPoint.z());
	prim_prob[isel]->Fill(max(-20.0, log10(prob)));
	
	if (isel==1 && nprim < nprimx-1) {
	  primpr[nprim] = prob;
	  ntkpm[nprim] = 1000*(1000*min(int(vert->tracksSize()),999) + min(ngoodtrk,999)) + min(999, nseltrk);
	  nprim++;
	}
      }
    }
    
    prim_hist[0]->Fill(tmpvert);
    prim_sel[0]->Fill(nprim);
    
    prim_hist_rewt[0]->Fill(tmpvert, weighttrg);
    prim_sel_rewt[0]->Fill(nprim, weighttrg);

    if (irecoht>=0 && irecoht<nHLTmx) { 
      prim_hist[irecoht]->Fill(tmpvert);
      prim_sel[irecoht]->Fill(nprim);
      
      prim_hist_rewt[irecoht]->Fill(tmpvert, weighttrg);
      prim_sel_rewt[irecoht]->Fill(nprim, weighttrg);   
    }
    prim_correl->Fill(tmpvert, nprim);
  } 
#endif 
  //cout<<"2 aveleadingpt "<<aveleadingpt<< " ; "<<ihltfill<<" "<<irecoht<<endl; 
  
  vector<double> jetptx[njecmx];
  vector<double> jetscl[njecmx];
  vector<int> jetindx[njecmx];

#ifndef GENPART
  if (ak4PFJets.isValid()) { 
    for(unsigned ijet = 0; ijet != ak4PFJets->size(); ijet++) {
      double pt = (*ak4PFJets)[ijet].pt();
      
      //#ifndef JETENERGY
      //#ifdef JETRESO
      
#if defined(JETRESO)&&(!defined(JETENERGY))
      //resolution file 
      JME::JetResolution resolution;
      //resolution = JME::JetResolution("Summer20UL16APV_JRV3_DATA_PtResolution_AK4PFchs.txt");    // for DATA APV
      //resolution = JME::JetResolution("Summer20UL16_JRV3_DATA_PtResolution_AK4PFchs.txt");       // for DATA
      resolution = JME::JetResolution("Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.txt");      // for MC APV
      //resolution = JME::JetResolution("Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt");         // for MC
 
      //Scalefactor file
      JME::JetResolutionScaleFactor res_sf;
      //cout<<"Filename="<<scalefile<<endl;
      //res_sf = JME::JetResolutionScaleFactor("Summer20UL16APV_JRV3_DATA_SF_AK4PFchs.txt");       // for DATA APV
      //res_sf = JME::JetResolutionScaleFactor("Summer20UL16_JRV3_DATA_SF_AK4PFchs.txt");          // for DATA
      res_sf = JME::JetResolutionScaleFactor("Summer20UL16APV_JRV3_MC_SF_AK4PFchs.txt");         // for MC APV
      //res_sf = JME::JetResolutionScaleFactor("Summer20UL16_JRV3_MC_SF_AK4PFchs.txt");            // for MC
      
      edm::Handle<double> rho;
      iEvent.getByToken(m_rho_token, rho);
      //cout<< "  rho=" << *rho << endl;
      
      //cout << "Write test 3 = ok " << endl;
      double eta = (*ak4PFJets)[ijet].eta();
      double reso = 1;
      JME::JetParameters parameters_5 = {{JME::Binning::JetPt, pt}, {JME::Binning::JetEta, eta}, {JME::Binning::Rho, *rho}};
      float rp = resolution.getResolution(parameters_5);
      float sf = res_sf.getScaleFactor({{JME::Binning::JetEta, eta}});
      float sf_up= res_sf.getScaleFactor({{JME::Binning::JetEta, eta}}, Variation::UP);
      float sf_dn= res_sf.getScaleFactor({{JME::Binning::JetEta, eta}}, Variation::DOWN);
      //#endif
      //#endif
#endif 
      for (int isrc = 0; isrc < njecmx; isrc++) {
	double sup = 1;
#ifdef JETENERGY
	double eta = (*ak4PFJets)[ijet].eta();
	if (isrc>0 && isrc<=nsrc) {
	  JetCorrectionUncertainty *jecUnc = vsrc[isrc-1];
	  jecUnc->setJetEta(eta);
	  jecUnc->setJetPt(pt);
	  
	  sup += jecUnc->getUncertainty(true);
	} else if (isrc>nsrc) {
	  JetCorrectionUncertainty *jecUnc = vsrc[isrc-nsrc-1];
	  jecUnc->setJetEta(eta);
	  jecUnc->setJetPt(pt);
	  sup -= jecUnc->getUncertainty(false);
	}
#elif defined(JETRESO)
	if (isrc==0) {  
	  reso = sqrt(sf*sf - 1)*rp;
	} else if (isrc==1) {
	  reso = sqrt(sf_up*sf_up - 1)*rp;
	} else if (isrc==2) {
	  reso = sqrt(sf_dn*sf_dn - 1)*rp;
	}
	sup = gRandom->Gaus(1.0, reso);			
#endif
	jetptx[isrc].push_back(sup*pt);
	jetscl[isrc].push_back(sup);
	jetindx[isrc].push_back(ijet);
      }
    }
    //#if defined(JETENERGY)||defined(JETRESO)
    
    for (int isrc = 0; isrc < njecmx; isrc++) {
      for (unsigned int ij=0; ij<jetptx[isrc].size()-1; ij++) {
	for (unsigned int jk=ij+1; jk<jetptx[isrc].size(); jk++) {
	  if (jetptx[isrc][jk]>jetptx[isrc][ij]) {
	    double tmppt = jetptx[isrc][ij];
	    double tmpscl = jetscl[isrc][ij];
	    int tmpindx = jetindx[isrc][ij];
	    
	    jetptx[isrc][ij] = jetptx[isrc][jk];
	    jetscl[isrc][ij] = jetscl[isrc][jk];
	    jetindx[isrc][ij] = jetindx[isrc][jk];					
	    
	    jetptx[isrc][jk] = tmppt;
	    jetscl[isrc][jk] = tmpscl;
	    jetindx[isrc][jk] = tmpindx;
	  }
	}
      }
    }
    //#endif
    
    //cout<<"1 aveleadingpt "<<endl; //aveleadingpt<< " ; "<<ihltfill<<" "<<irecoht<<endl;
    //double aveleadingptjec[njecmx] ={0};
    for (int isrc = 0; isrc < njecmx; isrc++) {
      if (jetptx[isrc].size()>=2) {
	//aveleadingptjec[isrc] = 0.5*(jetptx[isrc][0] + jetptx[isrc][1]);
	leadingptjec[isrc] = jetptx[isrc][0];
	//irecohtjec[isrc] = getbinid(aveleadingptjec[isrc], nHLTmx, leadingPtThreshold);
        irecohtjec[isrc] = getbinid(leadingptjec[isrc], nHLTmx, leadingPtThreshold);
      } else {
	irecohtjec[isrc] = -1;
      }
    }

    //GMA Need the corection on aveleadingpt
    if (ak4PFJets.isValid() && ak4PFJets->size() >=2) { //  && aveleadingpt >leadingPtThreshold[0]) { //GMA look on this
      
      for (int iet=0; iet<njetetamn; iet++) {
	for (int isrc = 0; isrc < njecmx; isrc++) {
	  //if (aveleadingptjec[isrc] >leadingPtThreshold[0]) {
	  if (leadingptjec[isrc] >leadingPtThreshold[0]) { 
	    //int njets=0;
	    ncount=0;
	    //recterm=0;
	    //ithird=-1;
	    //double sup = 1;	
	    //px=0;
	    //py=0;
	    //ptxy=0;
	    tmpjt4v.clear();
	    //tmpcand4v.clear();
	    tmpgen4v.clear();
	    
	    //if (abs((*ak4PFJets)[0].eta())<etarange[iet] && abs((*ak4PFJets)[1].eta())<etarange[iet]) {
	    //for(unsigned ijet = 0; ijet != ak4PFJets->size(); ijet++) {
	    
	    for(unsigned ijet = 0; ijet != ak4PFJets->size(); ijet++) {
	      if (abs((*ak4PFJets)[jetindx[isrc][0]].eta())<etarange[iet] && abs((*ak4PFJets)[jetindx[isrc][1]].eta())<etarange[iet]) {  // eta cut only on leading two jets ??
	      //if (abs((*ak4PFJets)[jetindx[isrc][0]].eta())<etarange[iet]) {
	      	//int flavour = (*ak4PFJets)[ijet].partonFlavour();
        	//cout << "Flavour : "<<flavour<<endl;
		
		int ireorjt = jetindx[isrc][ijet];
		
		int jetflav = (*ak4PFJets)[ireorjt].partonFlavour();	
		//if(jetflav==5){cout << "Flavour : "<<jetflav<<endl;}	
		
		double pt = jetptx[isrc][ijet];
		double sup = jetscl[isrc][ijet];
		double abseta = abs((*ak4PFJets)[ireorjt].eta());
	        if (pt<30.0 || abseta >etarange[iet]) continue;	
		//if (iet==0 && isrc==0) cout <<"pteta "<<pt<<" "<<abseta<<endl;
		bool isEta = (abseta<2.5) ? true : false;
		
		//if (isEta && pt>30.0) { njets++;}
		if (abseta>5.0) continue;
		bool isPt = (pt>30.0) ? true : false;
		if (isEta && isPt) {ncount++;}
		
		//cout<< "ncount = " << ncount << endl;
		//Jet ID ================= Tight ID 2017 Recomendation  check for 2018
		double NHF = (*ak4PFJets)[ireorjt].neutralHadronEnergyFraction();
		double NEMF = (*ak4PFJets)[ireorjt].neutralEmEnergyFraction();
		double CHF = (*ak4PFJets)[ireorjt].chargedHadronEnergyFraction();
		//double MUF = (*ak4PFJets)[ireorjt].muonEnergyFraction();
		//double CEMF = (*ak4PFJets)[ireorjt].chargedEmEnergyFraction();
		int NumConst = (*ak4PFJets)[ireorjt].chargedMultiplicity()+(*ak4PFJets)[ireorjt].neutralMultiplicity();
		//int NumNeutralParticles =(*ak4PFJets)[ireorjt].neutralMultiplicity();
		int CHM = (*ak4PFJets)[ireorjt].chargedMultiplicity();
                //cout<<"NHF== "<< NHF << "; NEF== " << NEMF <<" ; CHF==" <<CHF <<" ;cef==" << CEMF <<"; no= " << NumConst <<" ; nch==" << CHM <<" ; NO of part==" << NumNeutralParticles <<endl;
		bool TightJetID =false;
                if(abs((*ak4PFJets)[ireorjt].eta())<=2.7){
      		if(NHF<0.90 && NEMF<0.90 && NumConst>1 && CHF>0 && CHM>0 && abs((*ak4PFJets)[ireorjt].eta())<=2.4 )  TightJetID =true;
      		if(NHF<0.90 && NEMF<0.99 && abs((*ak4PFJets)[ireorjt].eta())>2.4 )  TightJetID =true;}
      		else {TightJetID =false;}

                if (abs((*ak4PFJets)[ireorjt].eta())>2.7) {TightJetID = false;}
                if ((*ak4PFJets)[ireorjt].pt()<30.0) {TightJetID = false;}
		
		if( ireorjt<=1 && !TightJetID) break;
		if (!TightJetID) continue;
		
		/*	
		if (ncount <=2 && ncount !=ijet+1) {
		  for (int ix=0; ix<ntype; ix++) { 
		    recomom[isrc][ix][iet].clear(); 
		  }
		  break;
		}
		if (isrc==0 && iet==0) {
		cout <<"recomom[isrc][0][iet].size() "<<iet<<" "<<isrc<<" "<<ijet<<" "<<ncount<<" "<<recomom[isrc][0][iet].size()<<" "<<recomom[isrc][1][iet].size()<<" "<<ncount<<" "<<ijet<<endl;}
		*/
		
		HepLorentzVector tmp4v((*ak4PFJets)[ireorjt].px(), (*ak4PFJets)[ireorjt].py(), (*ak4PFJets)[ireorjt].pz(), (*ak4PFJets)[ireorjt].energy());
		tmp4v *=sup;
		
		/*double respfact=respfun(1.02, 0.000004799, 0.000000007044,tmp4v.perp()); 
		  cout <<"Response factor = " <<respfact << endl;
		  tmp4v /=respfact;
		*/
		//cout << "Pt before correction = "<< tmp4v.perp()<< endl;
		/*double  respfact=0.;
		  bool isCorrect=false;   					
		  for (int iresp=0; iresp<7; iresp++){
		  if(abs(tmp4v.eta())> resetarange[iresp] && abs(tmp4v.eta())<resetarange[iresp+1]){
		  respfact=respfun(par0[iresp], par1[iresp], par2[iresp], tmp4v.perp());
		  isCorrect =true;
		  }
		  //		cout << "iresp = "<< iresp << " Eta = " <<tmp4v.eta() <<endl;
		  if(isCorrect) break;
		  }
		  //cout <<"Response factor = " <<respfact << endl;
		  double invrespfact=0;
		  if (respfact!=0) invrespfact =1/respfact;
		  tmp4v*=invrespfact;*/
		  //cout <<"Response factor = " <<respfact << "Inv Response factor =" << invrespfact << " Corrected Pt = "<< tmp4v.perp() <<endl;
		  //cout <<"perp "<<iet<<" "<<isrc<<" "<<ijet<<" "<<sup<<" "<<tmp4v<<" "<<tmp4v.eta()<<" "<<tmp4v.perp()<<" "<<pt<<endl;
		
		if (isEta && isPt) { tmpjt4v.push_back(tmp4v);}
		//tmpjt4v.push_back(tmp4v);	  
		//if (isEta && isPt) {allrecojetmom.push_back(tmp4v);}
		//if (ncount<=2) {  //change for all jet 26th June
		  if (isEta && isPt) {
		     recomom[isrc][0][iet].push_back(tmp4v);
		  }
		  //}
		  //cout <<"ncount filled "<<ncount<<" "<<isrc<<" "<<iet<<" "<<recomom[isrc][0][iet].size()<<endl;
		  //px +=tmp4v.px();
		  //py +=tmp4v.py();
		  //ptxy +=tmp4v.perp();
		  if (isrc==0) { 
		    if ((isInEtaRange[iet])) {recojt_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0) {recojt_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isEta && isPt) {recojt_phi->Fill(tmp4v.phi(), weighttrg);}
		    if (isEta && ncount==1) {recoht2_pt[iet]->Fill(aveleadingpt,weighttrg);}
			// quarks & antiquarks
			if(abs(jetflav==21 || jetflav==9)){
                        if (isInEtaRange[iet]) {recojtg_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtg_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtg_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
			if(abs(jetflav==5)){
                        if (isInEtaRange[iet]) {recojtb_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtb_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtb_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(abs(jetflav==4)){
                        if (isInEtaRange[iet]) {recojtc_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtc_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtc_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(abs(jetflav==3)){
                        if (isInEtaRange[iet]) {recojts_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojts_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojts_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(abs(jetflav==2)){
                        if (isInEtaRange[iet]) {recojtu_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtu_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtu_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(abs(jetflav==1)){
                        if (isInEtaRange[iet]) {recojtd_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtd_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtd_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
			// quarks
			if(jetflav==21 || jetflav==9){
                        if (isInEtaRange[iet]) {recojtqg_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtqg_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtqg_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==5){
                        if (isInEtaRange[iet]) {recojtqb_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtqb_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtqb_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==4){
                        if (isInEtaRange[iet]) {recojtqc_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtqc_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtqc_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==3){
                        if (isInEtaRange[iet]) {recojtqs_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtqs_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtqs_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==2){
                        if (isInEtaRange[iet]) {recojtqu_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtqu_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtqu_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==1){
                        if (isInEtaRange[iet]) {recojtqd_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtqd_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtqd_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
			// antiquark
			if(jetflav==-21 || jetflav==-9){
                        if (isInEtaRange[iet]) {recojtag_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtag_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtag_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-5){
                        if (isInEtaRange[iet]) {recojtab_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtab_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtab_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-4){
                        if (isInEtaRange[iet]) {recojtac_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtac_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtac_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-3){
                        if (isInEtaRange[iet]) {recojtas_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtas_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtas_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-2){
                        if (isInEtaRange[iet]) {recojtau_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtau_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtau_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-1){
                        if (isInEtaRange[iet]) {recojtad_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
                        if (isPt && iet==0) {recojtad_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojtad_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
		  }//if (isrc==0) {
		//} else {
		/*  if (isrc==0) { 
		    if ((isInEtaRange[iet])) {recojt_oth_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0) {recojt_oth_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isEta && isPt) {recojt_oth_phi->Fill(tmp4v.phi(), weighttrg);}
		  }*/
		/*  if (isEta && isPt) {
		    double tmppx = px + tmp4v.px();
		    double tmppy = py + tmp4v.py();
		    double tmppt = ptxy + tmp4v.perp();
		    double tmprec = sqrt(pow(tmppx, 2)+pow(tmppy, 2))/tmppt;
		    
		    if (tmprec>recterm) {
		      recterm = tmprec;
		      ithird = ireorjt;
		      //cout <<"ithird Data : "<< ijet<<endl;
		    }
		  }*/
		//}
		
		if (isrc==0) { 
		  if(ijet==0) { 
		    if (isInEtaRange[iet]) {recojt1_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1_pt = tmp4v.perp();}
		    if (isPt && iet==0) {recojt1_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isEta && isPt) {recojt1_phi->Fill(tmp4v.phi(), weighttrg);}
			// quarks & antiquarks
			if(abs(jetflav==21 || jetflav==9)){
			if (isInEtaRange[iet]) {recojt1g_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1g_pt = tmp4v.perp();}
			if (isPt && iet==0) {recojt1g_eta->Fill(tmp4v.eta(), weighttrg);}
			if (isEta && isPt) {recojt1g_phi->Fill(tmp4v.phi(), weighttrg);}
			}
			if(abs(jetflav==5)){
                        if (isInEtaRange[iet]) {recojt1b_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1b_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1b_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1b_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
			if(abs(jetflav==4)){
                        if (isInEtaRange[iet]) {recojt1c_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1c_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1c_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1c_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
			if(abs(jetflav==3)){
                        if (isInEtaRange[iet]) {recojt1s_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1s_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1s_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1s_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
			if(abs(jetflav==2)){
                        if (isInEtaRange[iet]) {recojt1u_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1u_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1u_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1u_phi->Fill(tmp4v.phi(), weighttrg);}
                        }	
			if(abs(jetflav==1)){
                        if (isInEtaRange[iet]) {recojt1d_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1d_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1d_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1d_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
			// quarks
			if(jetflav==21 || jetflav==9){
                        if (isInEtaRange[iet]) {recojt1qg_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1qg_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1qg_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1qg_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==5){
                        if (isInEtaRange[iet]) {recojt1qb_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1qb_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1qb_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1qb_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==4){
                        if (isInEtaRange[iet]) {recojt1qc_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1qc_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1qc_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1qc_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==3){
                        if (isInEtaRange[iet]) {recojt1qs_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1qs_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1qs_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1qs_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==2){
                        if (isInEtaRange[iet]) {recojt1qu_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1qu_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1qu_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1qu_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==1){
                        if (isInEtaRange[iet]) {recojt1qd_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1qd_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1qd_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1qd_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
			// antiquarks
			if(jetflav==-21 || jetflav==-9){
                        if (isInEtaRange[iet]) {recojt1ag_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1ag_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1ag_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1ag_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-5){
                        if (isInEtaRange[iet]) {recojt1ab_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1ab_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1ab_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1ab_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-4){
                        if (isInEtaRange[iet]) {recojt1ac_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1ac_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1ac_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1ac_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-3){
                        if (isInEtaRange[iet]) {recojt1as_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1as_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1as_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1as_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-2){
                        if (isInEtaRange[iet]) {recojt1au_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1au_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1au_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1au_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-1){
                        if (isInEtaRange[iet]) {recojt1ad_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet1ad_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt1ad_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt1ad_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
	
		  } else if(ijet==1){
		    if (isInEtaRange[iet]) {recojt2_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2_pt = tmp4v.perp();}
		    if (isPt && iet==0) {recojt2_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isInEtaRange[iet] && isPt) {recojt2_phi->Fill(tmp4v.phi(), weighttrg);}
			// quarks & antiquarks
			if(abs(jetflav==21 || jetflav==9)){
                        if (isInEtaRange[iet]) {recojt2g_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2g_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2g_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2g_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(abs(jetflav==5)){
                        if (isInEtaRange[iet]) {recojt2b_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2b_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2b_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2b_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(abs(jetflav==4)){
                        if (isInEtaRange[iet]) {recojt2c_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2c_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2c_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2c_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(abs(jetflav==3)){
                        if (isInEtaRange[iet]) {recojt2s_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2s_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2s_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2s_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(abs(jetflav==2)){
                        if (isInEtaRange[iet]) {recojt2u_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2u_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2u_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2u_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(abs(jetflav==1)){
                        if (isInEtaRange[iet]) {recojt2d_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2d_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2d_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2d_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
			// quarks
			if(jetflav==21 || jetflav==9){
                        if (isInEtaRange[iet]) {recojt2qg_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2qg_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2qg_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2qg_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==5){
                        if (isInEtaRange[iet]) {recojt2qb_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2qb_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2qb_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2qb_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==4){
                        if (isInEtaRange[iet]) {recojt2qc_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2qc_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2qc_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2qc_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==3){
                        if (isInEtaRange[iet]) {recojt2qs_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2qs_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2qs_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2qs_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==2){
                        if (isInEtaRange[iet]) {recojt2qu_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2qu_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2qu_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2qu_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==1){
                        if (isInEtaRange[iet]) {recojt2qd_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2qd_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2qd_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2qd_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
			// antiquarks
			if(jetflav==-21 || jetflav==-9){
                        if (isInEtaRange[iet]) {recojt2ag_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2ag_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2ag_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2ag_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-5){
                        if (isInEtaRange[iet]) {recojt2ab_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2ab_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2ab_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2ab_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-4){
                        if (isInEtaRange[iet]) {recojt2ac_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2ac_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2ac_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2ac_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-3){
                        if (isInEtaRange[iet]) {recojt2as_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2as_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2as_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2as_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-2){
                        if (isInEtaRange[iet]) {recojt2au_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2au_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2au_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2au_phi->Fill(tmp4v.phi(), weighttrg);}
                        }
                        if(jetflav==-1){
                        if (isInEtaRange[iet]) {recojt2ad_pt[iet]->Fill(tmp4v.perp(), weighttrg); recojet2ad_pt = tmp4v.perp();}
                        if (isPt && iet==0) {recojt2ad_eta->Fill(tmp4v.eta(), weighttrg);}
                        if (isEta && isPt) {recojt2ad_phi->Fill(tmp4v.phi(), weighttrg);}
                        }

		    if (isInEtaRange[iet] && ncount==2) { 
		      if (irecoht>=0 && irecoht<nHLTmx) { 
			recojtave_pt[iet][irecoht]->Fill(aveleadingpt, weighttrg);
			recojtavewt1_pt[iet][irecoht]->Fill(aveleadingpt);
		      }
		      
		      recojtallavewt1_pt[iet]->Fill(aveleadingpt);
		      recojtallave_pt[iet]->Fill(aveleadingpt, weighttrg);
		    }
		    
		  } else if(ijet==2) {
		    if (isInEtaRange[iet]) {recojt3_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0 ) {recojt3_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isInEtaRange[iet] && isPt) {recojt3_phi->Fill(tmp4v.phi(), weighttrg);} 
		  }
		  
		  if (tmpjt4v.size()==2 && isInEtaRange[iet]) { 
		    double dphi = dPhi(tmpjt4v[0].phi(), tmpjt4v[1].phi());
		    double dpt = tmpjt4v[0].perp() - tmpjt4v[1].perp();
		    double dperp = fabs(tmpjt4v[1].perp()*sin(dphi))/tmpjt4v[0].perp();
		    hjetdphi[iet]->Fill(dphi, weighttrg);
		    hjetdpt[iet]->Fill(dpt, weighttrg);
		    hjetptbypl[iet]->Fill(dperp, weighttrg);
		    hjetpt2bypt1[iet]->Fill(tmpjt4v[1].perp()/tmpjt4v[0].perp(), weighttrg);
		  }
		  
		  if (tmpjt4v.size()==3) {hjetpt3bypt2[iet]->Fill(tmpjt4v[2].perp()/tmpjt4v[1].perp(), weighttrg);}
		  } //if (isrc==0) {
		
	/* 
		if (ijet==0) {
			if (isInEtaRange[iet]) {recojet0_pt = tmp4v.perp();
				cout << "check0 :"<<recojet0_pt<<endl;}
			}
		if (ijet==1){ 
			if (isInEtaRange[iet]) {recojet1_pt = tmp4v.perp();}
			}
	*/

		//int nchg=0;
		nchg = 0;
		std::vector<reco::CandidatePtr> daus((*ak4PFJets)[ireorjt].daughterPtrVector());           
		std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); });                                                                                                  
		for (unsigned int i2 = 0; i2< daus.size(); ++i2) {   
		  const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
		  int charge = pfcand.charge();
		  HepLorentzVector cand4v(pfcand.px(), pfcand.py(), pfcand.pz(), pfcand.energy());
		  //tmpcand4v.push_back(cand4v);	
                  //nchg++;
		  if (charge !=0){
                  	//h_nchg[iet]->Fill(nchg, weighttrg);   // need to check
			nchg++;
                  	}
		  //if (cand4v.perp()<0.5) continue;
		  //if (ncount<=2 && isEta && isPt) { 
		  //recomom[isrc][1][iet].push_back(cand4v);
		    
		    //if (charge !=0) {
		      recomom[isrc][1][iet].push_back(cand4v);
#ifdef TRACKSYS
		      if (gRandom->Uniform() < 0.96) {recomom[isrc][2][iet].push_back(cand4v); }
#endif
		    //}
		   /* if (charge==0) { //other option if need open
		      if (cand4v.perp()>1.0) {
			recomom[isrc][3][iet].push_back(cand4v);
		      }  
		    } else {
		      if (cand4v.perp()>0.5) {
			recomom[isrc][3][iet].push_back(cand4v);
		      }
		    }*/
		    
		    //double dphi = dPhi(recomom[0][0][0].phi(), recomom[0][0][1].phi());
		    //double dpt = recomom[0][0][0].perp() - recomom[0][0][1].perp();
		    //double dperp = fabs(tmpcand4v[1].perp()*sin(dphi))/tmpjt4v[0].perp();
		    
		    //if (dpt<0) cout <<" "<< jk<<" "<<ij<<" "<<mn<<" "<<seljtvar4v[0]<<" "<<seljtvar4v[1]<<" "<<seljtvar4v[0].perp()<<" "<<seljtvar4v[1].perp()<<" "<<endl;
		    
		    //hjet1dphi->Fill(dphi, weighttrg);
		    //hjet1dpt->Fill(dpt, weighttrg);
		    //}
		  if (isrc==0) { 
		    //if (isEta && isPt) {
		      if (charge !=0) {
			recochg_phi->Fill(cand4v.phi(), weighttrg);
			recochg_pt->Fill(cand4v.perp(), weighttrg);
			recochg_eta->Fill(cand4v.eta(), weighttrg);
		      }
		    
		      if (ijet==0 && charge !=0) {
                        recochg1_phi->Fill(cand4v.phi(), weighttrg);
                        recochg1_pt->Fill(cand4v.perp(), weighttrg);
                        recochg1_eta->Fill(cand4v.eta(), weighttrg);
		      }
                     else if (ijet==1 && charge !=0) {
                        recochg2_phi->Fill(cand4v.phi(), weighttrg);
                        recochg2_pt->Fill(cand4v.perp(), weighttrg);
                        recochg2_eta->Fill(cand4v.eta(), weighttrg);
                      }
                     else if (ijet==2 && charge !=0) {
                        recochg3_phi->Fill(cand4v.phi(), weighttrg);
                        recochg3_pt->Fill(cand4v.perp(), weighttrg);
                        recochg3_eta->Fill(cand4v.eta(), weighttrg);
                      }                     

		    }//if (isrc==0) {
		     //}

		   // jet charge observables
		   if (isrc==0){
			for (int ik=0; ik<10; ik++){
				if (ijet==0 && cand4v.perp() > 1.0){
				//if(ijet==0){
					//recojet1_pt = tmp4v.perp();
					//cout<<"testing 1"<<endl;
                                        ijet1candsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);
					//ijet1candsmom.push_back(candsmom(charge, cand4v.perp(), kappa[ik]));

					ijet1_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
					ijet1_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

					ijet1_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
		                        ijet1_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
				if(abs(jetflav==21 || jetflav==9)){
					ijet1gcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);	

					ijet1g_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1g_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1g_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1g_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
					}
				if(abs(jetflav==5)){
					ijet1bcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1b_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1b_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1b_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1b_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
					}
				if(abs(jetflav==4)){
					ijet1ccandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1c_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1c_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1c_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1c_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
					}
				if(abs(jetflav==3)){
                                        ijet1scandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1s_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1s_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1s_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1s_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
				if(abs(jetflav==2)){
                                        ijet1ucandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1u_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1u_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1u_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1u_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
				if(abs(jetflav==1)){
                                        ijet1dcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1d_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1d_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1d_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1d_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
				// quarks
				if(jetflav==21 || jetflav==9){
                                        ijet1qgcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1qg_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1qg_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1qg_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1qg_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==5){
                                        ijet1qbcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1qb_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1qb_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1qb_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1qb_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==4){
                                        ijet1qccandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1qc_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1qc_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1qc_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1qc_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==3){
                                        ijet1qscandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1qs_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1qs_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1qs_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1qs_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==2){
                                        ijet1qucandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1qu_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1qu_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1qu_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1qu_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==1){
                                        ijet1qdcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1qd_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1qd_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1qd_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1qd_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
				// antiquark
				if(jetflav==-21 || jetflav==-9){
                                        ijet1agcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1ag_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1ag_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1ag_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1ag_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==-5){
                                        ijet1abcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1ab_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1ab_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1ab_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1ab_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==-4){
                                        ijet1accandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1ac_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1ac_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1ac_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1ac_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==-3){
                                        ijet1ascandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1as_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1as_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1as_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1as_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==-2){
                                        ijet1aucandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1au_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1au_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1au_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1au_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==-1){
                                        ijet1adcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet1ad_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet1ad_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet1ad_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet1ad_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
				}//if(ijet==0){
				if (ijet==1 && cand4v.perp() > 1.0){
				//else if(ijet==1){
					//recojet2_pt = tmp4v.perp();
					ijet2candsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

					ijet2_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
				if(abs(jetflav==21 || jetflav==9)){
                                        ijet2gcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2g_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2g_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2g_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2g_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(abs(jetflav==5)){
                                        ijet2bcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2b_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2b_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2b_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2b_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
				if(abs(jetflav==4)){
                                        ijet2ccandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2c_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2c_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2c_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2c_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(abs(jetflav==3)){
                                        ijet2scandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2s_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2s_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2s_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2s_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
				if(abs(jetflav==2)){
                                        ijet2ucandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2u_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2u_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2u_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2u_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(abs(jetflav==1)){
                                        ijet2dcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2d_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2d_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2d_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2d_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
				// quarks
				if(jetflav==21 || jetflav==9){
                                        ijet2qgcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2qg_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2qg_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2qg_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2qg_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==5){
                                        ijet2qbcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2qb_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2qb_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2qb_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2qb_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==4){
                                        ijet2qccandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2qc_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2qc_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2qc_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2qc_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==3){
                                        ijet2qscandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2qs_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2qs_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2qs_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2qs_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==2){
                                        ijet2qucandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2qu_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2qu_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2qu_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2qu_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==1){
                                        ijet2qdcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2qd_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2qd_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2qd_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2qd_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
				// antiquark
				 if(jetflav==-21 || jetflav==-9){
                                        ijet2agcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2ag_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2ag_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2ag_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2ag_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==-5){
                                        ijet2abcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2ab_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2ab_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2ab_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2ab_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==-4){
                                        ijet2accandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2ac_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2ac_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2ac_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2ac_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==-3){
                                        ijet2ascandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2as_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2as_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2as_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2as_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==-2){
                                        ijet2aucandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2au_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2au_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2au_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2au_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
                                if(jetflav==-1){
                                        ijet2adcandsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                        ijet2ad_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                        ijet2ad_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                        ijet2ad_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                        ijet2ad_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
					}				
				}//if(ijet==1){
			}//for (int ik=0; ik<10; ik++){
		}//if (isrc==0){
	} //for (unsigned int i2 = 0; i2< daus.size(); ++i2
            	h_nchg[iet]->Fill(nchg, weighttrg);    	
		//  if(isEta && isPt) {ncount++;}
	   	//  }//if (abs((*ak4PFJets)[jetindx[isrc][0]].eta())<etarange[iet] && abs((*ak4PFJets)[jetindx[isrc][1]].eta())<etarange[iet])
	  	//  }//for(unsigned ijet = 0; ijet != ak4PFJets->size(); ijet++)
/*	    if (ithird>=0) {
	      
	      recomom[isrc][0][iet].push_back(tmp4v);
	      //					cout <<"recomom[isrc][0][iet] "<< isrc<<" "<<iet<<" "<<recomom[isrc][0][iet].size()<<endl;
	      // tmpjt4v.push_back(tmp4v);   
	      
	      std::vector<reco::CandidatePtr> daus((*ak4PFJets)[ithird].daughterPtrVector());
	      std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2 ->pt(); });
	      for (unsigned int i2 = 0; i2< daus.size(); ++i2) {
		const pat::PackedCandidate &pfcand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
		int charge = pfcand.charge();
		HepLorentzVector cand4v(pfcand.px(), pfcand.py(), pfcand.pz(), pfcand.energy());
		//      if (cand4v.perp()<0.5) continue;                                                             
		recomom[isrc][1][iet].push_back(cand4v);
		
		if (charge !=0) {
		  recomom[isrc][2][iet].push_back(cand4v);
#ifdef TRACKSYS
		  if (gRandom->Uniform() < 0.96) {recomom[isrc][4][iet].push_back(cand4v); }
#endif
		  
		}
		if (charge==0){
		  if (cand4v.perp()>1.0) {
		    recomom[isrc][3][iet].push_back(cand4v);
		  }
		} else{
		  if (cand4v.perp()>0.5) {
		    recomom[isrc][3][iet].push_back(cand4v);
		  }
		}
	      }
	    }*/ //if (ithird>=0) 
	    //if (isrc==0) {h_njets[iet]->Fill(ncount, weighttrg);}
	    h_njets[iet]->Fill(ncount, weighttrg);
              } //if (abs((*ak4PFJets)[jetindx[isrc][0]].eta())<etarange[iet] && abs((*ak4PFJets)[jetindx[isrc][1]].eta())<etarange[iet])
            } // for(unsigned ijet = 0; ijet != ak4PFJets->size(); ijet++)
	  } //if (aveleadingptjec[isrc] >leadingPtThreshold[0])
	} // 	for (int isrc = 0; isrc < njecmx; isrc++)
      } //for (int iet=0; iet<njetetamn; iet++)	   
    } // if (ak4PFJets.isValid() && ak4PFJets->size()>=2 && (*ak4PFJets)[0].pt()>leadingPtThreshold[0])
  } // if (ak4PFJets.isValid())
#endif
  
  //  cout << "Write test 31 = ok " << endl;
  //===================********Trigger****============================================================
  
  //t2=clock();
  //float diff ((float)t2-(float)t1);
  //if(diff>30000) return;
  //cout << "Time T2 = " << t2 << " ;Time Diff to Run ="<< diff << endl;
/*  
#ifdef TRIGGER
  if(!isMC){
#ifndef DIJETAVE
    //  vector<triggervar> alltrgobj;
    if (trigRes.isValid() && isReconstruct  &&
	(tmpjt4v.size() ==2 || (tmpjt4v.size()>=3 && tmpjt4v[2].perp()<30.0)) &&
	abs(dPhi(tmpjt4v[0].phi(), tmpjt4v[1].phi()))>2.0){
      
      //  if (trigRes.isValid() && isReconstruct  &&
      //   (tmpjt4v.size() ==2) && abs(dPhi(tmpjt4v[0].phi(), tmpjt4v[1].phi()))>2.0){
      
      int ijet = int(2*gRandom->Uniform())%2;
      int ijet2 = (ijet==0) ? 1 : 0;
      //cout <<"gRandom= "<<gRandom->Uniform() << " ijet" <<ijet<<endl; 
      HepLorentzVector tagjet4v = tmpjt4v[ijet];
      HepLorentzVector probjet4v = tmpjt4v[ijet2];

      vector<triggervar> alltrgobj;
      alltrgobj.clear(); 
      const edm::TriggerNames &names = iEvent.triggerNames(*trigRes);
      for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
	obj.unpackPathNames(names);
	std::vector<std::string> pathNamesAll  = obj.pathNames(false);
	for (unsigned ih = 0, n = pathNamesAll.size(); ih < n; ++ih) {
	  variab2 = pathNamesAll[ih].c_str(); 
	  for (int jk=0; jk<nHLTmx; jk++) {
	    if (strstr(variab2,jethlt_name[jk]) && strlen(variab2)-strlen(jethlt_name[jk])<5){
	      triggervar tmpvec;
	      if( obj.pt()<jethlt_thr[jk] ) continue;
	      tmpvec.both = obj.hasPathName( pathNamesAll[ih], true, true );
	      if(obj.pt()>10){
		tmpvec.both = obj.hasPathName( pathNamesAll[ih], true, true );
		tmpvec.highl  = obj.hasPathName( pathNamesAll[ih], false, true );
		tmpvec.level1 = obj.hasPathName( pathNamesAll[ih], true, false );
		tmpvec.trg4v = HepLorentzVector(obj.px(), obj.py(), obj.pz(), obj.energy());
		tmpvec.prescl = 1;
		tmpvec.ihlt = jk;
		alltrgobj.push_back(tmpvec);
	      }
	    }
	  }
	}
      }
      
      for (unsigned ij=0; ij<alltrgobj.size(); ij++) {
	HepLorentzVector trigger4v = alltrgobj[ij].trg4v;
	int ihlt = -1;
	int tmphlt = alltrgobj[ij].ihlt;
	if( trigger4v.perp()<jethlt_thr[tmphlt] ) continue;
	//      bool isBoth=alltrgobj[ij].both;
	bool isLF =alltrgobj[ij].level1;
	bool isL3 =alltrgobj[ij].highl;
	double angle = deltaR(tagjet4v, trigger4v);
	if (isLF) { 
	  trgjet_angle[tmphlt][0]->Fill(angle);
	  trgjet_2dangle[tmphlt][0]->Fill(trigger4v.perp(), angle);	
	}
	if (isL3) { 
	  trgjet_angle[tmphlt][1]->Fill(angle);
	  trgjet_2dangle[tmphlt][1]->Fill(trigger4v.perp(), angle);		
	}	
	// bool tag=false;
	if (deltaR(tagjet4v, trigger4v)<0.2) {
	  // tag=true;
	  ihlt = alltrgobj[ij].ihlt;
	  if (isLF)  {
	    //        if (isLF && !isBoth)  {
	    trgjet_pt[ihlt][0]->Fill(probjet4v.perp());
	    trgjet_eta[ihlt][0]->Fill(probjet4v.eta());
	    trgjet_phi[ihlt][0]->Fill(probjet4v.phi());
	  }
	  if (isL3) {
	    // if (isLF && !isBoth)  {
	    trgjet_pt[ihlt][1]->Fill(probjet4v.perp());
	    trgjet_eta[ihlt][1]->Fill(probjet4v.eta());
	    trgjet_phi[ihlt][1]->Fill(probjet4v.phi());
	  }
	  
	  for (unsigned jk=0; jk<alltrgobj.size(); jk++) {
	    if (ij==jk || alltrgobj[jk].ihlt !=ihlt) continue;
	    //        if( trigprbjet4v.perp()<jethlt_thr[ihlt] ) continue;
	    HepLorentzVector trigprbjet4v = alltrgobj[jk].trg4v;
	    if( trigprbjet4v.perp()<jethlt_thr[ihlt] ) continue;
	    double angle1 = deltaR(probjet4v, trigprbjet4v);
	    if (isLF && angle1<0.5 ) {
	      //          if (isLF && !isBoth && angle1<0.5 ) {
	      prbjet_pt[ihlt][0]->Fill(probjet4v.perp());
	      prbjet_eta[ihlt][0]->Fill(probjet4v.eta());
	      prbjet_phi[ihlt][0]->Fill(probjet4v.phi());
	      isLF = false;
	    }
	    
	    if (isL3 && angle1<0.2) {
	      //          if (isL3 && !isBoth && angle1<0.2)
	      prbjet_pt[ihlt][1]->Fill(probjet4v.perp());
	      prbjet_eta[ihlt][1]->Fill(probjet4v.eta());
	      prbjet_phi[ihlt][1]->Fill(probjet4v.phi());
	      isL3 = false;
	    }
	    
	    if ((!isL3) && (!isLF)) continue;
	    
	  } //for (unsigned jk=0; jk<alltrgobj.size(); jk++) 
	} //if (deltaR(tagjet4v, trigger4v)<0.2)
	if (ihlt>=0) continue;
      } //for (int ij=0; ij<alltrgobj.size(); ij++)
    } // if (trigRes.isValid() && m_trigeff && isReconstruct  && tmpjt4v.size() ==2) &&
    //      abs(dPhi(tmpjt4v[0].phi(), tmpjt4v[1].phi()))>2.0
#endif
  }
#endif*/
  //======******Trigger Efficiency Normal=======================
  
  //cout << "Write test 1 = ok " << endl;
  //==================================***GenJets*****=================================
  //cout<<"0 aveleadingpt "<<aveleadingpt<< " ; "<<ihltfill<<" "<<irecoht<<endl;
  if(isMC) {
    
    edm::Handle<reco::GenJetCollection> genjets;
    iEvent.getByToken(genjetToken_,genjets);
    
    double avegenpt =0;
    double leadgenpt =0;
    //cout <<"HGebjet "<<endl;
    if (genjets.isValid() &&  genjets->size()>=2) {
#ifdef DIJETAVE
      for (int iet=0; iet<njetetamn; iet++) {
	isInEtaRange[iet] = true;
      }
      
      for (int ij=0; ij<2; ij++) {
	for (int iet=0; iet<njetetamn; iet++) {
	  if (abs((*genjets)[ij].eta())>etarange[iet]) { isInEtaRange[iet] = false;}
	}
	
	if (abs((*genjets)[ij].eta())<2.5 && (*genjets)[ij].pt()>30.0 ) { 
	  avegenpt +=(*genjets)[ij].pt();
          leadgenpt = (*genjets)[0].pt();
	} else {avegenpt -=100000;
		leadgenpt -=100000;}
      }
      avegenpt /=2.0;
      
#else 

#endif
    }
    
    //igenht = getbinid(avegenpt, njetptmn, leadingPtThreshold);
    igenht = getbinid(leadgenpt, njetptmn, leadingPtThreshold);
 
    //cout << "Write test 2 = ok " << endl;
    //cout << "Write test 321 = ok " << endl;
    vector<double> genjetptx[nGenReso];
    vector<double> genjetscl[nGenReso]; 
    vector<int> genjetindx[nGenReso];
    
    for(unsigned ijet = 0; ijet != genjets->size(); ijet++) {
      double pt = (*genjets)[ijet].pt();
      //#ifdef JETRESO		
      //			double eta = (*genjets)[ijet].eta();
      //			double reso = 1;
      //			JME::JetParameters parameters_5 = {{JME::Binning::JetPt, pt}, {JME::Binning::JetEta, eta}, {JME::Binning::Rho, *rho}};
      //			float rp = resolution.getResolution(parameters_5);
      //			float sf = res_sf.getScaleFactor({{JME::Binning::JetEta, eta}});
      //			float sf_up= res_sf.getScaleFactor({{JME::Binning::JetEta, eta}}, Variation::UP);
      //			float sf_dn= res_sf.getScaleFactor({{JME::Binning::JetEta, eta}}, Variation::DOWN);
      
      //#endif		
      for (int isrc = 0; isrc < nGenReso; isrc++) {
	double sup = 1.0;
	//#ifdef JETRESO
	//				if (isrc==0) {  
	//					reso = sqrt(sf*sf - 1)*rp;
	//				} else if (isrc==1) {
	//					reso = sqrt(sf_up*sf_up - 1)*rp;
	//				} else if (isrc==2) {
	//					reso = sqrt(sf_dn*sf_dn - 1)*rp;
	//				}
	
	//				sup = gRandom->Gaus(1.0, reso);
	//				//cout<<"isrc "<< ijet<<" "<< pt<<" "<<eta<<" "<<isrc<<" rp "<<rp<<" "<<sf<<" "<<sf_dn<<" "<<sf_up<<" "<<reso<<" "<<sup<<endl;
	//#endif
	
	genjetptx[isrc].push_back(sup*pt);
	genjetscl[isrc].push_back(sup);
	genjetindx[isrc].push_back(ijet);
      }
    }
    
       //cout << "Write test 3 = ok " << endl;
       //cout << "Write test 322 = ok "<<nGenReso << endl;
//////#ifdef JETRESO
    for (int isrc = 0; isrc < nGenReso; isrc++) {
      //cout << "Write test 31 = ok "<<isrc << " ; " << genjetptx[isrc].size() <<endl;
      if(genjetptx[isrc].size()==0) break;
      for (unsigned int ij=0; ij<genjetptx[isrc].size()-1; ij++) {
	//cout << "Write test 32 = ok "<<nGenReso << endl;
	for (unsigned int jk=ij+1; jk<genjetptx[isrc].size(); jk++) {  
	  //if(jk<genjetptx[isrc].size()) return;
	  //cout << "Write test 33 = ok "<<nGenReso << endl;
	  if (genjetptx[isrc][jk]>genjetptx[isrc][ij]) {
	    //cout << "Write test 34 = ok "<<nGenReso << endl;
	    double tmppt = genjetptx[isrc][ij];
	    double tmpscl = genjetscl[isrc][ij];
	    int tmpindx = genjetindx[isrc][ij];
	    
	    genjetptx[isrc][ij] = genjetptx[isrc][jk];
	    genjetscl[isrc][ij] = genjetscl[isrc][jk];
	    genjetindx[isrc][ij] = genjetindx[isrc][jk];					
	    
	    genjetptx[isrc][jk] = tmppt;
	    genjetscl[isrc][jk] = tmpscl;
	    genjetindx[isrc][jk] = tmpindx;
	    //cout << "Write test 35 = ok "<<nGenReso << endl;
	  }
	}
      }
    }
    //////#endif
    //cout << "Write test 4 = ok " << endl;
    //double avegenptres[nGenReso]={0};
    
    for (int isrc = 0; isrc < nGenReso; isrc++) {
      if (genjetptx[isrc].size()>=2) {
	//avegenptres[isrc] = 0.5*(genjetptx[isrc][0] + genjetptx[isrc][1]);
        leadgenptres[isrc] = (genjetptx[isrc][0]);
	//igenhtres[isrc] = getbinid(avegenptres[isrc], njetptmn, leadingPtThreshold);
        igenhtres[isrc] = getbinid(leadgenptres[isrc], njetptmn, leadingPtThreshold);
      } else {
	igenhtres[isrc] = -1;
      }
    }
    
    if(genjets.isValid() && genjets->size() >=2) { //  && avegenpt>leadingPtThreshold[0]) { 
      for (int iet=0; iet<njetetamn; iet++) {
	for (int isrc=0; isrc<nGenReso; isrc++) { 
	  //if (avegenptres[isrc] > leadingPtThreshold[0]) {
          if (leadgenptres[isrc] > leadingPtThreshold[0]) {
	    //double px =0;
	    //double py =0;
	    //double ptxy =0;
	    
	    ncount=0;
	    //int recterm=0;
	    //int ithird=-1;
	    
	    //for(unsigned ijet = 0; ijet < genjets->size(); ijet++) {
	    for(unsigned ijet = 0; ijet != genjets->size(); ijet++) {
	      //int igenjt = genjetindx[isrc][ijet];
	     /* if ((*genjets)[igenjt].pt()>25.0) {
		cout<<"ievt "<<ievt<<" "<<ijet<<" "<<igenjt<<" "<<genjetptx[isrc][ijet]<<" "<<(*genjets)[igenjt].pt()<<" "<<(*genjets)[igenjt].eta()<<" "<<(*genjets)[igenjt].phi()<<endl;
	      }*/

	      if (abs((*genjets)[genjetindx[isrc][0]].eta())<etarange[iet] && abs((*genjets)[genjetindx[isrc][1]].eta())<etarange[iet]) {     // need to check
	      //if (abs((*genjets)[genjetindx[isrc][0]].eta())<etarange[iet]) {	
		int igenjt = genjetindx[isrc][ijet];
	
		double pt = genjetptx[isrc][ijet];
		double sup = genjetscl[isrc][ijet];
		double abseta = abs((*genjets)[igenjt].eta());
		if (pt<30.0 || abseta >etarange[iet]) continue;
		bool isEta = (abseta<2.5) ? true : false;
		
		//if (iet==0 && isrc==0) 
		//cout <<"MC:pteta "<<ijet<<" "<<pt<<" "<<abseta<<endl;
		if (abseta>5.0) continue;
		bool isPt = (pt>30.0) ? true : false;
		if (isEta && isPt) {ncount++;}
		
		
		HepLorentzVector tmp4v((*genjets)[igenjt].px(), (*genjets)[igenjt].py(), (*genjets)[igenjt].pz(), (*genjets)[igenjt].energy());
		
		tmp4v *=sup;
		//bool isPt = (pt>30.0) ? true : false;
		//Response 
		/*if(isPt && isReconstruct) {
		  for(unsigned ijet = 0; ijet != ak4PFJets->size(); ijet++) {
		    HepLorentzVector tmp4vreco((*ak4PFJets)[ijet].px(), (*ak4PFJets)[ijet].py(), (*ak4PFJets)[ijet].pz(), (*ak4PFJets)[ijet].energy());
		    bool isResp=false;
		    for (int iresp=0; iresp<7; iresp++){
		      bool isEtaMatch=false; 
		      if(abs(tmp4v.eta())> resetarange[iresp] && abs(tmp4v.eta())<resetarange[iresp+1]){
			isEtaMatch=true;
			// cout << "Eta Match = " << tmp4v.eta() << " iresp = "<< iresp << endl; 
			double respangle=deltaR(tmp4v,tmp4vreco);
			if(respangle <0.2) {
			  resp_jet[iresp]->Fill(tmp4v.perp(), tmp4vreco.perp()/tmp4v.perp(), weighttrg);
			  resp_jet1[iresp]->Fill(abs((tmp4v.perp()-tmp4vreco.perp())/tmp4v.perp()), weighttrg);
			  //		cout << "Resolution = " << abs((tmp4v.perp()-tmp4vreco.perp())/tmp4v.perp()) << endl;
			  isResp=true;
			}
		      }		
		      if(isEtaMatch) break;
		    }
		    if(isResp) break;
		  }
		}*/
		//Response
		
		//								cout <<"isrc "<<iet<<" "<< isrc <<" "<<igenjt<<" "<<pt<<" " <<tmp4v.perp()<<" "<<tmp4v.eta()<<endl;
		//							pt = tmp4v.perp();
		
		//bool isPt = (pt>30.0) ? true : false;
		//if (isEta && isPt) {ncount++;}
		
		/*if (ncount <=2 && ncount !=ijet+1) {
		  for (int ix=0; ix<ntype; ix++) { 
		    genmom[isrc][ix][iet].clear(); 
		  }
		  break;
		}*/
		//							if (isrc==0 && iet==0) cout <<"	genmom[isrc][0][iet].size() "<<genmom[isrc][0][iet].size()<<" "<<ncount<<" "<<igenjt<<" "<<endl;
		
		if (isEta && isPt) { tmpgen4v.push_back(tmp4v);} 
		//if (ncount<=2) {
		  if (isEta && isPt) {
		    genmom[isrc][0][iet].push_back(tmp4v);
		   // cout <<"twoijx "<<isrc<<" "<<iet<<" "<< genmom[isrc][0][iet].size()<<" "<<tmp4v.perp()<<" "<<tmp4v.eta()<<" "<<tmp4v.phi()<<endl;
		  }
		  //px +=tmp4v.px();
		 // py +=tmp4v.py();
		 // ptxy +=tmp4v.perp();
		  if (isrc==0) { 
		    if (isInEtaRange[iet]) {genjt_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0) {genjt_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isEta && isPt) {genjt_phi->Fill(tmp4v.phi(), weighttrg);}  
		  }
		//} else {
		  /*if (isrc==0) { 
		    if (isInEtaRange[iet]) {genjt_oth_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0) {genjt_oth_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isEta && isPt) {genjt_oth_phi->Fill(tmp4v.phi(), weighttrg);}
		  }*/
		 /* if (isEta && isPt) {
		    double tmppx = px + tmp4v.px();
		    double tmppy = py + tmp4v.py();
		    double tmppt = ptxy + tmp4v.perp();
		    double tmprec = sqrt(pow(tmppx, 2)+pow(tmppy, 2))/tmppt;
		    
		    if (tmprec>recterm) {
		      recterm = tmprec;
		      ithird = igenjt;
		      // 		      cout <<"ithird MC : "<< igenjt<<" "<<tmp4v.perp()<<" "<<tmp4v.eta()<<endl;
		    }
		  }*/
		//}
		
		if (isrc==0) { 
		  if(ijet==0) {
		    //cout<<"Gen Pt= " << avegenpt <<endl;
		    if (isInEtaRange[iet]) {genjt1_pt[iet]->Fill(tmp4v.perp(), weighttrg); genrecojet1_pt = tmp4v.perp();}
		    if (isPt && iet==0) {genjt1_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isEta && isPt) {genjt1_phi->Fill(tmp4v.phi(), weighttrg);}
		  } else if(ijet==1){
		    //cout<<"okkkkkkkk" <<endl;
		    if (isInEtaRange[iet]) {genjt2_pt[iet]->Fill(tmp4v.perp(), weighttrg); genrecojet2_pt = tmp4v.perp();}
		    if (isPt && iet==0) {genjt2_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isInEtaRange[iet] && isPt) {genjt2_phi->Fill(tmp4v.phi(), weighttrg);}
		    if (isInEtaRange[iet] && ncount==2) {
		      //cout<<"Gen Pt 1= " << avegenpt <<endl;
		      genjtallave_pt[iet]->Fill(avegenpt, weighttrg);
		    }
		  } else if(ijet==2) {
		    if (isInEtaRange[iet]) {genjt3_pt[iet]->Fill(tmp4v.perp(), weighttrg);}
		    if (isPt && iet==0 ) {genjt3_eta->Fill(tmp4v.eta(), weighttrg);}
		    if (isInEtaRange[iet] && isPt) {genjt3_phi->Fill(tmp4v.phi(), weighttrg);}
		  }
		  if (tmpgen4v.size()==2 && isInEtaRange[iet]) {
		    double dphi = dPhi(tmpgen4v[0].phi(), tmpgen4v[1].phi());
		    double dpt = tmpgen4v[0].perp() - tmpgen4v[1].perp();
		    double dperp = fabs(tmpgen4v[1].perp()*sin(dphi))/tmpgen4v[0].perp();
		    
		    genjetdphi[iet]->Fill(dphi, weighttrg);
		    genjetdpt[iet]->Fill(dpt, weighttrg);
		    genjetptbypl[iet]->Fill(dperp, weight);
		    genjetpt2bypt1[iet]->Fill(tmpgen4v[1].perp()/tmpgen4v[0].perp(), weight);
		  }
		  
		  if (tmpgen4v.size()==3) {genjetpt3bypt2[iet]->Fill(tmpgen4v[2].perp()/tmpgen4v[1].perp(), weight);}
		}
#ifdef GENPART

		std::vector <const GenParticle*> daus ((*genjets)[igenjt].getGenConstituents ());
		//std::sort(daus.begin(),daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); });  //need to check 
		
		for (unsigned int i2 =0; i2< daus.size(); ++i2) {
		  const GenParticle* pfcand = daus[i2];
		  int charge = pfcand->charge();
		  HepLorentzVector cand4v(pfcand->px(), pfcand->py(), pfcand->pz(), pfcand->energy());
		  //int pdgid = pfcand->pdgId();
		 //std::cout<<"GENPART loop"<<endl; 
#else		
					
		  std::vector<reco::CandidatePtr> daus((*genjets)[igenjt].daughterPtrVector());
		  std::sort(daus.begin(),daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); });                               
		  
		  for (unsigned int i2 = 0; i2< daus.size(); ++i2) {
		    const pat::PackedCandidate &pfcand = static_cast<const pat::PackedCandidate &>(*daus[i2]);
		    int charge = pfcand.charge();
		    //std::cout<<"ALL loop"<<endl;
		    HepLorentzVector cand4v(pfcand.px(), pfcand.py(), pfcand.pz(), pfcand.energy());
#endif
		    //	    if (cand4v.perp()<0.5) continue;
		    
		    //if (ncount<=2 && isEta && isPt) {
		      //genmom[isrc][1][iet].push_back(cand4v);
		      //if (charge !=0) {
			genmom[isrc][1][iet].push_back(cand4v);
#ifdef TRACKSYS
			if (gRandom->Uniform() < 0.96) {genmom[isrc][2][iet].push_back(cand4v); }
#endif
		     // }
		      
		      //   if (charge ==0) {genmom[isrc][2][iet].push_back(cand4v);}
		      
		     /* if(charge ==0) {
			if (cand4v.perp()>1.0) {
			  genmom[isrc][3][iet].push_back(cand4v);
			}
		      } else {
			if (cand4v.perp()>0.5) {
			  genmom[isrc][3][iet].push_back(cand4v);
			}
		      }*/
		    //}
		    /*
		    if (isrc==0) { 
			if (charge !=0) {
		      if (isEta && isPt) {
			  genchg_phi->Fill(cand4v.phi(), weighttrg);
			} 
                      if (isEta) {
                         genchg_pt->Fill(cand4v.perp(), weighttrg);
                         }
                      if (isPt) {
                         genchg_eta->Fill(cand4v.eta(), weighttrg);
                        }

                       if(ijet==0) {
                      if (isEta && isPt) {
                          genchg1_phi->Fill(cand4v.phi(), weighttrg);
                        }
                      if (isEta) {
                         genchg1_pt->Fill(cand4v.perp(), weighttrg);
                         }
                      if (isPt) {
                         genchg1_eta->Fill(cand4v.eta(), weighttrg);
                        }                      
                      }

                      if(ijet==1) {
                      if (isEta && isPt) {
                          genchg2_phi->Fill(cand4v.phi(), weighttrg);
                        }
                      if (isEta) {
                         genchg2_pt->Fill(cand4v.perp(), weighttrg);
                         }
                      if (isPt) {
                         genchg2_eta->Fill(cand4v.eta(), weighttrg);
                        }
                      }
                  
                      if(ijet==2) {
                      if (isEta && isPt) {
                          genchg3_phi->Fill(cand4v.phi(), weighttrg);
                        }
                      if (isEta) {
                         genchg3_pt->Fill(cand4v.perp(), weighttrg);
                         }
                      if (isPt) {
                         genchg3_eta->Fill(cand4v.eta(), weighttrg);
                        }
                      }
			*/
	if (isrc==0) { 
               	if (charge !=0) {
			genchg_phi->Fill(cand4v.phi(), weighttrg);
			genchg_pt->Fill(cand4v.perp(), weighttrg);
			genchg_eta->Fill(cand4v.eta(), weighttrg);
		  	}
		if (ijet==0 && charge !=0) {
			genchg1_phi->Fill(cand4v.phi(), weighttrg);
                        genchg1_pt->Fill(cand4v.perp(), weighttrg);
                        genchg1_eta->Fill(cand4v.eta(), weighttrg);
			}
		else if (ijet==1 && charge !=0) {
			genchg2_phi->Fill(cand4v.phi(), weighttrg);
                        genchg2_pt->Fill(cand4v.perp(), weighttrg);
                        genchg2_eta->Fill(cand4v.eta(), weighttrg);
			}
		else if (ijet==2 && charge !=0) {
                        genchg3_phi->Fill(cand4v.phi(), weighttrg);
                        genchg3_pt->Fill(cand4v.perp(), weighttrg);
                        genchg3_eta->Fill(cand4v.eta(), weighttrg);
                        }
		}//if (isrc==0)
	
	if (isrc==0){
		for (int ik=0; ik<10; ik++){
			if (ijet==0 && cand4v.perp()>1.0){
				//genrecojet1_pt = tmp4v.perp();
				//cout << " Test1 : "<<genrecojet1_pt<<endl;
				igenjet1candsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);
				igenjet1_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                igenjet1_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                igenjet1_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                igenjet1_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
				}
			if (ijet==1 && cand4v.perp()>1.0){
				//genrecojet2_pt = tmp4v.perp();
                                igenjet2candsmom[ik] += candsmom(charge, cand4v.perp(), kappa[ik]);

                                igenjet2_long_num[ik] += (charge*(dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik])));
                                igenjet2_long_den[ik] += (dotproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));

                                igenjet2_tran_num[ik] += (charge*(crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(),kappa[ik])));
                                igenjet2_tran_den[ik] += (crossproduct(cand4v.px(), cand4v.py(), cand4v.pz(), tmp4v.px(), tmp4v.py(), tmp4v.pz(), tmp4v.perp(), kappa[ik]));
                                        }
			}//for (int ik=0; ik<10; ik++){
		}//if (isrc==0){
                    /*if (isEta) {
                        if (charge !=0) {
                          genchg_pt->Fill(tmp4v.perp(), weighttrg);
                        }
                    if (isPt) {
                        if (charge !=0) {
                          genchg_eta->Fill(tmp4v.eta(), weighttrg);
                        }                      

else {
			if (charge !=0) {
			  genchg_oth_phi->Fill(tmp4v.phi(), weighttrg);
			} else {
			  genneu_oth_phi->Fill(tmp4v.phi(), weighttrg);
			}
		      }
		      
		      if (isEta) {
			if (charge !=0) {
			  genchg_pt->Fill(tmp4v.perp(), weighttrg);
			} else {
			  genneu_pt->Fill(tmp4v.perp(), weighttrg);
			}
		      } else {
			if (charge !=0) {
			  genchg_oth_pt->Fill(tmp4v.perp(), weighttrg);
			} else {
			  genneu_oth_pt->Fill(tmp4v.perp(), weighttrg);
			}
		      }
		      if (isPt) {
			if (charge !=0) {
			  genchg_eta->Fill(tmp4v.eta(), weighttrg);
			} else {
			  genneu_eta->Fill(tmp4v.eta(), weighttrg);
			}
		      } else {
			if (charge !=0) {
			  genchg_oth_eta->Fill(tmp4v.eta(), weighttrg);
			} else {
			  genneu_oth_eta->Fill(tmp4v.eta(), weighttrg);
			}
		      }*/
		    //} //if (isrc==0)
		  } //for (unsigned int i2 = 0; i2< daus.size(); ++i2)

		  //  if (isEta && isPt) {ncount++;}
		} // if (abs((*genjets)[genjetindx[isrc][0]].eta())<etarange[iet] && 
		//								abs((*genjets)[genjetindx[isrc][1]].eta())<etarange[iet])
	      } //	for(unsigned ijet = 0; ijet != genjets->size(); ijet++) 
	      //cout << "Write test 324 = ok " << endl;
	     /* if (ithird>=0) {
		//							cout <<"ithird "<<isrc<<" "<< iet<<" "<< ithird<<endl;
		
		HepLorentzVector tmp4v((*genjets)[ithird].px(), (*genjets)[ithird].py(), (*genjets)[ithird].pz(), (*genjets)[ithird].energy());
		genmom[isrc][0][iet].push_back(tmp4v);
		//cout <<"thirdijxxx "<<isrc<<" "<<iet<<" "<< genmom[isrc][0][iet].size()<<" "<<genjets->size()<<" "<<ithird<<" "<<tmp4v.perp()<<" "<<tmp4v.eta()<<" "<<tmp4v.phi()<<" "<<setprecision(14)<<weighttrg<<endl;
#ifdef GENPART
		std::vector <const GenParticle*> daus ((*genjets)[ithird].getGenConstituents ());
		//								std::sort(daus.begin(),daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); 
		
		for (unsigned int i2 =0; i2< daus.size(); ++i2) {
		  const GenParticle* pfcand = daus[i2];
		  int charge = pfcand->charge();
		  HepLorentzVector cand4v(pfcand->px(), pfcand->py(), pfcand->pz(), pfcand->energy());
		  //								int pdgid = pfcand->pdgId();
		  
#else
		  std::vector<reco::CandidatePtr> daus((*genjets)[ithird].daughterPtrVector());
		  std::sort(daus.begin(),daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); });    
		  
		  for (unsigned int i2 = 0; i2< daus.size(); ++i2) {
		    const pat::PackedCandidate &pfcand = static_cast<const pat::PackedCandidate &>(*daus[i2]);
		    
		    int charge = pfcand.charge();
		    HepLorentzVector cand4v(pfcand.px(), pfcand.py(), pfcand.pz(), pfcand.energy());
#endif
		    
		    
		    //      if (cand4v.perp()<0.5) continue;                                                                                                                     
		    genmom[isrc][1][iet].push_back(cand4v);
		    if (charge !=0) {
		      genmom[isrc][2][iet].push_back(cand4v);
#ifdef TRACKSYS
		      if (gRandom->Uniform() < 0.96) {genmom[isrc][4][iet].push_back(cand4v); }
#endif
		    }
		    
		    if(charge ==0) {
		      if (cand4v.perp()>1.0) {
			genmom[isrc][3][iet].push_back(cand4v);
		      }
		    } else {
		      if (cand4v.perp()>0.5) {
			genmom[isrc][3][iet].push_back(cand4v);
		      }
		    }
		  } //for (unsigned int i2 = 0; i2< daus.size(); ++i2) 
		}*/// if (ithird>=0)
		gen_njets[iet]->Fill(ncount,weighttrg); 
	      } // if (avegenptres[isrc] > leadingPtThreshold[0])
	    } //	for (int isrc=0; isrc<nGenReso; isrc++)
	  } //for (int iet=0; iet<njetetamn; iet++)
	} // if(genjets.isValid() && genjets->size()>=2 && (*genjets)[0].pt()>leadingPtThreshold[0])
	// } //if (genjets.isValid() &&  genjets->size()>=2) 
	h_2ht->Fill(aveleadingpt,avegenpt, weighttrg);
	///////Response
      } //isMC
      //	cout<<"22 aveleadingpt "<<aveleadingpt<< " ; "<<ihltfill<<" "<<irecoht<<endl;
      // if(isMC) h_2ht->Fill(aveleadingpt,avegenpt, weighttrg);
      //cout << "Write test 325 = ok " << endl;
      //for(int rnum=0; rnum<10; rnum++) {
      /*double rand=gRandom->Uniform();
      int k = rand/0.1;
      //cout << "Rand Number " << k << endl;*/
 
//-----------------------------------------------Calculate And Fill Jet Charge Obserables------------------------------------
	//Profile histograms
	//hprof->Fill(recojet0_pt,ijet0_candsmom_10/(pow(recojet0_pt,1.0)),weighttrg);
	hchpt->Fill(recojet1_pt,nchg,weighttrg);
/*
bool isValue=false;
double v_recojc_D1J1[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double v_recojc_D1J2[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double v_recojc_D2J1[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double v_recojc_D2J2[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double v_recojc_D3J1[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double v_recojc_D3J2[nkappa] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

	for(int ik=0; ik<10; ik++){
		isValue = true;
		v_recojc_D1J1[ik] = (pow(recojet0_pt,kappa[ik]) > 0) ? (ijet1candsmom[ik]/(pow(recojet0_pt,kappa[ik]))) : 0 ;
                v_recojc_D1J2[ik] = (pow(recojet1_pt,kappa[ik]) > 0) ? (ijet2candsmom[ik]/(pow(recojet1_pt,kappa[ik]))) : 0 ;
                v_recojc_D2J1[ik] = (ijet1_long_den[ik] > 0) ? (ijet1_long_num[ik]/ijet1_long_den[ik]) : 0 ;
                v_recojc_D2J2[ik] = (ijet2_long_den[ik] > 0) ? (ijet2_long_num[ik]/ijet2_long_den[ik]) : 0 ;
                v_recojc_D3J1[ik] = (ijet1_tran_den[ik] > 0) ? (ijet1_tran_num[ik]/ijet1_tran_den[ik]) : 0 ;
                v_recojc_D3J2[ik] = (ijet2_tran_den[ik] > 0) ? (ijet2_tran_num[ik]/ijet2_tran_den[ik]) : 0 ;
		cout << "Test1 " <<" kappa : "<< ik << " value : "<<v_recojc_D1J1[ik]<<endl;	
		}
*/
for(int itp=0; itp<ntype; itp++){
	for (int iet=0; iet<njetetamn; iet++){
		for (int ik=0; ik<10; ik++){
                	if (isReconstruct) {
			recovar1.clear();
                        	for (int isrc=0; isrc<njecmx; isrc++) {
				recovar.clear();
				if (isrc==0) {isRECO[itp][iet]=false;} 
					//if(isrc==0){isRECOJC = false;}
					//isRECOJC = false;
                                	//if (irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn) {
					if(irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1){
					EventShape_vector  recoevtshape(recomom[isrc][itp][iet], 2.4, 0, 2, 1);
			                recovar =  recoevtshape.getEventShapes();
                			if(isrc==0){recovar1 =  recoevtshape.getEventShapes();}
					if (recovar[nvar]>=2) {
					if (isrc==0) {isRECO[itp][iet] = true;}
			                for (int ij=0; ij<nvar; ij++) {
                    			if (isItUsed(ij)) {
						//isRECOJC = true;
                                		if (isrc==0) {
						if (int(recovar[nvar])>=2) {	
							//isRECOJC = true;
							//cout<<"testing 2"<<endl;
							//test1[ik]=(ijet1candsmom[ik]/(pow(recojet0_pt,kappa[ik])));
							/*
							double v_recojc_D1J1 = (pow(recojet1_pt,kappa[ik]) > 0) ? (ijet1candsmom[ik]/(pow(recojet1_pt,kappa[ik]))) : 0 ;
							double v_recojc_D1J2 = (pow(recojet2_pt,kappa[ik]) > 0) ? (ijet2candsmom[ik]/(pow(recojet2_pt,kappa[ik]))) : 0 ;
							double v_recojc_D2J1 = (ijet1_long_den[ik] > 0) ? (ijet1_long_num[ik]/ijet1_long_den[ik]) : 0 ;
							double v_recojc_D2J2 = (ijet2_long_den[ik] > 0) ? (ijet2_long_num[ik]/ijet2_long_den[ik]) : 0 ;
							double v_recojc_D3J1 = (ijet1_tran_den[ik] > 0) ? (ijet1_tran_num[ik]/ijet1_tran_den[ik]) : 0 ;
							double v_recojc_D3J2 = (ijet2_tran_den[ik] > 0) ? (ijet2_tran_num[ik]/ijet2_tran_den[ik]) : 0 ;
								
							// reco jetcharge 1D	
							//if(isValue){
							h_recojc_D1J1[ik][irecohtjec[isrc]][iet]->Fill(v_recojc_D1J1,weighttrg);
							h_recojc_D1J2[ik][irecohtjec[isrc]][iet]->Fill(v_recojc_D1J2,weighttrg);
							h_recojc_D2J1[ik][irecohtjec[isrc]][iet]->Fill(v_recojc_D2J1,weighttrg);
							h_recojc_D2J2[ik][irecohtjec[isrc]][iet]->Fill(v_recojc_D2J2,weighttrg);
							h_recojc_D3J1[ik][irecohtjec[isrc]][iet]->Fill(v_recojc_D3J1,weighttrg);
							h_recojc_D3J2[ik][irecohtjec[isrc]][iet]->Fill(v_recojc_D3J2,weighttrg);
							//cout << "Test2 " <<" kappa : "<< ik << " value : "<<v_recojc_D1J1<<endl;
							*/
							h_recojc_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1candsmom[ik]/(pow(recojet1_pt,kappa[ik])),weighttrg);
							h_recojc_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2candsmom[ik]/(pow(recojet2_pt,kappa[ik])),weighttrg);
							h_recojc_D2J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1_long_num[ik]/ijet1_long_den[ik]),weighttrg);
                                                        h_recojc_D2J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2_long_num[ik]/ijet2_long_den[ik]),weighttrg);
							h_recojc_D3J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1_tran_num[ik]/ijet1_tran_den[ik]),weighttrg);
                                                        h_recojc_D3J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2_tran_num[ik]/ijet2_tran_den[ik]),weighttrg);

					// reco jetcharge 2D
					int irecbin_D1J1 = RecoBinning2D_D1J1[ik][iet]->GetGlobalBinNumber(ijet1candsmom[ik]/(pow(recojet1_pt,kappa[ik])),leadingptjec[isrc]);
					//int irecbin_D1J1 = RecoBinning2D_D1J1[ik][iet]->GetGlobalBinNumber(v_recojc_D1J1,leadingptjec[isrc]);
					//cout <<" Debug : "<<" kappa : "<<ik<< " value : "<<irecbin_D1J1<<endl;
					h_recovar_2D_D1J1[ik][iet]->Fill(irecbin_D1J1, weighttrg);
					
					int irecbin_D1J2 = RecoBinning2D_D1J2[ik][iet]->GetGlobalBinNumber(ijet2candsmom[ik]/(pow(recojet2_pt,kappa[ik])),leadingptjec[isrc]);	
					//int irecbin_D1J2 = RecoBinning2D_D1J2[ik][iet]->GetGlobalBinNumber(v_recojc_D1J2,leadingptjec[isrc]);
					h_recovar_2D_D1J2[ik][iet]->Fill(irecbin_D1J2, weighttrg);
							
					int irecbin_D2J1 = RecoBinning2D_D2J1[ik][iet]->GetGlobalBinNumber(ijet1_long_num[ik]/ijet1_long_den[ik],leadingptjec[isrc]);
					//int irecbin_D2J1 = RecoBinning2D_D2J1[ik][iet]->GetGlobalBinNumber(v_recojc_D2J1,leadingptjec[isrc]);
					h_recovar_2D_D2J1[ik][iet]->Fill(irecbin_D2J1, weighttrg);

					int irecbin_D2J2 = RecoBinning2D_D2J2[ik][iet]->GetGlobalBinNumber(ijet2_long_num[ik]/ijet2_long_den[ik],leadingptjec[isrc]);
					//int irecbin_D2J2 = RecoBinning2D_D2J2[ik][iet]->GetGlobalBinNumber(v_recojc_D2J2,leadingptjec[isrc]);
                                        h_recovar_2D_D2J2[ik][iet]->Fill(irecbin_D2J2, weighttrg);

					int irecbin_D3J1 = RecoBinning2D_D3J1[ik][iet]->GetGlobalBinNumber(ijet1_tran_num[ik]/ijet1_tran_den[ik],leadingptjec[isrc]);
					//int irecbin_D3J1 = RecoBinning2D_D3J1[ik][iet]->GetGlobalBinNumber(v_recojc_D3J1,leadingptjec[isrc]);
                                        h_recovar_2D_D3J1[ik][iet]->Fill(irecbin_D3J1, weighttrg);

                                        int irecbin_D3J2 = RecoBinning2D_D3J2[ik][iet]->GetGlobalBinNumber(ijet2_tran_num[ik]/ijet2_tran_den[ik],leadingptjec[isrc]);
					//int irecbin_D3J2 = RecoBinning2D_D3J2[ik][iet]->GetGlobalBinNumber(v_recojc_D3J2,leadingptjec[isrc]);
                                        h_recovar_2D_D3J2[ik][iet]->Fill(irecbin_D3J2, weighttrg);

	
							// gluon jets
							/*
							double v_recojc_gjt_D1J1 = (pow(recojet1g_pt,kappa[ik]) > 0) ? (ijet1gcandsmom[ik]/(pow(recojet1g_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_gjt_D1J2 = (pow(recojet2g_pt,kappa[ik]) > 0) ? (ijet2gcandsmom[ik]/(pow(recojet2g_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_gjt_D2J1 = (ijet1g_long_den[ik] > 0) ? (ijet1g_long_num[ik]/ijet1g_long_den[ik]) : 0 ;
                                                        double v_recojc_gjt_D2J2 = (ijet2g_long_den[ik] > 0) ? (ijet2g_long_num[ik]/ijet2g_long_den[ik]) : 0 ;
                                                        double v_recojc_gjt_D3J1 = (ijet1g_tran_den[ik] > 0) ? (ijet1g_tran_num[ik]/ijet1g_tran_den[ik]) : 0 ;
                                                        double v_recojc_gjt_D3J2 = (ijet2g_tran_den[ik] > 0) ? (ijet2g_tran_num[ik]/ijet2g_tran_den[ik]) : 0 ;

                                                        double v_recojc_qgjt_D1J1 = (pow(recojet1qg_pt,kappa[ik]) > 0) ? (ijet1qgcandsmom[ik]/(pow(recojet1qg_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_qgjt_D1J2 = (pow(recojet2qg_pt,kappa[ik]) > 0) ? (ijet2qgcandsmom[ik]/(pow(recojet2qg_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_qgjt_D2J1 = (ijet1qg_long_den[ik] > 0) ? (ijet1qg_long_num[ik]/ijet1qg_long_den[ik]) : 0 ;
                                                        double v_recojc_qgjt_D2J2 = (ijet2qg_long_den[ik] > 0) ? (ijet2qg_long_num[ik]/ijet2qg_long_den[ik]) : 0 ;
                                                        double v_recojc_qgjt_D3J1 = (ijet1qg_tran_den[ik] > 0) ? (ijet1qg_tran_num[ik]/ijet1qg_tran_den[ik]) : 0 ;
                                                        double v_recojc_qgjt_D3J2 = (ijet2qg_tran_den[ik] > 0) ? (ijet2qg_tran_num[ik]/ijet2qg_tran_den[ik]) : 0 ;

                                                        double v_recojc_agjt_D1J1 = (pow(recojet1ag_pt,kappa[ik]) > 0) ? (ijet1agcandsmom[ik]/(pow(recojet1ag_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_agjt_D1J2 = (pow(recojet2ag_pt,kappa[ik]) > 0) ? (ijet2agcandsmom[ik]/(pow(recojet2ag_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_agjt_D2J1 = (ijet1ag_long_den[ik] > 0) ? (ijet1ag_long_num[ik]/ijet1ag_long_den[ik]) : 0 ;
                                                        double v_recojc_agjt_D2J2 = (ijet2ag_long_den[ik] > 0) ? (ijet2ag_long_num[ik]/ijet2ag_long_den[ik]) : 0 ;
                                                        double v_recojc_agjt_D3J1 = (ijet1ag_tran_den[ik] > 0) ? (ijet1ag_tran_num[ik]/ijet1ag_tran_den[ik]) : 0 ;
                                                        double v_recojc_agjt_D3J2 = (ijet2ag_tran_den[ik] > 0) ? (ijet2ag_tran_num[ik]/ijet2ag_tran_den[ik]) : 0 ;
							*/
							
							// 1D	
							h_recojc_gjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1gcandsmom[ik]/(pow(recojet1g_pt,kappa[ik])),weighttrg);
							h_recojc_gjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2gcandsmom[ik]/(pow(recojet2g_pt,kappa[ik])),weighttrg);
							h_recojc_gjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1g_long_num[ik]/ijet1g_long_den[ik]),weighttrg);
							h_recojc_gjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2g_long_num[ik]/ijet2g_long_den[ik]),weighttrg);
							h_recojc_gjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1g_tran_num[ik]/ijet1g_tran_den[ik]),weighttrg);
							h_recojc_gjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2g_tran_num[ik]/ijet2g_tran_den[ik]),weighttrg);

							h_recojc_qgjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qgcandsmom[ik]/(pow(recojet1qg_pt,kappa[ik])),weighttrg);
                                                        h_recojc_qgjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qgcandsmom[ik]/(pow(recojet2qg_pt,kappa[ik])),weighttrg);
                                                        h_recojc_qgjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1qg_long_num[ik]/ijet1qg_long_den[ik]),weighttrg);
                                                        h_recojc_qgjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2qg_long_num[ik]/ijet2qg_long_den[ik]),weighttrg);
                                                        h_recojc_qgjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1qg_tran_num[ik]/ijet1qg_tran_den[ik]),weighttrg);
                                                        h_recojc_qgjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2qg_tran_num[ik]/ijet2qg_tran_den[ik]),weighttrg);
							
							h_recojc_agjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1agcandsmom[ik]/(pow(recojet1ag_pt,kappa[ik])),weighttrg);
                                                        h_recojc_agjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2agcandsmom[ik]/(pow(recojet2ag_pt,kappa[ik])),weighttrg);
                                                        h_recojc_agjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1ag_long_num[ik]/ijet1ag_long_den[ik]),weighttrg);
                                                        h_recojc_agjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2ag_long_num[ik]/ijet2ag_long_den[ik]),weighttrg);
                                                        h_recojc_agjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1ag_tran_num[ik]/ijet1ag_tran_den[ik]),weighttrg);
                                                        h_recojc_agjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2ag_tran_num[ik]/ijet2ag_tran_den[ik]),weighttrg);
							

					int irecbin_gjt_D1J1 = RecoBinning2D_gjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1gcandsmom[ik]/(pow(recojet1g_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_gjt_2D_D1J1[ik][iet]->Fill(irecbin_gjt_D1J1, weighttrg);

					int irecbin_gjt_D1J2 = RecoBinning2D_gjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2gcandsmom[ik]/(pow(recojet2g_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_gjt_2D_D1J2[ik][iet]->Fill(irecbin_gjt_D1J2, weighttrg);

                                        int irecbin_gjt_D2J1 = RecoBinning2D_gjt_D2J1[ik][iet]->GetGlobalBinNumber((ijet1g_long_num[ik]/ijet1g_long_den[ik]),leadingptjec[isrc]);
                                        h_recovar_gjt_2D_D2J1[ik][iet]->Fill(irecbin_gjt_D2J1, weighttrg);

                                        int irecbin_gjt_D2J2 = RecoBinning2D_gjt_D2J2[ik][iet]->GetGlobalBinNumber((ijet2g_long_num[ik]/ijet2g_long_den[ik]),leadingptjec[isrc]);
                                        h_recovar_gjt_2D_D2J2[ik][iet]->Fill(irecbin_gjt_D2J2, weighttrg);

                                        int irecbin_gjt_D3J1 = RecoBinning2D_gjt_D3J1[ik][iet]->GetGlobalBinNumber((ijet1g_tran_num[ik]/ijet1g_tran_den[ik]),leadingptjec[isrc]);
                                        h_recovar_gjt_2D_D3J1[ik][iet]->Fill(irecbin_gjt_D3J1, weighttrg);

                                        int irecbin_gjt_D3J2 = RecoBinning2D_gjt_D3J2[ik][iet]->GetGlobalBinNumber((ijet2g_tran_num[ik]/ijet2g_tran_den[ik]),leadingptjec[isrc]);
                                        h_recovar_gjt_2D_D3J2[ik][iet]->Fill(irecbin_gjt_D3J2, weighttrg);

					
					int irecbin_qgjt_D1J1 = RecoBinning2D_qgjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1qgcandsmom[ik]/(pow(recojet1qg_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_qgjt_2D_D1J1[ik][iet]->Fill(irecbin_qgjt_D1J1, weighttrg);

                                        int irecbin_qgjt_D1J2 = RecoBinning2D_qgjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2qgcandsmom[ik]/(pow(recojet2qg_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_qgjt_2D_D1J2[ik][iet]->Fill(irecbin_qgjt_D1J2, weighttrg);

                                        int irecbin_qgjt_D2J1 = RecoBinning2D_qgjt_D2J1[ik][iet]->GetGlobalBinNumber((ijet1qg_long_num[ik]/ijet1qg_long_den[ik]),leadingptjec[isrc]);
                                        h_recovar_qgjt_2D_D2J1[ik][iet]->Fill(irecbin_qgjt_D2J1, weighttrg);

                                        int irecbin_qgjt_D2J2 = RecoBinning2D_qgjt_D2J2[ik][iet]->GetGlobalBinNumber((ijet2qg_long_num[ik]/ijet2qg_long_den[ik]),leadingptjec[isrc]);
                                        h_recovar_qgjt_2D_D2J2[ik][iet]->Fill(irecbin_qgjt_D2J2, weighttrg);

                                        int irecbin_qgjt_D3J1 = RecoBinning2D_qgjt_D3J1[ik][iet]->GetGlobalBinNumber((ijet1qg_tran_num[ik]/ijet1qg_tran_den[ik]),leadingptjec[isrc]);
                                        h_recovar_qgjt_2D_D3J1[ik][iet]->Fill(irecbin_qgjt_D3J1, weighttrg);

                                        int irecbin_qgjt_D3J2 = RecoBinning2D_qgjt_D3J2[ik][iet]->GetGlobalBinNumber((ijet2qg_tran_num[ik]/ijet2qg_tran_den[ik]),leadingptjec[isrc]);
                                        h_recovar_qgjt_2D_D3J2[ik][iet]->Fill(irecbin_qgjt_D3J2, weighttrg);

					
					int irecbin_agjt_D1J1 = RecoBinning2D_agjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1agcandsmom[ik]/(pow(recojet1ag_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_agjt_2D_D1J1[ik][iet]->Fill(irecbin_agjt_D1J1, weighttrg);

                                        int irecbin_agjt_D1J2 = RecoBinning2D_agjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2agcandsmom[ik]/(pow(recojet2ag_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_agjt_2D_D1J2[ik][iet]->Fill(irecbin_agjt_D1J2, weighttrg);

                                        int irecbin_agjt_D2J1 = RecoBinning2D_agjt_D2J1[ik][iet]->GetGlobalBinNumber((ijet1ag_long_num[ik]/ijet1ag_long_den[ik]),leadingptjec[isrc]);
                                        h_recovar_agjt_2D_D2J1[ik][iet]->Fill(irecbin_agjt_D2J1, weighttrg);

                                        int irecbin_agjt_D2J2 = RecoBinning2D_agjt_D2J2[ik][iet]->GetGlobalBinNumber((ijet2ag_long_num[ik]/ijet2ag_long_den[ik]),leadingptjec[isrc]);
                                        h_recovar_agjt_2D_D2J2[ik][iet]->Fill(irecbin_agjt_D2J2, weighttrg);

                                        int irecbin_agjt_D3J1 = RecoBinning2D_agjt_D3J1[ik][iet]->GetGlobalBinNumber((ijet1ag_tran_num[ik]/ijet1ag_tran_den[ik]),leadingptjec[isrc]);
                                        h_recovar_agjt_2D_D3J1[ik][iet]->Fill(irecbin_agjt_D3J1, weighttrg);

                                        int irecbin_agjt_D3J2 = RecoBinning2D_agjt_D3J2[ik][iet]->GetGlobalBinNumber((ijet2ag_tran_num[ik]/ijet2ag_tran_den[ik]),leadingptjec[isrc]);
                                        h_recovar_agjt_2D_D3J2[ik][iet]->Fill(irecbin_agjt_D3J2, weighttrg);

				
							// b jets	
							/*
							double v_recojc_bjt_D1J1 = (pow(recojet1b_pt,kappa[ik]) > 0) ? (ijet1bcandsmom[ik]/(pow(recojet1b_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_bjt_D1J2 = (pow(recojet2b_pt,kappa[ik]) > 0) ? (ijet2bcandsmom[ik]/(pow(recojet2b_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_bjt_D2J1 = (ijet1b_long_den[ik] > 0) ? (ijet1b_long_num[ik]/ijet1b_long_den[ik]) : 0 ;
                                                        double v_recojc_bjt_D2J2 = (ijet2b_long_den[ik] > 0) ? (ijet2b_long_num[ik]/ijet2b_long_den[ik]) : 0 ;
                                                        double v_recojc_bjt_D3J1 = (ijet1b_tran_den[ik] > 0) ? (ijet1b_tran_num[ik]/ijet1b_tran_den[ik]) : 0 ;
                                                        double v_recojc_bjt_D3J2 = (ijet2b_tran_den[ik] > 0) ? (ijet2b_tran_num[ik]/ijet2b_tran_den[ik]) : 0 ;

							double v_recojc_qbjt_D1J1 = (pow(recojet1qb_pt,kappa[ik]) > 0) ? (ijet1qbcandsmom[ik]/(pow(recojet1qb_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_qbjt_D1J2 = (pow(recojet2qb_pt,kappa[ik]) > 0) ? (ijet2qbcandsmom[ik]/(pow(recojet2qb_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_qbjt_D2J1 = (ijet1qb_long_den[ik] > 0) ? (ijet1qb_long_num[ik]/ijet1qb_long_den[ik]) : 0 ;
                                                        double v_recojc_qbjt_D2J2 = (ijet2qb_long_den[ik] > 0) ? (ijet2qb_long_num[ik]/ijet2qb_long_den[ik]) : 0 ;
                                                        double v_recojc_qbjt_D3J1 = (ijet1qb_tran_den[ik] > 0) ? (ijet1qb_tran_num[ik]/ijet1qb_tran_den[ik]) : 0 ;
                                                        double v_recojc_qbjt_D3J2 = (ijet2qb_tran_den[ik] > 0) ? (ijet2qb_tran_num[ik]/ijet2qb_tran_den[ik]) : 0 ;

							double v_recojc_abjt_D1J1 = (pow(recojet1ab_pt,kappa[ik]) > 0) ? (ijet1abcandsmom[ik]/(pow(recojet1ab_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_abjt_D1J2 = (pow(recojet2ab_pt,kappa[ik]) > 0) ? (ijet2abcandsmom[ik]/(pow(recojet2ab_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_abjt_D2J1 = (ijet1ab_long_den[ik] > 0) ? (ijet1ab_long_num[ik]/ijet1ab_long_den[ik]) : 0 ;
                                                        double v_recojc_abjt_D2J2 = (ijet2ab_long_den[ik] > 0) ? (ijet2ab_long_num[ik]/ijet2ab_long_den[ik]) : 0 ;
                                                        double v_recojc_abjt_D3J1 = (ijet1ab_tran_den[ik] > 0) ? (ijet1ab_tran_num[ik]/ijet1ab_tran_den[ik]) : 0 ;
                                                        double v_recojc_abjt_D3J2 = (ijet2ab_tran_den[ik] > 0) ? (ijet2ab_tran_num[ik]/ijet2ab_tran_den[ik]) : 0 ;
							*/

							// 1D
							h_recojc_bjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1bcandsmom[ik]/(pow(recojet1b_pt,kappa[ik])),weighttrg);
							h_recojc_bjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2bcandsmom[ik]/(pow(recojet2b_pt,kappa[ik])),weighttrg);
							h_recojc_bjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1b_long_num[ik]/ijet1b_long_den[ik],weighttrg);
							h_recojc_bjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2b_long_num[ik]/ijet2b_long_den[ik],weighttrg);
							h_recojc_bjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1b_tran_num[ik]/ijet1b_tran_den[ik],weighttrg);
							h_recojc_bjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2b_tran_num[ik]/ijet2b_tran_den[ik],weighttrg);

							h_recojc_qbjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qbcandsmom[ik]/(pow(recojet1qb_pt,kappa[ik])),weighttrg);
							h_recojc_qbjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qbcandsmom[ik]/(pow(recojet2qb_pt,kappa[ik])),weighttrg);
							h_recojc_qbjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qb_long_num[ik]/ijet1qb_long_den[ik],weighttrg);
							h_recojc_qbjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qb_long_num[ik]/ijet2qb_long_den[ik],weighttrg);
							h_recojc_qbjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qb_tran_num[ik]/ijet1qb_tran_den[ik],weighttrg);
							h_recojc_qbjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qb_tran_num[ik]/ijet2qb_tran_den[ik],weighttrg);

							h_recojc_abjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1abcandsmom[ik]/(pow(recojet1ab_pt,kappa[ik])),weighttrg);
							h_recojc_abjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2abcandsmom[ik]/(pow(recojet2ab_pt,kappa[ik])),weighttrg);
							h_recojc_abjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1ab_long_num[ik]/ijet1ab_long_den[ik],weighttrg);
							h_recojc_abjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2ab_long_num[ik]/ijet2ab_long_den[ik],weighttrg);
							h_recojc_abjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1ab_tran_num[ik]/ijet1ab_tran_den[ik],weighttrg);
							h_recojc_abjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2ab_tran_num[ik]/ijet2ab_tran_den[ik],weighttrg);
							
	
					// 2D 
					int irecbin_bjt_D1J1 = RecoBinning2D_bjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1bcandsmom[ik]/(pow(recojet1b_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_bjt_2D_D1J1[ik][iet]->Fill(irecbin_bjt_D1J1, weighttrg);

                                        int irecbin_bjt_D1J2 = RecoBinning2D_bjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2bcandsmom[ik]/(pow(recojet2b_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_bjt_2D_D1J2[ik][iet]->Fill(irecbin_bjt_D1J2, weighttrg);

                                        int irecbin_bjt_D2J1 = RecoBinning2D_bjt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1b_long_num[ik]/ijet1b_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_bjt_2D_D2J1[ik][iet]->Fill(irecbin_bjt_D2J1, weighttrg);

                                        int irecbin_bjt_D2J2 = RecoBinning2D_bjt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2b_long_num[ik]/ijet2b_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_bjt_2D_D2J2[ik][iet]->Fill(irecbin_bjt_D2J2, weighttrg);

                                        int irecbin_bjt_D3J1 = RecoBinning2D_bjt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1b_tran_num[ik]/ijet1b_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_bjt_2D_D3J1[ik][iet]->Fill(irecbin_bjt_D3J1, weighttrg);

                                        int irecbin_bjt_D3J2 = RecoBinning2D_bjt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2b_tran_num[ik]/ijet2b_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_bjt_2D_D3J2[ik][iet]->Fill(irecbin_bjt_D3J2, weighttrg);

		
					int irecbin_qbjt_D1J1 = RecoBinning2D_qbjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1qbcandsmom[ik]/(pow(recojet1qb_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_qbjt_2D_D1J1[ik][iet]->Fill(irecbin_qbjt_D1J1, weighttrg);

                                        int irecbin_qbjt_D1J2 = RecoBinning2D_qbjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2qbcandsmom[ik]/(pow(recojet2qb_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_qbjt_2D_D1J2[ik][iet]->Fill(irecbin_qbjt_D1J2, weighttrg);

                                        int irecbin_qbjt_D2J1 = RecoBinning2D_qbjt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1qb_long_num[ik]/ijet1qb_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_qbjt_2D_D2J1[ik][iet]->Fill(irecbin_qbjt_D2J1, weighttrg);

                                        int irecbin_qbjt_D2J2 = RecoBinning2D_qbjt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2qb_long_num[ik]/ijet2qb_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_qbjt_2D_D2J2[ik][iet]->Fill(irecbin_qbjt_D2J2, weighttrg);

                                        int irecbin_qbjt_D3J1 = RecoBinning2D_qbjt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1qb_tran_num[ik]/ijet1qb_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_qbjt_2D_D3J1[ik][iet]->Fill(irecbin_qbjt_D3J1, weighttrg);

                                        int irecbin_qbjt_D3J2 = RecoBinning2D_qbjt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2qb_tran_num[ik]/ijet2qb_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_qbjt_2D_D3J2[ik][iet]->Fill(irecbin_qbjt_D3J2, weighttrg);

					
					int irecbin_abjt_D1J1 = RecoBinning2D_abjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1abcandsmom[ik]/(pow(recojet1ab_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_abjt_2D_D1J1[ik][iet]->Fill(irecbin_abjt_D1J1, weighttrg);

                                        int irecbin_abjt_D1J2 = RecoBinning2D_abjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2abcandsmom[ik]/(pow(recojet2ab_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_abjt_2D_D1J2[ik][iet]->Fill(irecbin_abjt_D1J2, weighttrg);

                                        int irecbin_abjt_D2J1 = RecoBinning2D_abjt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1ab_long_num[ik]/ijet1ab_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_abjt_2D_D2J1[ik][iet]->Fill(irecbin_abjt_D2J1, weighttrg);

                                        int irecbin_abjt_D2J2 = RecoBinning2D_abjt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2ab_long_num[ik]/ijet2ab_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_abjt_2D_D2J2[ik][iet]->Fill(irecbin_abjt_D2J2, weighttrg);

                                        int irecbin_abjt_D3J1 = RecoBinning2D_abjt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1ab_tran_num[ik]/ijet1ab_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_abjt_2D_D3J1[ik][iet]->Fill(irecbin_abjt_D3J1, weighttrg);

                                        int irecbin_abjt_D3J2 = RecoBinning2D_abjt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2ab_tran_num[ik]/ijet2ab_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_abjt_2D_D3J2[ik][iet]->Fill(irecbin_abjt_D3J2, weighttrg);


							// c jets
							/*
							double v_recojc_cjt_D1J1 = (pow(recojet1c_pt,kappa[ik]) > 0) ? (ijet1ccandsmom[ik]/(pow(recojet1c_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_cjt_D1J2 = (pow(recojet2c_pt,kappa[ik]) > 0) ? (ijet2ccandsmom[ik]/(pow(recojet2c_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_cjt_D2J1 = (ijet1c_long_den[ik] > 0) ? (ijet1c_long_num[ik]/ijet1c_long_den[ik]) : 0 ;
                                                        double v_recojc_cjt_D2J2 = (ijet2c_long_den[ik] > 0) ? (ijet2c_long_num[ik]/ijet2c_long_den[ik]) : 0 ;
                                                        double v_recojc_cjt_D3J1 = (ijet1c_tran_den[ik] > 0) ? (ijet1c_tran_num[ik]/ijet1c_tran_den[ik]) : 0 ;
                                                        double v_recojc_cjt_D3J2 = (ijet2c_tran_den[ik] > 0) ? (ijet2c_tran_num[ik]/ijet2c_tran_den[ik]) : 0 ;

							double v_recojc_qcjt_D1J1 = (pow(recojet1qc_pt,kappa[ik]) > 0) ? (ijet1qccandsmom[ik]/(pow(recojet1qc_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_qcjt_D1J2 = (pow(recojet2qc_pt,kappa[ik]) > 0) ? (ijet2qccandsmom[ik]/(pow(recojet2qc_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_qcjt_D2J1 = (ijet1qc_long_den[ik] > 0) ? (ijet1qc_long_num[ik]/ijet1qc_long_den[ik]) : 0 ;
                                                        double v_recojc_qcjt_D2J2 = (ijet2qc_long_den[ik] > 0) ? (ijet2qc_long_num[ik]/ijet2qc_long_den[ik]) : 0 ;
                                                        double v_recojc_qcjt_D3J1 = (ijet1qc_tran_den[ik] > 0) ? (ijet1qc_tran_num[ik]/ijet1qc_tran_den[ik]) : 0 ;
                                                        double v_recojc_qcjt_D3J2 = (ijet2qc_tran_den[ik] > 0) ? (ijet2qc_tran_num[ik]/ijet2qc_tran_den[ik]) : 0 ;

							double v_recojc_acjt_D1J1 = (pow(recojet1ac_pt,kappa[ik]) > 0) ? (ijet1accandsmom[ik]/(pow(recojet1ac_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_acjt_D1J2 = (pow(recojet2ac_pt,kappa[ik]) > 0) ? (ijet2accandsmom[ik]/(pow(recojet2ac_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_acjt_D2J1 = (ijet1ac_long_den[ik] > 0) ? (ijet1ac_long_num[ik]/ijet1ac_long_den[ik]) : 0 ;
                                                        double v_recojc_acjt_D2J2 = (ijet2ac_long_den[ik] > 0) ? (ijet2ac_long_num[ik]/ijet2ac_long_den[ik]) : 0 ;
                                                        double v_recojc_acjt_D3J1 = (ijet1ac_tran_den[ik] > 0) ? (ijet1ac_tran_num[ik]/ijet1ac_tran_den[ik]) : 0 ;
                                                        double v_recojc_acjt_D3J2 = (ijet2ac_tran_den[ik] > 0) ? (ijet2ac_tran_num[ik]/ijet2ac_tran_den[ik]) : 0 ;
							*/

							// 1D
                                                        h_recojc_cjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1ccandsmom[ik]/(pow(recojet1c_pt,kappa[ik])),weighttrg);
                                                        h_recojc_cjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2ccandsmom[ik]/(pow(recojet2c_pt,kappa[ik])),weighttrg);
                                                        h_recojc_cjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1c_long_num[ik]/ijet1c_long_den[ik],weighttrg);
                                                        h_recojc_cjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2c_long_num[ik]/ijet2c_long_den[ik],weighttrg);
                                                        h_recojc_cjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1c_tran_num[ik]/ijet1c_tran_den[ik],weighttrg);
                                                        h_recojc_cjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2c_tran_num[ik]/ijet2c_tran_den[ik],weighttrg);

							h_recojc_qcjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qccandsmom[ik]/(pow(recojet1qc_pt,kappa[ik])),weighttrg);
                                                        h_recojc_qcjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qccandsmom[ik]/(pow(recojet2qc_pt,kappa[ik])),weighttrg);
                                                        h_recojc_qcjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qc_long_num[ik]/ijet1qc_long_den[ik],weighttrg);
                                                        h_recojc_qcjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qc_long_num[ik]/ijet2qc_long_den[ik],weighttrg);
                                                        h_recojc_qcjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qc_tran_num[ik]/ijet1qc_tran_den[ik],weighttrg);
                                                        h_recojc_qcjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qc_tran_num[ik]/ijet2qc_tran_den[ik],weighttrg);

							h_recojc_acjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1accandsmom[ik]/(pow(recojet1ac_pt,kappa[ik])),weighttrg);
                                                        h_recojc_acjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2accandsmom[ik]/(pow(recojet2ac_pt,kappa[ik])),weighttrg);
                                                        h_recojc_acjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1ac_long_num[ik]/ijet1ac_long_den[ik],weighttrg);
                                                        h_recojc_acjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2ac_long_num[ik]/ijet2ac_long_den[ik],weighttrg);
                                                        h_recojc_acjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1ac_tran_num[ik]/ijet1ac_tran_den[ik],weighttrg);
                                                        h_recojc_acjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2ac_tran_num[ik]/ijet2ac_tran_den[ik],weighttrg);

					// 2D
					int irecbin_cjt_D1J1 = RecoBinning2D_cjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1ccandsmom[ik]/(pow(recojet1c_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_cjt_2D_D1J1[ik][iet]->Fill(irecbin_cjt_D1J1, weighttrg);

                                        int irecbin_cjt_D1J2 = RecoBinning2D_cjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2ccandsmom[ik]/(pow(recojet2c_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_cjt_2D_D1J2[ik][iet]->Fill(irecbin_cjt_D1J2, weighttrg);

                                        int irecbin_cjt_D2J1 = RecoBinning2D_cjt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1c_long_num[ik]/ijet1c_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_cjt_2D_D2J1[ik][iet]->Fill(irecbin_cjt_D2J1, weighttrg);

                                        int irecbin_cjt_D2J2 = RecoBinning2D_cjt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2c_long_num[ik]/ijet2c_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_cjt_2D_D2J2[ik][iet]->Fill(irecbin_cjt_D2J2, weighttrg);

                                        int irecbin_cjt_D3J1 = RecoBinning2D_cjt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1c_tran_num[ik]/ijet1c_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_cjt_2D_D3J1[ik][iet]->Fill(irecbin_cjt_D3J1, weighttrg);

                                        int irecbin_cjt_D3J2 = RecoBinning2D_cjt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2c_tran_num[ik]/ijet2c_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_cjt_2D_D3J2[ik][iet]->Fill(irecbin_cjt_D3J2, weighttrg);


					int irecbin_qcjt_D1J1 = RecoBinning2D_qcjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1qccandsmom[ik]/(pow(recojet1qc_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_qcjt_2D_D1J1[ik][iet]->Fill(irecbin_qcjt_D1J1, weighttrg);

                                        int irecbin_qcjt_D1J2 = RecoBinning2D_qcjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2qccandsmom[ik]/(pow(recojet2qc_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_qcjt_2D_D1J2[ik][iet]->Fill(irecbin_qcjt_D1J2, weighttrg);

                                        int irecbin_qcjt_D2J1 = RecoBinning2D_qcjt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1qc_long_num[ik]/ijet1qc_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_qcjt_2D_D2J1[ik][iet]->Fill(irecbin_qcjt_D2J1, weighttrg);

                                        int irecbin_qcjt_D2J2 = RecoBinning2D_qcjt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2qc_long_num[ik]/ijet2qc_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_qcjt_2D_D2J2[ik][iet]->Fill(irecbin_qcjt_D2J2, weighttrg);

                                        int irecbin_qcjt_D3J1 = RecoBinning2D_qcjt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1qc_tran_num[ik]/ijet1qc_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_qcjt_2D_D3J1[ik][iet]->Fill(irecbin_qcjt_D3J1, weighttrg);

                                        int irecbin_qcjt_D3J2 = RecoBinning2D_qcjt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2qc_tran_num[ik]/ijet2qc_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_qcjt_2D_D3J2[ik][iet]->Fill(irecbin_qcjt_D3J2, weighttrg);

					
					int irecbin_acjt_D1J1 = RecoBinning2D_acjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1accandsmom[ik]/(pow(recojet1ac_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_acjt_2D_D1J1[ik][iet]->Fill(irecbin_acjt_D1J1, weighttrg);

                                        int irecbin_acjt_D1J2 = RecoBinning2D_acjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2accandsmom[ik]/(pow(recojet2ac_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_acjt_2D_D1J2[ik][iet]->Fill(irecbin_acjt_D1J2, weighttrg);

                                        int irecbin_acjt_D2J1 = RecoBinning2D_acjt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1ac_long_num[ik]/ijet1ac_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_acjt_2D_D2J1[ik][iet]->Fill(irecbin_acjt_D2J1, weighttrg);

                                        int irecbin_acjt_D2J2 = RecoBinning2D_acjt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2ac_long_num[ik]/ijet2ac_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_acjt_2D_D2J2[ik][iet]->Fill(irecbin_acjt_D2J2, weighttrg);

                                        int irecbin_acjt_D3J1 = RecoBinning2D_acjt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1ac_tran_num[ik]/ijet1ac_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_acjt_2D_D3J1[ik][iet]->Fill(irecbin_acjt_D3J1, weighttrg);

                                        int irecbin_acjt_D3J2 = RecoBinning2D_acjt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2ac_tran_num[ik]/ijet2ac_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_acjt_2D_D3J2[ik][iet]->Fill(irecbin_acjt_D3J2, weighttrg);

	
							// s jets
							/*
							double v_recojc_sjt_D1J1 = (pow(recojet1s_pt,kappa[ik]) > 0) ? (ijet1scandsmom[ik]/(pow(recojet1s_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_sjt_D1J2 = (pow(recojet2s_pt,kappa[ik]) > 0) ? (ijet2scandsmom[ik]/(pow(recojet2s_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_sjt_D2J1 = (ijet1s_long_den[ik] > 0) ? (ijet1s_long_num[ik]/ijet1s_long_den[ik]) : 0 ;
                                                        double v_recojc_sjt_D2J2 = (ijet2s_long_den[ik] > 0) ? (ijet2s_long_num[ik]/ijet2s_long_den[ik]) : 0 ;
                                                        double v_recojc_sjt_D3J1 = (ijet1s_tran_den[ik] > 0) ? (ijet1s_tran_num[ik]/ijet1s_tran_den[ik]) : 0 ;
                                                        double v_recojc_sjt_D3J2 = (ijet2s_tran_den[ik] > 0) ? (ijet2s_tran_num[ik]/ijet2s_tran_den[ik]) : 0 ;

							double v_recojc_qsjt_D1J1 = (pow(recojet1qs_pt,kappa[ik]) > 0) ? (ijet1qscandsmom[ik]/(pow(recojet1qs_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_qsjt_D1J2 = (pow(recojet2qs_pt,kappa[ik]) > 0) ? (ijet2qscandsmom[ik]/(pow(recojet2qs_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_qsjt_D2J1 = (ijet1qs_long_den[ik] > 0) ? (ijet1qs_long_num[ik]/ijet1qs_long_den[ik]) : 0 ;
                                                        double v_recojc_qsjt_D2J2 = (ijet2qs_long_den[ik] > 0) ? (ijet2qs_long_num[ik]/ijet2qs_long_den[ik]) : 0 ;
                                                        double v_recojc_qsjt_D3J1 = (ijet1qs_tran_den[ik] > 0) ? (ijet1qs_tran_num[ik]/ijet1qs_tran_den[ik]) : 0 ;
                                                        double v_recojc_qsjt_D3J2 = (ijet2qs_tran_den[ik] > 0) ? (ijet2qs_tran_num[ik]/ijet2qs_tran_den[ik]) : 0 ;

							double v_recojc_asjt_D1J1 = (pow(recojet1as_pt,kappa[ik]) > 0) ? (ijet1ascandsmom[ik]/(pow(recojet1as_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_asjt_D1J2 = (pow(recojet2as_pt,kappa[ik]) > 0) ? (ijet2ascandsmom[ik]/(pow(recojet2as_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_asjt_D2J1 = (ijet1as_long_den[ik] > 0) ? (ijet1as_long_num[ik]/ijet1as_long_den[ik]) : 0 ;
                                                        double v_recojc_asjt_D2J2 = (ijet2as_long_den[ik] > 0) ? (ijet2as_long_num[ik]/ijet2as_long_den[ik]) : 0 ;
                                                        double v_recojc_asjt_D3J1 = (ijet1as_tran_den[ik] > 0) ? (ijet1as_tran_num[ik]/ijet1as_tran_den[ik]) : 0 ;
                                                        double v_recojc_asjt_D3J2 = (ijet2as_tran_den[ik] > 0) ? (ijet2as_tran_num[ik]/ijet2as_tran_den[ik]) : 0 ;
							*/

							// 1D
                                                        h_recojc_sjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1scandsmom[ik]/(pow(recojet1s_pt,kappa[ik])),weighttrg);
                                                        h_recojc_sjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2scandsmom[ik]/(pow(recojet2s_pt,kappa[ik])),weighttrg);
                                                        h_recojc_sjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1s_long_num[ik]/ijet1s_long_den[ik],weighttrg);
                                                        h_recojc_sjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2s_long_num[ik]/ijet2s_long_den[ik],weighttrg);
                                                        h_recojc_sjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1s_tran_num[ik]/ijet1s_tran_den[ik],weighttrg);
                                                        h_recojc_sjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2s_tran_num[ik]/ijet2s_tran_den[ik],weighttrg);

							h_recojc_qsjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qscandsmom[ik]/(pow(recojet1qs_pt,kappa[ik])),weighttrg);
                                                        h_recojc_qsjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qscandsmom[ik]/(pow(recojet2qs_pt,kappa[ik])),weighttrg);
                                                        h_recojc_qsjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qs_long_num[ik]/ijet1qs_long_den[ik],weighttrg);
                                                        h_recojc_qsjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qs_long_num[ik]/ijet2qs_long_den[ik],weighttrg);
                                                        h_recojc_qsjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qs_tran_num[ik]/ijet1qs_tran_den[ik],weighttrg);
                                                        h_recojc_qsjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qs_tran_num[ik]/ijet2qs_tran_den[ik],weighttrg);

							h_recojc_asjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1ascandsmom[ik]/(pow(recojet1as_pt,kappa[ik])),weighttrg);
                                                        h_recojc_asjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2ascandsmom[ik]/(pow(recojet2as_pt,kappa[ik])),weighttrg);
                                                        h_recojc_asjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1as_long_num[ik]/ijet1as_long_den[ik],weighttrg);
                                                        h_recojc_asjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2as_long_num[ik]/ijet2as_long_den[ik],weighttrg);
                                                        h_recojc_asjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1as_tran_num[ik]/ijet1as_tran_den[ik],weighttrg);
                                                        h_recojc_asjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2as_tran_num[ik]/ijet2as_tran_den[ik],weighttrg);


					// 2D
					int irecbin_sjt_D1J1 = RecoBinning2D_sjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1scandsmom[ik]/(pow(recojet1s_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_sjt_2D_D1J1[ik][iet]->Fill(irecbin_sjt_D1J1, weighttrg);

                                        int irecbin_sjt_D1J2 = RecoBinning2D_sjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2scandsmom[ik]/(pow(recojet2s_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_sjt_2D_D1J2[ik][iet]->Fill(irecbin_sjt_D1J2, weighttrg);

                                        int irecbin_sjt_D2J1 = RecoBinning2D_sjt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1s_long_num[ik]/ijet1s_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_sjt_2D_D2J1[ik][iet]->Fill(irecbin_sjt_D2J1, weighttrg);

                                        int irecbin_sjt_D2J2 = RecoBinning2D_sjt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2s_long_num[ik]/ijet2s_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_sjt_2D_D2J2[ik][iet]->Fill(irecbin_sjt_D2J2, weighttrg);

                                        int irecbin_sjt_D3J1 = RecoBinning2D_sjt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1s_tran_num[ik]/ijet1s_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_sjt_2D_D3J1[ik][iet]->Fill(irecbin_sjt_D3J1, weighttrg);

                                        int irecbin_sjt_D3J2 = RecoBinning2D_sjt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2s_tran_num[ik]/ijet2s_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_sjt_2D_D3J2[ik][iet]->Fill(irecbin_sjt_D3J2, weighttrg);


					int irecbin_qsjt_D1J1 = RecoBinning2D_qsjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1qscandsmom[ik]/(pow(recojet1qs_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_qsjt_2D_D1J1[ik][iet]->Fill(irecbin_qsjt_D1J1, weighttrg);

                                        int irecbin_qsjt_D1J2 = RecoBinning2D_qsjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2qscandsmom[ik]/(pow(recojet2qs_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_qsjt_2D_D1J2[ik][iet]->Fill(irecbin_qsjt_D1J2, weighttrg);

                                        int irecbin_qsjt_D2J1 = RecoBinning2D_qsjt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1qs_long_num[ik]/ijet1qs_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_qsjt_2D_D2J1[ik][iet]->Fill(irecbin_qsjt_D2J1, weighttrg);

                                        int irecbin_qsjt_D2J2 = RecoBinning2D_qsjt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2qs_long_num[ik]/ijet2qs_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_qsjt_2D_D2J2[ik][iet]->Fill(irecbin_qsjt_D2J2, weighttrg);

                                        int irecbin_qsjt_D3J1 = RecoBinning2D_qsjt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1qs_tran_num[ik]/ijet1qs_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_qsjt_2D_D3J1[ik][iet]->Fill(irecbin_qsjt_D3J1, weighttrg);

                                        int irecbin_qsjt_D3J2 = RecoBinning2D_qsjt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2qs_tran_num[ik]/ijet2qs_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_qsjt_2D_D3J2[ik][iet]->Fill(irecbin_qsjt_D3J2, weighttrg);

					
					int irecbin_asjt_D1J1 = RecoBinning2D_asjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1ascandsmom[ik]/(pow(recojet1as_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_asjt_2D_D1J1[ik][iet]->Fill(irecbin_asjt_D1J1, weighttrg);

                                        int irecbin_asjt_D1J2 = RecoBinning2D_asjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2ascandsmom[ik]/(pow(recojet2as_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_asjt_2D_D1J2[ik][iet]->Fill(irecbin_asjt_D1J2, weighttrg);

                                        int irecbin_asjt_D2J1 = RecoBinning2D_asjt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1as_long_num[ik]/ijet1as_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_asjt_2D_D2J1[ik][iet]->Fill(irecbin_asjt_D2J1, weighttrg);

                                        int irecbin_asjt_D2J2 = RecoBinning2D_asjt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2as_long_num[ik]/ijet2as_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_asjt_2D_D2J2[ik][iet]->Fill(irecbin_asjt_D2J2, weighttrg);

                                        int irecbin_asjt_D3J1 = RecoBinning2D_asjt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1as_tran_num[ik]/ijet1as_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_asjt_2D_D3J1[ik][iet]->Fill(irecbin_asjt_D3J1, weighttrg);

                                        int irecbin_asjt_D3J2 = RecoBinning2D_asjt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2as_tran_num[ik]/ijet2as_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_asjt_2D_D3J2[ik][iet]->Fill(irecbin_asjt_D3J2, weighttrg);



							// u jets
							/*
							double v_recojc_ujt_D1J1 = (pow(recojet1u_pt,kappa[ik]) > 0) ? (ijet1ucandsmom[ik]/(pow(recojet1u_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_ujt_D1J2 = (pow(recojet2u_pt,kappa[ik]) > 0) ? (ijet2ucandsmom[ik]/(pow(recojet2u_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_ujt_D2J1 = (ijet1u_long_den[ik] > 0) ? (ijet1u_long_num[ik]/ijet1u_long_den[ik]) : 0 ;
                                                        double v_recojc_ujt_D2J2 = (ijet2u_long_den[ik] > 0) ? (ijet2u_long_num[ik]/ijet2u_long_den[ik]) : 0 ;
                                                        double v_recojc_ujt_D3J1 = (ijet1u_tran_den[ik] > 0) ? (ijet1u_tran_num[ik]/ijet1u_tran_den[ik]) : 0 ;
                                                        double v_recojc_ujt_D3J2 = (ijet2u_tran_den[ik] > 0) ? (ijet2u_tran_num[ik]/ijet2u_tran_den[ik]) : 0 ;

							double v_recojc_qujt_D1J1 = (pow(recojet1qu_pt,kappa[ik]) > 0) ? (ijet1qucandsmom[ik]/(pow(recojet1qu_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_qujt_D1J2 = (pow(recojet2qu_pt,kappa[ik]) > 0) ? (ijet2qucandsmom[ik]/(pow(recojet2qu_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_qujt_D2J1 = (ijet1qu_long_den[ik] > 0) ? (ijet1qu_long_num[ik]/ijet1u_long_den[ik]) : 0 ;
                                                        double v_recojc_qujt_D2J2 = (ijet2qu_long_den[ik] > 0) ? (ijet2qu_long_num[ik]/ijet2u_long_den[ik]) : 0 ;
                                                        double v_recojc_qujt_D3J1 = (ijet1qu_tran_den[ik] > 0) ? (ijet1qu_tran_num[ik]/ijet1u_tran_den[ik]) : 0 ;
                                                        double v_recojc_qujt_D3J2 = (ijet2qu_tran_den[ik] > 0) ? (ijet2qu_tran_num[ik]/ijet2u_tran_den[ik]) : 0 ;

							double v_recojc_aujt_D1J1 = (pow(recojet1au_pt,kappa[ik]) > 0) ? (ijet1aucandsmom[ik]/(pow(recojet1au_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_aujt_D1J2 = (pow(recojet2au_pt,kappa[ik]) > 0) ? (ijet2aucandsmom[ik]/(pow(recojet2au_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_aujt_D2J1 = (ijet1au_long_den[ik] > 0) ? (ijet1au_long_num[ik]/ijet1u_long_den[ik]) : 0 ;
                                                        double v_recojc_aujt_D2J2 = (ijet2au_long_den[ik] > 0) ? (ijet2au_long_num[ik]/ijet2u_long_den[ik]) : 0 ;
                                                        double v_recojc_aujt_D3J1 = (ijet1au_tran_den[ik] > 0) ? (ijet1au_tran_num[ik]/ijet1u_tran_den[ik]) : 0 ;
                                                        double v_recojc_aujt_D3J2 = (ijet2au_tran_den[ik] > 0) ? (ijet2au_tran_num[ik]/ijet2u_tran_den[ik]) : 0 ;
							*/

							// 1D
                                                        h_recojc_ujt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1ucandsmom[ik]/(pow(recojet1u_pt,kappa[ik])),weighttrg);
                                                        h_recojc_ujt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2ucandsmom[ik]/(pow(recojet2u_pt,kappa[ik])),weighttrg);
                                                        h_recojc_ujt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1u_long_num[ik]/ijet1u_long_den[ik],weighttrg);
                                                        h_recojc_ujt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2u_long_num[ik]/ijet2u_long_den[ik],weighttrg);
                                                        h_recojc_ujt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1u_tran_num[ik]/ijet1u_tran_den[ik],weighttrg);
                                                        h_recojc_ujt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2u_tran_num[ik]/ijet2u_tran_den[ik],weighttrg);

							h_recojc_qujt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qucandsmom[ik]/(pow(recojet1qu_pt,kappa[ik])),weighttrg);
                                                        h_recojc_qujt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qucandsmom[ik]/(pow(recojet2qu_pt,kappa[ik])),weighttrg);
                                                        h_recojc_qujt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qu_long_num[ik]/ijet1u_long_den[ik],weighttrg);
                                                        h_recojc_qujt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qu_long_num[ik]/ijet2u_long_den[ik],weighttrg);
                                                        h_recojc_qujt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qu_tran_num[ik]/ijet1u_tran_den[ik],weighttrg);
                                                        h_recojc_qujt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qu_tran_num[ik]/ijet2u_tran_den[ik],weighttrg);

							h_recojc_aujt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1aucandsmom[ik]/(pow(recojet1au_pt,kappa[ik])),weighttrg);
                                                        h_recojc_aujt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2aucandsmom[ik]/(pow(recojet2au_pt,kappa[ik])),weighttrg);
                                                        h_recojc_aujt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1au_long_num[ik]/ijet1u_long_den[ik],weighttrg);
                                                        h_recojc_aujt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2au_long_num[ik]/ijet2u_long_den[ik],weighttrg);
                                                        h_recojc_aujt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1au_tran_num[ik]/ijet1u_tran_den[ik],weighttrg);
                                                        h_recojc_aujt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2au_tran_num[ik]/ijet2u_tran_den[ik],weighttrg);



					// 2D
					int irecbin_ujt_D1J1 = RecoBinning2D_ujt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1ucandsmom[ik]/(pow(recojet1u_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_ujt_2D_D1J1[ik][iet]->Fill(irecbin_ujt_D1J1, weighttrg);

                                        int irecbin_ujt_D1J2 = RecoBinning2D_ujt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2ucandsmom[ik]/(pow(recojet2u_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_ujt_2D_D1J2[ik][iet]->Fill(irecbin_ujt_D1J2, weighttrg);

                                        int irecbin_ujt_D2J1 = RecoBinning2D_ujt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1u_long_num[ik]/ijet1u_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_ujt_2D_D2J1[ik][iet]->Fill(irecbin_ujt_D2J1, weighttrg);

                                        int irecbin_ujt_D2J2 = RecoBinning2D_ujt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2u_long_num[ik]/ijet2u_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_ujt_2D_D2J2[ik][iet]->Fill(irecbin_ujt_D2J2, weighttrg);

                                        int irecbin_ujt_D3J1 = RecoBinning2D_ujt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1u_tran_num[ik]/ijet1u_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_ujt_2D_D3J1[ik][iet]->Fill(irecbin_ujt_D3J1, weighttrg);

                                        int irecbin_ujt_D3J2 = RecoBinning2D_ujt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2u_tran_num[ik]/ijet2u_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_ujt_2D_D3J2[ik][iet]->Fill(irecbin_ujt_D3J2, weighttrg);


					int irecbin_qujt_D1J1 = RecoBinning2D_qujt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1qucandsmom[ik]/(pow(recojet1qu_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_qujt_2D_D1J1[ik][iet]->Fill(irecbin_qujt_D1J1, weighttrg);

                                        int irecbin_qujt_D1J2 = RecoBinning2D_qujt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2qucandsmom[ik]/(pow(recojet2qu_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_qujt_2D_D1J2[ik][iet]->Fill(irecbin_qujt_D1J2, weighttrg);

                                        int irecbin_qujt_D2J1 = RecoBinning2D_qujt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1qu_long_num[ik]/ijet1u_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_qujt_2D_D2J1[ik][iet]->Fill(irecbin_qujt_D2J1, weighttrg);

                                        int irecbin_qujt_D2J2 = RecoBinning2D_qujt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2qu_long_num[ik]/ijet2u_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_qujt_2D_D2J2[ik][iet]->Fill(irecbin_qujt_D2J2, weighttrg);

                                        int irecbin_qujt_D3J1 = RecoBinning2D_qujt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1qu_tran_num[ik]/ijet1u_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_qujt_2D_D3J1[ik][iet]->Fill(irecbin_qujt_D3J1, weighttrg);

                                        int irecbin_qujt_D3J2 = RecoBinning2D_qujt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2qu_tran_num[ik]/ijet2u_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_qujt_2D_D3J2[ik][iet]->Fill(irecbin_qujt_D3J2, weighttrg);


					int irecbin_aujt_D1J1 = RecoBinning2D_aujt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1aucandsmom[ik]/(pow(recojet1au_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_aujt_2D_D1J1[ik][iet]->Fill(irecbin_aujt_D1J1, weighttrg);

                                        int irecbin_aujt_D1J2 = RecoBinning2D_aujt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2aucandsmom[ik]/(pow(recojet2au_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_aujt_2D_D1J2[ik][iet]->Fill(irecbin_aujt_D1J2, weighttrg);

                                        int irecbin_aujt_D2J1 = RecoBinning2D_aujt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1au_long_num[ik]/ijet1u_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_aujt_2D_D2J1[ik][iet]->Fill(irecbin_aujt_D2J1, weighttrg);

                                        int irecbin_aujt_D2J2 = RecoBinning2D_aujt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2au_long_num[ik]/ijet2u_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_aujt_2D_D2J2[ik][iet]->Fill(irecbin_aujt_D2J2, weighttrg);

                                        int irecbin_aujt_D3J1 = RecoBinning2D_aujt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1au_tran_num[ik]/ijet1u_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_aujt_2D_D3J1[ik][iet]->Fill(irecbin_aujt_D3J1, weighttrg);

                                        int irecbin_aujt_D3J2 = RecoBinning2D_aujt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2au_tran_num[ik]/ijet2u_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_aujt_2D_D3J2[ik][iet]->Fill(irecbin_aujt_D3J2, weighttrg);

							// d jets
							/*
							double v_recojc_djt_D1J1 = (pow(recojet1d_pt,kappa[ik]) > 0) ? (ijet1dcandsmom[ik]/(pow(recojet1d_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_djt_D1J2 = (pow(recojet2d_pt,kappa[ik]) > 0) ? (ijet2dcandsmom[ik]/(pow(recojet2d_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_djt_D2J1 = (ijet1d_long_den[ik] > 0) ? (ijet1d_long_num[ik]/ijet1d_long_den[ik]) : 0 ;
                                                        double v_recojc_djt_D2J2 = (ijet2d_long_den[ik] > 0) ? (ijet2d_long_num[ik]/ijet2d_long_den[ik]) : 0 ;
                                                        double v_recojc_djt_D3J1 = (ijet1d_tran_den[ik] > 0) ? (ijet1d_tran_num[ik]/ijet1d_tran_den[ik]) : 0 ;
                                                        double v_recojc_djt_D3J2 = (ijet2d_tran_den[ik] > 0) ? (ijet2d_tran_num[ik]/ijet2d_tran_den[ik]) : 0 ;

							double v_recojc_qdjt_D1J1 = (pow(recojet1qd_pt,kappa[ik]) > 0) ? (ijet1qdcandsmom[ik]/(pow(recojet1qd_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_qdjt_D1J2 = (pow(recojet2qd_pt,kappa[ik]) > 0) ? (ijet2qdcandsmom[ik]/(pow(recojet2qd_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_qdjt_D2J1 = (ijet1qd_long_den[ik] > 0) ? (ijet1qd_long_num[ik]/ijet1qd_long_den[ik]) : 0 ;
                                                        double v_recojc_qdjt_D2J2 = (ijet2qd_long_den[ik] > 0) ? (ijet2qd_long_num[ik]/ijet2qd_long_den[ik]) : 0 ;
                                                        double v_recojc_qdjt_D3J1 = (ijet1qd_tran_den[ik] > 0) ? (ijet1qd_tran_num[ik]/ijet1qd_tran_den[ik]) : 0 ;
                                                        double v_recojc_qdjt_D3J2 = (ijet2qd_tran_den[ik] > 0) ? (ijet2qd_tran_num[ik]/ijet2qd_tran_den[ik]) : 0 ;

							double v_recojc_adjt_D1J1 = (pow(recojet1ad_pt,kappa[ik]) > 0) ? (ijet1adcandsmom[ik]/(pow(recojet1ad_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_adjt_D1J2 = (pow(recojet2ad_pt,kappa[ik]) > 0) ? (ijet2adcandsmom[ik]/(pow(recojet2ad_pt,kappa[ik]))) : 0 ;
                                                        double v_recojc_adjt_D2J1 = (ijet1ad_long_den[ik] > 0) ? (ijet1ad_long_num[ik]/ijet1ad_long_den[ik]) : 0 ;
                                                        double v_recojc_adjt_D2J2 = (ijet2ad_long_den[ik] > 0) ? (ijet2ad_long_num[ik]/ijet2ad_long_den[ik]) : 0 ;
                                                        double v_recojc_adjt_D3J1 = (ijet1ad_tran_den[ik] > 0) ? (ijet1ad_tran_num[ik]/ijet1ad_tran_den[ik]) : 0 ;
                                                        double v_recojc_adjt_D3J2 = (ijet2ad_tran_den[ik] > 0) ? (ijet2ad_tran_num[ik]/ijet2ad_tran_den[ik]) : 0 ;
							*/

							// 1D
                                                        h_recojc_djt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1dcandsmom[ik]/(pow(recojet1d_pt,kappa[ik])),weighttrg);
                                                        h_recojc_djt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2dcandsmom[ik]/(pow(recojet2d_pt,kappa[ik])),weighttrg);
                                                        h_recojc_djt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1d_long_num[ik]/ijet1d_long_den[ik],weighttrg);
                                                        h_recojc_djt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2d_long_num[ik]/ijet2d_long_den[ik],weighttrg);
                                                        h_recojc_djt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1d_tran_num[ik]/ijet1d_tran_den[ik],weighttrg);
                                                        h_recojc_djt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2d_tran_num[ik]/ijet2d_tran_den[ik],weighttrg);

							h_recojc_qdjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qdcandsmom[ik]/(pow(recojet1qd_pt,kappa[ik])),weighttrg);
                                                        h_recojc_qdjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qdcandsmom[ik]/(pow(recojet2qd_pt,kappa[ik])),weighttrg);
                                                        h_recojc_qdjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qd_long_num[ik]/ijet1qd_long_den[ik],weighttrg);
                                                        h_recojc_qdjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qd_long_num[ik]/ijet2qd_long_den[ik],weighttrg);
                                                        h_recojc_qdjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1qd_tran_num[ik]/ijet1qd_tran_den[ik],weighttrg);
                                                        h_recojc_qdjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2qd_tran_num[ik]/ijet2qd_tran_den[ik],weighttrg);

							h_recojc_adjt_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1adcandsmom[ik]/(pow(recojet1ad_pt,kappa[ik])),weighttrg);
                                                        h_recojc_adjt_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2adcandsmom[ik]/(pow(recojet2ad_pt,kappa[ik])),weighttrg);
                                                        h_recojc_adjt_D2J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1ad_long_num[ik]/ijet1ad_long_den[ik],weighttrg);
                                                        h_recojc_adjt_D2J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2ad_long_num[ik]/ijet2ad_long_den[ik],weighttrg);
                                                        h_recojc_adjt_D3J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1ad_tran_num[ik]/ijet1ad_tran_den[ik],weighttrg);
                                                        h_recojc_adjt_D3J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2ad_tran_num[ik]/ijet2ad_tran_den[ik],weighttrg);

					// 2D
					int irecbin_djt_D1J1 = RecoBinning2D_djt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1dcandsmom[ik]/(pow(recojet1d_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_djt_2D_D1J1[ik][iet]->Fill(irecbin_djt_D1J1, weighttrg);

                                        int irecbin_djt_D1J2 = RecoBinning2D_djt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2dcandsmom[ik]/(pow(recojet2d_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_djt_2D_D1J2[ik][iet]->Fill(irecbin_djt_D1J2, weighttrg);

                                        int irecbin_djt_D2J1 = RecoBinning2D_djt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1d_long_num[ik]/ijet1d_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_djt_2D_D2J1[ik][iet]->Fill(irecbin_djt_D2J1, weighttrg);

                                        int irecbin_djt_D2J2 = RecoBinning2D_djt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2d_long_num[ik]/ijet2d_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_djt_2D_D2J2[ik][iet]->Fill(irecbin_djt_D2J2, weighttrg);

                                        int irecbin_djt_D3J1 = RecoBinning2D_djt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1d_tran_num[ik]/ijet1d_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_djt_2D_D3J1[ik][iet]->Fill(irecbin_djt_D3J1, weighttrg);

                                        int irecbin_djt_D3J2 = RecoBinning2D_djt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2d_tran_num[ik]/ijet2d_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_djt_2D_D3J2[ik][iet]->Fill(irecbin_djt_D3J2, weighttrg);


					int irecbin_qdjt_D1J1 = RecoBinning2D_qdjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1qdcandsmom[ik]/(pow(recojet1qd_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_qdjt_2D_D1J1[ik][iet]->Fill(irecbin_qdjt_D1J1, weighttrg);

                                        int irecbin_qdjt_D1J2 = RecoBinning2D_qdjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2qdcandsmom[ik]/(pow(recojet2qd_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_qdjt_2D_D1J2[ik][iet]->Fill(irecbin_qdjt_D1J2, weighttrg);

                                        int irecbin_qdjt_D2J1 = RecoBinning2D_qdjt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1qd_long_num[ik]/ijet1qd_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_qdjt_2D_D2J1[ik][iet]->Fill(irecbin_qdjt_D2J1, weighttrg);

                                        int irecbin_qdjt_D2J2 = RecoBinning2D_qdjt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2qd_long_num[ik]/ijet2qd_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_qdjt_2D_D2J2[ik][iet]->Fill(irecbin_qdjt_D2J2, weighttrg);

                                        int irecbin_qdjt_D3J1 = RecoBinning2D_qdjt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1qd_tran_num[ik]/ijet1qd_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_qdjt_2D_D3J1[ik][iet]->Fill(irecbin_qdjt_D3J1, weighttrg);

                                        int irecbin_qdjt_D3J2 = RecoBinning2D_qdjt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2qd_tran_num[ik]/ijet2qd_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_qdjt_2D_D3J2[ik][iet]->Fill(irecbin_qdjt_D3J2, weighttrg);


					int irecbin_adjt_D1J1 = RecoBinning2D_adjt_D1J1[ik][iet]->GetGlobalBinNumber(ijet1adcandsmom[ik]/(pow(recojet1ad_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_adjt_2D_D1J1[ik][iet]->Fill(irecbin_adjt_D1J1, weighttrg);

                                        int irecbin_adjt_D1J2 = RecoBinning2D_adjt_D1J2[ik][iet]->GetGlobalBinNumber(ijet2adcandsmom[ik]/(pow(recojet2ad_pt,kappa[ik])),leadingptjec[isrc]);
                                        h_recovar_adjt_2D_D1J2[ik][iet]->Fill(irecbin_adjt_D1J2, weighttrg);

                                        int irecbin_adjt_D2J1 = RecoBinning2D_adjt_D2J1[ik][iet]->GetGlobalBinNumber(ijet1ad_long_num[ik]/ijet1ad_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_adjt_2D_D2J1[ik][iet]->Fill(irecbin_adjt_D2J1, weighttrg);

                                        int irecbin_adjt_D2J2 = RecoBinning2D_adjt_D2J2[ik][iet]->GetGlobalBinNumber(ijet2ad_long_num[ik]/ijet2ad_long_den[ik],leadingptjec[isrc]);
                                        h_recovar_adjt_2D_D2J2[ik][iet]->Fill(irecbin_adjt_D2J2, weighttrg);

                                        int irecbin_adjt_D3J1 = RecoBinning2D_adjt_D3J1[ik][iet]->GetGlobalBinNumber(ijet1ad_tran_num[ik]/ijet1ad_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_adjt_2D_D3J1[ik][iet]->Fill(irecbin_adjt_D3J1, weighttrg);

                                        int irecbin_adjt_D3J2 = RecoBinning2D_adjt_D3J2[ik][iet]->GetGlobalBinNumber(ijet2ad_tran_num[ik]/ijet2ad_tran_den[ik],leadingptjec[isrc]);
                                        h_recovar_adjt_2D_D3J2[ik][iet]->Fill(irecbin_adjt_D3J2, weighttrg);
		
						}//if (int(recovar[nvar])>=2) {
							//}//if(isValue){
							//cout <<"optimization 2"<<endl;   
							}//if (isrc==0)
							}//if (isItUsed(ij)) {
							}//for (int ij=0; ij<nvar; ij++) {
							}//if (recovar[nvar]>=2) {
							
						}//if (irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn) {
						
					}//for (int isrc=0; isrc<njecmx; isrc++) {
					
				}//if (isReconstruct) {
			                
			if(isMC){
                        	for (int isrc=0; isrc<nGenReso; isrc++) {
				genvar.clear();
              			bool isGEN=false;
					//if(isrc==0){isRECOJC = false;}
					//isGENJC = false;
					//if (isMC && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn){
                                	if (isMC && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1) {
					EventShape_vector  genevtshape(genmom[isrc][itp][iet], 2.4, 0, 2, 1);

		                	genvar =  genevtshape.getEventShapes();
                			if (genvar[nvar]>=2) {
                  			isGEN = true;
                  			for (int ij=0; ij<nvar; ij++) {
                    			if (isItUsed(ij)) {
						//isGENJC = true;
                                		if(isrc==0) {
						if (int(genvar[nvar])>=2) {
							//isGENJC = true;
							/*
							double v_genjc_D1J1 = (pow(genrecojet1_pt,kappa[ik]) > 0) ? (igenjet1candsmom[ik]/(pow(genrecojet1_pt,kappa[ik]))) : 0 ;
							double v_genjc_D1J2 = (pow(genrecojet2_pt,kappa[ik]) > 0) ? (igenjet2candsmom[ik]/(pow(genrecojet2_pt,kappa[ik]))) : 0 ;
							double v_genjc_D2J1 = (igenjet1_long_den[ik] > 0) ? (igenjet1_long_num[ik]/igenjet1_long_den[ik]) : 0 ;
							double v_genjc_D2J2 = (igenjet2_long_den[ik] > 0) ? (igenjet2_long_num[ik]/igenjet2_long_den[ik]) : 0 ;
							double v_genjc_D3J1 = (igenjet1_tran_den[ik] > 0) ? (igenjet1_tran_num[ik]/igenjet1_tran_den[ik]) : 0 ;
							double v_genjc_D3J2 = (igenjet2_tran_den[ik] > 0) ? (igenjet2_tran_num[ik]/igenjet2_tran_den[ik]) : 0 ;
							//cout << " Test2 : "<<igenjet1candsmom[ik] << endl;		
							//cout << " Test2 : "<<genrecojet1_pt<<endl;	
							h_genjc_D1J1[ik][igenhtres[isrc]][iet]->Fill(v_genjc_D1J1,weighttrg);
							h_genjc_D1J2[ik][igenhtres[isrc]][iet]->Fill(v_genjc_D1J2,weighttrg);
							h_genjc_D2J1[ik][igenhtres[isrc]][iet]->Fill(v_genjc_D2J1,weighttrg);
							h_genjc_D2J2[ik][igenhtres[isrc]][iet]->Fill(v_genjc_D2J2,weighttrg);
							h_genjc_D3J1[ik][igenhtres[isrc]][iet]->Fill(v_genjc_D3J1,weighttrg);
							h_genjc_D3J2[ik][igenhtres[isrc]][iet]->Fill(v_genjc_D3J2,weighttrg);															*/
							h_genjc_D1J1[ik][igenhtres[isrc]][iet]->Fill(igenjet1candsmom[ik]/(pow(genrecojet1_pt,kappa[ik])),weighttrg);
							h_genjc_D1J2[ik][igenhtres[isrc]][iet]->Fill(igenjet2candsmom[ik]/(pow(genrecojet2_pt,kappa[ik])),weighttrg);
							h_genjc_D2J1[ik][igenhtres[isrc]][iet]->Fill((igenjet1_long_num[ik]/igenjet1_long_den[ik]),weighttrg);
                                                        h_genjc_D2J2[ik][igenhtres[isrc]][iet]->Fill((igenjet2_long_num[ik]/igenjet2_long_den[ik]),weighttrg);
							h_genjc_D3J1[ik][igenhtres[isrc]][iet]->Fill((igenjet1_tran_num[ik]/igenjet1_tran_den[ik]),weighttrg);
                                                        h_genjc_D3J2[ik][igenhtres[isrc]][iet]->Fill((igenjet2_tran_num[ik]/igenjet2_tran_den[ik]),weighttrg);

					int igenbin_D1J1 = GenBinning2D_D1J1[ik][iet]->GetGlobalBinNumber(igenjet1candsmom[ik]/(pow(genrecojet1_pt,kappa[ik])),leadgenptres[isrc]);
					//int igenbin_D1J1 = GenBinning2D_D1J1[ik][iet]->GetGlobalBinNumber(v_genjc_D1J1,leadgenptres[isrc]);
					h_genvar_2D_D1J1[ik][iet]->Fill(igenbin_D1J1, weighttrg);

					int igenbin_D1J2 = GenBinning2D_D1J2[ik][iet]->GetGlobalBinNumber(igenjet2candsmom[ik]/(pow(genrecojet2_pt,kappa[ik])),leadgenptres[isrc]);
					//int igenbin_D1J2 = GenBinning2D_D1J2[ik][iet]->GetGlobalBinNumber(v_genjc_D1J2,leadgenptres[isrc]);
                                        h_genvar_2D_D1J2[ik][iet]->Fill(igenbin_D1J2, weighttrg);

					int igenbin_D2J1 = GenBinning2D_D2J1[ik][iet]->GetGlobalBinNumber((igenjet1_long_num[ik]/igenjet1_long_den[ik]),leadgenptres[isrc]);
					//int igenbin_D2J1 = GenBinning2D_D2J1[ik][iet]->GetGlobalBinNumber(v_genjc_D2J1,leadgenptres[isrc]);
                                        h_genvar_2D_D2J1[ik][iet]->Fill(igenbin_D2J1, weighttrg);

					int igenbin_D2J2 = GenBinning2D_D2J2[ik][iet]->GetGlobalBinNumber((igenjet2_long_num[ik]/igenjet2_long_den[ik]),leadgenptres[isrc]);
					//int igenbin_D2J2 = GenBinning2D_D2J2[ik][iet]->GetGlobalBinNumber(v_genjc_D2J2,leadgenptres[isrc]);
                                        h_genvar_2D_D2J2[ik][iet]->Fill(igenbin_D2J2, weighttrg);

					int igenbin_D3J1 = GenBinning2D_D3J1[ik][iet]->GetGlobalBinNumber((igenjet1_tran_num[ik]/igenjet1_tran_den[ik]),leadgenptres[isrc]);
					//int igenbin_D3J1 = GenBinning2D_D3J1[ik][iet]->GetGlobalBinNumber(v_genjc_D3J1,leadgenptres[isrc]);
                                        h_genvar_2D_D3J1[ik][iet]->Fill(igenbin_D3J1, weighttrg);

                                        int igenbin_D3J2 = GenBinning2D_D3J2[ik][iet]->GetGlobalBinNumber((igenjet2_tran_num[ik]/igenjet2_tran_den[ik]),leadgenptres[isrc]);
					//int igenbin_D3J2 = GenBinning2D_D3J2[ik][iet]->GetGlobalBinNumber(v_genjc_D3J2,leadgenptres[isrc]);
                                        h_genvar_2D_D3J2[ik][iet]->Fill(igenbin_D3J2, weighttrg);
							}//if (int(genvar[nvar])>=2) {
							
							}//if(isrc==0)
							}//if (isItUsed(ij)) {
							}//for (int ij=0; ij<nvar; ij++) {
							}//if (genvar[nvar]>=2) {
							
						}//if (isMC && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn) {
				if (isrc==0 && isReconstruct){
				//if(igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn){
				//if(igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1 && recomom[isrc][itp][iet].size()>1){
				for(int ij=0; ij<nvar; ij++) {
	                    	if (isItUsed(ij)) {
                      		if(isRECO[itp][iet] && isGEN && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1 && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1) {

				if(recovar1[nvar]>=2 &&  genvar[nvar]>=2){
				//if(irecohtjec[isrc]==igenhtres[isrc] && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn){
				//if(isRECOJC && isGENJC && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1 && recomom[isrc][itp][iet].size()>1){	
				//if(isRECOJC && isGENJC && irecohtjec[isrc]==igenhtres[isrc] && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1 && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1){
                                	/*			
						double v1_recojc_D1J1 = (pow(recojet1_pt,kappa[ik]) > 0) ? (ijet1candsmom[ik]/(pow(recojet1_pt,kappa[ik]))) : 0 ;
                                                double v1_recojc_D1J2 = (pow(recojet2_pt,kappa[ik]) > 0) ? (ijet2candsmom[ik]/(pow(recojet2_pt,kappa[ik]))) : 0 ;
                                                double v1_recojc_D2J1 = (ijet1_long_den[ik] > 0) ? (ijet1_long_num[ik]/ijet1_long_den[ik]) : 0 ;
                                                double v1_recojc_D2J2 = (ijet2_long_den[ik] > 0) ? (ijet2_long_num[ik]/ijet2_long_den[ik]) : 0 ;
                                                double v1_recojc_D3J1 = (ijet1_tran_den[ik] > 0) ? (ijet1_tran_num[ik]/ijet1_tran_den[ik]) : 0 ;
                                                double v1_recojc_D3J2 = (ijet2_tran_den[ik] > 0) ? (ijet2_tran_num[ik]/ijet2_tran_den[ik]) : 0 ;

						double v1_genjc_D1J1 = (pow(genrecojet1_pt,kappa[ik]) > 0) ? (igenjet1candsmom[ik]/(pow(genrecojet1_pt,kappa[ik]))) : 0 ;
                                                double v1_genjc_D1J2 = (pow(genrecojet2_pt,kappa[ik]) > 0) ? (igenjet2candsmom[ik]/(pow(genrecojet2_pt,kappa[ik]))) : 0 ;
                                                double v1_genjc_D2J1 = (igenjet1_long_den[ik] > 0) ? (igenjet1_long_num[ik]/igenjet1_long_den[ik]) : 0 ;
                                                double v1_genjc_D2J2 = (igenjet2_long_den[ik] > 0) ? (igenjet2_long_num[ik]/igenjet2_long_den[ik]) : 0 ;
                                                double v1_genjc_D3J1 = (igenjet1_tran_den[ik] > 0) ? (igenjet1_tran_num[ik]/igenjet1_tran_den[ik]) : 0 ;
                                                double v1_genjc_D3J2 = (igenjet2_tran_den[ik] > 0) ? (igenjet2_tran_num[ik]/igenjet2_tran_den[ik]) : 0 ;

						h_RM_D1J1[ik][irecohtjec[isrc]][iet]->Fill(v1_recojc_D1J1,v1_genjc_D1J1,weighttrg);
						h_RM_D1J2[ik][irecohtjec[isrc]][iet]->Fill(v1_recojc_D1J2,v1_genjc_D1J2,weighttrg);
						h_RM_D2J1[ik][irecohtjec[isrc]][iet]->Fill(v1_recojc_D2J1,v1_genjc_D2J1,weighttrg);
						h_RM_D2J2[ik][irecohtjec[isrc]][iet]->Fill(v1_recojc_D2J2,v1_genjc_D2J2,weighttrg);
						h_RM_D3J1[ik][irecohtjec[isrc]][iet]->Fill(v1_recojc_D3J1,v1_genjc_D3J1,weighttrg);
						h_RM_D3J2[ik][irecohtjec[isrc]][iet]->Fill(v1_recojc_D3J2,v1_genjc_D3J2,weighttrg);
					*/
						h_RM_D1J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1candsmom[ik]/(pow(recojet1_pt,kappa[ik]))),(igenjet1candsmom[ik]/(pow(genrecojet1_pt,kappa[ik]))),weighttrg);
						h_RM_D1J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2candsmom[ik]/(pow(recojet2_pt,kappa[ik]))),(igenjet2candsmom[ik]/(pow(genrecojet2_pt,kappa[ik]))),weighttrg);
						
						h_RM_D2J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1_long_num[ik]/ijet1_long_den[ik]),(igenjet1_long_num[ik]/igenjet1_long_den[ik]),weighttrg);
                                                h_RM_D2J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2_long_num[ik]/ijet2_long_den[ik]),(igenjet2_long_num[ik]/igenjet2_long_den[ik]),weighttrg);
			
						h_RM_D3J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1_tran_num[ik]/ijet1_tran_den[ik]),(igenjet1_tran_num[ik]/igenjet1_tran_den[ik]),weighttrg);
                                                h_RM_D3J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2_tran_num[ik]/ijet2_tran_den[ik]),(igenjet2_tran_num[ik]/igenjet2_tran_den[ik]),weighttrg);


				        int irecbin_D1J1 = RecoBinning2D_D1J1[ik][iet]->GetGlobalBinNumber(ijet1candsmom[ik]/(pow(recojet1_pt,kappa[ik])),leadingptjec[isrc]);
       					int igenbin_D1J1 = GenBinning2D_D1J1[ik][iet]->GetGlobalBinNumber(igenjet1candsmom[ik]/(pow(genrecojet1_pt,kappa[ik])),leadgenptres[isrc]);
					//int irecbin_D1J1 = RecoBinning2D_D1J1[ik][iet]->GetGlobalBinNumber(v1_recojc_D1J1,leadingptjec[isrc]);
                                        //int igenbin_D1J1 = GenBinning2D_D1J1[ik][iet]->GetGlobalBinNumber(v1_genjc_D1J1,leadgenptres[isrc]);
						RM_2D_D1J1[ik][iet]->Fill(irecbin_D1J1, igenbin_D1J1, weighttrg);

                                        int irecbin_D1J2 = RecoBinning2D_D1J2[ik][iet]->GetGlobalBinNumber(ijet2candsmom[ik]/(pow(recojet2_pt,kappa[ik])),leadingptjec[isrc]);
                                	int igenbin_D1J2 = GenBinning2D_D1J2[ik][iet]->GetGlobalBinNumber(igenjet2candsmom[ik]/(pow(genrecojet2_pt,kappa[ik])),leadgenptres[isrc]); 
					//int irecbin_D1J2 = RecoBinning2D_D1J2[ik][iet]->GetGlobalBinNumber(v1_recojc_D1J2,leadingptjec[isrc]);
                                        //int igenbin_D1J2 = GenBinning2D_D1J2[ik][iet]->GetGlobalBinNumber(v1_genjc_D1J2,leadgenptres[isrc]);    
						RM_2D_D1J2[ik][iet]->Fill(irecbin_D1J2, igenbin_D1J2, weighttrg);

                                        int irecbin_D2J1 = RecoBinning2D_D2J1[ik][iet]->GetGlobalBinNumber(ijet1_long_num[ik]/ijet1_long_den[ik],leadingptjec[isrc]);
                                	int igenbin_D2J1 = GenBinning2D_D2J1[ik][iet]->GetGlobalBinNumber((igenjet1_long_num[ik]/igenjet1_long_den[ik]),leadgenptres[isrc]);
					//int irecbin_D2J1 = RecoBinning2D_D2J1[ik][iet]->GetGlobalBinNumber(v1_recojc_D2J1,leadingptjec[isrc]);
                                        //int igenbin_D2J1 = GenBinning2D_D2J1[ik][iet]->GetGlobalBinNumber(v1_genjc_D2J1,leadgenptres[isrc]);
						RM_2D_D2J1[ik][iet]->Fill(irecbin_D2J1, igenbin_D2J1, weighttrg);

                                        int irecbin_D2J2 = RecoBinning2D_D2J2[ik][iet]->GetGlobalBinNumber(ijet2_long_num[ik]/ijet2_long_den[ik],leadingptjec[isrc]);
                                	int igenbin_D2J2 = GenBinning2D_D2J2[ik][iet]->GetGlobalBinNumber((igenjet2_long_num[ik]/igenjet2_long_den[ik]),leadgenptres[isrc]);
					//int irecbin_D2J2 = RecoBinning2D_D2J2[ik][iet]->GetGlobalBinNumber(v1_recojc_D2J2,leadingptjec[isrc]);
                                        //int igenbin_D2J2 = GenBinning2D_D2J2[ik][iet]->GetGlobalBinNumber(v1_genjc_D2J2,leadgenptres[isrc]);
						RM_2D_D2J2[ik][iet]->Fill(irecbin_D2J2, igenbin_D2J2, weighttrg);

                                        int irecbin_D3J1 = RecoBinning2D_D3J1[ik][iet]->GetGlobalBinNumber(ijet1_tran_num[ik]/ijet1_tran_den[ik],leadingptjec[isrc]);
                                    	int igenbin_D3J1 = GenBinning2D_D3J1[ik][iet]->GetGlobalBinNumber((igenjet1_tran_num[ik]/igenjet1_tran_den[ik]),leadgenptres[isrc]);
					//int irecbin_D3J1 = RecoBinning2D_D3J1[ik][iet]->GetGlobalBinNumber(v1_recojc_D3J1,leadingptjec[isrc]);
                                        //int igenbin_D3J1 = GenBinning2D_D3J1[ik][iet]->GetGlobalBinNumber(v1_genjc_D3J1,leadgenptres[isrc]);
						RM_2D_D3J1[ik][iet]->Fill(irecbin_D3J1, igenbin_D3J1, weighttrg);

                                        int irecbin_D3J2 = RecoBinning2D_D3J2[ik][iet]->GetGlobalBinNumber(ijet2_tran_num[ik]/ijet2_tran_den[ik],leadingptjec[isrc]);
  					int igenbin_D3J2 = GenBinning2D_D3J2[ik][iet]->GetGlobalBinNumber((igenjet2_tran_num[ik]/igenjet2_tran_den[ik]),leadgenptres[isrc]); 
					//int irecbin_D3J2 = RecoBinning2D_D3J2[ik][iet]->GetGlobalBinNumber(v1_recojc_D3J2,leadingptjec[isrc]);
                                        //int igenbin_D3J2 = GenBinning2D_D3J2[ik][iet]->GetGlobalBinNumber(v1_genjc_D3J2,leadgenptres[isrc]);
                				RM_2D_D3J2[ik][iet]->Fill(irecbin_D3J2, igenbin_D3J2, weighttrg);                   
					}//if(recovar1[nvar]>=2 &&  genvar[nvar]>=2){	
								}//if(irecohtjec[isrc]==igenhtres[isrc] && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn){
				else {
					//if( irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn){
					//if( isRECOJC && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1){
					//if(irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1){
					if (isRECO[itp][iet] && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1 && recovar1[nvar]>=2) {
				/*	
						double v2_recojc_D1J1 = (pow(recojet1_pt,kappa[ik]) > 0) ? (ijet1candsmom[ik]/(pow(recojet1_pt,kappa[ik]))) : 0 ;
                                                double v2_recojc_D1J2 = (pow(recojet2_pt,kappa[ik]) > 0) ? (ijet2candsmom[ik]/(pow(recojet2_pt,kappa[ik]))) : 0 ;
                                                double v2_recojc_D2J1 = (ijet1_long_den[ik] > 0) ? (ijet1_long_num[ik]/ijet1_long_den[ik]) : 0 ;
                                                double v2_recojc_D2J2 = (ijet2_long_den[ik] > 0) ? (ijet2_long_num[ik]/ijet2_long_den[ik]) : 0 ;
                                                double v2_recojc_D3J1 = (ijet1_tran_den[ik] > 0) ? (ijet1_tran_num[ik]/ijet1_tran_den[ik]) : 0 ;
                                                double v2_recojc_D3J2 = (ijet2_tran_den[ik] > 0) ? (ijet2_tran_num[ik]/ijet2_tran_den[ik]) : 0 ;
							
						h_recofake_D1J1[ik][irecohtjec[isrc]][iet]->Fill(v2_recojc_D1J1,weighttrg);
						h_recofake_D1J2[ik][irecohtjec[isrc]][iet]->Fill(v2_recojc_D1J2,weighttrg);
						h_recofake_D2J1[ik][irecohtjec[isrc]][iet]->Fill(v2_recojc_D2J1,weighttrg);
						h_recofake_D2J2[ik][irecohtjec[isrc]][iet]->Fill(v2_recojc_D2J2,weighttrg);
						h_recofake_D3J1[ik][irecohtjec[isrc]][iet]->Fill(v2_recojc_D3J1,weighttrg);
						h_recofake_D3J2[ik][irecohtjec[isrc]][iet]->Fill(v2_recojc_D3J2,weighttrg);				
		*/
						h_recofake_D1J1[ik][irecohtjec[isrc]][iet]->Fill(ijet1candsmom[ik]/(pow(recojet1_pt,kappa[ik])),weighttrg);
						h_recofake_D1J2[ik][irecohtjec[isrc]][iet]->Fill(ijet2candsmom[ik]/(pow(recojet2_pt,kappa[ik])),weighttrg);
                                                h_recofake_D2J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1_long_num[ik]/ijet1_long_den[ik]),weighttrg);
                                                h_recofake_D2J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2_long_num[ik]/ijet2_long_den[ik]),weighttrg);
                                                h_recofake_D3J1[ik][irecohtjec[isrc]][iet]->Fill((ijet1_tran_num[ik]/ijet1_tran_den[ik]),weighttrg);
                                                h_recofake_D3J2[ik][irecohtjec[isrc]][iet]->Fill((ijet2_tran_num[ik]/ijet2_tran_den[ik]),weighttrg);
			
					int irecbin_D1J1 = RecoBinning2D_D1J1[ik][iet]->GetGlobalBinNumber(ijet1candsmom[ik]/(pow(recojet1_pt,kappa[ik])),leadingptjec[isrc]);
					//int irecbin_D1J1 = RecoBinning2D_D1J1[ik][iet]->GetGlobalBinNumber(v2_recojc_D1J1,leadingptjec[isrc]);
                                        h_recofake_2D_D1J1[ik][iet]->Fill(irecbin_D1J1, weighttrg);

                                        int irecbin_D1J2 = RecoBinning2D_D1J2[ik][iet]->GetGlobalBinNumber(ijet2candsmom[ik]/(pow(recojet2_pt,kappa[ik])),leadingptjec[isrc]);
					//int irecbin_D1J2 = RecoBinning2D_D1J2[ik][iet]->GetGlobalBinNumber(v2_recojc_D1J2,leadingptjec[isrc]);
                                        h_recofake_2D_D1J2[ik][iet]->Fill(irecbin_D1J2, weighttrg);

                                        int irecbin_D2J1 = RecoBinning2D_D2J1[ik][iet]->GetGlobalBinNumber(ijet1_long_num[ik]/ijet1_long_den[ik],leadingptjec[isrc]);
					//int irecbin_D2J1 = RecoBinning2D_D2J1[ik][iet]->GetGlobalBinNumber(v2_recojc_D2J1,leadingptjec[isrc]);
                                        h_recofake_2D_D2J1[ik][iet]->Fill(irecbin_D2J1, weighttrg);

                                        int irecbin_D2J2 = RecoBinning2D_D2J2[ik][iet]->GetGlobalBinNumber(ijet2_long_num[ik]/ijet2_long_den[ik],leadingptjec[isrc]);
					//int irecbin_D2J2 = RecoBinning2D_D2J2[ik][iet]->GetGlobalBinNumber(v2_recojc_D2J2,leadingptjec[isrc]);
                                        h_recofake_2D_D2J2[ik][iet]->Fill(irecbin_D2J2, weighttrg);

                                        int irecbin_D3J1 = RecoBinning2D_D3J1[ik][iet]->GetGlobalBinNumber(ijet1_tran_num[ik]/ijet1_tran_den[ik],leadingptjec[isrc]);
					//int irecbin_D3J1 = RecoBinning2D_D3J1[ik][iet]->GetGlobalBinNumber(v2_recojc_D3J1,leadingptjec[isrc]);
                                        h_recofake_2D_D3J1[ik][iet]->Fill(irecbin_D3J1, weighttrg);

                                        int irecbin_D3J2 = RecoBinning2D_D3J2[ik][iet]->GetGlobalBinNumber(ijet2_tran_num[ik]/ijet2_tran_den[ik],leadingptjec[isrc]);
					//int irecbin_D3J2 = RecoBinning2D_D3J2[ik][iet]->GetGlobalBinNumber(v2_recojc_D3J2,leadingptjec[isrc]);
                                        h_recofake_2D_D3J2[ik][iet]->Fill(irecbin_D3J2, weighttrg);

					}
					//if(igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn){
					//if(isGENJC && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1){
					//if(igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1){
					if (isGEN && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1 && genvar[nvar]>=2) {
		/*
						double v2_genjc_D1J1 = (pow(genrecojet1_pt,kappa[ik]) > 0) ? (igenjet1candsmom[ik]/(pow(genrecojet1_pt,kappa[ik]))) : 0 ;
                                                double v2_genjc_D1J2 = (pow(genrecojet2_pt,kappa[ik]) > 0) ? (igenjet2candsmom[ik]/(pow(genrecojet2_pt,kappa[ik]))) : 0 ;
                                                double v2_genjc_D2J1 = (igenjet1_long_den[ik] > 0) ? (igenjet1_long_num[ik]/igenjet1_long_den[ik]) : 0 ;
                                                double v2_genjc_D2J2 = (igenjet2_long_den[ik] > 0) ? (igenjet2_long_num[ik]/igenjet2_long_den[ik]) : 0 ;
                                                double v2_genjc_D3J1 = (igenjet1_tran_den[ik] > 0) ? (igenjet1_tran_num[ik]/igenjet1_tran_den[ik]) : 0 ;
                                                double v2_genjc_D3J2 = (igenjet2_tran_den[ik] > 0) ? (igenjet2_tran_num[ik]/igenjet2_tran_den[ik]) : 0 ;

						h_genmiss_D1J1[ik][igenhtres[isrc]][iet]->Fill(v2_genjc_D1J1,weighttrg);
						h_genmiss_D1J2[ik][igenhtres[isrc]][iet]->Fill(v2_genjc_D1J2,weighttrg);
						h_genmiss_D2J1[ik][igenhtres[isrc]][iet]->Fill(v2_genjc_D2J1,weighttrg);
						h_genmiss_D2J2[ik][igenhtres[isrc]][iet]->Fill(v2_genjc_D2J2,weighttrg);
						h_genmiss_D3J1[ik][igenhtres[isrc]][iet]->Fill(v2_genjc_D3J1,weighttrg);
						h_genmiss_D3J2[ik][igenhtres[isrc]][iet]->Fill(v2_genjc_D3J2,weighttrg);
*/
						h_genmiss_D1J1[ik][igenhtres[isrc]][iet]->Fill(igenjet1candsmom[ik]/(pow(genrecojet1_pt,kappa[ik])),weighttrg);
						h_genmiss_D1J2[ik][igenhtres[isrc]][iet]->Fill(igenjet2candsmom[ik]/(pow(genrecojet2_pt,kappa[ik])),weighttrg);
                                                h_genmiss_D2J1[ik][igenhtres[isrc]][iet]->Fill((igenjet1_long_num[ik]/igenjet1_long_den[ik]),weighttrg);
                                                h_genmiss_D2J2[ik][igenhtres[isrc]][iet]->Fill((igenjet2_long_num[ik]/igenjet2_long_den[ik]),weighttrg);
                                                h_genmiss_D3J1[ik][igenhtres[isrc]][iet]->Fill((igenjet1_tran_num[ik]/igenjet1_tran_den[ik]),weighttrg);
                                                h_genmiss_D3J2[ik][igenhtres[isrc]][iet]->Fill((igenjet2_tran_num[ik]/igenjet2_tran_den[ik]),weighttrg);

					int igenbin_D1J1 = GenBinning2D_D1J1[ik][iet]->GetGlobalBinNumber(igenjet1candsmom[ik]/(pow(genrecojet1_pt,kappa[ik])),leadgenptres[isrc]);
					//int igenbin_D1J1 = GenBinning2D_D1J1[ik][iet]->GetGlobalBinNumber(v2_genjc_D1J1,leadgenptres[isrc]);
                                        h_genmiss_2D_D1J1[ik][iet]->Fill(igenbin_D1J1, weighttrg);

                                        int igenbin_D1J2 = GenBinning2D_D1J2[ik][iet]->GetGlobalBinNumber(igenjet2candsmom[ik]/(pow(genrecojet2_pt,kappa[ik])),leadgenptres[isrc]);
					//int igenbin_D1J2 = GenBinning2D_D1J2[ik][iet]->GetGlobalBinNumber(v2_genjc_D1J2,leadgenptres[isrc]);
                                        h_genmiss_2D_D1J2[ik][iet]->Fill(igenbin_D1J2, weighttrg);

                                        int igenbin_D2J1 = GenBinning2D_D2J1[ik][iet]->GetGlobalBinNumber((igenjet1_long_num[ik]/igenjet1_long_den[ik]),leadgenptres[isrc]);
					//int igenbin_D2J1 = GenBinning2D_D2J1[ik][iet]->GetGlobalBinNumber(v2_genjc_D2J1,leadgenptres[isrc]);
                                        h_genmiss_2D_D2J1[ik][iet]->Fill(igenbin_D2J1, weighttrg);

                                        int igenbin_D2J2 = GenBinning2D_D2J2[ik][iet]->GetGlobalBinNumber((igenjet2_long_num[ik]/igenjet2_long_den[ik]),leadgenptres[isrc]);
					//int igenbin_D2J2 = GenBinning2D_D2J2[ik][iet]->GetGlobalBinNumber(v2_genjc_D2J2,leadgenptres[isrc]);
                                        h_genmiss_2D_D2J2[ik][iet]->Fill(igenbin_D2J2, weighttrg);

                                        int igenbin_D3J1 = GenBinning2D_D3J1[ik][iet]->GetGlobalBinNumber((igenjet1_tran_num[ik]/igenjet1_tran_den[ik]),leadgenptres[isrc]);
					//int igenbin_D3J1 = GenBinning2D_D3J1[ik][iet]->GetGlobalBinNumber(v2_genjc_D3J1,leadgenptres[isrc]);
                                        h_genmiss_2D_D3J1[ik][iet]->Fill(igenbin_D3J1, weighttrg);

                                        int igenbin_D3J2 = GenBinning2D_D3J2[ik][iet]->GetGlobalBinNumber((igenjet2_tran_num[ik]/igenjet2_tran_den[ik]),leadgenptres[isrc]);
					//int igenbin_D3J2 = GenBinning2D_D3J2[ik][iet]->GetGlobalBinNumber(v2_genjc_D3J2,leadgenptres[isrc]);
                                        h_genmiss_2D_D3J2[ik][iet]->Fill(igenbin_D3J2, weighttrg);
								}//if(igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1){
							}
							}//if (isItUsed(ij))
							}// for(int ij=0; ij<nvar; ij++)	
						}//if (isrc==0 && isReconstruct){
					}//for (int isrc=0; isrc<nGenReso; isrc++) {
				}//if(isMC){
			}//for (int ik=0; ik<10; ik++){
		}//for (int iet=0; iet<njetetamn; iet++){
	}//for(int itp=0; itp<ntype; itp++){
//-----------------------------------------------Calculate And Fill the EventShape Variables--------------------------------

      for(int itp=0; itp<ntype; itp++) {
	for (int iet=0; iet<njetetamn; iet++) {
	  if (isReconstruct) { 
	      recovar1.clear();
	  for (int isrc=0; isrc<njecmx; isrc++) { 
	  //for (int isrc=0; isrc<1; isrc++) { 
	      recovar.clear();
	      //recovar1.clear();
	      if (isrc==0) {isRECO[itp][iet]=false;}
	      
                if (irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1) {
		EventShape_vector  recoevtshape(recomom[isrc][itp][iet], 2.4, 0, 2, 1);
		recovar =  recoevtshape.getEventShapes();
		if(isrc==0){recovar1 =  recoevtshape.getEventShapes();}
		if (recovar[nvar]>=2) {
		  if (isrc==0) {isRECO[itp][iet] = true;}
		  for (int ij=0; ij<nvar; ij++) {
		    if (isItUsed(ij)) { 
		      if (isrc==0) { 
			if (int(recovar[nvar])>=2) {
			nreco++;
                        //if(ij==3 && itp==0 ){cout<<"reco: "<<ievt<<" "<<"Ty:" << itp  << " Nvar : "<<recovar[nvar]<<" "<<recomom[isrc][itp][iet].size() << " Ht2 Bins :" <<irecohtjec[isrc];}
                        //if(itp==0){ cout <<" Var :  " << ij <<" : "<< recovar[ij];}
		        h_recoevtvar[itp][irecohtjec[isrc]][iet][ij]->Fill(recovar[ij], weighttrg); 
                        int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar[ij],aveleadingptjec[isrc]);
			//int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar[ij],leadingptjec[isrc]);
                        h_recovar_2D[itp][iet][ij]->Fill(irecbin, weighttrg);
		        }
			//for (int irand=0; irand<10; irand++) {
			  //if(irand !=k ) h_recoevtvar[irand][itp][irecohtjec[isrc]][iet][ij]->Fill(recovar[ij], weighttrg); 
//#ifdef LHAPDF
			  //for (int ix=1; ix<nnnmx; ix++) {
			  //h_recoevtvarpdf[itp][irecohtjec[isrc]][iet][ij][ix]->Fill(recovar[ij], weighttrg*pdfwt[ix]); 
			  //			 	}
//#endif
			  //}
		      } else {
#ifdef JETENERGY
			if (int(recovar[nvar])>=2) {h_recoevtvarjec[itp][irecohtjec[isrc]][iet][ij][isrc]->Fill(recovar[ij], weighttrg);
                           int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar[ij],aveleadingptjec[isrc]);
			   //int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar[ij],leadingptjec[isrc]);
                           h_recoevtvarjec_2D[itp][iet][ij][isrc]->Fill(irecbin, weighttrg); }
#elif defined(JETRESO)
			if (int(recovar[nvar])>=2) {h_recoevtvarres[itp][irecohtjec[isrc]][iet][ij][isrc]->Fill(recovar[ij], weighttrg);
                           int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar[ij],aveleadingptjec[isrc]);
			   //int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar[ij],leadingptjec[isrc]);
                           h_recoevtvarres_2D[itp][iet][ij][isrc]->Fill(irecbin, weighttrg);}
#endif
		      }
		    }
		  }
		}
	      }
	    }
	  } // if (isReconstruct)
//	  cout << endl;
	  if(isMC) {
	    for (int isrc=0; isrc<nGenReso; isrc++) {
	      //for (int isrc=0; isrc<1; isrc++) { 
	      genvar.clear();
	      bool isGEN=false;
	      if (isMC && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1) { 
		EventShape_vector  genevtshape(genmom[isrc][itp][iet], 2.4, 0, 2, 1);
		
		genvar =  genevtshape.getEventShapes();
		if (genvar[nvar]>=2) {
		  isGEN = true;
		  for (int ij=0; ij<nvar; ij++) {
		    if (isItUsed(ij)) { 
		      if (isrc==0) { 
			if (int(genvar[nvar])>=2) {
			h_genevtvar[itp][igenhtres[isrc]][iet][ij]->Fill(genvar[ij], weighttrg);
			int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
                        //int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], leadgenptres[isrc]);
                        h_genvar_2D[itp][iet][ij]->Fill(igenbin, weighttrg);
                        //if(ij==3 && itp==0 ){cout<<"Gen: "<<ievt<<" "<<"Ty:" << itp  << " Nvar : "<<genvar[nvar]<<" "<<genmom[isrc][itp][iet].size() << " Ht2 Bins :" <<igenhtres[isrc];}
                        //if(itp==0){ cout <<" Var :  " << ij <<" : "<< genvar[ij];}
			} //else {
			  //h_genevtvar2[itp][igenhtres[isrc]][iet][ij]->Fill(genvar[ij], weighttrg);
			//}
#ifdef JETRESO
			//	} else {
			//    	h_genevtvarres[itp][igenhtres[isrc]][iet][ij][isrc]->Fill(genvar[ij], weighttrg);	
#endif
#ifdef LHAPDF
			for (int ix=1; ix<nnnmx; ix++) {
			if (int(genvar[nvar])>=2) {h_genevtvarpdf[itp][igenhtres[isrc]][iet][ij][ix]->Fill(genvar[ij], weighttrg*pdfwt[ix]);
                        int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
			//int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], leadgenptres[isrc]);
                        h_genevtvarpdf_2D[itp][iet][ij][ix]->Fill(igenbin, weighttrg*pdfwt[ix]);   
                         }
			}
#endif
		      }
		    }
		  }
		}
	      }

///cout <<endl;	
  	      if(isrc==0 && isReconstruct){ 
		  for(int ij=0; ij<nvar; ij++) {
		    if (isItUsed(ij)) { 	
		      if(isRECO[itp][iet] && isGEN && irecohtjec[isrc]==igenhtres[isrc] && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1 && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1) { 
			naa++;
	                if(recovar1[nvar]>=2 &&  genvar[nvar]>=2){
 
			 h_2devtvar[itp][irecohtjec[isrc]][iet][ij]->Fill(recovar1[ij], genvar[ij], weighttrg);
			 int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
			 //int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], leadgenptres[isrc]);
                         int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar1[ij],aveleadingptjec[isrc]);
			 //int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar1[ij],leadingptjec[isrc]);
                         RM_2D[itp][iet][ij]->Fill(irecbin, igenbin, weighttrg);
                        }else if (recovar1[nvar]>=2) {
			 
                         //h_2devtvar[itp][igenht][iet][ij]->Fill(recovar[ij],-10.0, weighttrg);	
			  h_recoevtfake[itp][irecohtjec[isrc]][iet][ij]->Fill(recovar1[ij], weighttrg);
                          int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar1[ij],aveleadingptjec[isrc]);
			  //int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar1[ij],leadingptjec[isrc]);
                          h_recofake_2D[itp][iet][ij]->Fill(irecbin, weighttrg);
                        }else if (genvar[nvar]>=2) {
			//h_2devtvar[itp][igenht][iet][ij]->Fill(-10.0, genvar[ij], weighttrg);	//Fill in Reco Underflow
			  h_genevtmiss[itp][igenhtres[isrc]][iet][ij]->Fill(genvar[ij], weighttrg);	
                          int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
			  //int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], leadgenptres[isrc]);
                          h_genmiss_2D[itp][iet][ij]->Fill(igenbin, weighttrg);
                            }
			//  h_2devtvar[itp][0][iet][ij]->Fill(recovar[ij], genvar[ij], weighttrg);
		      } else {
			if (isRECO[itp][iet] && irecohtjec[isrc]>=0 && irecohtjec[isrc]<njetptmn && recomom[isrc][itp][iet].size()>1 && recovar1[nvar]>=2) {
			  nbb++;
			    //h_2devtvar[itp][igenht][iet][ij]->Fill(recovar[ij],-10.0, weighttrg); //Fill Fake in Gen Underflow
			    h_recoevtfake[itp][irecohtjec[isrc]][iet][ij]->Fill(recovar1[ij], weighttrg);
                            int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar1[ij],aveleadingptjec[isrc]);
			    //int irecbin =  RecoBinning2D[itp][iet][ij]->GetGlobalBinNumber(recovar1[ij],leadingptjec[isrc]);
                            h_recofake_2D[itp][iet][ij]->Fill(irecbin, weighttrg);
			}
			if (isGEN && igenhtres[isrc]>=0 && igenhtres[isrc]<njetptmn && genmom[isrc][itp][iet].size()>1 && genvar[nvar]>=2) {
			  ncc++;
                             
			   h_2devtvar[itp][igenht][iet][ij]->Fill(-10.0, genvar[ij], weighttrg);	//Fill Miss in Reco Underflow
			   h_genevtmiss[itp][igenhtres[isrc]][iet][ij]->Fill(genvar[ij], weighttrg);	
			   int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], avegenptres[isrc]);
			   //int igenbin = GenBinning2D[itp][iet][ij]->GetGlobalBinNumber(genvar[ij], leadgenptres[isrc]);
                           h_genmiss_2D[itp][iet][ij]->Fill(igenbin, weighttrg);
			}
		      }
		    } //if (isItUsed(ij)) 
		  } // for(int ij=0; ij<nvar; ij++)	
		} // if (isrc==0 && isReconstruct)
	      //} // if (igenht>=0 && igenht<njetptmn && genmom[isrc][itp][iet].size()>1)
	    } // for (int isrc=0; isrc<nGenReso; isrc++)
	  }//isMC
	} // for (int iet=0; iet<njetetamn; iet++)
      } //for (int itp=0; itp<ntype; itp++) 
      
      //if (nevt%1000==1) { std::cout <<"nevt "<<nevt<<" naa "<<naa<<" nbb "<<nbb<<" ncc "<<ncc<< std::endl;}
      //if(nevt==100){cout <<igenht <<endl;}

      //cout <<"END EVENT"<< endl;
 
      }

// ------------ method called once each job just before starting event loop  ------------
void 
QCDEventShape::beginJob() {
//t1=clock();
  nevt = 0;
  if (isMC) { 
    double dattot[nHLTmx]={0};
    double mctot=0;
    for (int ij=0; ij<npileupmx; ij++) {
      for (int jk=0; jk<nHLTmx; jk++) {dattot[jk] +=datpileup[jk][ij];}
         mctot +=mcpileup[ij];
	}
    
      for (int ij=0; ij<npileupmx; ij++) {
        mcpileup[ij] /=max(1.e-6,mctot);
      	   for (int jk=0; jk<nHLTmx; jk++) {
		datpileup[jk][ij] /=max(1.e-6,dattot[jk]);
		rat_pileup[jk][ij] =  datpileup[jk][ij]/mcpileup[ij];
      		}
    	}
  }

#ifdef JETENERGY
  for (int isrc = 0; isrc < nsrc; isrc++) {
    const char *name = srcnames[isrc];                                         
    //JetCorrectorParameters *p = new JetCorrectorParameters("Summer19UL16APV_RunBCD_V7_DATA_UncertaintySources_AK4PFchs.txt", name);  // data chnage eras for different era
    //JetCorrectorParameters *p = new JetCorrectorParameters("Summer19UL16APV_RunEF_V7_DATA_UncertaintySources_AK4PFchs.txt", name);   // data chnage eras for different era
    //JetCorrectorParameters *p = new JetCorrectorParameters("Summer19UL16_RunFGH_V7_DATA_UncertaintySources_AK4PFchs.txt", name);     // data chnage eras for different era
    JetCorrectorParameters *p = new JetCorrectorParameters("Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.txt", name);           // for MC APV
    //JetCorrectorParameters *p = new JetCorrectorParameters("Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.txt", name);              // for MC
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
    //vsrc[isrc] = unc;vsrc.push_back(unc);
    vsrc.push_back(unc);
    }
#endif  

//cout << "Write test 34 = ok " << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
QCDEventShape::endJob() 
{

     TUnfoldBinng2D->cd();

     // TUnfoldBinning ESV
     for (int ityp=0; ityp<ntype; ityp++) {
      for (int iet=0; iet<njetetamn; iet++) {
        for (int ij=0; ij<nvar; ij++) {
          if (isItUsed(ij)) {

        h_recovar_2D[ityp][iet][ij]->Write();
        h_recofake_2D[ityp][iet][ij]->Write();
        h_genvar_2D[ityp][iet][ij]->Write();
        h_genmiss_2D[ityp][iet][ij]->Write();
        RM_2D[ityp][iet][ij]->Write();

#ifdef  LHAPDF
            for (int ix=1; ix<nnnmx; ix++) {h_genevtvarpdf_2D[ityp][iet][ij][ix]->Write();}
#endif
#ifdef  JETENERGY
            for (int ix=1; ix<njecmx; ix++) {h_recoevtvarjec_2D[ityp][iet][ij][ix]->Write();}
#elif defined(JETRESO)
            for (int ix=1; ix<njecmx; ix++ ) {h_recoevtvarres_2D[ityp][iet][ij][ix]->Write();}
#endif

           }
         }       
       }
     }
  //theFile->cd();
  //theFile->Write();
  //theFile->Close();
  //myfile1->Close();
  //fs->Write();
  //fs->Close();
  
  // TUnfoldBinning Jet Charge
  for (int iet=0; iet<njetetamn; iet++) {
	for(int ik=0; ik<10; ik++){
		h_recovar_2D_D1J1[ik][iet]->Write();
		h_recovar_2D_D1J2[ik][iet]->Write();
		h_recovar_2D_D2J1[ik][iet]->Write();
                h_recovar_2D_D2J2[ik][iet]->Write();
		h_recovar_2D_D3J1[ik][iet]->Write();
                h_recovar_2D_D3J2[ik][iet]->Write();

                h_recovar_gjt_2D_D1J1[ik][iet]->Write();
                h_recovar_gjt_2D_D1J2[ik][iet]->Write();
                h_recovar_gjt_2D_D2J1[ik][iet]->Write();
                h_recovar_gjt_2D_D2J2[ik][iet]->Write();
                h_recovar_gjt_2D_D3J1[ik][iet]->Write();
                h_recovar_gjt_2D_D3J2[ik][iet]->Write();

		h_recovar_qgjt_2D_D1J1[ik][iet]->Write();
                h_recovar_qgjt_2D_D1J2[ik][iet]->Write();
                h_recovar_qgjt_2D_D2J1[ik][iet]->Write();
                h_recovar_qgjt_2D_D2J2[ik][iet]->Write();
                h_recovar_qgjt_2D_D3J1[ik][iet]->Write();
                h_recovar_qgjt_2D_D3J2[ik][iet]->Write();

		h_recovar_agjt_2D_D1J1[ik][iet]->Write();
                h_recovar_agjt_2D_D1J2[ik][iet]->Write();
                h_recovar_agjt_2D_D2J1[ik][iet]->Write();
                h_recovar_agjt_2D_D2J2[ik][iet]->Write();
                h_recovar_agjt_2D_D3J1[ik][iet]->Write();
                h_recovar_agjt_2D_D3J2[ik][iet]->Write();


		h_recovar_bjt_2D_D1J1[ik][iet]->Write();
                h_recovar_bjt_2D_D1J2[ik][iet]->Write();
                h_recovar_bjt_2D_D2J1[ik][iet]->Write();
                h_recovar_bjt_2D_D2J2[ik][iet]->Write();
                h_recovar_bjt_2D_D3J1[ik][iet]->Write();
                h_recovar_bjt_2D_D3J2[ik][iet]->Write();

		h_recovar_qbjt_2D_D1J1[ik][iet]->Write();
                h_recovar_qbjt_2D_D1J2[ik][iet]->Write();
                h_recovar_qbjt_2D_D2J1[ik][iet]->Write();
                h_recovar_qbjt_2D_D2J2[ik][iet]->Write();
                h_recovar_qbjt_2D_D3J1[ik][iet]->Write();
                h_recovar_qbjt_2D_D3J2[ik][iet]->Write();

		h_recovar_abjt_2D_D1J1[ik][iet]->Write();
                h_recovar_abjt_2D_D1J2[ik][iet]->Write();
                h_recovar_abjt_2D_D2J1[ik][iet]->Write();
                h_recovar_abjt_2D_D2J2[ik][iet]->Write();
                h_recovar_abjt_2D_D3J1[ik][iet]->Write();
                h_recovar_abjt_2D_D3J2[ik][iet]->Write();

 	
		h_recovar_cjt_2D_D1J1[ik][iet]->Write();
                h_recovar_cjt_2D_D1J2[ik][iet]->Write();
                h_recovar_cjt_2D_D2J1[ik][iet]->Write();
                h_recovar_cjt_2D_D2J2[ik][iet]->Write();
                h_recovar_cjt_2D_D3J1[ik][iet]->Write();
                h_recovar_cjt_2D_D3J2[ik][iet]->Write();

		h_recovar_qcjt_2D_D1J1[ik][iet]->Write();
                h_recovar_qcjt_2D_D1J2[ik][iet]->Write();
                h_recovar_qcjt_2D_D2J1[ik][iet]->Write();
                h_recovar_qcjt_2D_D2J2[ik][iet]->Write();
                h_recovar_qcjt_2D_D3J1[ik][iet]->Write();
                h_recovar_qcjt_2D_D3J2[ik][iet]->Write();

		h_recovar_acjt_2D_D1J1[ik][iet]->Write();
                h_recovar_acjt_2D_D1J2[ik][iet]->Write();
                h_recovar_acjt_2D_D2J1[ik][iet]->Write();
                h_recovar_acjt_2D_D2J2[ik][iet]->Write();
                h_recovar_acjt_2D_D3J1[ik][iet]->Write();
                h_recovar_acjt_2D_D3J2[ik][iet]->Write();


		h_recovar_sjt_2D_D1J1[ik][iet]->Write();
                h_recovar_sjt_2D_D1J2[ik][iet]->Write();
                h_recovar_sjt_2D_D2J1[ik][iet]->Write();
                h_recovar_sjt_2D_D2J2[ik][iet]->Write();
                h_recovar_sjt_2D_D3J1[ik][iet]->Write();
                h_recovar_sjt_2D_D3J2[ik][iet]->Write();

		h_recovar_qsjt_2D_D1J1[ik][iet]->Write();
                h_recovar_qsjt_2D_D1J2[ik][iet]->Write();
                h_recovar_qsjt_2D_D2J1[ik][iet]->Write();
                h_recovar_qsjt_2D_D2J2[ik][iet]->Write();
                h_recovar_qsjt_2D_D3J1[ik][iet]->Write();
                h_recovar_qsjt_2D_D3J2[ik][iet]->Write();
		
		h_recovar_asjt_2D_D1J1[ik][iet]->Write();
                h_recovar_asjt_2D_D1J2[ik][iet]->Write();
                h_recovar_asjt_2D_D2J1[ik][iet]->Write();
                h_recovar_asjt_2D_D2J2[ik][iet]->Write();
                h_recovar_asjt_2D_D3J1[ik][iet]->Write();
                h_recovar_asjt_2D_D3J2[ik][iet]->Write();


		h_recovar_djt_2D_D1J1[ik][iet]->Write();
                h_recovar_djt_2D_D1J2[ik][iet]->Write();
                h_recovar_djt_2D_D2J1[ik][iet]->Write();
                h_recovar_djt_2D_D2J2[ik][iet]->Write();
                h_recovar_djt_2D_D3J1[ik][iet]->Write();
                h_recovar_djt_2D_D3J2[ik][iet]->Write();
		
		h_recovar_qdjt_2D_D1J1[ik][iet]->Write();
                h_recovar_qdjt_2D_D1J2[ik][iet]->Write();
                h_recovar_qdjt_2D_D2J1[ik][iet]->Write();
                h_recovar_qdjt_2D_D2J2[ik][iet]->Write();
                h_recovar_qdjt_2D_D3J1[ik][iet]->Write();
                h_recovar_qdjt_2D_D3J2[ik][iet]->Write();

		h_recovar_adjt_2D_D1J1[ik][iet]->Write();
                h_recovar_adjt_2D_D1J2[ik][iet]->Write();
                h_recovar_adjt_2D_D2J1[ik][iet]->Write();
                h_recovar_adjt_2D_D2J2[ik][iet]->Write();
                h_recovar_adjt_2D_D3J1[ik][iet]->Write();
                h_recovar_adjt_2D_D3J2[ik][iet]->Write();

		h_recovar_ujt_2D_D1J1[ik][iet]->Write();
                h_recovar_ujt_2D_D1J2[ik][iet]->Write();
                h_recovar_ujt_2D_D2J1[ik][iet]->Write();
                h_recovar_ujt_2D_D2J2[ik][iet]->Write();
                h_recovar_ujt_2D_D3J1[ik][iet]->Write();
                h_recovar_ujt_2D_D3J2[ik][iet]->Write();
		
		h_recovar_qujt_2D_D1J1[ik][iet]->Write();
                h_recovar_qujt_2D_D1J2[ik][iet]->Write();
                h_recovar_qujt_2D_D2J1[ik][iet]->Write();
                h_recovar_qujt_2D_D2J2[ik][iet]->Write();
                h_recovar_qujt_2D_D3J1[ik][iet]->Write();
		h_recovar_qujt_2D_D3J2[ik][iet]->Write();

		h_recovar_aujt_2D_D1J1[ik][iet]->Write();
                h_recovar_aujt_2D_D1J2[ik][iet]->Write();
                h_recovar_aujt_2D_D2J1[ik][iet]->Write();
                h_recovar_aujt_2D_D2J2[ik][iet]->Write();
                h_recovar_aujt_2D_D3J1[ik][iet]->Write();
		h_recovar_aujt_2D_D3J2[ik][iet]->Write();

		h_genvar_2D_D1J1[ik][iet]->Write();
                h_genvar_2D_D1J2[ik][iet]->Write();
                h_genvar_2D_D2J1[ik][iet]->Write();
                h_genvar_2D_D2J2[ik][iet]->Write();
                h_genvar_2D_D3J1[ik][iet]->Write();
                h_genvar_2D_D3J2[ik][iet]->Write();

		h_recofake_2D_D1J1[ik][iet]->Write();
		h_recofake_2D_D1J2[ik][iet]->Write();

		h_recofake_2D_D2J1[ik][iet]->Write();
                h_recofake_2D_D2J2[ik][iet]->Write();

		h_recofake_2D_D3J1[ik][iet]->Write();
                h_recofake_2D_D3J2[ik][iet]->Write();

 		h_genmiss_2D_D1J1[ik][iet]->Write();
		h_genmiss_2D_D1J2[ik][iet]->Write();

		h_genmiss_2D_D2J1[ik][iet]->Write();
                h_genmiss_2D_D2J2[ik][iet]->Write();

		h_genmiss_2D_D3J1[ik][iet]->Write();
                h_genmiss_2D_D3J2[ik][iet]->Write();


		RM_2D_D1J1[ik][iet]->Write();
		RM_2D_D1J2[ik][iet]->Write();

		RM_2D_D2J1[ik][iet]->Write();
                RM_2D_D2J2[ik][iet]->Write();
			
		RM_2D_D3J1[ik][iet]->Write();
                RM_2D_D3J2[ik][iet]->Write();
		
		}//for(int ik=0; ik<10; ik++){
	}//for (int iet=0; iet<njetetamn; iet++) {
/*	
	binning2D->cd();
	for (int iet=0; iet<njetetamn; iet++) {
	        for(int ik=0; ik<10; ik++){
			binsRec2D_D1J1[ik][iet]->Write();
			binsRec2D_D1J2[ik][iet]->Write();

			binsRec2D_D2J1[ik][iet]->Write();
                        binsRec2D_D2J2[ik][iet]->Write();

			binsRec2D_D3J1[ik][iet]->Write();
                        binsRec2D_D3J2[ik][iet]->Write();

			binsGen2D_D1J1[ik][iet]->Write();
			binsGen2D_D1J2[ik][iet]->Write();			

			binsGen2D_D2J1[ik][iet]->Write();
                        binsGen2D_D2J2[ik][iet]->Write();

			binsGen2D_D3J1[ik][iet]->Write();
                        binsGen2D_D3J2[ik][iet]->Write();
			}
		}
*/	
}

// ------------ method called when starting to processes a run  ------------

void 
QCDEventShape::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
// Initialize hltConfig

#ifdef TRIGGER

// cout << "Write test 4 = ok " << endl;
  bool changed(true);
  if (hltPrescaleProvider_.init(iRun,iSetup,theHLTTag.c_str(),changed)) {
  HLTConfigProvider const&  hltConfig = hltPrescaleProvider_.hltConfigProvider();
  hltConfig.dump("Triggers");
  hltConfig.dump("PrescaleTable");

    for (unsigned int ij=0; ij<nHLTmx; ij++) {
      l1pres[ij] = hltpres[ij]=-7;
      }
	} 
	else {
         	}

/*   bool changedConfig;
   if (!hltConfig_.init(iRun, iSetup, theHLTTag.c_str(), changedConfig)) {
     LogError("HLTMuonVal") << "Initialization of HLTConfigProvider failed!!"; 
     return;
 
    //for (unsigned int ij=0; ij<nHLTmx; ij++) {
      //l1pres[ij] = hltpres[ij]=-7;
    //}

  }*/
 
/* 
bool changed(true);
  if (hltConfig_.init(iRun,iSetup,theHLTTag.c_str(),changed)) {
//    cout <<"Trigger tables "<<endl;
    hltConfig_.dump("Triggers");
//    cout <<"Prescale tables "<<endl;
    hltConfig_.dump("PrescaleTable");

    for (unsigned int ij=0; ij<nHLTmx; ij++) {
      l1pres[ij] = hltpres[ij]=-7;
    }

    //trig_init=0; //GMA
    //    // ..
       } else {
              // ..
         }

*/

/* bool changed(true);
   if (hltConfig_.init(iRun,iSetup,theHLTTag.c_str(),changed)) {
     if (changed) {
      // check if trigger name in (new) config
       if (triggerName_!="@") { // "@" means: analyze all triggers in config
     const unsigned int n(hltConfig_.size());
     const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
     if (triggerIndex>=n) {
       LogVerbatim("HLTEventAnalyzerAOD") << "HLTEventAnalyzerAOD::analyze:"
            << " TriggerName " << triggerName_ 
            << " not available in (new) config!" << endl;
       LogVerbatim("HLTEventAnalyzerAOD") << "Available TriggerNames are: " << endl;
       hltConfig_.dump("Triggers");
     }
       }
       hltConfig_.dump("ProcessName");
       hltConfig_.dump("GlobalTag");
       hltConfig_.dump("TableName");
       hltConfig_.dump("Streams");
       hltConfig_.dump("Datasets");
       hltConfig_.dump("PrescaleTable");
       hltConfig_.dump("ProcessPSet");
     }
   } else {
     LogVerbatim("HLTEventAnalyzerAOD") << "HLTEventAnalyzerAOD::analyze:"
      << " config extraction failure with process name "
      << processName_ << endl;
   }
 
*/
#endif
 
  std::cout<<" End of QCDEventShape::beginRun"<<std::endl; //"nevt "<<nevt<<" naa "<<naa<<" nbb "<<nbb<<" ncc "<<ncc<< std::endl;
}


// ------------ method called when ending the processing of a run  ------------
/*
void 
QCDEventShape::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
std::cout<<" End of QCDEventShape::beginRun"<<std::endl;
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
QCDEventShape::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{


}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
QCDEventShape::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QCDEventShape::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  //Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

double PhiInRange(const double& phi) {
      double phiout = phi;
      
      if( phiout > 2*M_PI || phiout < -2*M_PI) {
	phiout = fmod( phiout, 2*M_PI);
      }
      if (phiout <= -M_PI) phiout += 2*M_PI;
      else if (phiout >  M_PI) phiout -= 2*M_PI;
      
      return phiout;
}

template <class T, class U>
double deltaR(const T& t, const U& u) {
return sqrt(pow(t.eta()-u.eta(),2) +pow(PhiInRange(t.phi()-u.phi()),2));
}

// Default jet charge observable NOT USED
double JetCharge1(int charge, double candspt, double jpt, double k) {
	double Q1 = 0.0;
	//double j1 = 0.0;
	Q1 +=1.0*( charge*(pow(candspt,k)));
	//j1 = Q1/pow(jpt,k);
	return Q1/pow(jpt,k);
	//return j1;
}

// For default definition Q
double candsmom(int charge, double candspt, double k){	
	double q = 0.0;
	q =1.0*( charge*(pow(candspt,k)));
	return q;
}

// For longitudinal definition Q(L)
double dotproduct(double candspx, double candspy, double candspz, double jpx, double jpy, double jpz, double jpt, double k) {
	double dot = 0.0;
	//double dotresult =0.0;
	dot = (pow((((candspx*jpx) + (candspy*jpy) + (candspz*jpz))/jpt),k));
	//dotresult = (charge*(pow(dot,k)));
	return dot;
	//return dotresult;
}

// For transverse definition Q(T)
double crossproduct(double candspx, double candspy, double candspz, double jpx, double jpy, double jpz, double jpt, double k){
	double cross = 0.0;
	//double crossresult =0.0;
	cross = (pow(((sqrt((pow(((candspy*jpz) - (candspz*jpy)),2)) + (pow(((candspz*jpx) - (candspx*jpz)),2)) + (pow(((candspx*jpy) - (candspy*jpx)),2))))/jpt),k));
	//crossreult = (charge*(pow(cross,k)));
	return cross;
	//return crossresult;
}

//define this as a plug-in
DEFINE_FWK_MODULE(QCDEventShape);

/*
L1_ZeroBias 	29989	5327	5327	5327	5327	1601	801	801	801	801	801	801	801
17	L1_SingleJet52 	3000	10000	6000	4000	3000	1500	800	500	400	300	150	100	262139
71	HLT_DiPFJetAve80_v2 (2013432) 	7 	7 	7 	7 	7 	7 	7 	7 	7 	7 	7 	7 	7 	L1_SingleJet52
58	HLT_DiPFJetAve140_v2 (2013433) 	2 	2 	2 	2 	1 	1 	1 	1 	1 	1 	1 	1 	1 	L1_SingleJet92
60	HLT_DiPFJetAve200_v2 (2013434) 	250 	250 	125 	85 	60 	35 	16 	12 	9 	6 	4 	2 	1 	L1_SingleJet128
62	HLT_DiPFJetAve260_v2 (2013435) 	85 	85 	85 	60 	42 	24 	12 	8 	6 	4 	2 	1 	1 	L1_SingleJet128 OR L1_SingleJet176
64	HLT_DiPFJetAve320_v2 (2013436) 	15 	15 	15 	10 	6 	4 	2 	1 	1 	1 	1 	1 	1 	L1_SingleJet128 OR L1_SingleJet176
65	HLT_DiPFJetAve400_v2 (2013437) 	5 	5 	5 	3 	2 	1 	1 	1 	1 	1 	1 	1 	1 	L1_SingleJet128 OR L1_SingleJet176
67	HLT_DiPFJetAve500_v2 (2013438) 	1 	1 	1 	1 	1 	1 	1 	1 	1 	1 	1 	1 	1 	L1_SingleJet128 OR L1_SingleJet176
*/
