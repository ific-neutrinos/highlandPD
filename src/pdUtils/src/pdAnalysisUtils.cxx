#include "pdAnalysisUtils.hxx"
#include "TSpline.h"
#include "CategoryManager.hxx"
#include "standardPDTree.hxx"
#include <TH3F.h>
#include <TH2F.h>
#include <TF1.h>
#include <Math/VavilovAccurate.h>

//data for range-momentum conversion, muons
//http://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf divided by LAr density for cm

float Range_grampercm_new[29] = {
  9.833E-1/1.396, 1.786E0/1.396, 3.321E0/1.396, 6.598E0/1.396, 1.058E1/1.396, 3.084E1/1.396, 4.250E1/1.396, 6.732E1/1.396,
  1.063E2/1.396,  1.725E2/1.396, 2.385E2/1.396, 4.934E2/1.396, 6.163E2/1.396, 8.552E2/1.396, 1.202E3/1.396, 1.758E3/1.396,
  2.297E3/1.396,  4.359E3/1.396, 5.354E3/1.396, 7.298E3/1.396, 1.013E4/1.396, 1.469E4/1.396, 1.910E4/1.396, 3.558E4/1.396,
  4.326E4/1.396,  5.768E4/1.396, 7.734E4/1.396, 1.060E5/1.396, 1.307E5/1.396};

float KE_MeV_new[29]= {
  10,    14,    20,    30,    40,     80,     100,    140,    200,   300,
  400,   800,   1000,  1400,  2000,   3000,   4000,   8000,   10000, 14000,
  20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000};

//data for range-momentum conversion, protons

double Range_gpercm_P_Nist[31]={
  1.887E-1/1.396, 3.823E-1/1.396, 6.335E-1/1.396, 1.296/1.396,   2.159/1.396,   7.375/1.396,   1.092E1/1.396, 2.215E1/1.396,
  3.627E1/1.396,  5.282E1/1.396,  7.144E1/1.396,  9.184E1/1.396, 1.138E2/1.396, 1.370E2/1.396, 1.614E2/1.396, 1.869E2/1.396, 
  2.132E2/1.396,  2.403E2/1.396,  2.681E2/1.396,  2.965E2/1.396, 3.254E2/1.396, 3.548E2/1.396, 3.846E2/1.396, 4.148E2/1.396, 
  4.454E2/1.396,  7.626E2/1.396,  1.090E3/1.396,  1.418E3/1.396, 1.745E3/1.396, 2.391E3/1.396, 3.022E3/1.396};

double KE_MeV_P_Nist[31]={
  10,   15,   20,   30,   40,   80,   100, 150,
  200,  250,  300,  350,  400,  450,  500, 550, 
  600,  650,  700,  750,  800,  850,  900, 950, 
  1000, 1500, 2000, 2500, 3000, 4000, 5000};

TGraph const KEvsR(29, Range_grampercm_new, KE_MeV_new);
TSpline3 const KEvsR_spline3("KEvsRS", &KEvsR);

TGraph const RvsKE(29, KE_MeV_new, Range_grampercm_new);
TSpline3 const RvsKE_spline3("RvsKES", &RvsKE);

TGraph const RvsKE_P(31, KE_MeV_P_Nist, Range_gpercm_P_Nist);
TSpline3 const RvsKE_P_spline3("RvsKE_P_S", &RvsKE_P);


TFile* dEdX_template_file = new TFile( (std::string(getenv("PDUTILSROOT"))+"/data/dEdxrestemplates.root").c_str(), "OPEN" );
std::map< int, TProfile* > templates;

TProfile* ProtonTemplate = (TProfile*)dEdX_template_file->Get( "dedx_range_pro" );
TProfile* MuonTemplate   = (TProfile*)dEdX_template_file->Get( "dedx_range_mu" );
TProfile* KaonTemplate   = (TProfile*)dEdX_template_file->Get( "dedx_range_ka" );
/*
templates[ 211 ]  = (TProfile*)dEdX_template_file->Get( "dedx_range_pi"  );
templates[ 321 ]  = (TProfile*)dEdX_template_file->Get( "dedx_range_ka"  );
templates[ 13 ]   = (TProfile*)dEdX_template_file->Get( "dedx_range_mu"  );
templates[ 2212 ] = (TProfile*)dEdX_template_file->Get( "dedx_range_pro" );
*/

//*****************************************************************************
Float_t pdAnaUtils::ComputeRangeMomentum(double trkrange, int pdg){
//*****************************************************************************
  
  /* Muon range-momentum tables from CSDA (Argon density = 1.4 g/cm^3)
     website:
     http://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf
     
     CSDA table values:
     float Range_grampercm[30] = {9.833E-1, 1.786E0, 3.321E0,
     6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1, 1.063E2, 1.725E2,
     2.385E2, 4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3, 2.297E3,
     4.359E3, 5.354E3, 7.298E3, 1.013E4, 1.469E4, 1.910E4, 3.558E4,
     4.326E4, 5.768E4, 7.734E4, 1.060E5, 1.307E5}; float KE_MeV[30] = {10, 14,
     20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000, 1400, 2000, 3000,
     4000, 8000, 10000, 14000, 20000, 30000, 40000, 80000, 100000, 140000,
     200000, 300000, 400000};
     
     Functions below are obtained by fitting polynomial fits to KE_MeV vs
     Range (cm) graph. A better fit was obtained by splitting the graph into
     two: Below range<=200cm,a polynomial of power 4 was a good fit; above
     200cm, a polynomial of power 6 was a good fit
     
     Fit errors for future purposes:
     Below 200cm, Forpoly4 fit: p0 err=1.38533;p1 err=0.209626; p2
     err=0.00650077; p3 err=6.42207E-5; p4 err=1.94893E-7; Above 200cm,
     Forpoly6 fit: p0 err=5.24743;p1 err=0.0176229; p2 err=1.6263E-5; p3
     err=5.9155E-9; p4 err=9.71709E-13; p5 err=7.22381E-17;p6
     err=1.9709E-21;*/
  
  //*********For muon, the calculations are valid up to 1.91E4 cm range
  //corresponding to a Muon KE of 40 GeV**********//
  
  /*Proton range-momentum tables from CSDA (Argon density = 1.4 g/cm^3):
    website: https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html
    
    CSDA values:
    double KE_MeV_P_Nist[31]={10, 15, 20, 30, 40, 80, 100, 150, 200, 250, 300,
    350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
    1500, 2000, 2500, 3000, 4000, 5000};
    
    double Range_gpercm_P_Nist[31]={1.887E-1,3.823E-1, 6.335E-1, 1.296,
    2.159, 7.375, 1.092E1, 2.215E1, 3.627E1, 5.282E1, 7.144E1,
    9.184E1, 1.138E2, 1.370E2, 1.614E2, 1.869E2, 2.132E2, 2.403E2,
    2.681E2, 2.965E2, 3.254E2, 3.548E2, 3.846E2, 4.148E2, 4.454E2,
    7.626E2, 1.090E3, 1.418E3, 1.745E3, 2.391E3, 3.022E3};
    
    Functions below are obtained by fitting power and polynomial fits to
    KE_MeV vs Range (cm) graph. A better fit was obtained by splitting the
    graph into two: Below range<=80cm,a a*(x^b) was a good fit; above 80cm, a
    polynomial of power 6 was a good fit
    
    Fit errors for future purposes:
    For power function fit: a=0.388873; and b=0.00347075
    Forpoly6 fit: p0 err=3.49729;p1 err=0.0487859; p2 err=0.000225834; p3
    err=4.45542E-7; p4 err=4.16428E-10; p5 err=1.81679E-13;p6
    err=2.96958E-17;*/
  
  //*********For proton, the calculations are valid up to 3.022E3 cm range
  //corresponding to a Muon KE of 5 GeV**********//

  if (trkrange < 0 || std::isnan(trkrange)) {
    //    std::cout << "TrackMomentumCalculator   " 
    //              << "Invalid track range " << trkrange << " return -1" << std::endl;
    return -1.;
  }

    
  double KE, Momentum, M;
  constexpr double Muon_M = 105.7, Proton_M = 938.272;
  
  if (abs(pdg) == 13) {
    M = Muon_M;
    KE = KEvsR_spline3.Eval(trkrange);
  } else if (abs(pdg) == 2212) {
    M = Proton_M;
    if (trkrange > 0 && trkrange <= 80)
      KE = 29.9317 * std::pow(trkrange, 0.586304);
    else if (trkrange > 80 && trkrange <= 3.022E3)
      KE =
        149.904 + (3.34146 * trkrange) + (-0.00318856 * trkrange * trkrange) +
        (4.34587E-6 * trkrange * trkrange * trkrange) +
        (-3.18146E-9 * trkrange * trkrange * trkrange * trkrange) +
        (1.17854E-12 * trkrange * trkrange * trkrange * trkrange * trkrange) +
        (-1.71763E-16 * trkrange * trkrange * trkrange * trkrange * trkrange *
         trkrange);
    else
      KE = -999;
  } else
    KE = -999;
  
  if (KE < 0)
    Momentum = -999;
  else
    Momentum = std::sqrt((KE * KE) + (2 * M * KE));
  
  Momentum = Momentum / 1000;
  
  return Momentum;
}

//*****************************************************************************
Float_t pdAnaUtils::ComputeCSDARange(double beammom, int pdg){
//*****************************************************************************
    
  if (beammom < 0 || std::isnan(beammom)) {
    //    std::cout << "CSDARangeCalculator   " 
    //              << "Invalid beam mom " << beammom << " return -1" << std::endl;
    return -1.;
  }

  double KE, M, CSDARange;
  constexpr double Muon_M = 105.7, Proton_M = 938.272;
  
  if (abs(pdg) == 13) {
    M = Muon_M;
    KE = sqrt(beammom * beammom + M * M) - M;
    CSDARange = RvsKE_spline3.Eval(KE);
  } else if (abs(pdg) == 2212) {
    M = Proton_M;
    KE = sqrt(beammom * beammom + M * M) - M;
    CSDARange = RvsKE_P_spline3.Eval(KE);
  }
  else 
    CSDARange=-1;
  
  return CSDARange;
}

//********************************************************************
Float_t pdAnaUtils::ComputeKineticEnergy(const AnaParticlePD &part) {
//********************************************************************

  int plane = 2;

  if(part.Hits[plane].size() < 1)return -1;
  
  double kinetic = 0;
  double res     = 0;

  for(size_t i = 0; i < part.Hits[plane].size(); i++){
    if(part.Hits[plane][i].dEdx_calib > 1000. || part.Hits[plane][i].dEdx_calib==-999)continue;
    double dedxi = part.Hits[plane][i].dEdx_calib;
    double Residualrangei = part.Hits[plane][i].ResidualRange;
    kinetic = kinetic + dedxi * fabs(Residualrangei - res);
    res = Residualrangei;
  }

  // convert to GeV
  return kinetic/units::GeV;
}

//********************************************************************
void pdAnaUtils::ComputeDistanceToVertex(AnaParticlePD* part, std::vector<Float_t>& distance){
//********************************************************************

  distance.clear();
  if (!part) return;
  if (part->Daughters.empty()) return;

  for(size_t i = 0; i < part->Daughters.size(); i++){
    AnaParticlePD* dau = static_cast<AnaParticlePD*>(part->Daughters[i]);

    //compute distance between part and daughter. Assume dau can be reconstructed backwards
    Float_t dis1 = sqrt(pow(part->PositionEnd[0]-dau->PositionStart[0],2)+pow(part->PositionEnd[1]-dau->PositionStart[1],2)+pow(part->PositionEnd[2]-dau->PositionStart[2],2));
    Float_t dis2 = sqrt(pow(part->PositionEnd[0]-dau->PositionEnd[0],2)+pow(part->PositionEnd[1]-dau->PositionEnd[1],2)+pow(part->PositionEnd[2]-dau->PositionEnd[2],2));

    //return lowest value
    if(dis1<dis2)distance.push_back(dis1);
    else distance.push_back(dis2);
  }
}

//********************************************************************
AnaTrueParticle* pdAnaUtils::GetBeamTrueParticle(const AnaSpillB& spill){  
//********************************************************************

  AnaTrueParticle* beampart = NULL;
  
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill.Beam);
  AnaParticleMomB* beamPart = beam->BeamParticle;

  int true_id = 0;
  if(beamPart)
    if(beamPart->TrueObject)
      true_id = static_cast<AnaTrueParticleB*>(beamPart->TrueObject)->ID;
  
  if(spill.TrueParticles.size() > 0){
    for(int i = 0; i < (int)spill.TrueParticles.size(); i++){
      if(true_id == spill.TrueParticles[i]->ID){
        beampart = static_cast<AnaTrueParticle*>(spill.TrueParticles[i]);
        break;
      }
    }
  }

  return beampart;
}

//*****************************************************************************
AnaTrueParticlePD* pdAnaUtils::GetTrueParticle(AnaEventB* event, Int_t ID){
//*****************************************************************************

  // Get all reconstructed tracks in the event
  AnaTrueParticleB** trueParticles = event->TrueParticles;
  Int_t nTrueParts                 = event->nTrueParticles;

  for (Int_t i=0;i<nTrueParts;i++){
    if(trueParticles[i]->ID == ID){
	return static_cast<AnaTrueParticlePD*>(trueParticles[i]);
      }
    }
  return NULL;
}

//*****************************************************************************
AnaTrueParticlePD* pdAnaUtils::GetTrueParticle(const std::vector<AnaTrueParticleB*>& trueParticles, Int_t ID){
//*****************************************************************************
 
  // Get all reconstructed tracks in the event
  for (UInt_t i=0;i<trueParticles.size();i++){
    if (trueParticles[i]->ID == ID){
      return static_cast<AnaTrueParticlePD*>(trueParticles[i]);
    }
  }

  return NULL;
}

//*****************************************************************************
AnaParticlePD* pdAnaUtils::GetRecoParticleWithAssociatedTrueID(const std::vector<AnaParticleB*> particles, Int_t true_ID){
//*****************************************************************************
 
  // loop over reconstructed tracks
  for(UInt_t i = 0; i < particles.size(); i++){
    AnaTrueParticlePD* truepart = static_cast<AnaTrueParticlePD*>(particles[i]->TrueObject);
    if(!truepart)continue;
    if(truepart->ID == true_ID)
      return static_cast<AnaParticlePD*>(particles[i]);
  }

  return NULL;
}

//*****************************************************************************
std::pair< double, int > pdAnaUtils::Chi2PID(const AnaParticlePD& part, const int pdg ){
//*****************************************************************************	

  double pid_chi2 = 0.; 
  int npt = 0;

  Int_t plane=2;

  TProfile* profile;
  
  if(pdg == 2212)profile = ProtonTemplate;
  else if(pdg == 13)profile = MuonTemplate;
  else if(pdg == 321)profile = KaonTemplate;
  else{
    std::cout << "no profile for pdg " << pdg << std::endl;
    return std::make_pair(9999., -1);
  }
  
  if( part.Hits[plane].size() < 1 )
    return std::make_pair(9999., -1);

  //Ignore first and last point
  for( UInt_t i = 1; i < part.Hits[plane].size()-1; ++i ){
    //Skip large pulse heights
    if( part.Hits[plane][i].dEdx > 1000. || part.Hits[plane][i].dEdx==-999)
      continue;

    int bin = profile->FindBin( part.Hits[plane][i].ResidualRange );

    if( bin >= 1 && bin <= profile->GetNbinsX() ){
      
      double template_dedx = profile->GetBinContent( bin );
      if( template_dedx < 1.e-6 ){
        template_dedx = ( profile->GetBinContent( bin - 1 ) + profile->GetBinContent( bin + 1 ) ) / 2.;        
      }
      
      
      double template_dedx_err = profile->GetBinError( bin );
      if( template_dedx_err < 1.e-6 ){
        template_dedx_err = ( profile->GetBinError( bin - 1 ) + profile->GetBinError( bin + 1 ) ) / 2.;        
      }

      double dedx_res = 0.04231 + 0.0001783 * part.Hits[plane][i].dEdx * part.Hits[plane][i].dEdx;      
      dedx_res *= part.Hits[plane][i].dEdx; 
      
      
      //Chi2 += ( track_dedx - template_dedx )^2  / ( (template_dedx_err)^2 + (dedx_res)^2 )      
      pid_chi2 += ( pow( (part.Hits[plane][i].dEdx - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) ); 
            
      ++npt;      
    }	
  }
		
  if( npt == 0 )	
    return std::make_pair(9999., -1);
	  	
  return std::make_pair(pid_chi2, npt); 	
}

//*****************************************************************************
std::pair< double, int > pdAnaUtils::Chi2PID_UpToRR(const AnaParticlePD& part, const int pdg, const double RR){
//*****************************************************************************	

  double pid_chi2 = 0.; 
  int npt = 0;

  Int_t plane=2;

  TProfile* profile;
  
  if(pdg == 2212)profile = ProtonTemplate;
  else if(pdg == 13)profile = MuonTemplate;
  else if(pdg == 321)profile = KaonTemplate;
  else{
    std::cout << "no profile for pdg " << pdg << std::endl;
    return std::make_pair(9999., -1);
  }
  
  if( part.Hits[plane].size() < 1 )
    return std::make_pair(9999., -1);
  
  //check particles' length is at least as long as maximum RR
  // if(part.Length<RR) 
  //   return std::make_pair(9999., -1);
  //Ignore first and last point
  for( UInt_t i = 1; i < part.Hits[plane].size()-1; ++i ){
    //Skip large pulse heights
    if( part.Hits[plane][i].dEdx > 1000. || part.Hits[plane][i].dEdx==-999)
      continue;

    //break whenever above RR upper limit, 26 cm is maximum
    if(part.Hits[plane][i].ResidualRange > RR)continue;
    
    int bin = profile->FindBin( part.Hits[plane][i].ResidualRange );

    if( bin >= 1 && bin <= profile->GetNbinsX() ){
      
      double template_dedx = profile->GetBinContent( bin );
      if( template_dedx < 1.e-6 ){
        template_dedx = ( profile->GetBinContent( bin - 1 ) + profile->GetBinContent( bin + 1 ) ) / 2.;        
      }
      
      
      double template_dedx_err = profile->GetBinError( bin );
      if( template_dedx_err < 1.e-6 ){
        template_dedx_err = ( profile->GetBinError( bin - 1 ) + profile->GetBinError( bin + 1 ) ) / 2.;        
      }

      double dedx_res = 0.04231 + 0.0001783 * part.Hits[plane][i].dEdx * part.Hits[plane][i].dEdx;      
      dedx_res *= part.Hits[plane][i].dEdx; 
      
      
      //Chi2 += ( track_dedx - template_dedx )^2  / ( (template_dedx_err)^2 + (dedx_res)^2 )      
      pid_chi2 += ( pow( (part.Hits[plane][i].dEdx - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) ); 
            
      ++npt;      
    }	
  }

  if( npt == 0 )	
    return std::make_pair(9999., -1);
	  	
  return std::make_pair(pid_chi2, npt); 	
}

//*****************************************************************************
bool pdAnaUtils::isBeamLike(AnaParticlePD* part, AnaBeamPD* beam ){
//*****************************************************************************	

  if (!part) return false;
  
  //From Owen Goodwins studies
  Float_t mccuts[7]  ={-3.,  7., -8.,  7., 27.5, 32.5, 0.93};
  Float_t datacuts[7]={ 0., 10., -5., 10., 30.,  35.0, 0.93};

  //cast the beam particle
  AnaParticlePD* beampart = static_cast<AnaParticlePD*>(beam->BeamParticle);
  if (!beampart) return false;
  
  Float_t beampos[3],beamdir[3], dist[3], dcos=0, cuts[7];      

  // different way of obtaining the beam position and angle for DATA and MC
  // Use the true beam particle to discriminate between data and MC
  AnaTrueParticle* trueBeamPart = static_cast<AnaTrueParticle*>(beam->BeamParticle->TrueObject);

  // For MC
  if (trueBeamPart){
    for (int i=0;i<3;i++){
      beampos[i] = trueBeamPart->Position[i]-trueBeamPart->Position[2]*(trueBeamPart->Direction[i]/trueBeamPart->Direction[2]);
      beamdir[i] = trueBeamPart->Direction[i];
    }
    for (int i=0;i<7;i++) cuts[i] = mccuts[i];
  }
  else{
    // For Data
    if(beam->nMomenta != 1 || beam->nTracks != 1)return false;
    for (int i=0;i<3;i++){
      beampos[i] = beampart->PositionEnd[i];
      beamdir[i] = beampart->DirectionEnd[i];
    }
    for (int i=0;i<7;i++) cuts[i] = datacuts[i];
  }
  
  // compute the difference in position and cos(angle) considering that particle could have been reconstructed backwards
  if(part->PositionEnd[2]<part->PositionStart[2] && part->PositionEnd[2]!=-999){
    for (int i=0;i<3;i++){
      dist[i] = part->PositionEnd[i] - beampos[i];
      dcos   += -1*(part->DirectionEnd[i])*beamdir[i];
    }
  }
  else{
    for (int i=0;i<3;i++){
      dist[i] = part->PositionStart[i] - beampos[i];
      dcos   += part->DirectionStart[i]*beamdir[i];
    }
  }

  if(dist[0] < cuts[0] || dist[0] > cuts[1] || dist[1] < cuts[2] || dist[1] > cuts[3] || dist[2] < cuts[4] || dist[2] > cuts[5] || dcos < cuts[6]) return false;
  else return true;
}

//***************************************************************
AnaParticlePD* pdAnaUtils::GetBeamParticle(const AnaEventC& event){
//***************************************************************

  //Get the beam 
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<const AnaEventB*>(&event)->Beam);
  if (!beam->BeamParticle) return NULL;

  // Get the beam particle
  return static_cast<AnaParticlePD*>(beam->BeamParticle);
}

//***************************************************************
AnaParticlePD* pdAnaUtils::GetBeamParticleInTPC(const AnaEventPD* event){
//***************************************************************

  // Get the array of parts from the event
  AnaParticleB** parts = event->Particles;
  int nParts           = event->nParticles;

  AnaParticlePD* result = NULL;
  
  //look over the particles in the event
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    //select the one tagged as beam particle by Pandora
    if(part->isPandora){
      result = part;
      break;
    }
  }
  return result;
}

//***************************************************************
AnaTrueParticlePD* pdAnaUtils::GetTrueBeamParticle(const AnaEventC& event){
//***************************************************************

  // Get the beam particle
  AnaParticlePD* beampart = GetBeamParticle(event);
  if (!beampart) return NULL;

  // Get the true beam particle
  return static_cast<AnaTrueParticlePD*>(beampart->TrueObject);
}

//***************************************************************
Float_t pdAnaUtils::ComputeTrackLengthFromHitPosition(const AnaParticlePD* part){
//***************************************************************

  Float_t length = 0.;

  //check if hit vector is empty
  if(part->Hits[2].empty()){
    //std::cout << "HitPosition vector empty! Returning -1" << std::endl;
    return -1;
  }

  // Initial position
  TVector3 disp(part->Hits[2][0].Position.X(),part->Hits[2][0].Position.Y(),part->Hits[2][0].Position.Z());

  // Add subsequent hits
  for(size_t i = 1; i < part->Hits[2].size(); ++i){
    if (part->Hits[2][i].Position.X() == -999) break;
    TVector3 pos(part->Hits[2][i].Position.X(),part->Hits[2][i].Position.Y(),part->Hits[2][i].Position.Z());
    disp -= pos;
    length += disp.Mag();
    disp = pos;
  }
  
  return length;
}

//***************************************************************
Float_t pdAnaUtils::ComputeTrackLengthFromTrajectoryPoints(AnaParticlePD* part){
//***************************************************************
  
  if(!part)return -999;
  
  Float_t length = 0.;
  
  int ntps = part->TrjPoints.size();
  int i0 = 0;
  double x0 = 0,y0 = 0,z0 = 0,dx = 0,dy = 0,dz = 0;
  for(int itp = 0; itp < ntps; itp++){
    if(part->TrjPoints[itp].IsValid()){
      x0 = part->TrjPoints[itp].Position.X();
      y0 = part->TrjPoints[itp].Position.Y();
      z0 = part->TrjPoints[itp].Position.Z();
      i0 = itp;
      break;
    }
  }

  for(int itp = i0; itp < ntps; itp++){
    if(!part->TrjPoints[itp].IsValid())continue;
    dx = part->TrjPoints[itp].Position.X()-x0;
    dy = part->TrjPoints[itp].Position.Y()-y0;
    dz = part->TrjPoints[itp].Position.Z()-z0;
    
    length += sqrt(dx*dx+dy*dy+dz*dz);
    
    x0 = part->TrjPoints[itp].Position.X();
    y0 = part->TrjPoints[itp].Position.Y();
    z0 = part->TrjPoints[itp].Position.Z();
  }
  
  return length;
}

//***************************************************************
void pdAnaUtils::ComputeParticlePositionAndDirection(AnaParticlePD* part){
//***************************************************************
  
  if(!part)return;

  int ntps = part->TrjPoints.size();
  if(ntps<1)return;

  int ifirst = -1;
  for(int itp = 0; itp < ntps; itp++){
    if(part->TrjPoints[itp].IsValid()){
      ifirst = itp;
      break;
    }
  }
  if(ifirst != -1){
    part->PositionStart[0] = part->TrjPoints[ifirst].Position.X();
    part->PositionStart[1] = part->TrjPoints[ifirst].Position.Y();
    part->PositionStart[2] = part->TrjPoints[ifirst].Position.Z();
    part->DirectionStart[0] = part->TrjPoints[ifirst].Direction.X();
    part->DirectionStart[1] = part->TrjPoints[ifirst].Direction.Y();
    part->DirectionStart[2] = part->TrjPoints[ifirst].Direction.Z();
  }
  
  int ilast  = -1;
  for(int itp = 1; itp < ntps; itp++){
    if(part->TrjPoints[ntps-itp].IsValid()){
      ilast = ntps-itp;
      break;
    }
  }
  if(ilast != -1){
    part->PositionEnd[0] = part->TrjPoints[ilast].Position.X();
    part->PositionEnd[1] = part->TrjPoints[ilast].Position.Y();
    part->PositionEnd[2] = part->TrjPoints[ilast].Position.Z();
    part->DirectionEnd[0] = part->TrjPoints[ilast].Direction.X();
    part->DirectionEnd[1] = part->TrjPoints[ilast].Direction.Y();
    part->DirectionEnd[2] = part->TrjPoints[ilast].Direction.Z();
  } 
}

//***************************************************************
Float_t pdAnaUtils::ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<double> dEdx){
//***************************************************************

  //check levels are ok
  truncate_high = 1 - truncate_high;
  if((truncate_low < 0 || truncate_low > 1) || (truncate_high < 0 || truncate_high > 1) || truncate_low > truncate_high){
    //std::cout << "invalid truncation levels, returning -999" << std::endl;
    return -999;
  }

  //check vector is not empty
  if(dEdx.empty()){
    //std::cout << "empty dEdx vector, returning -999" << std::endl;
    return -999;
  }

  //compute limits
  int size   = dEdx.size();
  int i_low  = rint(truncate_low*size);
  int i_high = rint(truncate_high*size);

  //compute mean value
  Float_t accumulated = 0;
  int counter = 0;

  for(int i = i_low; i < i_high; i++){
    accumulated = accumulated + dEdx[i];
    counter ++;
  }
  
  return accumulated/counter;
}

//***************************************************************
Float_t pdAnaUtils::ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<AnaHitPD> hits){
//***************************************************************

  //check levels are ok
  truncate_high = 1 - truncate_high;
  if((truncate_low < 0 || truncate_low > 1) || (truncate_high < 0 || truncate_high > 1) || truncate_low > truncate_high){
    //std::cout << "invalid truncation levels, returning -999" << std::endl;
    return -999;
  }

  //check vector is not empty
  if(hits.empty()){
    //std::cout << "empty dEdx vector, returning -999" << std::endl;
    return -999;
  }

  //compute limits
  int size   = hits.size();
  int i_low  = rint(truncate_low*size);
  int i_high = rint(truncate_high*size);

  //compute mean value
  Float_t accumulated = 0;
  int counter = 0;

  for(int i = i_low; i < i_high; i++){
    accumulated = accumulated + hits.at(i).dEdx;
    counter ++;
  }
  
  return accumulated/counter;
}

//***************************************************************
Float_t pdAnaUtils::ComputeDistanceMotherDaughter(AnaParticlePD* mother, AnaParticlePD* daughter){
//***************************************************************

  if(!mother || !daughter){
    std::cout << "ComputeDistanceMotherDaughter: one of the particles does not exist" << std::endl;
    std::cout << "Returning -999" << std::endl;
    return -999.;
  }
  
  double dis = 0;
  for(int i = 0; i < 3; i++)dis = dis + pow(mother->PositionEnd[i] - daughter->PositionStart[i],2);
  dis = sqrt(dis);
  return dis;
}
  
//***************************************************************
Float_t pdAnaUtils::ComputeCosMotherDaughter(AnaParticlePD* mother, AnaParticlePD* daughter){
//***************************************************************

  if(!mother || !daughter){
    std::cout << "ComputeCosMotherDaughter: one of the particles does not exist" << std::endl;
    std::cout << "Returning -999" << std::endl;
    return -999.;
  }
  
  double cos = 0;
  for(int i = 0; i < 3; i++)cos = cos + mother->DirectionEnd[i] * daughter->DirectionStart[i];
  return cos;
}

//***************************************************************
Float_t pdAnaUtils::ComputeAveragedEdxOverResRange(AnaParticlePD* part, double maxresrange){
//***************************************************************

  if(!part){
    std::cout << "ComputeAveragedEdxOverResRange: particle does not exist" << std::endl;
    std::cout << "Returning -999" << std::endl;
    return -999.;
  }
  
  if(part->Hits[2].empty()){
    //std::cout << "ComputeAveragedEdxOverResRange: has no hits" << std::endl;
    //std::cout << "Returning -999" << std::endl;
    return -999.;
  }

  double sumdedx = 0;
  int nhits      = 0;
  for(size_t i = 0; i < part->Hits[2].size(); i++){
    if(part->Hits[2][i].ResidualRange < maxresrange){
      if(part->Hits[2][i].dEdx != -999 && part->Hits[2][i].dEdx < 1000){
        sumdedx += part->Hits[2][i].dEdx;
        nhits++;
      }
    }
  }

  return sumdedx/nhits;
}

//***************************************************************
bool pdAnaUtils::IsStoppingInFV(AnaParticlePD *part){
//***************************************************************

  if(!part)return false;

  bool ItIs = true;
  
  if((TMath::Abs(part->PositionStart[0])>350 || 
      part->PositionStart[1]<50 || part->PositionStart[1]>550 || 
      part->PositionStart[2]<50 || part->PositionStart[2]>645)
      &&
      (TMath::Abs(part->PositionEnd[0])>350 || 
       part->PositionEnd[1]<50 || part->PositionEnd[1]>550 || 
       part->PositionEnd[2]<50 || part->PositionEnd[2]>645))
    ItIs = false;

  return ItIs;
}

//***************************************************************
int pdAnaUtils::GetHitTPCid(AnaHitPD& hit){
//***************************************************************

  return GetPosTPCid(hit.Position);
}

//***************************************************************
int pdAnaUtils::GetPosTPCid(TVector3 pos){
//***************************************************************

  int TPCid = -1;
  
  if(pos.X() < 0){
    if(pos.Z() > 0 && pos.Z() < 230)       TPCid = 1;
    else if(pos.Z() > 230 && pos.Z() < 460)TPCid = 5;
    else if(pos.Z() > 460 && pos.Z() < 690)TPCid = 9;
  }
  else{
    if(pos.Z() > 0 && pos.Z() < 230)       TPCid = 2;
    else if(pos.Z() > 230 && pos.Z() < 460)TPCid = 6;
    else if(pos.Z() > 460 && pos.Z() < 690)TPCid = 10;
  }
  
  return TPCid;
}

//***************************************************************
void pdAnaUtils::EstimateHitsDirection(AnaParticlePD* part){
//***************************************************************
  
  if(!part)return;
  
  //loop over hits
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
    if(ihit < (int)part->Hits[2].size()-1){
      TVector3 dir = part->Hits[2][ihit+1].Position-part->Hits[2][ihit].Position;
      if(dir.Mag()!=0)dir.SetMag(1);
      part->Hits[2][ihit].Direction.SetXYZ(dir.X(),dir.Y(),dir.Z());
      part->Hits[2][ihit].Direction_NoSCE.SetXYZ(dir.X(),dir.Y(),dir.Z());
    }
    else{
      part->Hits[2][ihit].Direction.SetXYZ(part->Hits[2][ihit-1].Direction.X(),part->Hits[2][ihit-1].Direction.Y(),part->Hits[2][ihit-1].Direction.Z());
      part->Hits[2][ihit].Direction_NoSCE.SetXYZ(part->Hits[2][ihit-1].Direction.X(),part->Hits[2][ihit-1].Direction.Y(),part->Hits[2][ihit-1].Direction.Z());
    }
  } 
}

//***************************************************************
void pdAnaUtils::ComputeResidualRange(AnaParticlePD* part){
//***************************************************************
  
  if(!part)return;

  std::vector<double> delta; delta.clear();
  //loop over hits
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
    TVector3 diff = part->Hits[2][ihit+1].Position-part->Hits[2][ihit].Position;
    delta.push_back(diff.Mag());
  }

  //compute residual range
  std::vector<double> new_rr; new_rr.clear();
  new_rr.push_back(delta[0]/2);
  for(int i = 1; i < (int)delta.size(); i++)
    new_rr.push_back(new_rr[i-1]+delta[i-1]);

  //associate new rr to each hit
  for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++)
    part->Hits[2][ihit].ResidualRange = new_rr[ihit];
}

//***************************************************************
Double_t pdAnaUtils::ComputeDepositedEnergy(AnaParticlePD* part){
//***************************************************************
  
  Double_t E = -999;
  if(!part)return E;
  
  int nhits = part->Hits[2].size();
  if(nhits <= 0)return E;

  E = 0;
  for(int ihit = 1; ihit < nhits-1; ihit++){
    if(part->Hits[2][ihit].dEdx > 1000. || part->Hits[2][ihit].dEdx==-999 || part->Hits[2][ihit].Pitch < 0)
      continue;
    E += part->Hits[2][ihit].dEdx * part->Hits[2][ihit].Pitch;
  }
  
  return E;
}

//***************************************************************
Double_t pdAnaUtils::EstimateTrueMomAtAPABorder(AnaParticlePD* part){
//***************************************************************

  Double_t momf = -999;
  if(!part)return momf;

  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
  if(!truePart)return momf;
  
  if(part->PositionStart[2]>220)return momf;

  int nhits = part->Hits[2].size();
  if(nhits <= 0)return momf;

  double depE = 0;
  for(int ihit = 0; ihit < nhits; ihit++){
    if(part->Hits[2][ihit].Position.Z()>220)break;
    depE += part->Hits[2][ihit].dEdx * part->Hits[2][ihit].Pitch;
  }

  double mass;
  if(abs(truePart->PDG) == 211)mass       = 139.57;
  else if(abs(truePart->PDG) == 321)mass  = 493.7;
  else if(abs(truePart->PDG) == 2212)mass = 938.27;
  else return momf;

  double momi = truePart->Momentum*1000;
  momf = sqrt(pow(sqrt(momi*momi+mass*mass)-depE,2)-mass*mass);

  return momf;
}

//***************************************************************
Double_t pdAnaUtils::ComputeDistanceToClosestParticle(AnaParticlePD* part, AnaParticleB** parts, const int nparts){
//***************************************************************

  double distance = 9999;
  for(int ipart = 0; ipart < nparts; ipart++){
    AnaParticlePD* other = static_cast<AnaParticlePD*>(parts[ipart]);
    if(part->UniqueID == other->UniqueID)continue;
    double dis1 = 0;
    double dis2 = 0;
    for(int idis = 0; idis < 3; idis++){
      dis1 += pow(part->PositionEnd[idis]-other->PositionStart[idis],2);
      dis2 += pow(part->PositionEnd[idis]-other->PositionEnd[idis],2);
    }
    dis1 = sqrt(dis1);
    dis2 = sqrt(dis2);
    if(dis1 < distance)distance = dis1;
    if(dis2 < distance)distance = dis2;
  }

  return distance;
}

//***************************************************************
bool pdAnaUtils::GetBeamQuality(AnaEventPD* event){
//***************************************************************

  //get cuts depending on MC/data and nominal momentum
  double mean_x, mean_y, mean_z, sigma_x, sigma_y, sigma_z, min_cos;
  pdAnaUtils::GetBeamQualityCuts(event,
				 mean_x,mean_y,mean_z,
				 sigma_x,sigma_y,sigma_z,
				 min_cos);

  // Main track must exist
  AnaParticlePD* pandora = GetBeamParticleInTPC(event);
  if(!pandora)return false;
  
  // get beam particle
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(event->Beam);
  AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
  if(!beamPart)return false;

  //computations
  double normalized_x = 0;
  double normalized_y = 0;
  double normalized_z = 0;
  double cos = 0;

  normalized_x = abs(pandora->PositionStart[0]-mean_x)/sigma_x;
  normalized_y = abs(pandora->PositionStart[1]-mean_y)/sigma_y;
  normalized_z = abs(pandora->PositionStart[2]-mean_z)/sigma_z;
  for(int i = 0; i < 3; i++)cos = cos + pandora->DirectionStart[i]*beamPart->DirectionEnd[i];

  //check cut
  if(normalized_x < 3 && normalized_y < 3 && normalized_z < 3 && cos > min_cos)return true;
  else return false;
}


//***************************************************************
void pdAnaUtils::GetBeamQualityCuts(AnaEventPD* event, 
				    double &mean_x, double &mean_y, double &mean_z,
				    double &sigma_x, double &sigma_y, double &sigma_z,
				    double &cos){
//***************************************************************

  //get nominal beam momentum. If none, set it to 1.
  AnaEventInfoPD* EventInfo = static_cast<AnaEventInfoPD*>(event->EventInfo);
  int NomBeamMom = (int)EventInfo->NominalBeamMom;
  if(NomBeamMom < 0 || NomBeamMom > 3)NomBeamMom = 1; //we still have no values for 6 and 7 GeV

  //get BQC parameters depending on beam mom and MC/data
  std::stringstream ssmom;
  ssmom << NomBeamMom;
  if(event->GetIsMC()){
    std::string parameter = "pdUtils.AnalysisUtils.BeamQualityCuts.MC."+ssmom.str()+".";
    mean_x = ND::params().GetParameterD((parameter+"meanx").c_str());
    mean_y = ND::params().GetParameterD((parameter+"meany").c_str());
    mean_z = ND::params().GetParameterD((parameter+"meanz").c_str());
    sigma_x = ND::params().GetParameterD((parameter+"sigmax").c_str());
    sigma_y = ND::params().GetParameterD((parameter+"sigmay").c_str());
    sigma_z = ND::params().GetParameterD((parameter+"sigmaz").c_str());
    cos = ND::params().GetParameterD((parameter+"cos").c_str());
  }
  else{
    std::string parameter = "pdUtils.AnalysisUtils.BeamQualityCuts.Data."+ssmom.str()+".";
    mean_x = ND::params().GetParameterD((parameter+"meanx").c_str());
    mean_y = ND::params().GetParameterD((parameter+"meany").c_str());
    mean_z = ND::params().GetParameterD((parameter+"meanz").c_str());
    sigma_x = ND::params().GetParameterD((parameter+"sigmax").c_str());
    sigma_y = ND::params().GetParameterD((parameter+"sigmay").c_str());
    sigma_z = ND::params().GetParameterD((parameter+"sigmaz").c_str());
    cos = ND::params().GetParameterD((parameter+"cos").c_str());
  }
}

//***************************************************************
double pdAnaUtils::GetDensityCorrection(double beta, double gamma){
//***************************************************************

  //Parameters for the density correction
  const double density_C  = 5.2146; 
  const double density_y0 = 0.2; 
  const double density_y1 = 3.0; 
  const double density_a  = 0.19559;
  const double density_k  = 3.0; 
  
  //Estimate the density correction
  double density_y = TMath::Log10(beta * gamma); 
  double ln10 = TMath::Log(10);
  double this_delta = 0.;
  if(density_y > density_y1){
    this_delta = 2.0 * ln10 * density_y - density_C; 
  }
  else if (density_y < density_y0){
    this_delta = 0.; 
  }
  else{
    this_delta = 2.0 * ln10 * density_y - density_C + density_a * pow(density_y1 - density_y, density_k);
  }
  
  return this_delta; 
} 

//***************************************************************
double pdAnaUtils::GetdEdxBetheBloch(double KE, double mass){ 
//***************************************************************

  //Bethe-Bloch parameters, https://indico.fnal.gov/event/14933/contributions/28526/attachments/17961/22583/Final_SIST_Paper.pdf 
  const double rho = 1.39; // [g/cm3], density of LAr
  const double K   = 0.307075; // [MeV cm2 / mol]
  const double A   = 39.948; // [g / mol], atomic mass of Ar 
  const double I   = 188.0e-6; // [MeV], mean excitation energy
  const double Me  = 0.511; // [Mev], mass of electron

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));
  double delta = GetDensityCorrection(beta, gamma);
  
  // == dE/dx with the density correction
  double f = rho * K * (18.0 / A) * pow(1. / beta, 2); 
  double a0 = 0.5 * TMath::Log(2.0 * Me * pow(beta * gamma, 2) * Wmax / (I * I));
  double this_dEdx = f * ( a0 - pow(beta, 2) - delta / 2.0); // [MeV/cm] 
  
  return this_dEdx;
}

//***************************************************************
double pdAnaUtils::GetWmax(double KE, double mass){ 
//***************************************************************

  const double Me  = 0.511; // [Mev], mass of electron
  
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));
  
  return Wmax; 
} 

//***************************************************************
double pdAnaUtils::GetLandauXi(double KE, double dx, double mass){
//***************************************************************
 
  const double rho = 1.39; // [g/cm3], density of LAr
  const double K   = 0.307075; // [MeV cm2 / mol]
  const double A   = 39.948; // [g / mol], atomic mass of Ar 
  
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double xi = rho * dx * 0.5 * K * (18.0 / A) * pow(1. / beta, 2); 
  return xi; 
}    

//***************************************************************
double pdAnaUtils::dEdxPDF(double *x, double *par){ 
//***************************************************************  

  ROOT::Math::VavilovAccurate vav;

  double a = par[2] / par[4];
  double b = (0.422784 + par[1] + log(par[0])) * par[2] / par[4] + par[3]; 
  double y = (x[0] - b) / a; 
  
  double this_vav = 0.;
  
  if(par[0] < 0.01){ // == Landau
    this_vav = TMath::Landau(y); 
    this_vav =this_vav / a;
  }
  else if(par[0] > 10.){ // == Gaussian
    double mu = vav.Mean(par[0], par[1]);
    double sigma = sqrt(vav.Variance(par[0], par[1])); 
    this_vav =TMath::Gaus(y, mu, sigma); 
  }
  else{ // == Vavilov
    this_vav = vav.Pdf(y, par[0], par[1]);
    this_vav = this_vav / a;
  }
  
  return this_vav;
}

//***************************************************************
Float_t pdAnaUtils::dEdxLikelihood(TGraph* tg, TGraph* tg_ke, 
				   Float_t mass){
//***************************************************************

  double width = 0.65;
  TF1* pdf = new TF1("pdf", dEdxPDF, -10., 20., 5);
  double likelihood = 0;
  for(int i = 0; i < tg->GetN(); i++){
    double range = tg->GetPointX(i);//+L;
    double dEdx = tg->GetPointY(i);
    double ke = tg_ke->Eval(range);
    double gamma = (ke/mass)+1.0;
    double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
    double xi = GetLandauXi(ke, width, mass);
    double Wmax = GetWmax(ke, mass);
    double kappa = xi / Wmax;
    double dEdx_BB = GetdEdxBetheBloch(ke, mass);
    double par[5] = {kappa, beta * beta, xi, dEdx_BB, width};
    pdf->SetParameters(par);
    if(pdf->Eval(dEdx) == 0)continue;
    else
      likelihood += log(pdf->Eval(dEdx));
    std::cout << "MASS: " << mass << " RR: " << range << " DEDX: " << dEdx << " KE: " << ke << " KAPPA: " << kappa << " WMAX: " << Wmax << " XI: " << xi << " LKL: " << likelihood << std::endl; 
  }
  delete pdf;
  return likelihood;
}

//***************************************************************
Float_t pdAnaUtils::GetdEdxLikelihood(AnaParticlePD* part, Int_t PDG){
//***************************************************************

  //basic checks
  if(part->Hits[2].empty())return -999.;
  if(PDG != 13 && PDG != 211 && PDG != 321 && PDG != 2212)return -999.;

  //get necessary information
  std::string ssparticle;
  Float_t mass;
  if(PDG == 13){
    ssparticle = "muon";
    mass = 105.66;
  }
  else if(PDG == 211){
    ssparticle = "pion";
    mass = 139.57;
  }
  else if(PDG == 321){
    ssparticle = "kaon";
    mass = 493.677;
  }
  else{
    ssparticle = "proton";
    mass = 938.272;
  }
  TFile* file_ke = TFile::Open((std::string(getenv("PDUTILSROOT"))+"/data/ke_vs_range.root").c_str(),"OPEN");
  TGraph* tg_ke  = (TGraph*)file_ke->Get(ssparticle.c_str());

  //get dedx vs rr graph for particle
  std::vector<double> dedx,rr;
  dedx.clear();
  rr.clear();
  for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){ //ignore first and last hit
    dedx.push_back(part->Hits[2][ihit].dEdx);
    rr.push_back(part->Hits[2][ihit].ResidualRange);
  }
  TGraph* tg = new TGraph(dedx.size(),&rr[0],&dedx[0]);
//   double dedx[189] = {2.59096431732178,3.3262083530426,2.99936604499817,2.69621539115906,2.45239639282227,2.65889501571655,2.8721923828125,2.96113324165344,2.5631046295166,2.13499164581299,2.84897756576538,2.28435754776001,4.85947513580322,2.57258558273315,2.53099274635315,2.31136155128479,2.19671654701233,2.39540076255798,2.44027614593506,2.40085601806641,2.33159995079041,2.39275121688843,2.94427633285522,2.28337812423706,2.57541728019714,2.33782124519348,3.14544320106506,2.46237850189209,2.18722152709961,2.31510615348816,2.5407919883728,2.4052369594574,2.60192775726318,2.1620500087738,1.97512233257294,2.00497245788574,3.67113018035889,2.78016662597656,2.16116213798523,3.2123761177063,2.15956282615662,2.53732967376709,2.12903571128845,2.31120705604553,3.32276225090027,2.54401159286499,2.45249199867249,2.63295531272888,2.40191197395325,2.15561771392822,2.46324372291565,2.19734621047974,6.28847408294678,2.84497785568237,4.08778953552246,2.87325549125671,3.86758995056152,2.26380562782288,2.15557932853699,3.91243505477905,2.46982455253601,1.93973767757416,1.78229963779449,2.80157232284546,2.56042718887329,2.23831176757812,2.21581387519836,2.08650326728821,2.61188101768494,2.26256489753723,2.48247075080872,2.20683670043945,2.16188549995422,2.31100535392761,2.11773681640625,1.74697887897491,1.89291179180145,3.67076539993286,2.72570419311523,2.11984324455261,2.94880390167236,2.44372057914734,2.29839706420898,2.6672351360321,2.67290925979614,2.7120316028595,2.70571875572205,2.80990552902222,2.4881649017334,2.24165368080139,2.17152237892151,2.0438506603241,2.6602041721344,2.1686863899231,2.14035820960999,2.39983797073364,2.32761549949646,2.25096344947815,2.2850079536438,1.85133767127991,1.90941596031189,2.02054238319397,2.68489789962769,2.25075006484985,2.2542781829834,2.00240349769592,2.31677460670471,2.22072243690491,2.3314356803894,2.25312113761902,2.0907871723175,2.25341844558716,2.48657608032227,1.97184157371521,2.06611680984497,2.08378982543945,2.42258906364441,4.92590999603272,3.39344120025635,8.54634475708008,2.18713617324829,2.55733180046082,2.22490668296814,2.10700130462646,3.42171597480774,2.06779813766479,2.91834139823914,2.5253643989563,3.14987707138062,3.16592216491699,2.8147292137146,2.87871074676514,1.98284554481506,2.51280331611633,1.95752394199371,2.10781645774841,2.20413517951965,2.26310658454895,2.11322045326233,2.10560870170593,2.10125184059143,2.04906320571899,1.87636649608612,1.9877393245697,2.12836146354675,2.33427929878235,2.23855590820312,2.06583428382874,1.99633371829987,1.81923615932465,2.3602192401886,2.4761016368866,2.90728330612183,2.19378471374512,2.47629284858704,1.87241816520691,1.91761481761932,1.92383933067322,2.03943872451782,2.158371925354,2.13005781173706,2.14477610588074,3.60203289985657,1.97866332530975,1.96905040740967,2.13348412513733,5.16601037979126,1.9453991651535,2.23490238189697,3.21706128120422,2.13215160369873,1.95216715335846,2.44398474693298,1.82735991477966,1.81957876682282,1.95420336723328,1.89899015426636,2.00407576560974,2.29312038421631,1.92866826057434,2.10341835021973,2.09143447875977,2.774493932724,3.07616472244263,1.80989301204681,2.04485034942627,2.58910512924194,1.81456351280212,2.05662274360657
// };
//   double rr[189] = {0.222535088658333,0.66760528087616,1.1058691740036,1.54694879055023,1.99009585380554,2.42369365692139,2.86316871643066,3.30938839912415,3.74931573867798,4.18554067611694,4.63362073898315,5.05915403366089,5.52104043960571,5.9430456161499,6.38564872741699,6.82544183731079,7.25781774520874,7.69618654251099,8.13615703582764,8.57535171508789,9.01869964599609,9.45083999633789,9.89880657196045,10.3279542922974,10.7715797424316,11.2091064453125,11.6591958999634,12.086986541748,12.5238018035889,12.9615325927734,13.4078454971313,13.8413610458374,14.2800521850586,14.7141857147217,15.1500186920166,15.584002494812,16.0340061187744,16.4695453643799,16.9027881622314,17.3575820922852,17.7901878356934,18.2339515686035,18.6728515625,19.1138095855713,19.5603103637695,19.9906692504883,20.4264736175537,20.86549949646,21.3052558898926,21.7423305511475,22.1878128051758,22.6257591247559,23.080753326416,23.5092277526855,23.9724578857422,24.4023838043213,24.853235244751,25.2788200378418,25.7238521575928,26.1729393005371,26.6052684783936,27.0447673797607,27.4812450408936,27.9357490539551,28.3746566772461,28.8232860565186,29.2710113525391,29.711311340332,30.1601734161377,30.5947074890137,31.0491218566895,31.4856472015381,31.9258213043213,32.3680572509766,32.8065490722656,33.2526550292969,33.6904220581055,34.1508293151855,34.5819396972656,35.0196495056152,35.4701232910156,35.9138832092285,36.3561859130859,36.8015632629395,37.2419548034668,37.6784858703613,38.1164779663086,38.5722732543945,39.0155029296875,39.4599800109863,39.9048156738281,40.352668762207,40.8044776916504,41.2427253723145,41.6892585754395,42.1353263854981,42.5770606994629,43.020435333252,43.4637145996094,43.9044456481934,44.3581008911133,44.805835723877,45.2565460205078,45.6969451904297,46.149227142334,46.5925178527832,47.0412406921387,47.4877738952637,47.9330062866211,48.3752365112305,48.8154602050781,49.2630462646484,49.7139129638672,50.1576042175293,50.6105499267578,51.0555419921875,51.5015754699707,51.9541740417481,52.3883171081543,52.8491706848145,53.2814712524414,53.7383422851562,54.1857452392578,54.6323585510254,55.0837707519531,55.5232048034668,55.983325958252,56.430778503418,56.8926887512207,57.3442077636719,57.7977561950684,58.2565460205078,58.7037239074707,59.1593437194824,59.6045227050781,60.0606460571289,60.5189323425293,60.9719696044922,61.4251022338867,61.8861846923828,62.3372268676758,62.7908744812012,63.2462768554688,63.7039794921875,64.1611175537109,64.6171569824219,65.0737075805664,65.5265960693359,65.9785690307617,66.4308776855469,66.8908538818359,67.3483963012695,67.8054809570313,68.2578506469727,68.7184143066406,69.1698303222656,69.6272277832031,70.0845260620117,70.5432815551758,71.0007858276367,71.4587326049805,71.9147567749023,72.3719177246094,72.8220977783203,73.2834014892578,73.7457733154297,74.2179489135742,74.6691360473633,75.1364822387695,75.5971832275391,76.0552368164062,76.5172500610352,76.9868927001953,77.4440689086914,77.9096374511719,78.3684997558594,78.8243637084961,79.2866058349609,79.7575988769531,80.2283477783203,80.6964569091797,81.1637725830078,81.6315002441406,82.0994262695313,82.5601654052734,83.0321960449219,83.503044128418,83.9679946899414,84.4388580322266
// };
//   TGraph* tg = new TGraph(189,rr,dedx);


  std::cout << "MIGUEEEEEEEEEEEEEEEEEEEE" << std::endl;
  Float_t result = dEdxLikelihood(tg,tg_ke,mass);
  std::cout << "MIGUEEEEEEEEEEEEEEEEEEEE END" << std::endl;

  delete tg;
  file_ke->Close("R");

  return result;
}

//***************************************************************
Float_t pdAnaUtils::GetdEdxLikelihood_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR){
//***************************************************************

  //basic checks
  if(part->Hits[2].empty())return -999.;
  if(PDG != 13 && PDG != 211 && PDG != 321 && PDG != 2212)return -999.;

  //get necessary information
  std::string ssparticle;
  Float_t mass;
  if(PDG == 13){
    ssparticle = "muon";
    mass = 105.66;
  }
  else if(PDG == 211){
    ssparticle = "pion";
    mass = 139.57;
  }
  else if(PDG == 321){
    ssparticle = "kaon";
    mass = 493.677;
  }
  else{
    ssparticle = "proton";
    mass = 938.272;
  }
  TFile* file_ke = TFile::Open((std::string(getenv("PDUTILSROOT"))+"/data/ke_vs_range.root").c_str(),"OPEN");
  TGraph* tg_ke  = (TGraph*)file_ke->Get(ssparticle.c_str());

  //get dedx vs rr graph for particle
  std::vector<double> dedx,rr;
  dedx.clear();
  rr.clear();
  for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){ //ignore first and last hit
    if(part->Hits[2][ihit].ResidualRange > maxRR)continue;
    dedx.push_back(part->Hits[2][ihit].dEdx);
    rr.push_back(part->Hits[2][ihit].ResidualRange);
  }
  TGraph* tg = new TGraph(dedx.size(),&rr[0],&dedx[0]);

  Float_t result = dEdxLikelihood(tg,tg_ke,mass);

  delete tg;
  file_ke->Close("R");

  return result;
}

//***************************************************************
std::pair<Float_t,Float_t> pdAnaUtils::dEdxLikelihoodFreeRange(TGraph* tg, TGraph* tg_ke, 
					    Float_t mass){
//***************************************************************

  double width = 0.65;
  TF1* pdf = new TF1("pdf", dEdxPDF, -10., 20., 5);
  double L  = 0;
  double Lf = 10;
  double step = 0.1;
  std::vector<double> L_v,Likelihood_v;
  L_v.clear();
  Likelihood_v.clear();
  while(L<Lf){
    double likelihood = 0;
    for(int i = 0; i < tg->GetN(); i++){
      double range = tg->GetPointX(i)+L;
      double dEdx = tg->GetPointY(i);
      double ke = tg_ke->Eval(range);
      double gamma = (ke/mass)+1.0;
      double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
      double xi = GetLandauXi(ke, width, mass);
      double Wmax = GetWmax(ke, mass);
      double kappa = xi / Wmax;
      double dEdx_BB = GetdEdxBetheBloch(ke, mass);
      double par[5] = {kappa, beta * beta, xi, dEdx_BB, width};
      pdf->SetParameters(par);
      if(pdf->Eval(dEdx) == 0)continue;
      else
	likelihood += log(pdf->Eval(dEdx));
    }
    L_v.push_back(L);
    Likelihood_v.push_back(likelihood);
    L += step;
  }
  auto it = std::max_element(Likelihood_v.begin(), Likelihood_v.end());
  int index = std::distance(Likelihood_v.begin(), it);
  delete pdf;
  return std::make_pair(Likelihood_v[index],L_v[index]);
}

//***************************************************************
std::pair<Float_t,Float_t> pdAnaUtils::GetdEdxLikelihoodFreeRange(AnaParticlePD* part, Int_t PDG){
//***************************************************************

  //basic checks
  if(part->Hits[2].empty())return std::make_pair(-999.,-999.);
  if(PDG != 13 && PDG != 211 && PDG != 321 && PDG != 2212)return std::make_pair(-999.,-999.);

  //get necessary information
  std::string ssparticle;
  Float_t mass;
  if(PDG == 13){
    ssparticle = "muon";
    mass = 105.66;
  }
  else if(PDG == 211){
    ssparticle = "pion";
    mass = 139.57;
  }
  else if(PDG == 321){
    ssparticle = "kaon";
    mass = 493.677;
  }
  else{
    ssparticle = "proton";
    mass = 938.272;
  }
  TFile* file_ke = TFile::Open((std::string(getenv("PDUTILSROOT"))+"/data/ke_vs_range.root").c_str(),"OPEN");
  TGraph* tg_ke  = (TGraph*)file_ke->Get(ssparticle.c_str());

  //get dedx vs rr graph for particle
  std::vector<double> dedx,rr;
  dedx.clear();
  rr.clear();
  for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){ //ignore first and last hit
    if(part->Hits[2][ihit].dEdx<1 || part->Hits[2][ihit].dEdx>100)continue;
    dedx.push_back(part->Hits[2][ihit].dEdx);
    rr.push_back(part->Hits[2][ihit].ResidualRange);
  }
  TGraph* tg = new TGraph(dedx.size(),&rr[0],&dedx[0]);

  std::pair<Float_t,Float_t>result = dEdxLikelihoodFreeRange(tg,tg_ke,mass);
  delete tg;
  file_ke->Close("R");
  
  return result;
}

//***************************************************************
std::pair<Float_t,Float_t> pdAnaUtils::GetdEdxLikelihoodFreeRange_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR){
//***************************************************************

  //basic checks
  if(part->Hits[2].empty())return std::make_pair(-999.,-999.);
  if(PDG != 13 && PDG != 211 && PDG != 321 && PDG != 2212)return std::make_pair(-999.,-999.);

  //get necessary information
  std::string ssparticle;
  Float_t mass;
  if(PDG == 13){
    ssparticle = "muon";
    mass = 105.66;
  }
  else if(PDG == 211){
    ssparticle = "pion";
    mass = 139.57;
  }
  else if(PDG == 321){
    ssparticle = "kaon";
    mass = 493.677;
  }
  else{
    ssparticle = "proton";
    mass = 938.272;
  }
  TFile* file_ke = TFile::Open((std::string(getenv("PDUTILSROOT"))+"/data/ke_vs_range.root").c_str(),"OPEN");
  TGraph* tg_ke  = (TGraph*)file_ke->Get(ssparticle.c_str());

  //get dedx vs rr graph for particle
  std::vector<double> dedx,rr;
  dedx.clear();
  rr.clear();
  for(int ihit = 1; ihit < (int)part->Hits[2].size()-1; ihit++){ //ignore first and last hit
    if(part->Hits[2][ihit].ResidualRange > maxRR)continue;
    dedx.push_back(part->Hits[2][ihit].dEdx);
    rr.push_back(part->Hits[2][ihit].ResidualRange);
  }
  TGraph* tg = new TGraph(dedx.size(),&rr[0],&dedx[0]);

  std::pair<Float_t,Float_t>result = dEdxLikelihoodFreeRange(tg,tg_ke,mass);
  delete tg;
  file_ke->Close("R");
  
  return result;
}
