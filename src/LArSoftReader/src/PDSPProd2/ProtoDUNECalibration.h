#ifndef NewProtoDUNECalibration_h
#define NewProtoDUNECalibration_h

#include "ProtoDUNETrackUtils.h"
#include "LArSoftReaderProjectHeaders.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"


namespace protoana{


  class ProtoDUNECalibration{

    public:
      ProtoDUNECalibration();
      std::vector< float > GetCalibratedCalorimetry(  const recob::Track &track, const std::string trackModule, const std::string caloModule );
      double HitToEnergy(recob::Hit* hit, double X, double Y, double Z);

    private:
      float calc_dEdX(double dqdx, double betap, double Rho, double Efield, double Wion, double alpha);

      double tot_Ef( double, double, double );

      unsigned int planeID;
      double betap;
      double Rho;
      //double Efield;
      double Wion;
      double alpha;
      double norm_factor;
      double calib_factor;
      std::string X_correction_name;
      TFile * X_correction_file;

      std::string YZ_correction_name;
      TFile * YZ_correction_file;

      std::string E_field_correction_name;
      TFile * E_field_file;

      TH1F * X_correction_hist;
      TH2F * YZ_neg;
      TH2F * YZ_pos;

      TH3F * ex_neg;
      TH3F * ey_neg;
      TH3F * ez_neg;

      TH3F * ex_pos;
      TH3F * ey_pos;
      TH3F * ez_pos;

      ProtoDUNETrackUtils trackUtil;

  };

}


#endif
