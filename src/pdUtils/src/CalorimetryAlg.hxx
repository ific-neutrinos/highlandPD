// \file CalorimetryAlg.h
//
// \brief Functions to calculate dE/dx. Based on code in Calorimetry.cxx
//
// \author andrzej.szelc@yale.edu
//
#ifndef UTIL_CALORIMETRYALG_H
#define UTIL_CALORIMETRYALG_H
/*
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
*/
#include <vector>
#include "pdDataClasses.hxx"

namespace detinfo{
  class DetectorClocksData{

  };
  class DetectorPropertiesData{

  };
}

namespace calo {


  Int_t trigger_offset(detinfo::DetectorClocksData const& data);
  double sampling_rate(detinfo::DetectorClocksData const& data);

  class CalorimetryAlg {
  public:
    struct Config {
      /*
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Sequence<double> CalAmpConstants{
        Name("CalAmpConstants"),
        Comment("ADC to electrons constants for each plane.")};

      fhicl::Sequence<double> CalAreaConstants{
        Name("CalAreaConstants"),
        Comment("Area to electrons constants for each plane.")};

      fhicl::Atom<bool> CaloUseModBox{Name("CaloUseModBox"),
                                      Comment("Use modified box model if true, birks otherwise")};

      fhicl::Atom<int> CaloLifeTimeForm{Name("CaloLifeTimeForm"),
                                        Comment("0 = exponential, 1 = exponential + constant")};

      fhicl::Atom<bool> CaloDoLifeTimeCorrection{Name("CaloDoLifeTimeCorrection"),
                                                 Comment("Apply lifetime correction if true")};
      */
    };

    CalorimetryAlg(const std::string& pset):fUseModBox(false), fLifeTimeForm(0),fDoLifeTimeCorrection(false){}
    //      : CalorimetryAlg(fhicl::Table<Config>(pset, {})())
    //    {}

    //    CalorimetryAlg(const Config& config);
        CalorimetryAlg();

    double dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                    detinfo::DetectorPropertiesData const& det_prop,
                    AnaHitPD const& hit,
                    double pitch,
                    double T0 = 0) const;
    double dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                    detinfo::DetectorPropertiesData const& det_prop,
                    double dQ,
                    double time,
                    double pitch,
                    unsigned int plane,
                    double T0 = 0) const;
    double dEdx_AMP(detinfo::DetectorClocksData const& clock_data,
                    detinfo::DetectorPropertiesData const& det_prop,
                    double dQdx,
                    double time,
                    unsigned int plane,
                    double T0 = 0) const;

    // FIXME: How may of these are actually used?
    double dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                     detinfo::DetectorPropertiesData const& det_prop,
                     AnaHitPD const& hit,
                     double pitch,
                     double T0 = 0) const;
    double dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                     detinfo::DetectorPropertiesData const& det_prop,
                     double dQ,
                     double time,
                     double pitch,
                     unsigned int plane,
                     double T0 = 0) const;
    double dEdx_AREA(detinfo::DetectorClocksData const& clock_data,
                     detinfo::DetectorPropertiesData const& det_prop,
                     double dQdx,
                     double time,
                     unsigned int plane,
                     double T0 = 0) const;

    double
    ElectronsFromADCPeak(double adc, unsigned short plane) const
    {
      return adc / fCalAmpConstants[plane];
    }

    double
    ElectronsFromADCArea(double area, unsigned short plane) const
    {
      return area / fCalAreaConstants[plane];
    }

    double LifetimeCorrection(detinfo::DetectorClocksData const& clock_data,
                              detinfo::DetectorPropertiesData const& det_prop,
                              double time,
                              double T0 = 0) const;


    double BirksCorrection 	( 	double  	dQdX	) 	const;
    double ModBoxCorrection 	( 	double  	dQdX	) const;    
    //  private:
    //    art::ServiceHandle<geo::Geometry const> geom;

    double dEdx_from_dQdx_e(detinfo::DetectorClocksData const& clock_data,
                            detinfo::DetectorPropertiesData const& det_prop,
                            double dQdx_e,
                            double time,
                            double T0 = 0) const;

    std::vector<double> fCalAmpConstants;
    std::vector<double> fCalAreaConstants;
    bool const fUseModBox;
    int const fLifeTimeForm;
    bool const fDoLifeTimeCorrection;

  }; // class CalorimetryAlg
} // namespace calo
#endif // UTIL_CALORIMETRYALG_H