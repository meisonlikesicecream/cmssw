#ifndef _PCAFITTER_H_
#define _PCAFITTER_H_

#include "TrackFitter.h"
#include <math.h>

#include <iomanip>
#include <map>
#include <set>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>


/**
   \brief PCA fitter
**/
class PCATrackFitter:public TrackFitter
{

 private:
  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const
  {
    ar << boost::serialization::base_object<TrackFitter>(*this);
    ar << sec_phi;
  }
  
  template<class Archive> void load(Archive & ar, const unsigned int version)
  {
    ar >> boost::serialization::base_object<TrackFitter>(*this);
    ar >> sec_phi;
  }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

  Track *track_;
  bool paramscmpt_;
  std::string cfname_;

  int m_charge;

 public:

  /**
     \brief Default constructor   
  **/
  PCATrackFitter();
  /**
     \brief Constructor   
     \param nb Layers number
  **/  
  PCATrackFitter(int nb);
  ~PCATrackFitter();

  void initialize();
  void fit(vector<Hit*> hits);
  void fit();
  void setTrack(Track *intc); 

  void mergePatterns();
  void mergeTracks();

  void set_const_filename(const std::string & in)
  {
    cfname_ = in;
  }

  const std::string & get_const_filename () const
  {
    return cfname_;
  }



  TrackFitter* clone();
};
#endif