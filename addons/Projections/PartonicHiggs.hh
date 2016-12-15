// -*- C++ -*-
#ifndef RIVET_PartonicHiggs_HH
#define RIVET_PartonicHiggs_HH

#include "Rivet/Projections/ParticleFinder.hh"

namespace Rivet {


  /// @brief Convenience finder of partonic top quarks
  ///
  /// @warning Requires there to be tops in the event record. A fiducial pseudo-top
  /// analysis approach is strongly recommended instead of this.
  class PartonicHiggs : public ParticleFinder {
  public:


    /// @brief Enum for categorising top quark decay modes
    ///
    /// Consider bbbar and all decays for now only 
    enum DecayMode { BBBAR, ALL };


    /// @name Constructors
    //@{

    /// Constructor optionally taking cuts object
    PartonicHiggs(const Cut& c=Cuts::OPEN)
      : ParticleFinder(c), _decaymode(ALL)
    {  }

    /// Constructor taking decay mode details (and an optional cuts object)
    PartonicHiggs(DecayMode decaymode, const Cut& c=Cuts::OPEN)
      : ParticleFinder(c), _decaymode(decaymode)
    {  }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(PartonicHiggs);

    //@}


    /// Access to the found partonic tops
    const Particles& higgs() const { return _theParticles; }


    /// Clear the projection
    void clear() {
      _theParticles.clear();
    }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& event) {
      // Find partonic tops
      _theParticles = filter_select(event.allParticles(_cuts), lastParticleWith(isHiggs));
      // Filtering by decay mode 
      if (_decaymode == BBBAR) {
        const auto fn = [&](const Particle& h) {
	  if (filter_select(h.children(), isBottom).size()==2 ) return true;
	  return false;
	};
	ifilter_select(_theParticles, fn);
      }
    }
    

    /// Compare projections.
    int compare(const Projection& p) const {
      const PartonicHiggs& other = dynamic_cast<const PartonicHiggs&>(p);
      return cmp(_cuts, other._cuts) || cmp(_decaymode, other._decaymode) ;
    }


  private:

    DecayMode _decaymode;

  };


}


#endif
