#ifndef ROO_SUBTHRESHOLD
#define ROO_SUBTHRESHOLD

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooSubThreshold : public RooAbsPdf {

	public:

		inline RooSubThreshold() {};

		RooSubThreshold(const char *name,
				const char *title,
				RooAbsReal& _x,
				RooAbsReal& _mean,
				RooAbsReal& _width,
				RooAbsReal& _radius,
				RooAbsReal& _mass_a,
				RooAbsReal& _mass_b,
				RooAbsReal& _mass_c,
				RooAbsReal& _mass_parent,
				RooAbsReal& _spin);
		RooSubThreshold(const RooSubThreshold& other,
				const char* name=0) ;
		virtual TObject* clone(const char* newname) const {
			return new RooSubThreshold(*this,newname);
		}
		inline virtual ~RooSubThreshold() { }

	protected:

		RooRealProxy x ;      // observable (mass)
		RooRealProxy mean ;    // mean
		RooRealProxy width ;  // width
		RooRealProxy radius ; // meson radius parameter
		RooRealProxy mass_a ; // mass of first daughter from decay
		RooRealProxy mass_b ; // mass of second daughter from decay
		RooRealProxy mass_c ; // mass of third daughter from decay
		RooRealProxy mass_parent ; // mass of parent particle
		RooRealProxy spin ; // spin (= 0,1,2) //need to upgrade to baryons
    double mean_eff; // will be needed for computation of effective mass
		Double_t evaluate() const ;


	private:
		Double_t getWidth() const;
		Double_t getQterm(Double_t q, Double_t q0) const;
		Double_t getFF(Double_t q) const;
		Double_t getQ(Double_t mass) const;
	  Double_t compute_meff() const;


};

#endif
