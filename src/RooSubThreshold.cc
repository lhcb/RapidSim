/*****************************************************************************
 * Project:
 * Package:
 *    File: $Id: RooSubThreshold.cc,v 1.7 2020/12/10 08:53:43
 * Authors: Y.Amhis
 *
 *****************************************************************************/

// -- CLASS DESCRIPTION [PDF] --
// This is an implementation of Relativistic Breit Wigner for Spin 0,1,2 particles
// with Blatt-Weisskopf form factors and barrier functions fort RooRarFit.
// The meson radius is set to 3.1 GeV^-1 which is similar to the EvtGen number
// but may not be appropriate for every meson. See Phys Rev D 72, 052002 (2005).
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This is an implementation of Relativistic Breit Wigner for Spin 0,1,2 particles
// with Blatt-Weisskopf form factors and barrier functions for RooRarFit.
// The meson radius is set to 3.1 GeV^-1 which is similar to the EvtGen number
// but may not be appropriate for every meson. See Phys Rev D 72, 052002 (2005).
// END_HTML
//



#include "Riostream.h"

#include "RooSubThreshold.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
//#include "RooComplex.h"
#include <complex>

//ClassImp(RooSubThreshold)

//------------------------------------------------------------------
RooSubThreshold::RooSubThreshold(const char *name, const char *title,
		RooAbsReal& _x, RooAbsReal& _mean,
		RooAbsReal& _width, RooAbsReal& _radius,
		RooAbsReal& _mass_a, RooAbsReal& _mass_b, RooAbsReal& _mass_c, RooAbsReal& _mass_parent,
		RooAbsReal& _spin) :
	RooAbsPdf(name,title),
	x("x","Dependent",this,_x),
	mean("mean","Mean",this,_mean),
	width("width","Width",this,_width),
	radius("radius","Radius",this,_radius),  // meson radius e.g. 3.1 GeV^-1 //to be updated for baryons
	mass_a("mass_a","Mass of first daughter",this,_mass_a), // e.g. proton from Lambda*
	mass_b("mass_b","Mass of second daugher",this,_mass_b), // e.g. kaon from Lambda*
	mass_c("mass_c","Mass of third daughter", this, _mass_c),
	mass_parent("mass_parent","Mass of parent particle", this, _mass_parent), //e.g. Lb
	spin("spin","Spin",this,_spin)
{
}

//------------------------------------------------------------------
RooSubThreshold::RooSubThreshold(const RooSubThreshold& other,
		const char* name) :
	RooAbsPdf(other,name),
	x("x",this,other.x),
	mean("mean",this,other.mean),
	width("width",this,other.width),
	radius("radius",this,other.radius),
	mass_a("mass_a",this,other.mass_a),
	mass_b("mass_b",this,other.mass_b),
	mass_c("mass_c",this,other.mass_c),
	mass_parent("mass_parent",this,other.mass_parent),
	spin("spin",this,other.spin)
{
}

//------------------------------------------------------------------
Double_t RooSubThreshold::evaluate() const
{

	double top = getWidth();
	double dm2 = std::pow(x,2) - std::pow(mean,2);
	double bottom = std::pow(dm2,2) + pow(mean*getWidth(),2);
	return top/bottom;
}
//------------------------------------------------------------------
Double_t RooSubThreshold::compute_meff() const
{
double minMass = mass_a+mass_b;
double maxMass = mass_parent - mass_c ; // need to substract the mass of the 3rd particle
double tanhTerm = std::tanh( mean - ((minMass+maxMass)/2)/(maxMass-minMass));
double mean_eff = minMass + 1./2*(maxMass-minMass)*(1.+tanhTerm)	;
return mean_eff;
}
//------------------------------------------------------------------
Double_t RooSubThreshold::getWidth() const
{
	Double_t q  = getQ(x);

	Double_t q0 = getQ(mean_eff); //pass here m_eff for the subthrehold
	Double_t result(0.0);

	if (q>0 && q0>0 && x>0 && mean>0) {
		result = width * getQterm(q,q0) * (mean/x)
			* (getFF(q) / getFF(q0));
	}
	return (result);
}

//------------------------------------------------------------------
Double_t RooSubThreshold::getQterm(Double_t q, Double_t q0) const
{
	Double_t result = pow(q/q0,2*spin+1);
	return (result); // result = (q/q0)^(2*spin+1);
}

//------------------------------------------------------------------
Double_t RooSubThreshold::getFF(Double_t q) const
{

	Double_t z = q * radius;
	Double_t result(1.0);

	if (spin==1) {
		result = 1.0/(1+z*z);
	} else if (spin==2) {
		result = 1.0/(9 + 3*z*z + z*z*z*z);
	}

	return (result); //square of the Blatt-Weisskopf form factor
}

//-----------------------------------------------------
Double_t RooSubThreshold::getQ(Double_t mass) const
{

	if (mass < (mass_b+mass_a)) {return(0);}

	const Double_t mDaugSumSq  = (mass_b+mass_a)*(mass_b+mass_a);
	const Double_t mDaugDiffSq = (mass_b-mass_a)*(mass_b-mass_a);

	Double_t q = sqrt((mass*mass-mDaugSumSq)*(mass*mass-mDaugDiffSq))/(2*mass);
	return(q);
}
