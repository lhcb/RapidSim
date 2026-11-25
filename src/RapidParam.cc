#include "RapidParam.h"

#include "TRandom.h"

#include "RapidParticle.h"
#include "RapidParticleData.h"
#include "RapidPID.h"

double RapidParam::eval() {

	if(type_ == RapidParam::THETA || type_ == RapidParam::COSTHETA) return evalTheta();
	if(type_ == RapidParam::MCORR) return evalCorrectedMass();
	if(type_ == RapidParam::ProbNNmu || type_ == RapidParam::ProbNNpi || type_ == RapidParam::ProbNNe ||
	   type_ == RapidParam::ProbNNk || type_ == RapidParam::ProbNNp) return evalPID();
	if(type_ >= RapidParam::VtxX && type_ <= RapidParam::OrigZ) return evalPos();

	double fd      = particles_[0]->getFD(truth_);
	double ip      = particles_[0]->getIP();
	double minip   = particles_[0]->getMinIP();
	double ipSigma = particles_[0]->getSigmaIP();
	double minipSigma = particles_[0]->getSigmaMinIP();
	
	double decaytime = particles_[0]->t();
	int id = particles_[0]->id();
	bool hasOsc = particles_[0]->hasOsc();

	mom_.SetPxPyPzE(0.,0.,0.,0.);
	if(truth_) {
		for(unsigned int i=0; i<particles_.size(); ++i) {
			mom_ += particles_[i]->getP();
		}
	} else {
		ip = particles_[0]->getIPSmeared();
		minip = particles_[0]->getMinIPSmeared();
		for(unsigned int i=0; i<particles_.size(); ++i) {
			mom_ += particles_[i]->getPSmeared();
		}
	}

	switch(type_) {
		case RapidParam::M:
			return mom_.M();
		case RapidParam::M2:
			return mom_.M2();
		case RapidParam::MT:
			return mom_.Mt();
		case RapidParam::E:
			return mom_.E();
		case RapidParam::ET:
			return mom_.Et();
		case RapidParam::P:
			return mom_.P();
		case RapidParam::PX:
			return mom_.Px();
		case RapidParam::PY:
			return mom_.Py();
		case RapidParam::PZ:
			return mom_.Pz();
		case RapidParam::PT:
			return mom_.Pt();
		case RapidParam::ETA:
			return mom_.Eta();
		case RapidParam::PHI:
			return mom_.Phi();
		case RapidParam::RAPIDITY:
			return mom_.Rapidity();
		case RapidParam::GAMMA:
			return mom_.Gamma();
		case RapidParam::BETA:
			return mom_.Beta();
		case RapidParam::IP:
			return ip;
		case RapidParam::SIGMAIP:
			return ipSigma;
		case RapidParam::MINIP:
			return minip;
		case RapidParam::SIGMAMINIP:
			return minipSigma;
		case RapidParam::FD:
			return fd;
		case RapidParam::decaytime:
			return decaytime;
		case RapidParam::ID:
			return id;	
		case RapidParam::hasOsc:
			return hasOsc;						
		case RapidParam::ProbNNmu:
		case RapidParam::ProbNNpi:
		case RapidParam::ProbNNk:
		case RapidParam::ProbNNp:
		case RapidParam::ProbNNe:
		case RapidParam::VtxX:
		case RapidParam::VtxY:
		case RapidParam::VtxZ:
		case RapidParam::OrigX:
		case RapidParam::OrigY:
		case RapidParam::OrigZ:
		case RapidParam::THETA:
		case RapidParam::COSTHETA:
		case RapidParam::MCORR:
			std::cout << "WARNING in RapidParam::eval : parameter type " << type_ << " cannot be used with a single momentum." << std::endl
				  << "                              returning 0." << std::endl;
			break;
		default:
			std::cout << "WARNING in RapidParam::eval : unknown parameter type " << type_ << std::endl
				  << "                              returning 0." << std::endl;

	}

	return 0.;

}

bool RapidParam::canBeSmeared() {
	switch(type_) {
		case RapidParam::M:
			return true;
		case RapidParam::M2:
			return true;
		case RapidParam::MT:
			return true;
		case RapidParam::E:
			return true;
		case RapidParam::ET:
			return true;
		case RapidParam::P:
			return true;
		case RapidParam::PX:
			return true;
		case RapidParam::PY:
			return true;
		case RapidParam::PZ:
			return true;
		case RapidParam::PT:
			return true;
		case RapidParam::ETA:
			return true;
		case RapidParam::PHI:
			return true;
		case RapidParam::RAPIDITY:
			return true;
		case RapidParam::GAMMA:
			return true;
		case RapidParam::BETA:
			return true;
		case RapidParam::IP:
			return true;
		case RapidParam::SIGMAIP:
			return false;
		case RapidParam::MINIP:
			return true;
		case RapidParam::SIGMAMINIP:
			return false;
		case RapidParam::FD:
			return true;
		case RapidParam::ProbNNmu:
			return true;
		case RapidParam::ProbNNpi:
			return true;
		case RapidParam::ProbNNk:
			return true;
		case RapidParam::ProbNNp:
			return true;
		case RapidParam::ProbNNe:
			return true;
		case RapidParam::THETA:
			return true;
		case RapidParam::COSTHETA:
			return true;
		case RapidParam::MCORR:
			return true;
		case RapidParam::VtxX:
			return true;
		case RapidParam::VtxY:
			return true;
		case RapidParam::VtxZ:
			return true;
		case RapidParam::OrigX:
			return true;
		case RapidParam::OrigY:
			return true;
		case RapidParam::OrigZ:
			return true;
		case RapidParam::ID:
			return false;	
		case RapidParam::decaytime:
			return false;	
		case RapidParam::hasOsc:
			return false;								
		case RapidParam::UNKNOWN:
			return true;
	}

	return 0.;

}

bool RapidParam::canBeTrue() {
	switch(type_) {
		case RapidParam::M:
			return true;
		case RapidParam::M2:
			return true;
		case RapidParam::MT:
			return true;
		case RapidParam::E:
			return true;
		case RapidParam::ET:
			return true;
		case RapidParam::P:
			return true;
		case RapidParam::PX:
			return true;
		case RapidParam::PY:
			return true;
		case RapidParam::PZ:
			return true;
		case RapidParam::PT:
			return true;
		case RapidParam::ETA:
			return true;
		case RapidParam::PHI:
			return true;
		case RapidParam::RAPIDITY:
			return true;
		case RapidParam::GAMMA:
			return true;
		case RapidParam::BETA:
			return true;
		case RapidParam::IP:
			return true;
		case RapidParam::SIGMAIP:
			return true;
		case RapidParam::MINIP:
			return true;
		case RapidParam::SIGMAMINIP:
			return true;
		case RapidParam::FD:
			return true;
		case RapidParam::ProbNNmu:
			return false;
		case RapidParam::ProbNNpi:
			return false;
		case RapidParam::ProbNNk:
			return false;
		case RapidParam::ProbNNp:
			return false;
		case RapidParam::ProbNNe:
			return false;
		case RapidParam::THETA:
			return true;
		case RapidParam::COSTHETA:
			return true;
		case RapidParam::MCORR:
			return true;
		case RapidParam::VtxX:
			return true;
		case RapidParam::VtxY:
			return true;
		case RapidParam::VtxZ:
			return true;
		case RapidParam::OrigX:
			return true;
		case RapidParam::OrigY:
			return true;
		case RapidParam::OrigZ:
			return true;
		case RapidParam::ID:
			return false;	
		case RapidParam::decaytime:
			return false;		
		case RapidParam::hasOsc:
			return false;									
		case RapidParam::UNKNOWN:
			return true;
	}

	return 0.;

}

double RapidParam::evalPID() {
	double pid(0.);
	if (particles_[0]->stable() && particles_[0]->mass() > 0.) {
		RapidParticleData * particleData = RapidParticleData::getInstance();
		unsigned int id(0);
		if (particles_[0]->massHypothesisName() == "") id = TMath::Abs(particles_[0]->id());
		else id = TMath::Abs(particleData->pdgCode(particles_[0]->massHypothesisName()));

		if(!pidHist_) {
			std::cout << "WARNING in RapidParam::evalPID : PID histograms not set for " << particles_[0]->id() << std::endl;
			std::cout << "                                 returning 0" << std::endl;
			return 0.;
		}

		pid = pidHist_->getPID(id,particles_[0]->getP().P()*1000.,particles_[0]->getP().Eta());
	}
	return pid;
}

// double RapidParam::evalID() {
	
// 	RapidParticleData * particleData = RapidParticleData::getInstance();
// 	int id(0);
// 	if (particles_[0]->massHypothesisName() == "") id = particles_[0]->id();
// 	else id = particleData->pdgCode(particles_[0]->massHypothesisName());
	
// 	return id;
// }
// double RapidParam::evalDecayTime() {
	
// 	RapidParticleData * particleData = RapidParticleData::getInstance();
// 	int decaytime(0);
// 	decaytime = particles_[0]->t();
// 	return decaytime;
// }

double RapidParam::evalCorrectedMass() {

	TLorentzVector momS, momT;

	//load the true (inc. invisible) and smeared momenta
	for(unsigned int i=0; i<particles_.size(); ++i) {
		momT += particles_[i]->getP();
		momS += particles_[i]->getPSmeared();
	}

	//get the sine and cosine of the angle to the B flight direction with a Gaussian smearing applied
	double cosDir = momT.Vect().Dot(momS.Vect())/(momT.P()*momS.P());
	double dir = TMath::ACos(cosDir) + gRandom->Gaus(0.,0.001); //TODO smearing parameter should be configurable
	double sinDir = TMath::Sin(dir);
	cosDir = TMath::Cos(dir);

	//get the longitudinal and transverse momentum of the visible daughters wrt the B flight direction
	//double pLong  = TMath::Abs(cosDir * momS.P());
	double pTran = TMath::Abs(sinDir * momS.P());

	//invariant masses of the visible daughters and the parent as well as the missing mass
	double mVis2 = momS.M2();
	//double mPar2 = momT.M2();
	//double mMiss2 = mPar2 - mVis2;

	//the corrected mass
	double mCorr = TMath::Sqrt( mVis2 + pTran*pTran ) + pTran;

	return mCorr;
}

double RapidParam::evalTheta() {
	if(particles_.size()<2) {
		return -99.;
	}

	TLorentzVector pA, pB, pBoost;
	bool doBoost(false);

	if(truth_) {
		pA = particles_[0]->getP();
		pB = particles_[1]->getP();
		for(unsigned int i=2; i<particles_.size(); ++i) {
			pBoost += particles_[i]->getP();
			doBoost=true;
		}
	} else {
		pA = particles_[0]->getPSmeared();
		pB = particles_[1]->getPSmeared();
		for(unsigned int i=2; i<particles_.size(); ++i) {
			pBoost += particles_[i]->getPSmeared();
			doBoost=true;
		}
	}

	if(doBoost) {
		pA.Boost(-pBoost.BoostVector());
		pB.Boost(-pBoost.BoostVector());
	}

	if(type_ == RapidParam::THETA) {
		return pA.Angle(pB.Vect());
	} else {
		return pA.Vect().Dot(pB.Vect())/(pA.P()*pB.P());
	}
}

double RapidParam::evalPos() {
	RapidVertex* vtx(0);
	if(type_==RapidParam::OrigX ||
	   type_==RapidParam::OrigY ||
	   type_==RapidParam::OrigZ ) {
		vtx = particles_[0]->getOriginVertex();
	} else if(type_==RapidParam::VtxX ||
		  type_==RapidParam::VtxY ||
		  type_==RapidParam::VtxZ ) {
		vtx = particles_[0]->getDecayVertex();
		if(particles_[0]->stable()) {
			std::cout << "WARNING in RapidParam::evalPos : requesting " << type_ << " of a particle that has not decayed." << std::endl
				  << "Will use origin vertex instead." << std::endl;
		}
	}
	if(!vtx) {
		std::cout << "WARNING in RapidParam::evalPos : vertex for parameter " << type_ << " not found for particle " << particles_[0]->name() << "." << std::endl;
		return -999.;
	}
	ROOT::Math::XYZPoint point = vtx->getVertex(truth_);
	if(type_==RapidParam::OrigX || type_==RapidParam::VtxX) return point.X();
	if(type_==RapidParam::OrigY || type_==RapidParam::VtxY) return point.Y();
	if(type_==RapidParam::OrigZ || type_==RapidParam::VtxZ) return point.Z();

	std::cout << "WARNING in RapidParam::evalPos : unknown position parameter type " << type_ << "." << std::endl;
	return -999.;
}

TString RapidParam::name() {
	if(name_.Length()==0) {
		//construct default name
		std::vector<RapidParticle*>::iterator it = particles_.begin();
		for( ; it!= particles_.end(); ++it) {
			name_ += (*it)->name();
			name_ += "_";
		}

		name_ += typeName();

		if(truth_) {
			name_ += "_TRUE";
		}
	}
	return name_;
}

TString RapidParam::typeName() {
	switch(type_) {
		case RapidParam::M:
			return "M";
		case RapidParam::M2:
			return "M2";
		case RapidParam::MT:
			return "MT";
		case RapidParam::IP:
			return "IP";
		case RapidParam::SIGMAIP:
			return "SIGMAIP";
		case RapidParam::MINIP:
			return "MINIP";
		case RapidParam::SIGMAMINIP:
			return "SIGMAMINIP";
		case RapidParam::FD:
			return "FD";
		case RapidParam::E:
			return "E";
		case RapidParam::ET:
			return "ET";
		case RapidParam::P:
			return "P";
		case RapidParam::PX:
			return "PX";
		case RapidParam::PY:
			return "PY";
		case RapidParam::PZ:
			return "PZ";
		case RapidParam::PT:
			return "PT";
		case RapidParam::ETA:
			return "eta";
		case RapidParam::PHI:
			return "phi";
		case RapidParam::RAPIDITY:
			return "y";
		case RapidParam::GAMMA:
			return "gamma";
		case RapidParam::BETA:
			return "beta";
		case RapidParam::THETA:
			return "theta";
		case RapidParam::COSTHETA:
			return "costheta";
		case RapidParam::MCORR:
			return "Mcorr";
		case RapidParam::ProbNNmu:
			return "ProbNNmu";
		case RapidParam::ProbNNpi:
			return "ProbNNpi";
		case RapidParam::ProbNNk:
			return "ProbNNk";
		case RapidParam::ProbNNp:
			return "ProbNNp";
		case RapidParam::ProbNNe:
			return "ProbNNe";
		case RapidParam::VtxX:
			return "vtxX";
		case RapidParam::VtxY:
			return "vtxY";
		case RapidParam::VtxZ:
			return "vtxZ";
		case RapidParam::OrigX:
			return "origX";
		case RapidParam::OrigY:
			return "origY";
		case RapidParam::OrigZ:
			return "origZ";
		case RapidParam::ID:
			return "ID";	
		case RapidParam::decaytime:
			return "decaytime";		
		case RapidParam::hasOsc:
			return "hasOsc";									
		default:
			std::cout << "WARNING in RapidParam::typeName : unknown type " << type_ << "." << std::endl
				  << "                                  returning empty string." << std::endl;
			return "";
	}
}

RapidParam::ParamType RapidParam::typeFromString(TString str) {
	if(str=="M") {
		return RapidParam::M;
	} else if(str=="M2") {
		return RapidParam::M2;
	} else if(str=="MT") {
		return RapidParam::MT;
	} else if(str=="IP") {
		return RapidParam::IP;
	} else if(str=="SIGMAIP") {
		return RapidParam::SIGMAIP;
	} else if(str=="MINIP") {
		return RapidParam::MINIP;
	} else if(str=="SIGMAMINIP") {
		return RapidParam::SIGMAMINIP;
	} else if(str=="FD") {
		return RapidParam::FD;
	} else if(str=="E") {
		return RapidParam::E;
	} else if(str=="ET") {
		return RapidParam::ET;
	} else if(str=="P") {
		return RapidParam::P;
	} else if(str=="PX") {
		return RapidParam::PX;
	} else if(str=="PY") {
		return RapidParam::PY;
	} else if(str=="PZ") {
		return RapidParam::PZ;
	} else if(str=="PT") {
		return RapidParam::PT;
	} else if(str=="eta") {
		return RapidParam::ETA;
	} else if(str=="phi") {
		return RapidParam::PHI;
	} else if(str=="y") {
		return RapidParam::RAPIDITY;
	} else if(str=="gamma") {
		return RapidParam::GAMMA;
	} else if(str=="beta") {
		return RapidParam::BETA;
	} else if(str=="theta") {
		return RapidParam::THETA;
	} else if(str=="costheta") {
		return RapidParam::COSTHETA;
	} else if(str=="Mcorr") {
		return RapidParam::MCORR;
	} else if(str=="ProbNNmu") {
		return RapidParam::ProbNNmu;
	} else if(str=="ProbNNpi") {
		return RapidParam::ProbNNpi;
	} else if(str=="ProbNNk") {
		return RapidParam::ProbNNk;
	} else if(str=="ProbNNp") {
		return RapidParam::ProbNNp;
	} else if(str=="ProbNNe") {
		return RapidParam::ProbNNe;
	} else if(str=="vtxX") {
		return RapidParam::VtxX;
	} else if(str=="vtxY") {
		return RapidParam::VtxY;
	} else if(str=="vtxZ") {
		return RapidParam::VtxZ;
	} else if(str=="origX") {
		return RapidParam::OrigX;
	} else if(str=="origY") {
		return RapidParam::OrigY;
	} else if(str=="origZ") {
		return RapidParam::OrigZ;
	} else if(str=="ID") {
		return RapidParam::ID;
	} else if(str=="decaytime") {
		return RapidParam::decaytime;
	} else if(str=="hasOsc") {
		return RapidParam::hasOsc;		
	} else {
		std::cout << "WARNING in RapidParam::typeFromString : unknown type name " << str << "." << std::endl
			  << "                                        returning mass parameter type." << std::endl;
		return RapidParam::M;
	}
}

void RapidParam::getMinMax(const std::vector<RapidParticle*>& parts, double& min, double& max) {
	setDefaultMinMax(parts,min,max);
}

void RapidParam::getMinMax(double& min, double& max, bool recalculate) {
	if(!recalculate) {
		min = minVal_;
		max = maxVal_;
	} else {
		setDefaultMinMax(particles_,min,max);
	}
}

void RapidParam::getMinMax(RapidParticle* part, double& min, double& max) {
	std::vector<RapidParticle*> parts;
	parts.push_back(part);
	setDefaultMinMax(parts,min,max);
}

void RapidParam::getMinMax(RapidParticle* partA, RapidParticle* partB, double& min, double& max) {
	std::vector<RapidParticle*> parts;
	parts.push_back(partA);
	parts.push_back(partB);
	setDefaultMinMax(parts,min,max);
}

void RapidParam::getMinMax(RapidParticle* partA, RapidParticle* partB, RapidParticle* partC, double& min, double& max) {
	std::vector<RapidParticle*> parts;
	parts.push_back(partA);
	parts.push_back(partB);
	parts.push_back(partC);
	setDefaultMinMax(parts,min,max);
}

void RapidParam::setDefaultMinMax(const std::vector<RapidParticle*>& parts, double& min, double& max) {
	switch(type_) {
		case RapidParam::M:
			setMassMinMax(parts,min,max);
			break;
		case RapidParam::MT:
			setMassMinMax(parts,min,max);
			min = 0.;
			break;
		case RapidParam::MCORR:
			setMassMinMax(parts,min,max);
			min = -max;
			break;
		case RapidParam::M2:
			setMassMinMax(parts,min,max);
			min = min*min;
			max = max*max;
			break;
		case RapidParam::IP:
		case RapidParam::SIGMAIP:
		case RapidParam::MINIP:
		case RapidParam::SIGMAMINIP:
		case RapidParam::FD:
		case RapidParam::E:
		case RapidParam::ET:
		case RapidParam::P:
		case RapidParam::PX:
		case RapidParam::PY:
		case RapidParam::PZ:
		case RapidParam::PT:
			min = 0.;
			max = 100.;
			break;
		case RapidParam::ETA:
			min = -7.;
			max =  7.;
			break;
		case RapidParam::PHI:
			min = -3.5;
			max =  3.5;
			break;
		case RapidParam::RAPIDITY:
			min = 0.;
			max = 10.;
			break;
		case RapidParam::GAMMA:
			min = 1.;
			max = 50.;
			break;
		case RapidParam::BETA:
			min = 0.;
			max = 1.;
			break;
		case RapidParam::THETA:
			min = 0.0;
			max = 3.5;
			break;
		case RapidParam::COSTHETA:
			min = -1.0;
			max =  1.0;
			break;
		case RapidParam::ProbNNmu:
			min = 0.0;
			max = 1.0;
			break;
		case RapidParam::ProbNNpi:
			min = 0.0;
			max = 1.0;
			break;
		case RapidParam::ProbNNk:
			min = 0.0;
			max = 1.0;
			break;
		case RapidParam::ProbNNp:
			min = 0.0;
			max = 1.0;
			break;
		case RapidParam::ProbNNe:
			min = 0.0;
			max = 1.0;
			break;
		case RapidParam::VtxX:
		case RapidParam::VtxY:
		case RapidParam::OrigX:
		case RapidParam::OrigY:
			min = -100.0;
			max = 100.0;
			break;
		case RapidParam::VtxZ:
		case RapidParam::OrigZ:
			min = -1000.0;
			max = 1000.0;
			break;
		case RapidParam::UNKNOWN:
		default:
			return;
	}
}

void RapidParam::setMassMinMax(const std::vector<RapidParticle*>& parts, double& min, double& max) {
	if(parts.size()==0) {
		//no particles so mass is 0
		return;
	} else if(parts.size()==1) {
		//single particle is easy
		min = parts[0]->minMass() - 0.1;
		max = parts[0]->maxMass() + 0.1;

		//now check to see if we've used any alternative mass hypotheses and extend range if necessary
		double deltaDown(0.), deltaUp(0.);

		RapidParticleData::getInstance()->getMaxAltHypothesisMassShifts(parts,deltaDown,deltaUp);

		//TODO factor of 2. is arbitrary - there is probably a better way to calculate these
		min+=2.*deltaDown;
		max+=2.*deltaUp;

		if(min<0.) min=0.;
	} else {
		RapidParticleData* rpd = RapidParticleData::getInstance();

		//first check whether any of the particles we've been given have a hierarchical relationship
		if(rpd->checkHierarchy(parts)) {
			std::cout << "WARNING in RapidParam::setMassMinMax : some of the particles used in parameter " << name() << " are repeated or have a hierarchical relationship." << std::endl
				  << "                                       setting default mass range." << std::endl;
			min=0.;
			max=10.;
			return;
		}

		//find the most recent common ancestor of all our particles
		RapidParticle* commonAncestor = rpd->findCommonAncestor(parts);

		//this shouldn't happen
		if(!commonAncestor) {
			std::cout << "WARNING in RapidParam::setMassMinMax : some of the particles used in parameter " << name() << " have no common ancestors." << std::endl
				  << "                                       setting default mass range." << std::endl;
			min=0.;
			max=10.;
			return;
		}

		//get the list of stable particles descended from our particles
		std::vector<RapidParticle*> selStableParts;
		rpd->findStableDaughters(parts,selStableParts);

		//now compare against the list of stable particles from the ancestor to get the list of particles we did not include
		std::vector<RapidParticle*> otherParts;
		rpd->findOtherDaughters(commonAncestor, selStableParts, otherParts);

		//now combine any particles where all of the daughters are in this list as the breakup mass is also lost
		std::vector<RapidParticle*> otherPartsCombined;
		rpd->combineCompleteAncestors(otherParts, otherPartsCombined);

		//minimum mass is the sum of the masses of the included particles
		min = -0.1;
		std::vector<RapidParticle*>::const_iterator cit = parts.begin();
		for( ; cit!=parts.end(); ++cit) {
			min += (*cit)->minMass();
		}

		//maximum mass is the mass of the parent minus the sum of the masses of the missing daughters
		max = commonAncestor->maxMass() + 0.1;
		std::vector<RapidParticle*>::iterator it = otherPartsCombined.begin();
		for( ; it!=otherPartsCombined.end(); ++it) {
			max -= (*it)->minMass();
		}

		//now check to see if we've used any alternative mass hypotheses and extend range if necessary
		double deltaDown(0.), deltaUp(0.);

		RapidParticleData::getInstance()->getMaxAltHypothesisMassShifts(parts,deltaDown,deltaUp);

		//TODO factor of 2. is arbitrary - there is probably a better way to calculate these
		min+=2.*deltaDown;
		max+=2.*deltaUp;

		if(min<0.) min=0.;
	}
}