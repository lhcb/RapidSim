#include "RapidParticle.h"

#include <iostream>

#include "TMath.h"
#include "TRandom.h"

#include "RapidIPSmearGauss.h"
#include "RapidMomentumSmearEnergyGauss.h"
#include "RapidMomentumSmearGauss.h"
#include "RapidMomentumSmearHisto.h"
#include "RapidParticleData.h"

void RapidParticle::addDaughter(RapidParticle* part) {
	if(!daughters_.empty()) {
		daughters_[daughters_.size()-1]->next_ = part;
	}
	part->index_ = daughters_.size();
	daughters_.push_back(part);
	daughterMasses_.push_back(part->mass_);
}

void RapidParticle::addMassHypothesis(TString name, double mass) {
	altMasses_.push_back(mass);
	massHypothesisNames_.push_back(name);
}

void RapidParticle::setMassHypothesis(unsigned int i) {
	if(currentHypothesis_==i) {
		return;
	}
	if(i >= nMassHypotheses()) {
		std::cout << "WARNING in RapidParticle::setMassHypothesis : " << name() << " does not have mass hypothesis " << i << std::endl;
		return;
	}

	currentHypothesis_ = i;
	if(currentHypothesis_==0) {
		p_.SetXYZM(p_.Px(), p_.Py(), p_.Pz(), mass_);
		pSmeared_.SetXYZM(pSmeared_.Px(), pSmeared_.Py(), pSmeared_.Pz(), mass_);
	} else {
		p_.SetXYZM(p_.Px(), p_.Py(), p_.Pz(), altMasses_[i-1]);
		pSmeared_.SetXYZM(pSmeared_.Px(), pSmeared_.Py(), pSmeared_.Pz(), altMasses_[i-1]);
	}
	if(mother_) {
		mother_->updateMomenta();
	}
}

void RapidParticle::smearMomentum() {
	if(nDaughters() == 0) {
		//smear momentum
		if(invisible_) {
			pSmeared_.SetPxPyPzE(0.,0.,0.,0.);
		} else if(momSmear_) {
			pSmeared_ = momSmear_->smearMomentum(p_);
		} else {
			pSmeared_ = p_;
		}
	} else {
		//reconstruct mothers from their daughters
		pSmeared_.SetPxPyPzE(0.,0.,0.,0.);
		RapidParticle* daug = daughter(0);
		for( ; daug!=0; daug = daug->next()) {
			pSmeared_ += daug->pSmeared_;
		}
	}
}

double RapidParticle::getFD(bool truth) {
	if(nDaughters() == 0) { // If stable true flight distance is infinite, return -1 as default
		return -1.;
	} else {
		ROOT::Math::XYZPoint decVtx = decayVertex_->getVertex(truth);
		ROOT::Math::XYZPoint oriVtx = originVertex_->getVertex(truth);

		return sqrt(pow(decVtx.X()-oriVtx.X(),2) + \
				pow(decVtx.Y()-oriVtx.Y(),2) + \
				pow(decVtx.Z()-oriVtx.Z(),2) );
	}
}

void RapidParticle::smearIP() {
	if(nDaughters() != 0) {
		// Do not smear the IP of decaying particles
		// It is a derived quantity which can be computed from their
		// smeared origin & decay vertices and momentum
		ipSmeared_ = ip_;
		sigmaip_ = 0.;
		minipSmeared_ = minip_;
		sigmaminip_ = 0.;
	} else {
		if(invisible_) {
			ipSmeared_ = ip_;
			sigmaip_ = 0.;
			minipSmeared_ = minip_;
			sigmaminip_ = 0.;
		} else if (ipSmear_) {
			std::pair<double,double> smearedips = ipSmear_->smearIP(ip_,p_.Pt());
			ipSmeared_ = smearedips.first;
			sigmaip_   = smearedips.second;
			std::pair<double,double> smearedminips = ipSmear_->smearIP(minip_,p_.Pt());
			minipSmeared_ = smearedminips.first;
			sigmaminip_   = smearedminips.second;
		} else {
			ipSmeared_ = ip_;
			sigmaip_ = 0.;
			minipSmeared_ = minip_;
			sigmaminip_ = 0.;
		}
	}
}

double RapidParticle::deltaMass() {
	if(currentHypothesis_==0) {
		return 0.;
	} else {
		return altMasses_[currentHypothesis_-1] - mass_;
	}
}

TString RapidParticle::massHypothesisName() {
	if(currentHypothesis_==0) return "";
	else return massHypothesisNames_[currentHypothesis_-1];
}

bool RapidParticle::getAltMass(TString name, double& mass) {
	//look up an alternative mass hypothesis (added with altMass) by its sanitised name
	for(unsigned int i=0; i<massHypothesisNames_.size(); ++i) {
		if(massHypothesisNames_[i] == name) {
			mass = altMasses_[i];
			return true;
		}
	}
	return false;
}

bool RapidParticle::hasFlavour(int flavour) {
	int id = TMath::Abs(id_);

	if(id==flavour) return true;
	if(id<100) return false;

	if((id/10) %10 == flavour) return true;
	if((id/100) %10 == flavour) return true;
	if((id/1000) %10 == flavour) return true;

	return false;
}

bool RapidParticle::hasCharm() {
	return hasFlavour(4);
}

bool RapidParticle::hasBeauty() {
	return hasFlavour(5);
}

RapidParticle* RapidParticle::daughter(unsigned int i) {
	if(daughters_.size()>i) {
		return daughters_[i];
	}
	return 0;
}

void RapidParticle::setMass(double mass) {
	mass_=mass;
	if(mother_) {
		mother_->updateDaughterMass(index_);
	}
}

void RapidParticle::print(int index) {
	TString mname = "---";
	if(mother_) mname = mother_->name_;
	TString dname = "";
	RapidParticle* daug = daughter(0);
	if(!daug) dname = "---";
	for( ; daug!=0; daug=daug->next()) {
		dname += daug->name_;
		if(daug->next()) dname+=", ";
	}

	printf("%3d\t%-15s\t%6d\t\t%.6f\t%-15s\t%2d\t\t%-15s\n", index, name_.Data(), id_, mass_, mname.Data(), nDaughters(), dname.Data());
}

void RapidParticle::setMassShape(RooDataSet* ds, double minMass, double maxMass, TString varName) {
	massData_ = ds;
	minMass_ = minMass;
	maxMass_ = maxMass;
	varName_ = varName;
}

void RapidParticle::floatMass() {
	if(massData_) {
		int entry = static_cast<int>(massData_->numEntries() * gRandom->Uniform());
		const RooArgSet* row = massData_->get(entry);
		double value = row->getRealValue(varName_);
		setMass(value);
	}
}

void RapidParticle::updateDaughterMass(unsigned int index) {
	if(index<daughters_.size()) {
		daughterMasses_[index] = daughters_[index]->mass_;
	}
}

void RapidParticle::updateMomenta() {
	p_.SetPxPyPzE(0.,0.,0.,0.);
	pSmeared_.SetPxPyPzE(0.,0.,0.,0.);

	RapidParticle* daug = daughter(0);

	for( ; daug!=0; daug = daug->next()) {
		p_ += daug->p_;
		pSmeared_ += daug->pSmeared_;
	}
	if(mother_) {
		mother_->updateMomenta();
	}
}

void RapidParticle::setupVertices() {
	if(mother_) {
		originVertex_ = mother_->decayVertex_;
	} else {
		originVertex_ = new RapidVertex(0.,0.,0.);
	}
	if(ctau_>0) {
		decayVertex_ = new RapidVertex(0.,0.,0);
	} else decayVertex_ = originVertex_;
}
