#include "RapidParticleData.h"

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "RooRelBreitWigner.h"
#include "RooSubThreshold.h"
#include "RooGounarisSakurai.h"

#include "RapidParticle.h"

RapidParticleData* RapidParticleData::instance_=0;

RapidParticleData* RapidParticleData::getInstance() {
	if(!instance_) {
		instance_ = new RapidParticleData();
		TString path;
		path=getenv("RAPIDSIM_CONFIG");
		if(path!="") instance_->loadData(path+"/config/particles.dat");

		path=getenv("RAPIDSIM_ROOT");
		instance_->loadData(path+"/config/particles.dat");
	}
	return instance_;
}

void RapidParticleData::loadData(TString file) {
	std::cout << "INFO in RapidParticleData::loadData : loading particle data from " << file << std::endl;

	std::ifstream fin;
	fin.open(file, std::ifstream::in);
	TString buffer;
	buffer.ReadLine(fin);//ignore title line

	int id;
	TString part, anti;
	double mass, width, spin, charge, ctau;
	TString lineshape;

	while(fin.good()) {
		buffer.ReadToken(fin);
		id = buffer.Atoi();
		part.ReadToken(fin);
		anti.ReadToken(fin);
		buffer.ReadToken(fin);
		mass = buffer.Atof();
		buffer.ReadToken(fin);
		width = buffer.Atof();
		buffer.ReadToken(fin);
		charge = buffer.Atof();
		buffer.ReadToken(fin);
		spin = buffer.Atof();
		lineshape.ReadToken(fin);
		buffer.ReadToken(fin);
		ctau = buffer.Atof();

		if(!fin.good()) break;
		addEntry(id,part,mass,width,spin,charge,lineshape,ctau);
		if(anti != "---") addEntry(-1*id,anti,mass,width,spin,-1.*charge,lineshape,ctau);

	}
	fin.close();
}

double RapidParticleData::getMass(int id) {
	if(idToMass_.count(id)) {
		return idToMass_[id];
	} else {
		std::cout << "WARNING in RapidParticleData::getMass : Unknown particle ID " << id << std::endl;
		return 0.;
	}
}

double RapidParticleData::getCT(int id) {
	if(idToCT_.count(id)) {
		return idToCT_[id];
	} else {
		std::cout << "WARNING in RapidParticleData::getCT : Unknown particle ID " << id << std::endl;
		return 0.;
	}
}

double RapidParticleData::getWidth(int id) {
	if(idToWidth_.count(id)) {
		return idToWidth_[id];
	} else {
		std::cout << "WARNING in RapidParticleData::getWidth : Unknown particle ID " << id << std::endl;
		return 0.;
	}
}

double RapidParticleData::getSpin(int id) {
	if(idToSpin_.count(id)) {
		return idToSpin_[id];
	} else {
		std::cout << "WARNING in RapidParticleData::getSpin : Unknown particle ID " << id << std::endl;
		return 0.;
	}
}

double RapidParticleData::getCharge(int id) {
	if(idToCharge_.count(id)) {
		return idToCharge_[id];
	} else {
		std::cout << "WARNING in RapidParticleData::getCharge : Unknown particle ID " << id << std::endl;
		return 0.;
	}
}

TString RapidParticleData::getName(int id) {
	if(idToName_.count(id)) {
		return idToName_[id];
	} else {
		std::cout << "WARNING in RapidParticleData::getName : Unknown particle ID " << id << std::endl;
		return "";
	}
}

TString RapidParticleData::getSanitisedName(int id) {
	TString name = getName(id);
	return sanitiseName(name);
}

int RapidParticleData::pdgCode(TString name) {
	if(nameToId_.count(name)) {
		return nameToId_[name];
	} else if(sanitisedNameToId_.count(name)) {
		return sanitisedNameToId_[name];
	} else {
		std::cout << "WARNING in RapidParticleData::getId : Unknown particle name " << name << std::endl;
		return 0.;
	}
}

RapidParticle* RapidParticleData::makeParticle(int id, RapidParticle* mother) {

	double ctau = getCT(id);
	double mass = getMass(id);
	double charge = getCharge(id);
	TString name = getName(id);
	name = makeUniqName(name);

	return new RapidParticle(id, name, mass, charge, ctau, mother);
}

RapidParticle* RapidParticleData::makeParticle(TString name, RapidParticle* mother) {

	int id = pdgCode(name);
	double ctau = getCT(id);
	double mass = getMass(id);
	double charge = getCharge(id);
	name = makeUniqName(name);

	return new RapidParticle(id, name, mass, charge, ctau, mother);
}

void RapidParticleData::setupMass(RapidParticle* part) {
	unsigned int id = TMath::Abs(part->id());

	//check if this is a resonance that we know how to make
	if(!idToShape_.count(id)) return;
	ResLineShape shape = idToShape_[id];
	if(part->nDaughters() != 2) return;

	unsigned int idA = TMath::Abs(part->daughter(0)->id());
	unsigned int idB = TMath::Abs(part->daughter(1)->id());

	double mA = part->daughter(0)->mass();
	double mB = part->daughter(1)->mass();
	double mC = part->daughter(2)->mass();
	double mPar = part->mass();

	if(idA>idB) {
		int tmpId = idB;
		idB = idA;
		idA = tmpId;
		double tmpMass = mB;
		mB = mA;
		mA = tmpMass;
	}

	TString varName = "m";
	varName += part->name();

	TString name = part->name();
	double mass = getMass(id);
	double width = getWidth(id);
	double spin = getSpin(id);

	if(width<narrowWidth_) {
		std::cout << "INFO in RapidParticleData::setupMass : resonance " << name << " is narrow (" << width << " GeV) and will generated using a fixed mass." << std::endl;
		return;
	}

	double mmin = mass - 100.*width;
	double mmax = mass + 100.*width;

	RooRealVar m(varName,varName,mmin,mmax);
	RooAbsPdf* pdf(0);
	switch(shape) {
		case RapidParticleData::GS:
			pdf = makeGS(m, mass, width, spin, mA, mB, name);
			break;
		case RapidParticleData::SubTh:
			  pdf = makeSubTh(m, mass, width, spin, mA, mB, mC, mPar, name);
				break;
		default:
			std::cout << "WARNING in RapidParticleData::setupMass : unknown lineshape for " << name << "." << std::endl
				  << "                                        : using a relativistic Breit-Wigner." << std::endl;
		/* FALLTHRU */
		case RapidParticleData::RelBW:
			pdf = makeRelBW(m, mass, width, spin, mA, mB, name);

	}
	RooDataSet* massdata = pdf->generate(RooArgSet(m),100000);
	massdata->getRange(m,mmin,mmax);
	part->setMassShape(massdata,mmin,mmax,varName);
}

void RapidParticleData::addEntry(int id, TString name, double mass, double width, double spin, double charge, TString lineshape, double ctau) {
	if(idToName_.find(id) != idToName_.end()) {
		std::cout << "INFO in RapidParticleData::addEntry : particle with ID " << id << " already defined with name " << idToName_[id] << std::endl;
		std::cout << "                                      second definition will be ignored." << std::endl;
	}

	idToCT_[id] = ctau;
	idToMass_[id] = mass;
	idToWidth_[id] = width;
	idToSpin_[id] = spin;
	idToCharge_[id] = charge;
	idToName_[id] = name;
	nameToId_[name] = id;
	sanitisedNameToId_[sanitiseName(name)] = id;
	if(lineshape=="RBW") idToShape_[id] = RapidParticleData::RelBW;
	else if(lineshape =="SubTh")idToShape_[id] = RapidParticleData::SubTh;
	else if(lineshape=="GS") idToShape_[id] = RapidParticleData::GS;
}

TString RapidParticleData::sanitiseName(TString name) {
	name = name.ReplaceAll("+","p");
	name = name.ReplaceAll("-","m");
	name = name.ReplaceAll("*","st");
	name = name.ReplaceAll("(","_");
	name = name.ReplaceAll(")","_");
	name = name.ReplaceAll("[","_");
	name = name.ReplaceAll("]","_");
	name = name.ReplaceAll("<","_");
	name = name.ReplaceAll(">","_");
	name = name.ReplaceAll("{","_");
	name = name.ReplaceAll("}","_");
	name = name.ReplaceAll(" ","_");
	name = name.ReplaceAll("$","");
	name = name.ReplaceAll("%","");
	name = name.ReplaceAll("&","");
	name = name.ReplaceAll("/","");
	name = name.ReplaceAll(":","");
	name = name.ReplaceAll(";","");
	name = name.ReplaceAll("=","");
	name = name.ReplaceAll("\\","");
	name = name.ReplaceAll("^","");
	name = name.ReplaceAll("|","");
	name = name.ReplaceAll(",","");
	name = name.ReplaceAll(".","");
	name.Remove(TString::kBoth,'_');

	return name;

}

TString RapidParticleData::makeUniqName(TString name) {

	TString uniqName("");

	int i=-1;
	name = sanitiseName(name);

	do {
		++i;
		uniqName = name;
		uniqName+= "_";
		uniqName+= i;
	} while(usedNames_.count(uniqName)>0);

	usedNames_.insert(uniqName);

	return uniqName;

}

RooRelBreitWigner* RapidParticleData::makeRelBW(RooRealVar& m, double mean, double gamma, double thespin, double m1, double m2, TString name) {

	double barrierFactor=3.;//TODO

	RooRealVar* m0     = new RooRealVar(name+"m0",    name+"m0",     mean);
	RooRealVar* g0     = new RooRealVar(name+"g0",    name+"g0",    gamma);
	RooRealVar* spin   = new RooRealVar(name+"spin",  name+"spin",  thespin);
	RooRealVar* radius = new RooRealVar(name+"radius",name+"radius",barrierFactor); // not used
	RooRealVar* ma     = new RooRealVar(name+"ma",    name+"ma",     m1);
	RooRealVar* mb     = new RooRealVar(name+"mb",    name+"mb",     m2);

	return new RooRelBreitWigner(name,name, m,*m0,*g0,*radius,*ma,*mb,*spin);
}



RooSubThreshold* RapidParticleData::makeSubTh(RooRealVar& m, double mean, double gamma, double thespin, double m1, double m2, double m3, double mPar,TString name) {

	double barrierFactor=3.;//TODO

	RooRealVar* m0     = new RooRealVar(name+"m0",    name+"m0",     mean);
	RooRealVar* g0     = new RooRealVar(name+"g0",    name+"g0",    gamma);
	RooRealVar* spin   = new RooRealVar(name+"spin",  name+"spin",  thespin);
	RooRealVar* radius = new RooRealVar(name+"radius",name+"radius",barrierFactor); // not used
	RooRealVar* ma     = new RooRealVar(name+"ma",    name+"ma",     m1);
	RooRealVar* mb     = new RooRealVar(name+"mb",    name+"mb",     m2);
	RooRealVar* mc     = new RooRealVar(name+"mc",    name+"mc",     m3);
	RooRealVar* mParent= new RooRealVar(name+"mParent",    name+"mParent",     mPar);
	return new RooSubThreshold(name,name, m,*m0,*g0,*radius,*ma, *mb, *mc, *mParent, *spin);
}


RooGounarisSakurai* RapidParticleData::makeGS(RooRealVar& m, double mean, double gamma, double thespin, double m1, double m2, TString name) {

	double barrierFactor=3.;//TODO

	RooRealVar* m0     = new RooRealVar(name+"m0",    name+"m0",     mean);
	RooRealVar* g0     = new RooRealVar(name+"g0",    name+"g0",    gamma);
	RooRealVar* spin   = new RooRealVar(name+"spin",  name+"spin",  thespin);
	RooRealVar* radius = new RooRealVar(name+"radius",name+"radius",barrierFactor); // not used
	RooRealVar* ma     = new RooRealVar(name+"ma",    name+"ma",     m1);
	RooRealVar* mb     = new RooRealVar(name+"mb",    name+"mb",     m2);

	return new RooGounarisSakurai(name,name, m,*m0,*g0,*spin, *radius,*ma,*mb);
}

bool RapidParticleData::checkHierarchy(const std::vector<RapidParticle*>& parts) {
	std::vector<RapidParticle*>::const_iterator it1 = parts.begin();

	for( ; it1!=parts.end(); ++it1) {
		RapidParticle* part1 = *it1;

		std::vector<RapidParticle*>::const_iterator it2 = parts.begin();
		for( ; it2!=parts.end(); ++it2) {
			if(it1==it2) continue;
			RapidParticle* part2 = *it2;

			//check hierarchy for each pair of particles in both orders
			if(checkHierarchy(part1,part2)) return true;
		}
	}

	return false;
}

bool RapidParticleData::checkHierarchy(RapidParticle* part, RapidParticle* ancestor) {
	while(part) {
		//search the hierarchy of part until we either find ancestor or 0
		if(ancestor == part) return true;
		part = part->mother();
	}

	return false;
}

RapidParticle* RapidParticleData::findCommonAncestor(const std::vector<RapidParticle*>& parts) {
	std::vector<RapidParticle*>::const_iterator it = parts.begin();

	//common ancestor of the first particle is itself
	RapidParticle* ancestor = *it;

	for( ; it!=parts.end(); ++it) {
		RapidParticle* part = *it;
		//if the current ancestor is not an ancestor of the next particle then search upwards until we find one that is
		while(!checkHierarchy(part,ancestor)) {
			ancestor = ancestor->mother();
			if(!ancestor) {
				//if we arrive here then the particles do not belong to the same decay chain
				std::cout << "WARNING in RapidParticlleData::findCommonAncestor : some of the particles have no common ancestors." << std::endl;
				return 0;
			}
		}
	}

	return ancestor;

}

void RapidParticleData::findStableDaughters(const std::vector<RapidParticle*>& parts, std::vector<RapidParticle*>& daughters) {
	std::vector<RapidParticle*>::const_iterator it = parts.begin();

	for( ; it!=parts.end(); ++it) {
		RapidParticle* part = *it;
		//add the stable daughters of each particle into daughters
		findStableDaughters(part, daughters);
	}
}

void RapidParticleData::findStableDaughters(RapidParticle* part, std::vector<RapidParticle*>& daughters) {
	if(part->nDaughters()==0) {
		//if the particle is stable then keep it
		daughters.push_back(part);
	} else {
		//otherwise iterate over all daughters
		RapidParticle* daug = part->daughter(0);
		while(daug) {
			findStableDaughters(daug,daughters);
			daug = daug->next();
		}
	}
}

void RapidParticleData::findOtherDaughters(RapidParticle* ancestor, const std::vector<RapidParticle*>& parts, std::vector<RapidParticle*>& others) {
	std::vector<RapidParticle*> all;
	std::vector<RapidParticle*> sel;

	//find all daughters of the ancestor
	findStableDaughters(ancestor, all);

	//find all daughters of the particles
	findStableDaughters(parts, sel);

	//sort the vectors so we can diff them
	std::sort(all.begin(), all.end());
	std::sort(sel.begin(), sel.end());

	//fill others with all elements that are in all but not sel
	std::set_difference( all.begin(), all.end(), sel.begin(), sel.end(),std::back_inserter( others ) );
}

void RapidParticleData::combineCompleteAncestors(const std::vector<RapidParticle*>& parts, std::vector<RapidParticle*>& partsCombined) {

	bool updated(false);

	std::map<RapidParticle*,unsigned int> daughtersFound;

	//conut the number of daughters found for the mother of each particle
	std::vector<RapidParticle*>::const_iterator it = parts.begin();
	for( ; it!=parts.end(); ++it) {
		RapidParticle* mother = (*it)->mother();
		if(mother) {
			++daughtersFound[mother];
		}
	}

	//add in any mothers which are complete
	std::map<RapidParticle*,unsigned int>::iterator mit = daughtersFound.begin();
	for( ; mit!=daughtersFound.end(); ++mit) {
		RapidParticle* mother = (*mit).first;
		if(mother->nDaughters() == (*mit).second) {
			partsCombined.push_back(mother);
			updated = true;
		}
	}

	//now add in the particles whose mothers are not complete
	it = parts.begin();
	for( ; it!=parts.end(); ++it) {
		RapidParticle* part = (*it);
		RapidParticle* mother = part->mother();
		if(!mother || mother->nDaughters()!=daughtersFound[mother]) {
			partsCombined.push_back(part);
		}
	}

	//if we successfully combined any particles then we may also be able to combine their parents so iterate
	if(updated) {
		std::vector<RapidParticle*> partsUpdated;
		std::swap(partsCombined,partsUpdated);//TODO C++11 - can use std::move
		combineCompleteAncestors(partsUpdated,partsCombined);
	}
}

void RapidParticleData::getMaxAltHypothesisMassShifts(const std::vector<RapidParticle*>& parts, double& deltaDown, double& deltaUp) {
	std::vector<RapidParticle*>::const_iterator it = parts.begin();
	for( ; it!=parts.end(); ++it) {
		RapidParticle* part = (*it);
		getMaxAltHypothesisMassShifts(part,deltaDown,deltaUp);
	}
}

void RapidParticleData::getMaxAltHypothesisMassShifts(RapidParticle* part, double& deltaDown, double& deltaUp) {
	double deltaMass = part->deltaMass();

	if(deltaMass<deltaDown) deltaDown=deltaMass;
	else if(deltaMass>deltaUp) deltaUp=deltaMass;

	//now iterate over all daughters
	RapidParticle* daug = part->daughter(0);
	while(daug) {
		getMaxAltHypothesisMassShifts(daug,deltaDown,deltaUp);
		daug = daug->next();
	}
}
