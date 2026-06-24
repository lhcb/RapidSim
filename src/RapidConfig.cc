#include "RapidConfig.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <queue>

#include "TFile.h"
#include "TRandom.h"
#include "TSystem.h"

#include "RapidAcceptance.h"
#include "RapidAcceptanceLHCb.h"
#include "RapidCut.h"
#include "RapidDecay.h"
#include "RapidExternalEvtGen.h"
#include "RapidHistWriter.h"
#include "RapidIPSmearGauss.h"
#include "RapidMomentumSmearGauss.h"
#include "RapidMomentumSmearEnergyGauss.h"
#include "RapidMomentumSmearGaussPtEtaDep.h"
#include "RapidMomentumSmearHisto.h"
#include "RapidParam.h"
#include "RapidParticle.h"
#include "RapidParticleData.h"
#include "RapidPID.h"

RapidConfig::~RapidConfig() {
	std::map<TString, RapidMomentumSmear*>::iterator itr = momSmearCategories_.begin();
	while (itr != momSmearCategories_.end()) {
		delete itr->second;
		momSmearCategories_.erase(itr++);
	}
	std::map<TString, RapidIPSmear*>::iterator itr2 = ipSmearCategories_.begin();
	while (itr2 != ipSmearCategories_.end()) {
		delete itr2->second;
		ipSmearCategories_.erase(itr2++);
	}
	std::map<RapidParam::ParamType, RapidPID*>::iterator itr3 = pidHists_.begin();
	while (itr3 != pidHists_.end()) {
		delete itr3->second;
		pidHists_.erase(itr3++);
	}
	while(!parts_.empty()) {
		delete parts_[parts_.size()-1];
		parts_.pop_back();
	}

	while(!params_.empty()) {
		delete params_[params_.size()-1];
		params_.pop_back();
	}

	while(!paramsStable_.empty()) {
		delete paramsStable_[paramsStable_.size()-1];
		paramsStable_.pop_back();
	}

	while(!paramsDecaying_.empty()) {
		delete paramsDecaying_[paramsDecaying_.size()-1];
		paramsDecaying_.pop_back();
	}

	while(!paramsTwoBody_.empty()) {
		delete paramsTwoBody_[paramsTwoBody_.size()-1];
		paramsTwoBody_.pop_back();
	}

	while(!paramsThreeBody_.empty()) {
		delete paramsThreeBody_[paramsThreeBody_.size()-1];
		paramsThreeBody_.pop_back();
	}

	while(!cuts_.empty()) {
		delete cuts_[cuts_.size()-1];
		cuts_.pop_back();
	}

	if(accRejHisto_) delete accRejHisto_;

	if(acceptance_) delete acceptance_;
	if(decay_) delete decay_;
	if(writer_) delete writer_;
	if(external_) delete external_;
}

bool RapidConfig::load(TString fileName) {
	fileName_ = fileName;

	if(!loadDecay()) return false;
	if(!loadConfig()) return false;

	//automatically add the corrected mass variable if we have any invisible particles
	for(unsigned int i=0; i<parts_.size(); ++i) {
		if(parts_[i]->invisible()) {
			std::cout << "INFO in RapidConfig::load : invisible daughter found." << std::endl
				  << "                            Adding corrected mass variable for parent particle." << std::endl;
			std::vector<int> mother;
			mother.push_back(0);
			RapidParam* param = new RapidParam("MCorr", RapidParam::MCORR, parts_[0], false);
			params_.push_back(param);
			break;
		}
	}

	//if we're using an external generator then set it up
	if(external_) {
		std::cout << "INFO in RapidConfig::load : initialising external generator." << std::endl;
		//if the external generator is EvtGen we need to write the DEC file
		if(dynamic_cast<RapidExternalEvtGen*>(external_)) {
			dynamic_cast<RapidExternalEvtGen*>(external_)->writeDecFile(fileName_,parts_,usePhotos_);
		}
		external_->setup();
	}

	return true;
}

RapidDecay* RapidConfig::getDecay() {
	if(!decay_) {
		//check we have a particle to decay
		if(parts_.empty()) return 0;

		decay_ = new RapidDecay(parts_);

		//check decay
		if(!decay_->checkDecay()) {
			return 0;
		}

		//check that we have an FONLL model for the parent kinematics
		if(!loadParentKinematics()) {
			return 0;
		}
		decay_->setParentKinematics(ptHisto_,etaHisto_);

		if(!loadPVntracks()) {
			return 0;
		}
		decay_->setPVntracks(pvHisto_);


		decay_->setMaxGen(maxgen_);

		//load any PDF to generate
		if(accRejHisto_ && accRejParameterX_) {
			if(accRejParameterY_) {
				decay_->setAcceptRejectHist(accRejHisto_,accRejParameterX_,accRejParameterY_);
			}
			else {
				decay_->setAcceptRejectHist(accRejHisto_,accRejParameterX_);
			}
		}
		if(external_) {
			decay_->setExternal(external_);
		}
	}

	return decay_;
}

RapidAcceptance* RapidConfig::getAcceptance() {
	if(!acceptance_) {
		switch(detectorGeometry_) {
			case RapidAcceptance::LHCB:
				acceptance_ = new RapidAcceptanceLHCb(acceptanceType_, parts_, cuts_);
				break;
			case RapidAcceptance::FOURPI:
			default:
				acceptance_ = new RapidAcceptance(acceptanceType_, parts_, cuts_);
		}
	}
	return acceptance_;
}

RapidHistWriter* RapidConfig::getWriter(bool saveTree) {
	setupDefaultParams();

	if(!writer_) {
		//strip away path for name of histogram/tuple files - save in PWD
		TString histFileName(fileName_( fileName_.Last('/')+1, fileName_.Length()));
		if(!outputDir_.empty()) histFileName.Prepend((outputDir_+"/").data());
		writer_ = new RapidHistWriter(parts_, params_, paramsStable_, paramsDecaying_, paramsTwoBody_, paramsThreeBody_, histFileName, saveTree);
	}

	return writer_;
}

bool RapidConfig::loadDecay() {
	std::cout << "INFO in RapidConfig::loadDecay : loading decay descriptor from file: " << fileName_ << ".decay" << std::endl;
	TString decayStr;
	std::queue<TString> decays;
	std::queue<RapidParticle*> mothers;

	std::ifstream fin;
	fin.open(fileName_+".decay");
	if(!fin.good()) {
		std::cout << "ERROR in RapidConfig::loadDecay : file " << fileName_ << ".decay not found." << std::endl;
		return false;
	}
	decayStr.ReadLine(fin);
	fin.close();

	std::cout << "INFO in RapidConfig::loadDecay : Decay descriptor is:" << std::endl
		  << "                                 " << decayStr << std::endl;

	decays.push(decayStr);

	RapidParticleData* particleData = RapidParticleData::getInstance();

	while(!decays.empty()) {

		decayStr = decays.front();
		decays.pop();

		//first strip out any subdecays and add them to the queue
		while(decayStr.First('{')!=-1) {
			int start = decayStr.Index('{');
			int end = decayStr.Index('}',start);
			if(end < 0) {
				std::cout << "ERROR in RapidConfig::loadDecay : malformed decay descriptor." << std::endl
					  << "                                  Mismatched brackets in:" << decayStr << std::endl;
				return false;
			}

			//move sub decay into its own string
			TString subDecay = decayStr(start+1,end-start-1);
			subDecay = subDecay.Strip(TString::kBoth);

			//if the sub decay has subdecays of its own then the first '}' won't be the correct one
			while(subDecay.CountChar('{') != subDecay.CountChar('}')) {
				end = decayStr.Index("}",end+1);
				subDecay = decayStr(start+1,end-start-1);
				subDecay = subDecay.Strip(TString::kBoth);
			}
			//add to the queue and label the mother so we know which particle to decay
			decays.push(subDecay);
			std::cout << "INFO in RapidConfig::loadDecay : Found sub-decay:" << std::endl
				  << "                                 " << subDecay << std::endl;
			decayStr.Replace(start,end-start+1,"^"+subDecay(0,subDecay.First(" ")));
		}

		//the decaying particle is the first particle remaining in the mothers list
		RapidParticle* theMother(0);
		if(!mothers.empty()) {
			theMother = mothers.front();
			mothers.pop();
		}

		//now get info from decay string
		TString token;
		int from(0);

		//first is the parent
		decayStr.Tokenize(token, from, " ");
		//only need to add this if it is the top particle
		if(parts_.empty()) {
			theMother = particleData->makeParticle(token, 0);
			parts_.push_back(theMother);

			//get flavour of mother for FONLL
			if(theMother->hasBeauty()) {
				motherFlavour_ = "b";
				std::cout << "INFO in RapidConfig::loadDecay : Parent has beauty." << std::endl;
				std::cout << "                                 setting b-quark kinematics." << std::endl;
			} else if(theMother->hasCharm()) {
				motherFlavour_ = "c";
				std::cout << "INFO in RapidConfig::loadDecay : Parent has charm." << std::endl;
				std::cout << "                                 setting c-quark kinematics." << std::endl;
			} else {
				std::cout << "WARNING in RapidConfig::loadDecay : Parent has neither beauty nor charm." << std::endl;
				std::cout << "                                    defaulting to b-quark kinematics." << std::endl;
				motherFlavour_ = "b";
			}
		}

		//second should be ->
		decayStr.Tokenize(token, from, " ");
		while (decayStr.Tokenize(token, from, " ")) {
			bool stable = true;
			if(token[0] == '^') {
				token = token.Strip(TString::kBoth,'^');
				stable = false;//flags it to edit later
			}
			RapidParticle* part = particleData->makeParticle(token, theMother);
			if(!stable) mothers.push(part);
			parts_.push_back(part);
			theMother->addDaughter(part);

		}

	}

	return true;
}

bool RapidConfig::loadConfig() {
	std::cout << "INFO in RapidConfig::loadConfig : attempting to load configuration from file: " << fileName_+".config" << std::endl;

	gRandom->TRandom::SetSeed(0.);

	std::ifstream fin;
	fin.open(fileName_+".config", std::ifstream::in);
	if( ! fin.good()) {
		std::cout << "INFO in RapidConfig::loadConfig : failed to load configuration. Will write default config file to: " << fileName_+".config" << std::endl;
		fin.close();
		writeConfig();
		fin.open(fileName_+".config", std::ifstream::in);
	}

	TString buffer;
	unsigned int currentPart(parts_.size());
	while(fin.good()) {
		buffer.ReadLine(fin);
		switch(buffer[0]) {
			case '@': //particle config
				buffer.Remove(TString::kBoth,'@');
				currentPart = buffer.Atoi();
				break;
			case '!': //global config
				currentPart = parts_.size();
				break;
			case '#': //comment
				continue;
			default: //continue particle or global config
				int colon = buffer.Index(":");
				TString command = buffer(0,colon);
				command = command.Strip(TString::kBoth);
				TString value = buffer(colon+1, buffer.Length()-colon-1);
				value = value.Strip(TString::kBoth);

				if(currentPart==parts_.size()) {
					if(!configGlobal(command, value)) {
						fin.close();
						return false;
					}
				}
				else {
					if(!configParticle(currentPart, command, value)) {
						fin.close();
						return false;
					}
				}
				break;
		}
	}
	std::cout << "INFO in RapidConfig::loadConfig : finished loading configuration." << std::endl;
	fin.close();

	return true;
}

void RapidConfig::writeConfig() {
	std::ofstream fout;
	fout.open(fileName_+".config", std::ofstream::out);

	fout << "geometry : LHCb\n";
	fout << "paramsDecaying : M, P, PT\n";
	fout << "paramsStable : P, PT\n";

	std::vector<RapidParticle*>::iterator it = parts_.begin();
	for( ; it!=parts_.end(); ++it) {
		if((*it)->nDaughters() > 2) {
			fout << "paramsTwoBody : M2\n";
			break;
		}
	}

	it = parts_.begin();
	for( ; it!=parts_.end(); ++it) {
		if((*it)->nDaughters() > 3) {
			fout << "paramsThreeBody : M2\n";
			break;
		}
	}

	for(unsigned int i=0; i<parts_.size(); ++i) {
		RapidParticle* part = parts_[i];

		fout << "@" << i << "\n";
		fout << "\tname : " << part->name() << "\n";
		if(part->nDaughters()==0) {
			if(TMath::Abs(part->id()) == 11) {
				fout << "\tsmear : LHCbElectron\n";
			} else if(TMath::Abs(part->id()) == 22) {
				fout << "\tsmear : LHCbPhoton\n";
			} else {
				fout << "\tsmear : LHCbGeneric\n";
			}
			if(TMath::Abs(part->id()) == 12 ||
					TMath::Abs(part->id()) == 14 ||
					TMath::Abs(part->id()) == 16 ) {
				fout << "\tinvisible : true\n";
			}
		}

	}

	fout.close();
}

bool RapidConfig::configParticle(unsigned int part, TString command, TString value) {
	if(part > parts_.size()) {
		std::cout << "ERROR in RapidConfig::configParticle : no particle at index " << part << std::endl
			  << "                                       fix your configuration file." << std::endl;
		return false;
	}

	if(command=="name") {
		parts_[part]->setName(value);
	} else if(command=="smear") {
		setSmearing(part, value);
	} else if(command=="invisible") {
		if(value=="false") {
			parts_[part]->setInvisible(false);
		} else if(value=="true") {
			parts_[part]->setInvisible(true);
		}
	} else if(command=="altMass") {
		RapidParticleData* rpd = RapidParticleData::getInstance();

		int from(0);
		TString buffer;

		while(value.Tokenize(buffer,from," ")) {
			int altId = rpd->pdgCode(buffer);
			if(altId==0) {
				std::cout << "WARNING in RapidConfig::configParticle : unknown particle type, " << buffer << ", in alternative mass hypothesis for " << parts_[part]->name() << std::endl;
				std::cout << "                                         this hypothesis will be ignored." << std::endl;
				continue;
			}
			TString altName = rpd->getSanitisedName(altId);
			double altMass = rpd->getMass(altId);

			parts_[part]->addMassHypothesis(altName, altMass);
			std::cout << "INFO in RapidConfig::configParticle : added alternative mass hypothesis of " << buffer << " for particle " << parts_[part]->name() << std::endl;
		}
	} else if(command=="evtGenModel") {
		parts_[part]->setEvtGenDecayModel(value);
		std::cout << "INFO in RapidConfig::configParticle : set EvtGen decay model for particle " << parts_[part]->name() << std::endl
			  << "                                    : " << value << std::endl;
	} else if(command=="Define"){
		std::cout<<"Define: "<<value<<std::endl;
		parts_[part]->setEvtGenDefinitions(value);
	}

	return true;
}

bool RapidConfig::configGlobal(TString command, TString value) {
	if(command=="seed") {
		int seed = value.Atoi();
		gRandom->TRandom::SetSeed(seed);
		std::cout << "INFO in RapidConfig::configGlobal : setting seed for random number generation to " << seed << "." << std::endl
			  << "                                    seed is now " << gRandom->TRandom::GetSeed() << "." << std::endl;
	} else if(command=="acceptance") {
		std::cout << "INFO in RapidConfig::configGlobal : setting acceptance type to " << value << "." << std::endl;
		acceptanceType_ = RapidAcceptance::typeFromString(value);
	} else if(command=="geometry") {
		std::cout << "INFO in RapidConfig::configGlobal : setting detector geometry type to " << value << "." << std::endl;
		detectorGeometry_ = RapidAcceptance::detectorFromString(value);
	} else if(command=="ptRange") {
		std::cout << "INFO in RapidConfig::configGlobal : setting pT range to " << value << "." << std::endl;
		if(!loadRange("ptRange",value,ptMin_,ptMax_)) {
			std::cout << "ERROR in RapidConfig::configGlobal : failed to load pT range." << std::endl
				  << "                                     fix your configuration file." << std::endl;
			return false;
		}
	} else if(command=="etaRange") {
		std::cout << "INFO in RapidConfig::configGlobal : setting eta range to " << value << "." << std::endl;
		if(!loadRange("etaRange",value,etaMin_,etaMax_)) {
			std::cout << "ERROR in RapidConfig::configGlobal : failed to load eta range." << std::endl
				  << "                                     fix your configuration file." << std::endl;
			return false;
		}
	} else if(command=="energy") {
		ppEnergy_ = value.Atof();
		std::cout << "INFO in RapidConfig::configGlobal : pp CoM energy set to be " << ppEnergy_ << " TeV." << std::endl;
	} else if(command=="parent") {
		motherFlavour_ = value;
		std::cout << "INFO in RapidConfig::configGlobal : parent flavour forced to be " << motherFlavour_ << "." << std::endl;
	} else if(command=="minWidth") {
		RapidParticleData::getInstance()->setNarrowWidth(value.Atof());
		std::cout << "INFO in RapidConfig::configGlobal : minimum resonance width to be generated set to " << value.Atof() << " GeV." << std::endl;
	} else if(command=="maxAttempts") {
		maxgen_ = value.Atof();
		std::cout << "INFO in RapidConfig::configGlobal : maximum number of attempts to generated an event set to " << maxgen_ << "." << std::endl;
	} else if(command=="paramsStable") {
		paramStrStable_ = value;
		std::cout << "INFO in RapidConfig::configGlobal : will use the following parameters for all stable particles:" << std::endl
			  << "                                    " << paramStrStable_ << "." << std::endl;
	} else if(command=="paramsDecaying") {
		paramStrDecaying_ = value;
		std::cout << "INFO in RapidConfig::configGlobal : will use the following parameters for all decaying particles:" << std::endl
			  << "                                    " << paramStrDecaying_ << "." << std::endl;
	} else if(command=="paramsTwoBody") {
		paramStrTwoBody_ = value;
		std::cout << "INFO in RapidConfig::configGlobal : will use the following parameters for all two-body daughter combinations:" << std::endl
			  << "                                    " << paramStrTwoBody_ << "." << std::endl;
	} else if(command=="paramsThreeBody") {
		paramStrThreeBody_ = value;
		std::cout << "INFO in RapidConfig::configGlobal : will use the following parameters for all three-body daughter combinations:" << std::endl
			  << "                                    " << paramStrThreeBody_ << "." << std::endl;
	} else if(command=="param") {
		RapidParam* param = loadParam(value);
		if(!param) {
			std::cout << "ERROR in RapidConfig::configGlobal : failed to load parameter." << std::endl
				  << "                                     fix your configuration file." << std::endl;
			return false;
		} else {
			std::cout << "INFO in RapidConfig::configGlobal : adding parameter " << param->name() << std::endl;
			params_.push_back(param);
		}
	} else if(command=="cut") {
		RapidCut* cut = loadCut(value);
		if(!cut) {
			std::cout << "ERROR in RapidConfig::configGlobal : failed to load cut." << std::endl
				  << "                                     fix your configuration file." << std::endl;
			return false;
		} else {
			std::cout << "INFO in RapidConfig::configGlobal : adding cut " << cut->name() << std::endl;
			cuts_.push_back(cut);
		}
	} else if(command=="shape") {
		int from(0);
		TString histName, histFile;

		value.Tokenize(histFile,from," ");
		value.Tokenize(histName,from," ");

		histFile = histFile.Strip(TString::kBoth);
		histName = histName.Strip(TString::kBoth);

		TString paramNameX, paramNameY;
		RapidParam *paramX(0), *paramY(0);

		if(value.Tokenize(paramNameX,from," ")) {
			paramNameX = paramNameX.Strip(TString::kBoth);
			paramX = findParam(paramNameX);
			if(!paramX) {
				std::cout << "ERROR in RapidConfig::configGlobal : failed to setup shape PDF - unknown parameter." << std::endl;
				return false;
			}
		}
		if(value.Tokenize(paramNameY,from," ")) {
			paramNameY = paramNameY.Strip(TString::kBoth);
			paramY = findParam(paramNameY);
			if(!paramY) {
				std::cout << "ERROR in RapidConfig::configGlobal : failed to setup shape PDF - unknown parameter." << std::endl;
				return false;
			}
		}

		if(!loadAcceptRejectHist(histFile, histName, paramX, paramY)) return false;
	} else if(command=="useEvtGen") {
		if(!getenv("EVTGEN_ROOT")) {
			std::cout << "ERROR in RapidConfig::configGlobal : EVTGEN_ROOT environment variable must be set to use external EvtGen generator." << std::endl;
			return false;
		}
		external_ = new RapidExternalEvtGen();
		std::cout << "INFO in RapidConfig::configGlobal : will use external EvtGen generator to decay particles." << std::endl;
	} else if(command=="evtGenUsePHOTOS") {
		usePhotos_ = true;
		std::cout << "INFO in RapidConfig::configGlobal : external EvtGen generator will use PHOTOS." << std::endl;
	} else if (command=="pid") {
		std::cout << "INFO in RapidConfig::configGlobal : setting pid type to " << value << "." << std::endl;
		int from(0);
		TString histFile;
		value.Tokenize(histFile,from," ");
		histFile = histFile.Strip(TString::kBoth);

		pidLoaded_ = loadPID(histFile);
		if(!pidLoaded_) return false;
	}
	else if (command=="outputDirectory") {
		outputDir_ = gSystem->ExpandPathName(value.Data());
		std::cout << "INFO in RapidConfig::configGlobal : setting output directory to " << outputDir_ << "." << std::endl;
		auto vptr = static_cast<const char*>(gSystem->OpenDirectory(outputDir_.data()));
		if(vptr == nullptr){
			std::cout << "INFO in RapidConfig::configGlobal : Creating new directory: " << outputDir_ << "." << std::endl;
			auto status = gSystem->mkdir(outputDir_.data(),true);
			if(status == -1)
				throw std::runtime_error("RapidConfig::configGlobal : Directory "+outputDir_+" can not be created. Please check path");
		}
		delete vptr;
	}
	return true;
}

bool RapidConfig::loadRange(TString name, TString str, double& min, double& max) {
	int from(0);
	TString buffer;

	if(str.Tokenize(buffer,from," ")) {
		min = buffer.Atof();
	} else {
		std::cout << "ERROR in RapidConfig::loadRange : failed to setup range \"" << name << "\" - no minimum given." << std::endl;
		return false;
	}
	if(str.Tokenize(buffer,from," ")) {
		max = buffer.Atof();
	} else {
		std::cout << "ERROR in RapidConfig::loadRange : failed to setup range \"" << name << "\" - no maximum given." << std::endl;
		return false;
	}
	return true;
}

RapidParam* RapidConfig::loadParam(TString paramStr) {
	int from(0);
	TString name, buffer;
	RapidParam::ParamType type;
	bool truth = false;
	std::vector<RapidParticle*> partlist;
	paramStr.Tokenize(name,from," ");
	paramStr.Tokenize(buffer,from," ");
	type = RapidParam::typeFromString(buffer);

	while(paramStr.Tokenize(buffer,from," ")) {
		if(buffer.IsDec()) {
			unsigned int index = buffer.Atoi();
			if(index>=parts_.size()) {
				std::cout << "WARNING in RapidConfig::loadParam : particle " << index << "is out of range." << std::endl
					  << "                                    parameter " << name << " has not been added." << std::endl;
				return 0;
			} else {
				partlist.push_back(parts_[index]);
			}
		} else if(buffer=="TRUE") {
			truth = true;
		} else {
			std::cout << "WARNING in RapidConfig::loadParam : unknown argument " << buffer << "in parameter " << name << "." << std::endl
				  << "                                    argument will be ignored." << std::endl;
		}
	}

	if((type == RapidParam::THETA || type == RapidParam::COSTHETA) && partlist.size()<2) {
		std::cout << "WARNING in RapidConfig::loadParam : (cos)theta parameter expects at least two particles, " << partlist.size() << "were given." << std::endl
			  << "                                    parameter " << name << " has not been added." << std::endl;
		return 0;
	}

	RapidParam * param = new RapidParam(name, type, partlist, truth);
	return param;
}

RapidCut* RapidConfig::loadCut(TString cutStr) {
	int from(0);
	TString paramName, cutType, buffer;
	double min(0.), max(0.);

	cutStr.Tokenize(paramName,from," ");
	RapidParam* param = findParam(paramName);
	if(!param) {
		std::cout << "ERROR in RapidConfig::loadCut : failed to setup cut - unknown parameter." << std::endl;
		return 0;
	}

	if(!cutStr.Tokenize(cutType,from," ")) {
		std::cout << "ERROR in RapidConfig::loadCut : failed to setup cut - no cut type given." << std::endl;
		return 0;
	}

	if(cutType!="max") {
		if(cutStr.Tokenize(buffer,from," ")) {
			min = buffer.Atof();
		} else {
			std::cout << "ERROR in RapidConfig::loadCut : failed to setup cut - no minimum given." << std::endl;
			return 0;
		}
	}
	if(cutType!="min") {
		if(cutStr.Tokenize(buffer,from," ")) {
			max = buffer.Atof();
		} else {
			std::cout << "ERROR in RapidConfig::loadCut : failed to setup cut - no maximum given." << std::endl;
			return 0;
		}
	}

	RapidCut* cut(0);

	if(cutType=="min") {
		cut = new RapidCut(param,min,RapidCut::NOLIMIT);
	} else if(cutType=="max") {
		cut = new RapidCut(param,RapidCut::NOLIMIT,max);
	} else if(cutType=="range") {
		cut = new RapidCut(param,min,max);
	} else if(cutType=="veto") {
		cut = new RapidCut(param,min,max,true);
	} else {
		std::cout << "ERROR in RapidConfig::loadCut : failed to setup cut - unknown cut type " << cutType << "." << std::endl;
	}

	return cut;

}

RapidParam* RapidConfig::findParam(TString name) {
	std::vector<RapidParam*>::iterator it = params_.begin();

	for( ; it!=params_.end(); ++it) {
		RapidParam* param = *it;
		if(param->name() == name) return param;
	}

	std::cout << "WARNING in RapidConfig::findParam : parameter " << name << " not found." << std::endl;
	return 0;
}

bool RapidConfig::loadSmearing(TString category) {
	TString path;
	std::ifstream fin;
	bool found(false);

	path = getenv("RAPIDSIM_CONFIG");

	if(path!="") {
		fin.open(path+"/config/smear/"+category, std::ifstream::in);
		if( fin.good()) {
			std::cout << "INFO in RapidConfig::loadSmearing : found smearing category " << category << " in RAPIDSIM_CONFIG." << std::endl;
			std::cout << "                                    this version will be used." << std::endl;
			found = true;
		} else {
			std::cout << "INFO in RapidConfig::loadSmearing : smearing category " << category << " not found in RAPIDSIM_CONFIG." << std::endl;
			std::cout << "                                    checking RAPIDSIM_ROOT." << std::endl;
			fin.close();
		}
	}

	if(!found) {
		path = getenv("RAPIDSIM_ROOT");
		fin.open(path+"/config/smear/"+category, std::ifstream::in);
		if( ! fin.good()) {
			std::cout << "WARNING in RapidConfig::loadSmearing : failed to load smearing category" << category << std::endl;
			fin.close();
			return false;
		}
	}

	TString filename("");
	filename.ReadToken(fin);
	TFile* file = NULL;
	if (filename != "NULL") {
		file = TFile::Open(path+"/rootfiles/smear/"+filename);
		if(!file) {
			std::cout << "WARNING in RapidConfig::loadSmearing : failed to load root file " << filename << std::endl;
			fin.close();
			return false;
		}
	}

	TString type("");
	type.ReadToken(fin);
	if(type=="GAUSS") {
		TString histname("");
		histname.ReadToken(fin);
		TGraphErrors* graph = dynamic_cast<TGraphErrors*>(file->Get(histname));
		if(!graph) {
			std::cout << "WARNING in RapidConfig::loadSmearing : failed to load graph " << histname << std::endl;
			file->Close();
			fin.close();
			return false;
		}

		momSmearCategories_[category] = new RapidMomentumSmearGauss(graph);

	}else if(type=="GAUSSIP") {
		double intercept(0.);
		double slope(0.);
		fin >> intercept;
		fin >> slope;
		if ((intercept < 0) || (slope < 0) ) {
			std::cout << "WARNING in RapidConfig::loadSmearing : failed to load IP smearing" << std::endl;
			file->Close();
			fin.close();
			return false;
		}
		ipSmearCategories_[category] = new RapidIPSmearGauss(intercept,slope);
	}else if(type=="GAUSSPHOTON") {
		double stochastic(0.);
		double constant(0.);
		fin >> stochastic;
		fin >> constant;
		if ((stochastic < 0) || (constant < 0) ) {
			std::cout << "WARNING in RapidConfig::loadSmearing : failed to load photon smearing" << std::endl;
			file->Close();
			fin.close();
			return false;
		}
		std::cout << "INFO in RapidConfig::loadSmearing : load photon smearing for particle " << std::endl;
		momSmearCategories_[category] = new RapidMomentumSmearEnergyGauss(stochastic,constant);
	}else if(type=="GAUSSPTETA") {
		TString histname("");
		histname.ReadToken(fin);
		TH2* hist = dynamic_cast<TH2*>(file->Get(histname));
		if(!hist) {
			std::cout << "WARNING in RapidConfig::loadSmearing : failed to load histogram " << histname << std::endl;
			file->Close();
			fin.close();
			return false;
		}

		momSmearCategories_[category] = new RapidMomentumSmearGaussPtEtaDep(hist);

	} else if(type=="HISTS") {
		double threshold(0.);
		TString histname("");

		std::vector<TH1*> hists;
		std::vector<double> thresholds;

		while( true ) {
			fin >> threshold;
			histname.ReadToken(fin);
			if(fin.good()) {
				TH1* hist = dynamic_cast<TH1*>(file->Get(histname));

				if(!hist || !check1D(hist)) {
					std::cout << "WARNING in RapidConfig::loadSmearing : failed to load histogram " << histname << std::endl
						  << "                                       threshold will be ignored." << std::endl;
				}

				hists.push_back(hist);
				thresholds.push_back(threshold);

			} else {
				break;
			}
		}

		if(hists.size() == 0) {
			std::cout << "WARNING in RapidConfig::loadSmearing : failed to load any histograms for smearing category " << category << std::endl;
			file->Close();
			fin.close();
			return false;
		}

		momSmearCategories_[category] = new RapidMomentumSmearHisto(thresholds, hists);

	} else {
		std::cout << "WARNING in RapidConfig::loadSmearing : unknown smearing type. Category " << category << " not added." << std::endl;
		file->Close();
		fin.close();
		return false;
	}

	fin.close();
	return true;
}

void RapidConfig::setSmearing(unsigned int particle, TString category) {
	if(particle > parts_.size()) {
		std::cout << "WARNING in RapidConfig::setSmearing : particle " << particle << " does not exist - smearing functions not set." << std::endl;
		return;
	}
	if(parts_[particle]->nDaughters() != 0) {
		std::cout << "WARNING in RapidConfig::setSmearing : particle " << particle << " is composite - smearing functions not set." << std::endl;
		return;
	}

	std::cout << "INFO in RapidConfig::setSmearing : setting smearing functions for particle " << particle << " (category: " << category << ")" << std::endl;

	bool loadedmomsmear = true;
	if(!momSmearCategories_.count(category)) {
		if(!loadSmearing(category)) {
			std::cout << "WARNING in RapidConfig::setSmearing : failed to load momentum smearing category " << category << "." << std::endl
				  << "                                      smearing functions not set for particle " << particle << "." << std::endl;
			loadedmomsmear = false;
		} else if(!momSmearCategories_.count(category)) {
			loadedmomsmear = false;
		}
	}
	if (loadedmomsmear) {
		parts_[particle]->setSmearing(momSmearCategories_[category]);
	}

	bool loadedipsmear = true;
	if(!ipSmearCategories_.count(category)) {
		if(!loadSmearing(category)) {
			std::cout << "WARNING in RapidConfig::setSmearing : failed to load IP smearing category " << category << "." << std::endl
				  << "                                      smearing functions not set for particle " << particle << "." << std::endl;
			loadedipsmear = false;
		} else if(!ipSmearCategories_.count(category)) {
			loadedipsmear = false;
		}
	}
	if (loadedipsmear) {
		parts_[particle]->setSmearing(ipSmearCategories_[category]);
	}

	return;
}

bool RapidConfig::loadPID(TString category) {
	TString path;
	std::ifstream fin;
	bool found(false);

	path = getenv("RAPIDSIM_CONFIG");
	std::cout << "INFO in RapidConfig::loadPID " << std::endl;

	if(path!="") {
		fin.open(path+"/config/pid/"+category, std::ifstream::in);
		if( fin.good()) {
			std::cout << "INFO in RapidConfig::loadPID : found pid category " << category << " in RAPIDSIM_CONFIG." << std::endl;
			std::cout << "                                    this version will be used." << std::endl;
			found = true;
		} else {
			std::cout << "INFO in RapidConfig::loadPID : pid category " << category << " not found in RAPIDSIM_CONFIG." << std::endl;
			std::cout << "                                    checking RAPIDSIM_ROOT." << std::endl;
			fin.close();
		}
	}

	if(!found) {
		path = getenv("RAPIDSIM_ROOT");
		fin.open(path+"/config/pid/"+category, std::ifstream::in);
		if( ! fin.good()) {
			std::cout << "WARNING in RapidConfig::loadPID : failed to load pid category" << category << std::endl;
			fin.close();
			return false;
		}
	}

	TString line("");
	TString buffer("");
	while (fin.good()) {
		line.ReadLine(fin);
		int from(0);
		bool fileLoaded(false);
		bool idLoaded(false);
		unsigned int id(0);
		TFile* file = NULL;
		while (line.Tokenize(buffer, from)) {
			if (buffer.Contains(".root") && !fileLoaded) {
				std::cout << "INFO in RapidConfig::loadPID : loading root file " << buffer << std::endl;
				if (buffer.BeginsWith("/") || buffer.BeginsWith(".")) file = TFile::Open(buffer);
				else file = TFile::Open(path+"/rootfiles/pid/"+buffer);
				fileLoaded = true;
				if(!file) {
					std::cout << "WARNING in RapidConfig::loadPID : failed to load root file " << buffer << std::endl;
					fin.close();
					return false;
				}
				continue;
			}
			if ( fileLoaded && !idLoaded) {
				id = buffer.Atoi();
				idLoaded = true;
				continue;
			}
			if ( fileLoaded && idLoaded && buffer.Contains("Prob")) {
				std::cout << "INFO in RapidConfig::loadPID : loading histogram " << buffer << std::endl;
				TH3D * hist = dynamic_cast<TH3D*>(file->Get(buffer));
				if(!hist) {
					std::cout << "WARNING in RapidConfig::loadPID : failed to load histogram " << buffer << std::endl;
					continue;
				}
				if ( hist->GetMinimum() < 0 ) {
					for (int i = 0; i < hist->GetNcells(); ++i) {
						if (hist->GetBinContent(i) < 0.) hist->SetBinContent(i, 0.);
					}
				}
				RapidParam::ParamType type = RapidParam::typeFromString(buffer);
				if(pidHists_.find(type)==pidHists_.end()) {
					pidHists_[type] = new RapidPID(buffer);
				}
				pidHists_[type]->addPID(id, hist);
			}
			if(pidHists_.empty()) {
				std::cout << "WARNING in RapidConfig::loadPID : failed to load any histograms for PID category " << category << std::endl;
				file->Close();
				fin.close();
				return false;
			}
		}
	}
	fin.close();
	return true;
}

bool RapidConfig::loadAcceptRejectHist(TString histFile, TString histName, RapidParam* paramX, RapidParam* paramY) {
	if(accRejHisto_) {
		std::cout << "WARNING in RapidConfig::loadAcceptRejectHist : accept/reject histogram already set." << std::endl
			  << "                                               original histogram will be used." << std::endl;
		return false;
	}

	TFile* file = TFile::Open(histFile);
	if(!file) {
		std::cout << "WARNING in RapidConfig::loadAcceptRejectHist : could not open file " << histFile << "." << std::endl
			  << "                                               path should be absolute or relative to '" << getenv("PWD") << "'." << std::endl
			  << "                                               accept/reject histogram not set." << std::endl;
		return false;
	}
	TH1* hist = dynamic_cast<TH1*>(file->Get(histName));
	if(!hist) {
		std::cout << "WARNING in RapidConfig::loadAcceptRejectHist : could not load histogram " << histName << "." << std::endl
			  << "                                               accept/reject histogram not set." << std::endl;
		return false;
	}

	if(paramX) {
		if(paramY) {
			//2 parameters - check we have a 2D histogram
			if(!check2D(hist)) {
				std::cout << "ERROR in RapidConfig::loadAcceptRejectHist : two parameters provided for shape PDF but histogram is neither TH2F nor TH2D." << std::endl;
				return false;
			}
		} else {
			//1 parameter - check we have a 1D histogram
			if(!check1D(hist)) {
				std::cout << "ERROR in RapidConfig::loadAcceptRejectHist : one parameter provided for shape PDF but histogram is neither TH1F nor TH1D." << std::endl;
				return false;
			}
		}
	} else {
		std::cout << "ERROR in RapidConfig::loadAcceptRejectHist : no parameters provided for shape PDF." << std::endl;
		return false;
	}

	std::cout << "INFO in RapidConfig::loadAcceptRejectHist : setting the required kinematic distribution." << std::endl;
	accRejHisto_ = hist;
	accRejParameterX_ = paramX;
	accRejParameterY_ = paramY;

	return true;
}

bool RapidConfig::loadPVntracks() {
	TString path;
	TString fileName;
	TFile* file(0);
	bool found(false);

	path = getenv("RAPIDSIM_CONFIG");

	if(path!="") {
		fileName = path;
		fileName +="/rootfiles/smear/PVNTRACKS.root";
		file = TFile::Open(fileName);

		if(file) {
			std::cout << "INFO in RapidConfig::loadPVntracks : " << fileName << " in RAPIDSIM_CONFIG." << std::endl
				  << "                                            this version will be used." << std::endl;
			found = true;
		} else {
			std::cout << "INFO in RapidConfig::loadPVntracks : " << fileName << " not found in RAPIDSIM_CONFIG." << std::endl
				  << "                                            checking RAPIDSIM_ROOT." << std::endl;
		}
	}

	if(!found) {
		path = getenv("RAPIDSIM_ROOT");
		fileName = path;
		fileName +="/rootfiles/smear/PVNTRACKS.root";
		file = TFile::Open(fileName);

		if(!file) {
			std::cout << "ERROR in RapidConfig::loadPVntracks : " << fileName << " not found." << std::endl;
			return false;
		}
	}

	pvHisto_ = dynamic_cast<TH1*>(file->Get("h"));
	if(!pvHisto_) {
		std::cout << "ERROR in RapidConfig::loadPVntracks : PV ntracks histogram is neither TH1F nor TH1D." << std::endl;
		return false;
	}

	return true;
}


bool RapidConfig::loadParentKinematics() {
	TString path;
	TString fileName;
	TFile* file(0);
	bool found(false);

	path = getenv("RAPIDSIM_CONFIG");

	if(path!="") {
		fileName = path;
		fileName +="/rootfiles/fonll/LHC";
		fileName += motherFlavour_;
		fileName += ppEnergy_;
		fileName += ".root";
		file = TFile::Open(fileName);

		if(file) {
			std::cout << "INFO in RapidConfig::loadParentKinematics : found kinematics LHC" << motherFlavour_ << ppEnergy_ << " in RAPIDSIM_CONFIG." << std::endl
				  << "                                            this version will be used." << std::endl;
			found = true;
		} else {
			std::cout << "INFO in RapidConfig::loadParentKinematics : kinematics LHC" << motherFlavour_ << ppEnergy_ << " not found in RAPIDSIM_CONFIG." << std::endl
				  << "                                            checking RAPIDSIM_ROOT." << std::endl;
		}
	}

	if(!found) {
		path = getenv("RAPIDSIM_ROOT");
		fileName = path;
		fileName += "/rootfiles/fonll/LHC";
		fileName += motherFlavour_;
		fileName += ppEnergy_;
		fileName += ".root";
		file = TFile::Open(fileName);

		if(!file) {
			std::cout << "ERROR in RapidConfig::loadParentKinematics : unknown kinematics " << motherFlavour_ << "-quark from " << ppEnergy_ << " TeV pp collision." << std::endl
				  << "                                             file " << fileName << " not found." << std::endl;
			return false;
		}
	}

	TH1* ptHisto = dynamic_cast<TH1*>(file->Get("pT"));
	TH1* etaHisto = dynamic_cast<TH1*>(file->Get("eta"));

	if(!ptHisto || !check1D(ptHisto)) {
		std::cout << "ERROR in RapidConfig::loadParentKinematics : pT histogram is neither TH1F nor TH1D." << std::endl;
		return false;
	}

	if(!etaHisto || !check1D(etaHisto)) {
		std::cout << "ERROR in RapidConfig::loadParentKinematics : eta histogram is neither TH1F nor TH1D." << std::endl;
		return false;
	}

	if( ptMin_==-999. || ptMax_==-999. ) {
		std::cout << "INFO in RapidConfig::loadParentKinematics : pt range not defined by user." << std::endl;
		std::cout << "                                            Will take default for detector geometry." << std::endl;

		getAcceptance()->getDefaultPtRange(ptMin_,ptMax_);
	}
	if( etaMin_==-999. || etaMax_==-999. ) {
		std::cout << "INFO in RapidConfig::loadParentKinematics : pseudorapidity range not defined by user." << std::endl;
		std::cout << "                                            Will take default for detector geometry." << std::endl;

		getAcceptance()->getDefaultEtaRange(etaMin_,etaMax_);
	}

	ptHisto_ = reduceHistogram(ptHisto,ptMin_,ptMax_);
	etaHisto_ = reduceHistogram(etaHisto,etaMin_,etaMax_);

	return true;
}

void RapidConfig::setupDefaultParams() {
	int from(0);
	TString buffer;
	TString baseName;

	// we could factor some of this code out into separate functions
	// Stable particles
	while(TString(paramStrStable_).Tokenize(buffer,from," ")) {
		buffer = buffer.Strip(TString::kBoth,',');
		RapidParam::ParamType type = RapidParam::typeFromString(buffer);

		if ( buffer.Contains("ProbNN") && !pidLoaded_) loadPID("LHCbGenericPID");

		if(type==RapidParam::UNKNOWN) {
			std::cout << "WARNING in RapidConfig::setDefaultParams : Unknown parameter type " << buffer << "ignored." << std::endl;
			continue;
		} else {
			for(unsigned int i=0; i<parts_.size(); ++i) {
				RapidParticle* part = parts_[i];
				if(part->nDaughters() == 0) {
					RapidPID* pidHists=0;
					if (pidHists_.find(type)!=pidHists_.end() && pidHists_[type] ) {
						pidHists = pidHists_[type];
					}
					RapidParam* param = new RapidParam("", type, part, false, pidHists);
					if ( param->canBeSmeared() ) {
						param->name();
						paramsStable_.push_back(param);
					} else delete param;
					param = new RapidParam("", type, part, true, pidHists);
					if ( param->canBeTrue() ) {
						param->name();
						paramsStable_.push_back(param);
					} else delete param;
				}
			}
		}
	}

	from = 0;
	// Decaying particles
	while(TString(paramStrDecaying_).Tokenize(buffer,from," ")) {
		buffer = buffer.Strip(TString::kBoth,',');
		RapidParam::ParamType type = RapidParam::typeFromString(buffer);
		if(type==RapidParam::UNKNOWN) {
			std::cout << "WARNING in RapidConfig::setDefaultParams : Unknown parameter type " << buffer << "ignored." << std::endl;
			continue;
		} else {
			for(unsigned int i=0; i<parts_.size(); ++i) {
				RapidParticle* part = parts_[i];
				if(part->nDaughters() > 0) {
					RapidParam* param = new RapidParam("", type, part, false);
					if ( param->canBeSmeared() ) {
						param->name();
						paramsDecaying_.push_back(param);
					} else delete param;
					param = new RapidParam("", type, part, true);
					if ( param->canBeTrue() ) {
						param->name();
						paramsDecaying_.push_back(param);
					} else delete param;
				}
			}
		}
	}

	from = 0;
	//2-body IMs
	while(TString(paramStrTwoBody_).Tokenize(buffer,from," ")) {
		buffer = buffer.Strip(TString::kBoth,',');
		RapidParam::ParamType type = RapidParam::typeFromString(buffer);
		if(type==RapidParam::UNKNOWN) {
			std::cout << "WARNING in RapidConfig::setDefaultParams : Unknown parameter type " << buffer << "ignored." << std::endl;
			continue;
		} else {
			for(unsigned int i=0; i<parts_.size(); ++i) {
				RapidParticle* part = parts_[i];
				if(part->nDaughters() > 2) {
					RapidParticle* jDaug = part->daughter(0);

					for( ; jDaug!=0; jDaug=jDaug->next()) {
						for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
							std::vector<RapidParticle*> partlist;
							partlist.push_back(jDaug);
							partlist.push_back(kDaug);
							baseName  = jDaug->name()+"_";
							baseName += kDaug->name()+"_";
							RapidParam* param = new RapidParam("", type, partlist, false);
							if ( param->canBeSmeared() ) {
								param->name();
								paramsTwoBody_.push_back(param);
							} else delete param;
							param = new RapidParam("", type, partlist, true);
							if ( param->canBeTrue() ) {
								param->name();
								paramsTwoBody_.push_back(param);
							} else delete param;
						}
					}
				}
			}
		}
	}

	from = 0;
	//3-body IMs
	while(TString(paramStrThreeBody_).Tokenize(buffer,from," ")) {
		buffer = buffer.Strip(TString::kBoth,',');
		RapidParam::ParamType type = RapidParam::typeFromString(buffer);
		if(type==RapidParam::UNKNOWN) {
			std::cout << "WARNING in RapidConfig::setDefaultParams : Unknown parameter type " << buffer << "ignored." << std::endl;
			continue;
		} else {
			for(unsigned int i=0; i<parts_.size(); ++i) {
				RapidParticle* part = parts_[i];
				if(part->nDaughters() > 3) {
					RapidParticle* jDaug = part->daughter(0);

					for( ; jDaug!=0; jDaug=jDaug->next()) {
						for(RapidParticle* kDaug=jDaug->next(); kDaug!=0; kDaug=kDaug->next()) {
							for(RapidParticle* lDaug=kDaug->next(); lDaug!=0; lDaug=lDaug->next()) {
								std::vector<RapidParticle*> partlist;
								partlist.push_back(jDaug);
								partlist.push_back(kDaug);
								partlist.push_back(lDaug);
								baseName  = jDaug->name()+"_";
								baseName += kDaug->name()+"_";
								baseName += lDaug->name()+"_";
								RapidParam* param = new RapidParam("", type, partlist, false);
								if ( param->canBeSmeared() ) {
									param->name();
									paramsThreeBody_.push_back(param);
								} else delete param;
								param = new RapidParam("", type, partlist, true);
								if ( param->canBeTrue() ) {
									param->name();
									paramsThreeBody_.push_back(param);
								} else delete param;
							}
						}
					}
				}
			}
		}
	}
}

TH1* RapidConfig::reduceHistogram(TH1* histo, double min, double max) {
	TAxis* axis = histo->GetXaxis();

	//first check if histogram is defined over the required range
	if( min < axis->GetXmin() ) {
		std::cout << "WARNING in RapidConfig::reduceHistogram : Minimum set to " << min << " but histogram " << histo->GetName() << " has minimum " << axis->GetXmin() << std::endl;
		std::cout << "                                        : Minimum changed to " << axis->GetXmin() << std::endl;
		min = axis->GetXmin();
	}
	if( max > axis->GetXmax() ) {
		std::cout << "WARNING in RapidConfig::reduceHistogram : Maximum set to " << max << " but histogram " << histo->GetName() << " has maximum " << axis->GetXmax() << std::endl;
		std::cout << "                                        : Maximum changed to " << axis->GetXmax() << std::endl;
		max = axis->GetXmax();
	}

	//if we want the full histogram then just return it
	if(min == axis->GetXmin() && max == axis->GetXmax() ) {
		return histo;
	}

	//load the bin boundaries from the input histogram
	int nbins = axis->GetNbins();
	double* bounds = new double[nbins+1];
	axis->GetLowEdge(bounds);
	bounds[nbins] = axis->GetXmax();

	int firstBin=0;
	int lastBin=nbins;
	bool partFirstBin=false;
	bool partLastBin=false;

	//find the first bin to copy and check if it's partially removed
	for(int i=0; i<=nbins; ++i) {
		if(bounds[i]>=min) {
			firstBin=i;
			if(bounds[i]!=min) partFirstBin=true;
			break;
		}
	}
	//find the final bin to copy and check if it's partially removed
	for(int i=nbins; i>=0; --i) {
		if(bounds[i]<=max) {
			lastBin=i;
			if(bounds[i]!=max) partLastBin=true;
			break;
		}
	}

	//set up the bin boundaries for the output histogram
	int nbinsOut = lastBin - firstBin + partFirstBin + partLastBin;
	double* boundsOut = new double[nbinsOut+1];

	boundsOut[0] = min;
	boundsOut[nbinsOut] = max;

	for(int i=1; i<nbinsOut; ++i) {
		boundsOut[i] = bounds[firstBin+i-partFirstBin];
	}

	//create histogram - clone everything except the bins
	TH1* histoOut = dynamic_cast<TH1*>(histo->Clone());
	histoOut->SetBins(nbinsOut,boundsOut);

	//set bin contents - for partial bins, scale by the ratio of bin widths
	for(int i=1; i<=nbinsOut; ++i) {
		double value = histo->GetBinContent(firstBin+i-partFirstBin) * histoOut->GetBinWidth(i) / histo->GetBinWidth(firstBin+i-partFirstBin);
		histoOut->SetBinContent(i,value);
	}

	delete[] bounds;
	delete[] boundsOut;

	return histoOut;
}
