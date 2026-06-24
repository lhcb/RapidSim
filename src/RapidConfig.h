#ifndef RAPIDCONFIGPARSER_H
#define RAPIDCONFIGPARSER_H

#include <map>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "RapidAcceptance.h"
#include "RapidParam.h"

class RapidCut;
class RapidDecay;
class RapidExternalGenerator;
class RapidHistWriter;
class RapidMomentumSmear;
class RapidIPSmear;
class RapidVtxSmear;
class RapidParam;
class RapidParticle;
class RapidPID;

class RapidConfig {
	public:
		RapidConfig()
			: fileName_(""), accRejHisto_(0),
			  accRejParameterX_(0), accRejParameterY_(0),
			  acceptanceType_(RapidAcceptance::ANY),
			  detectorGeometry_(RapidAcceptance::FOURPI),
			  ppEnergy_(8.), motherFlavour_("b"),
			  ptHisto_(0), etaHisto_(0), pvHisto_(0), ptMin_(-999.), ptMax_(-999.), etaMin_(-999.), etaMax_(-999.),
			  maxgen_(1000), decay_(0), acceptance_(0), writer_(0), external_(0), usePhotos_(false)
		{}

		~RapidConfig();

		bool load(TString fileName);

		RapidDecay* getDecay();
		RapidAcceptance* getAcceptance();
		RapidHistWriter* getWriter(bool saveTree=false);

	private:
		bool loadDecay();
		bool loadConfig();
		void writeConfig();

		bool configParticle(unsigned int part, TString command, TString value);
		bool configGlobal(TString command, TString value);
		bool loadRange(TString name, TString str, double& min, double& max);
		RapidParam* loadParam(TString paramStr);
		RapidCut* loadCut(TString cutStr);

		RapidParam* findParam(TString name);

		void setSmearing(unsigned int particle, TString category);
		bool loadSmearing(TString category);
		bool loadPID(TString category);

		bool loadAcceptRejectHist(TString histFile, TString histName, RapidParam* paramX, RapidParam* paramY);
		bool loadParentKinematics();
		bool loadPVntracks();

		void setupDefaultParams();

		bool check1D(TH1* hist) { return (dynamic_cast<TH1F*>(hist) || dynamic_cast<TH1D*>(hist)); }
		bool check2D(TH1* hist) { return (dynamic_cast<TH2F*>(hist) || dynamic_cast<TH2D*>(hist)); }

		TH1* reduceHistogram(TH1* histo, double min, double max);

		TString fileName_;
    std::string outputDir_;

		std::vector<RapidParticle*> parts_;

		//particle specific parameters
		std::vector<RapidParam*> params_;

		//default parameter sets
		std::vector<RapidParam*> paramsStable_;
		std::vector<RapidParam*> paramsDecaying_;
		std::vector<RapidParam*> paramsTwoBody_;
		std::vector<RapidParam*> paramsThreeBody_;

		//parameters to load
		std::string paramStrStable_;
		std::string paramStrDecaying_;
		std::string paramStrTwoBody_;
		std::string paramStrThreeBody_;

		std::vector<RapidCut*> cuts_;

		//mometum smearing lookup for each smearing category
		std::map<TString, RapidMomentumSmear*> momSmearCategories_;
		//IP smearing lookup for each smearing category
		std::map<TString, RapidIPSmear*> ipSmearCategories_;
		//Vtx smearing lookup for each smearing category, placeholder for now
		//std::map<TString, RapidVtxSmear*> vtxSmearCategories_;

		// std::vector<std::map<TString, double>> definitions_;
		std::vector<TString> definitions_;
		//accept reject hist to sculpt kinematics
		TH1* accRejHisto_;
		RapidParam* accRejParameterX_;
		RapidParam* accRejParameterY_;

		// PID histogram file
		bool pidLoaded_;
		std::map<RapidParam::ParamType, RapidPID*> pidHists_;

		//type of geometric acceptance to apply
		RapidAcceptance::AcceptanceType acceptanceType_;

		//detector geometry to use
		RapidAcceptance::DetectorType detectorGeometry_;

		//parameters to determine parent kinematics
		double ppEnergy_;
		TString motherFlavour_;

		//parent kinematic distributions
		TH1* ptHisto_;
		TH1* etaHisto_;

		// PVNTRACKS
		TH1* pvHisto_;

		//ranges to take from histograms
		double ptMin_;
		double ptMax_;
		double etaMin_;
		double etaMax_;

		//max attempts to generate
		double maxgen_;

		RapidDecay* decay_;
		RapidAcceptance* acceptance_;
		RapidHistWriter* writer_;
		RapidExternalGenerator* external_;

		//flag to track whether an external EvtGen generator should use PHOTOS or not
		bool usePhotos_;
};

#endif
