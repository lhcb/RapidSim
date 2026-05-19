#ifndef RAPIDPARTICLEDATA_H
#define RAPIDPARTICLEDATA_H

#include <map>
#include <set>

#include "TString.h"

#include "RooRealVar.h"

class RooRelBreitWigner;
class RooSubThreshold;
class RooGounarisSakurai;
class RapidParticle;

class RapidParticleData {
	public:
		enum ResLineShape {
			RelBW,
			SubTh,
			GS
		};

		static RapidParticleData* getInstance();

		void loadData(TString file);

		double getCT(int id);
		double getMass(int id);
		double getWidth(int id);
		double getSpin(int id);
		double getCharge(int id);
		TString getName(int id);
		TString getSanitisedName(int id);

		int pdgCode(TString name);

		RapidParticle* makeParticle(int id, RapidParticle* mother);
		RapidParticle* makeParticle(TString name, RapidParticle* mother);

		void setupMass(RapidParticle* part);
		void setNarrowWidth(double narrowWidth) { narrowWidth_ = narrowWidth; }

		bool checkHierarchy(const std::vector<RapidParticle*>& parts);
		bool checkHierarchy(RapidParticle* part, RapidParticle* ancestor);

		RapidParticle* findCommonAncestor(const std::vector<RapidParticle*>& parts);

		void findStableDaughters(const std::vector<RapidParticle*>& parts, std::vector<RapidParticle*>& daughters);
		void findStableDaughters(RapidParticle* parts, std::vector<RapidParticle*>& daughters);

		void findOtherDaughters(RapidParticle* ancestor, const std::vector<RapidParticle*>& parts, std::vector<RapidParticle*>& others);

		void combineCompleteAncestors(const std::vector<RapidParticle*>& parts, std::vector<RapidParticle*>& partsCombined);

		void getMaxAltHypothesisMassShifts(const std::vector<RapidParticle*>& parts, double& deltaDown, double& deltaUp);
		void getMaxAltHypothesisMassShifts(RapidParticle* parts, double& deltaDown, double& deltaUp);

	private:
		static RapidParticleData* instance_;

		RapidParticleData()
		: narrowWidth_(0.001) {}

		~RapidParticleData() {}

		//copy constructor and copy assignment operator not implemented
		RapidParticleData( const RapidParticleData& other );
		RapidParticleData& operator=( const RapidParticleData& other );

		void addEntry(int id, TString name, double mass, double width, double spin, double charge, TString lineshape, double ctau);
		TString sanitiseName(TString name);
		TString makeUniqName(TString name);

		RooRelBreitWigner* makeRelBW(RooRealVar& m, double mean, double gamma, double thespin, double m1, double m2, TString name);
		RooSubThreshold* makeSubTh(RooRealVar& m, double mean, double gamma, double thespin, double m1, double m2, double m3, double mPar, TString name);
  	RooGounarisSakurai* makeGS(RooRealVar& m, double mean, double gamma, double thespin, double m1, double m2, TString name);

		std::map<int, double> idToCT_;
		std::map<int, double> idToMass_;
		std::map<int, double> idToWidth_;
		std::map<int, double> idToSpin_;
		std::map<int, double> idToCharge_;
		std::map<int, TString> idToName_;
		std::map<int, ResLineShape> idToShape_;
		std::map<TString, int> nameToId_;
		std::map<TString, int> sanitisedNameToId_;

		//register of particle names that have been used
		std::set<TString> usedNames_;

		double narrowWidth_;
};

#endif
