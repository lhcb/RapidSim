#include "RapidExternalEvtGen.h"

#include <cstdlib>
#include <fstream>
#include <queue>

#include "TRandom.h"
#include "TSystem.h"
#include "TObjString.h"

#ifdef RAPID_EVTGEN
#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtMTRandomEngine.hh"

#include "EvtGenExternal/EvtExternalGenList.hh"
#endif

int B0Id( 511 ), B0barId( -511 );
int BsId( 531 ), BsbarId( -531 );


bool RapidExternalEvtGen::decay(std::vector<RapidParticle*>& parts) {
#ifdef RAPID_EVTGEN

	if(parts.size() < 1) {
		std::cout << "WARNING in RapidExternalEvtGen::decay : There are no particles to decay." << std::endl;
		return false;
	}

	EvtSpinDensity* spinDensity = 0;
	TString parentName = getEvtGenName(parts[0]->id());
	EvtId theId = EvtPDL::evtIdFromLundKC(parts[0]->id());

	int PDGId = EvtPDL::getStdHep( theId );
	// std::cout<<"PDGId here:   "<<PDGId<<std::endl;
    if ( theId.getId() == -1 && theId.getAlias() == -1 ) {
        std::cout << "Error. Could not find valid EvtId for " << parentName << std::endl;
        return -1;
    }

	EvtSpinType::spintype baseSpin = EvtPDL::getSpinType( theId );

    if ( baseSpin == EvtSpinType::VECTOR ) {
        std::cout << "Setting spin density for vector particle " << parentName << std::endl;
        spinDensity = new EvtSpinDensity();
        spinDensity->setDiag( EvtSpinType::getSpinStates( EvtSpinType::VECTOR ) );
        spinDensity->set( 1, 1, EvtComplex( 0.0, 0.0 ) );
    }

	// RapidVertex* vtx(0);
	// vtx = parts[0]->getOriginVertex();
	
	// Parent particle XYZ-position
	// ROOT::Math::XYZPoint point = vtx->getVertex(true);
	// EvtVector4R origin(0.0, point.X(), point.Y(), point.Z());
	EvtVector4R origin(0.0, 0.0, 0.0, 0.0);

	// Parent particle 4-momentum
	TLorentzVector pIn = parts[0]->getP();
	EvtVector4R pInit(pIn.E(),pIn.Px(),pIn.Py(),pIn.Pz());
	// double mass = EvtPDL::getMeanMass( theId );
	// EvtVector4R pInit(mass,0.0, 0.0, 0.0);

    EvtHepMCEvent* theEvent =
        evtGen_->generateDecay( PDGId, pInit, origin, spinDensity );	
	// Retrieve the HepMC event information
    GenEvent* hepMCEvent = theEvent->getEvent();

	std::list<GenVertexPtr> allVertices;
	
	double flightTime;

	int iVtx(0);
	int iPart(1);
	bool hasOsc(false);

	TLorentzVector p4TLV;
	FourVector FourMom;
	FourVector DecayVtx;
	FourVector OrigVtx;

	std::queue<int> nExpectedChildren;
	nExpectedChildren.push(parts[0]->nDaughters());

	for ( auto theVertex : hepMCEvent->vertices() ) {
		if ( theVertex == 0 ) {
		    continue;
		}
		auto nin  = theVertex->particles_in_size();
		auto nout = theVertex->particles_out_size();

		while (nin == 1 && nout ==1){
			//B meson oscillation
			continue;
		}

		uint nChildren = nExpectedChildren.front();

		// For these, get the mother decay vertex position and the 4-momentum to calculate
		// the flight time.

		for ( auto inParticle : theVertex->particles_in() ) {

        		if ( inParticle == 0 ) {
        	    	continue;
        		}

				int inPDGId = inParticle->pdg_id();
				FourMom = inParticle->momentum();
				DecayVtx = theVertex->position();

				if ( inPDGId == B0Id || inPDGId == BsId || inPDGId == B0barId || inPDGId == BsbarId ) {

					if (PDGId!=inPDGId){
						hasOsc=true;
					}
					else{
						hasOsc=false;

					}

					parts[iVtx]->setHasOsc(hasOsc);
					
				}				

				flightTime = calcFlightTime( DecayVtx, FourMom );
				parts[iVtx]->setId(inPDGId);
				parts[iVtx]->setDecaytime(flightTime);
				parts[iVtx]->getDecayVertex()->setXYZ(DecayVtx.x(), DecayVtx.y(), DecayVtx.x());
							
		}  

		uint nRecorded(0);
		for ( auto outParticle : theVertex->particles_out() ) {
			
			if (nRecorded < nChildren) {

				FourMom = outParticle->momentum();
				OrigVtx = theVertex->position();

				p4TLV.SetPxPyPzE(FourMom.px(),FourMom.py(),FourMom.pz(),FourMom.e());

				// std::cout<<"ID:"<<outParticle->pdg_id()<<std::endl;
				// std::cout<<"px: "<<FourMom.px()<<std::endl;
				// std::cout<<"py: "<<FourMom.py()<<std::endl;
				// std::cout<<"pz: "<<FourMom.pz()<<std::endl;
				// std::cout<<""<<std::endl;

				parts[iPart]->setP(p4TLV);

				parts[iPart]->getOriginVertex()->setXYZ(OrigVtx.x(),OrigVtx.y(),OrigVtx.z());
				parts[iPart]->setId(outParticle->pdg_id());

				nExpectedChildren.push(parts[iPart]->nDaughters());

				++iPart;

			}
			++nRecorded;
			
		}
		nExpectedChildren.pop();
		++iVtx;
	}

	delete hepMCEvent;

	return true;
#else
	if(!suppressWarning_) {
		std::cout << "WARNING in RapidExternalEvtGen::decay : EvtGen extension not compiled. Will not use EvtGen to decay " << parts[0]->name() << "." << std::endl;
		suppressWarning_=true;
	}

	return false;
#endif
}

bool RapidExternalEvtGen::setup() {
#ifdef RAPID_EVTGEN
	std::cout << "INFO in RapidExternalEvtGen::setup : Setting decay for external EvtGen generator." << std::endl;
	if(!evtGen_) setupGenerator();
	evtGen_->readUDecay(decFileName_.Data());
	return true;
#else
	std::cout << "WARNING in RapidExternalEvtGen::setup : EvtGen extension not compiled." << std::endl;
	return false;
#endif
}

bool RapidExternalEvtGen::setupGenerator() {
#ifdef RAPID_EVTGEN
	std::cout << "INFO in RapidExternalEvtGen::setupGenerator : Setting up external EvtGen generator." << std::endl;
	EvtRandomEngine* randomEngine = 0;
	EvtAbsRadCorr* radCorrEngine = 0;
	std::list<EvtDecayBase*> extraModels;

	// Define the random number generator
	uint seed = gRandom->TRandom::GetSeed();
	randomEngine = new EvtMTRandomEngine(seed);

	bool useEvtGenRandom(false);
	EvtExternalGenList genList(true, "", "gamma", useEvtGenRandom);
	radCorrEngine = genList.getPhotosModel();
	extraModels = genList.getListOfModels();

	TString evtPDLPath;
	evtPDLPath += getenv("EVTGEN_ROOT");
	evtPDLPath += "/evt.pdl";

	bool foundDec=false;

	TString decPath;
	decPath += getenv("RAPIDSIM_CONFIG");
	if(decPath!="") {
		decPath += "/config/evtgen/DECAY.DEC";
		if(!gSystem->AccessPathName(decPath)) foundDec=true;
	}

	// We want to initialise EvtGen before we define our DEC file so we can use EvtPDL
	// To do this pass an empty DEC file as the main decay file and pass our file later as a user file
	if(!foundDec) {
		decPath += getenv("RAPIDSIM_ROOT");
		decPath += "/config/evtgen/DECAY.DEC";
	}

	int mixingType = EvtCPUtil::Incoherent;
	// int mixingType = EvtCPUtil::Coherent;

	evtGen_ = new EvtGen(decPath.Data(), evtPDLPath.Data(), randomEngine,
			radCorrEngine, &extraModels, mixingType);

	return true;
#else
	std::cout << "WARNING in RapidExternalEvtGen::setup : EvtGen extension not compiled." << std::endl;
	return false;
#endif
}

void RapidExternalEvtGen::writeDecFile(TString fname, std::vector<RapidParticle*>& parts, bool usePhotos) {
#ifdef RAPID_EVTGEN
	if(!evtGen_) setupGenerator();

	decFileName_ = fname+".DEC";
	std::cout << "INFO in RapidExternalEvtGen::writeDecFile : Writing EvtGen DEC file : " << decFileName_ << std::endl;

	std::ofstream fout;
	fout.open(decFileName_, std::ofstream::out);

	if(usePhotos) {
		fout << "yesFSR\n" << std::endl;
	} else {
		fout << "noFSR\n" << std::endl;
	}
	TString defString;
	for(unsigned int iPart=0; iPart<parts.size(); ++iPart) {
		if(parts[iPart]->evtGenDefinitions()){
			defString = parts[iPart]->evtGenDefinitions();
			TObjArray* tokens = defString.Tokenize(",");
			for (int i = 0; i < tokens->GetEntries(); i++) {
    			TString token = tokens->At(i)->GetName();
    			token = token.Strip(TString::kBoth);  // removes leading and trailing whitespace
    			std::cout << "INFO: Defining in DecFile: " << token << std::endl;
				fout << "Define " << token << std::endl;
			}
			delete tokens;
		}

	}

	fout << "\n" ;

	// Loop over all particles and write out Decay rule for each
	for(unsigned int iPart=0; iPart<parts.size(); ++iPart) {
		unsigned int nChildren = parts[iPart]->nDaughters();
		if(nChildren>0) {
			int id = parts[iPart]->id();
			fout << "Decay " << getEvtGenName(id) << "\n1.00\t";
			if ( !(parts[iPart]->evtGenDecayModel()).Contains("TAUOLA") ) {
				for(unsigned int iChild=0; iChild<nChildren; ++iChild) {
					fout << getEvtGenName(parts[iPart]->daughter(iChild)->id()) << "\t";
				}
			}
			fout << parts[iPart]->evtGenDecayModel() << ";" << std::endl;
			fout <<"Enddecay" << std::endl;

			// Workaround to deal with mixing of B0 and Bs
			if(TMath::Abs(id)==531||TMath::Abs(id)==511) fout <<"CDecay " << getEvtGenConjName(id) << std::endl << std::endl;
		}
	}
	fout <<"End\n" << std::endl;
	fout.close();
#else
	std::cout << "WARNING in RapidExternalEvtGen::writeDecFile : EvtGen extension not compiled. Cannot write DEC file "
		  << fname << " for " << parts.size() << "particles with usePhotos=" << usePhotos << "." << std::endl;
#endif
}

TString RapidExternalEvtGen::getEvtGenName(int id) {
#ifdef RAPID_EVTGEN
	EvtId evtId = EvtPDL::evtIdFromStdHep(id);
	TString name = EvtPDL::name(evtId);
	return name;
#else
	std::cout << "WARNING in RapidExternalEvtGen::getEvtGenName : EvtGen extension not compiled. Cannot lookup name for particle ID " << id << "." << std::endl;
	return "";
#endif
}

TString RapidExternalEvtGen::getEvtGenConjName(int id) {
#ifdef RAPID_EVTGEN
	EvtId evtId = EvtPDL::evtIdFromStdHep(id);
	EvtId evtConjId = EvtPDL::chargeConj(evtId);
	TString name = EvtPDL::name(evtConjId);
	return name;
#else
	std::cout << "WARNING in RapidExternalEvtGen::getEvtGenConjName : EvtGen extension not compiled. Cannot lookup conjugate name for particle ID " << id << "." << std::endl;
	return "";
#endif
}

double RapidExternalEvtGen::calcFlightTime( FourVector& DecayVtx, FourVector& P4mtm )
{
    double flightTime( 0.0 );

#ifdef EVTGEN_HEPMC3
    double distance = DecayVtx.length();    // mm
    double momentum = P4mtm.length();       // GeV/c
	double mass = P4mtm.m();
#else
    double distance = DecayVtx.rho();    // mm
    double momentum = P4mtm.rho();  
	double mass = P4mtm.m();
#endif

	double c0 = EvtConst::c ;//mm/s
    if ( momentum > 0.0 ) {
        flightTime = 1.0e12 * distance * mass /
                     ( momentum * c0 );    // picoseconds
    }

	return flightTime;
}