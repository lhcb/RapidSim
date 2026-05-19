#include <cstdlib>
#include <iostream>
#include <ctime>

#include "TString.h"

#include "RapidAcceptance.h"
#include "RapidConfig.h"
#include "RapidDecay.h"
#include "RapidHistWriter.h"

int rapidSim(const TString mode, const int nEvtToGen, signed int nEvtToSelect=-1, bool saveTree=false, int nToReDecay=0) {

	clock_t t0,t1,t2;

	t0=clock();

	if(!getenv("RAPIDSIM_ROOT")) {
		std::cout << "ERROR in RapidSim : environment variable RAPIDSIM_ROOT is not set" << std::endl
			  << "                    Terminating" << std::endl;
		return 1;
	}

	if(nEvtToGen==0) {
		std::cout << "ERROR in RapidSim : nEvtToGenerate cannot be 0" << std::endl
			  << "                     Terminating" << std::endl;	
		return 1;
	} else if (nEvtToGen<0) {
		std::cout << "ERROR in RapidSim : nEvtToGenerate cannot be negative (" << nEvtToGen << ")" << std::endl
			  << "                     Terminating" << std::endl;	
		return 1;
	}

	std::cout << "INFO in RapidSim : nEvtToGenerate set to " << nEvtToGen << std::endl;
	
	if(nEvtToSelect==-1) {
		std::cout << "INFO in RapidSim : nEvtToSelect set to -1 --> will select all possible events" << std::endl
			  << "                   Specify nEvtToSelect if you want to limit the number of selected events" << std::endl;
	} else if (nEvtToSelect>nEvtToGen) {
		std::cout << "ERROR in RapidSim : nEvtToSelect (" << nEvtToSelect << ") is larger than nEvtToGenerate (" << nEvtToGen << ")" << std::endl
			  << "                     If you want to select as many events as possible, set nEvtToSelect to -1" << std::endl
			  << "                     Terminating" << std::endl;
		return 1;
	} else if (nEvtToSelect < 0) {
		std::cout << "ERROR in RapidSim : nEvtToSelect set to negative value other than -1 (" << nEvtToSelect << ")" << std::endl
			  << "                     Terminating" << std::endl;
		return 1;
	} else if (nEvtToSelect==0) {
		std::cout << "ERROR in RapidSim : nEvtToSelect cannot be 0" << std::endl
			  << "                     Terminating" << std::endl;	
		return 1;
	} else {
		std::cout << "INFO in RapidSim : nEvtToSelect set to " << nEvtToSelect 
				  << " --> will stop event loop early if this number is reached" << std::endl;
	}

	TString configEnv=getenv("RAPIDSIM_CONFIG");
	if(configEnv!="") {
		std::cout << "INFO in RapidSim : environment variable RAPIDSIM_CONFIG is set" << std::endl
			  << "                   Settings in " << configEnv << " will be used" << std::endl;
	}

	RapidConfig config;
	if(!config.load(mode)) {
		std::cout << "ERROR in RapidSim : failed to load configuration for decay mode " << mode << std::endl
			  << "                    Terminating" << std::endl;
		return 1;
	}

	RapidDecay* decay = config.getDecay();
	if(!decay) {
		std::cout << "ERROR in RapidSim : failed to setup decay for decay mode " << mode << std::endl
			  << "                    Terminating" << std::endl;
		return 1;
	}

	if(nToReDecay>0) {
		std::cout << "INFO in RapidSim : re-decay mode is active" << std::endl
			  << "                   Each parent will be re-decayed " << nToReDecay << " times" << std::endl;
	}

	RapidAcceptance* acceptance = config.getAcceptance();

	RapidHistWriter* writer = config.getWriter(saveTree);

	t1=clock();

	int ngenerated = 0; int nselected = 0;
	bool stopGeneration = false;
	for (Int_t n=0; n<nEvtToGen; ++n) {
		writer->setNEvent(n);
		if (!decay->generate()) continue;
		++ngenerated;

		if(acceptance->isSelected()) {
			++nselected;
			writer->fill();
			if (nselected>=nEvtToSelect && nEvtToSelect != -1) {
				std::cout << "INFO in RapidSim : Selected " << nselected << " events. Stopping generation early." << std::endl;
				stopGeneration = true;
			}
		}
		if (stopGeneration) break;

		for (Int_t nrd=0; nrd<nToReDecay; ++nrd) {
			if (!decay->generate(false)) continue;
			++ngenerated;

			if(!acceptance->isSelected()) continue;
			++nselected;

			writer->fill();
			if (nselected>=nEvtToSelect && nEvtToSelect != -1) {
				std::cout << "INFO in RapidSim : Selected " << nselected << " events. Stopping generation early." << std::endl;
				stopGeneration = true;
			}
			if (stopGeneration) break;
		}
		if (stopGeneration) break;
	}
	
	writer->save();
	
	t2=clock();
	
	std::cout << "INFO in RapidSim : Generated " << ngenerated << std::endl;
	std::cout << "INFO in RapidSim : Selected " << nselected << std::endl;
	if (nselected < nEvtToSelect && nEvtToSelect != -1) {
		std::cout << "WARNING in RapidSim : Only selected " << nselected << " events out of requested " << nEvtToSelect << std::endl;
		std::cout << "                      Consider increasing nEvtToGenerate" << std::endl;
	}
	std::cout << "INFO in RapidSim : " << (float(t1) - float(t0)) / CLOCKS_PER_SEC << " seconds to initialise." << std::endl;
	std::cout << "INFO in RapidSim : " << (float(t2) - float(t1)) / CLOCKS_PER_SEC << " seconds to generate." << std::endl;

	return 0;
}

int main(int argc, char * argv[])
{
	if (argc < 3) {
		printf("Usage: %s mode numberToGenerate [saveTree=0] [numberToRedecay=0] [numberToSelect=-1 (select all possible events)]\n", argv[0]);
		return 1;
	}

	const TString mode = argv[1];
	const int numberToGenerate = static_cast<int>(atoll(argv[2]));
	signed int numberToSelect = -1;
	bool saveTree = false;
	int nToReDecay = 0;

	if(argc>3) {
		saveTree = atoi(argv[3]);
	}
	if(argc>4) {
		nToReDecay = atoi(argv[4]);
	}
	if(argc>5) {
		numberToSelect = static_cast<signed int>(atoll(argv[5]));
	}
	
	int status = rapidSim(mode, numberToGenerate, numberToSelect, saveTree, nToReDecay);

	return status;
}
