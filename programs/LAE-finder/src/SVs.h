/*
 * SVs.h
 *
 *  Created on: Mar 6, 2012
 *      Author: alexander.platzer
 */

#ifndef SVS_H_
#define SVS_H_

#include "general.h"


class SVs_class
{
	friend class SamFile;

	char chr1mask[CHR1SIZE + 1];
	char chr2mask[CHR2SIZE + 1];
	char chr3mask[CHR3SIZE + 1];
	char chr4mask[CHR4SIZE + 1];
	char chr5mask[CHR5SIZE + 1];

	void setRange(char * chr, int from, int to, int max);

	struct CoverageCollection01
	{
		int minimalCoverage;
		double meanCoverage;
	};

	CoverageCollection01 computeCoverageSingleReads(int chr, int from, int to, ReadPairShort * selectedReadPairs, int selectedReadPairsn);

	int computeCoverageReadpairs(int chr, int from, int to, ReadPairShort * selectedReadPairs, int selectedReadPairsn);

public:
//-------------------------------------------
//given SVs
	struct INV_RESULTs;
	struct TL_RESULTs;

	struct INV
	{
		int chr;
		int BP1_1;
		int BP1_2;
		int BP2_1;
		int BP2_2;

		INV_RESULTs * results;
	};
	struct TL
	{
		int chrBP;
		int BP_1;
		int BP_2;

		int chrRange;
		int from_1;
		int from_2;
		int to_1;
		int to_2;

		TL_RESULTs * results;
	};

	INV* invs;
	int invsn;

	TL* tls;
	int tlsn;
//-------------------------------------------
//results for SVs..

	struct INV_RESULTs
	{
		int index_INV;

		int readsSupporting;

		int continuousCoverageAsingle;
		int continuousCoverageBsingle;
		int continuousCoverageApair;
		int continuousCoverageBpair;

		double meanCoverageA;
		double meanCoverageB;

		INV_RESULTs()
		{
			memset(this, 0, sizeof(INV_RESULTs));
			index_INV = -1;
		}
	};

	INV_RESULTs* inv_results;


	struct TL_RESULTs
	{
		int index_TL;

		int readsSupportingFrom;
		int readsSupportingTo;
		int readsSupportingDeletion;

		int readsSupportingForward;
		int readsSupportingReverse;


		int amountOfStacks;
		int amountOfOverlappingSupportReads;

		double bestBPratioFrom;
		double bestBPratioTo;
		int bestBPC1;
		int bestBPC2;


		int continuousCoverageBPsingle;
		int continuousCoverageFROMsingle;
		int continuousCoverageTOsingle;
		int continuousCoverageBPpair;
		int continuousCoverageFROMpair;
		int continuousCoverageTOpair;


		double meanCoverageBP;
		double meanCoverageFROM;
		double meanCoverageTO;

		TL_RESULTs()
		{
			memset(this, 0, sizeof(INV_RESULTs));
			index_TL = -1;
		}
	};

	TL_RESULTs* tl_results;





//-------------------------------------------

	SVs_class();
	~SVs_class();

	static bool Comparison_INV(const SVs_class::INV &a, const SVs_class::INV &b);
	static bool Comparison_TL(const SVs_class::TL &a, const SVs_class::TL &b);
	int loadFile(string filename);


	void setMask(void);

	int analyseINVs(ReadPairShort * selectedReadPairs, int selectedReadPairsn);
	int analyseTLs(ReadPairShort * selectedReadPairs, int selectedReadPairsn);

	int saveINVresults(string filename);
	int saveTLresults(string filename);

};


#endif /* SVS_H_ */
