/*
 * samFile.h
 *
 *  Created on: Mar 7, 2012
 *      Author: alexander.platzer
 */

#ifndef SAMFILE_H_
#define SAMFILE_H_

#include "general.h"

#include "SVs.h"

class SamFile
{
	friend class SVs_class;

	ifstream file;
	vector<string> substrings;

	int singleReadCount;

	struct SamReadPair
	{
		string name;

		struct SamRead
		{
			//chr... 1 to #, 0 if error..
			int chr;
			int leftmost;
			int flag;
			int quality;
			int length;

			string cigar;
			string hitTag;

			//to search for something which can be on different positions..
			string fullLine;

		};

		SamRead firstRead;
		SamRead secondRead;
	};

	SamReadPair samReadPair;

	int internNextReadPair(void);
	int filterReadPair(void);

	int cigarFilterpattern1(string &s);

public:
	ReadPairShort currentReadPairShort;

	ReadPairShort * selectedReadPairs;
	int selectedReadPairsn;


	SamFile();
	~SamFile();


	int initFile(string filename);

	int toNextReadPair(void);

	int selectReadPairs(string filename, SVs_class & SVs);

	int saveSelectedReadPairs(string filename);
	int loadSelectedReadPairs(string filename);
};




#endif /* SAMFILE_H_ */
