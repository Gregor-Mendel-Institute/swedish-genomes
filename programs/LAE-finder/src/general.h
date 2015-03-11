/*
 * general.h
 *
 *  Created on: Mar 6, 2012
 *      Author: alexander.platzer
 */

#ifndef GENERAL_H_
#define GENERAL_H_

#include <stdio.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>


using namespace std;


//-----------------------------------------------------------------------------
//section: general parameters

//for filtering:
#define minimalQ 30

//thresholds for distances
#define maximalDistanceFromBP 700
#define minimalDistanceofBPs 3*maximalDistanceFromBP

//trimming of reads for coverage and for the building of stacks
#define TRIMREADS 5

//minimal length to search for stacks, should be about 2 * expected read length
#define MINIMALLENGTHFORSTACKSEARCH 75*2

//maximal insert (for counting as read-pair-coverage
#define MAXIMALINSERT 1000


//-----------------------------------------------------------------------------
//section: general defs

#define CHR1SIZE 30427671
#define CHR2SIZE 19698289
#define CHR3SIZE 23459830
#define CHR4SIZE 18585056
#define CHR5SIZE 26975502

//-----------------------------------------------------------------------------
//section: general enums

enum RelativePos {NOWHERE, BP1, BP2, BP3, BP1_INNER, BP1_OUTER, BP2_INNER, BP2_OUTER, BP3_INNER, BP3_OUTER};
enum Orientation {FORWARD, REVERSE, MISSING};

//-----------------------------------------------------------------------------
//section: general structures


struct ReadPairShort
{
	int chr1;
	int leftmost1;

	int chr2;
	int leftmost2;

	int length;

	Orientation orientation;
};


struct ReadShort
{
	int chr;
	int leftmost1;
	int length;

	Orientation direction;

	int pairSide;
};

static bool Comparison_ReadShort(const ReadShort &a, const ReadShort &b)
{
	if(a.chr > b.chr)
		return 0;
	if(a.chr < b.chr)
		return 1;

	if(a.leftmost1 > b.leftmost1)
		return 0;
	if(a.leftmost1 < b.leftmost1)
		return 1;

	return 0;
}

struct ReadShortExtra
{
	int overlappingOfOthertype;
	int mask;
};


//-----------------------------------------------------------------------------
//section: general functions

inline void split(string & text, char * separators, vector<string> & words)
{
	string sep;
	sep = separators;

	int n = (int)text.length();
	int start, stop;
	start = 0;
	while ((start >= 0) && (start < n))
	{
		stop = (int)text.find_first_of(sep, start);
		if ((stop < 0) || (stop > n))
			stop = n;
		words.push_back(text.substr(start, stop - start));
		start = stop + 1;
	}
}




#endif /* GENERAL_H_ */
