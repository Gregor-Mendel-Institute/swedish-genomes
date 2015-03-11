/*
 * SVs.cpp
 *
 *  Created on: Mar 6, 2012
 *      Author: alexander.platzer
 */

#include "SVs.h"

int debug24 = 6;

SVs_class::SVs_class()
{
	memset(this, 0, sizeof(this));
}

SVs_class::~SVs_class()
{
	if(invs != NULL)
		delete [] invs;

	if(tls != NULL)
		delete [] tls;

	if(inv_results != NULL)
		delete [] inv_results;
}

bool SVs_class::Comparison_INV(const SVs_class::INV &a, const SVs_class::INV &b)
{
	if(a.chr > b.chr)
		return 0;
	if(a.chr < b.chr)
		return 1;

	if(a.BP1_1 > b.BP1_1)
		return 0;
	if(a.BP1_1 < b.BP1_1)
		return 1;

	return 0;
}


bool SVs_class::Comparison_TL(const SVs_class::TL &a, const SVs_class::TL &b)
{
	if(a.chrBP > b.chrBP)
		return 0;
	if(a.chrBP < b.chrBP)
		return 1;

	if(a.BP_1 > b.BP_1)
		return 0;
	if(a.BP_1 < b.BP_1)
		return 1;

	return 0;
}


int SVs_class::loadFile(string filename)
{
	int i,j;
	int invi, tli;

	string s1,s2;

	stringstream sx1;

	vector<string> substrings;


	ifstream file(filename.c_str());

	if(!file.is_open())
	{
		printf("file error\n");
		return 0;
	}


	getline(file, s1);
	while(!file.eof())
	{
		getline(file, s1);

		if(s1.length() > 2)
		{
			if(s1.find("INV\t") == 0)
				invsn++;
			if(s1.find("TL\t") == 0)
				tlsn++;
		}
	}
	file.close();
	file.clear();

	invs = new INV[invsn];
	tls = new TL[tlsn];


	file.open(filename.c_str());

	invi = 0;
	tli = 0;
	getline(file, s1);
	while(!file.eof())
	{
		getline(file, s1);

		if(s1.length() > 2)
		{
			j = 0;

			substrings.clear();

			split(s1, "\t", substrings);

			if(s1.find("INV\t") == 0)
			{
				if(substrings.size() < 6)
					cout << "format error at INV " << j << endl;

				for(j=0;j<(signed)substrings.size();j++)
				{
					sx1.str("");
					sx1.clear();

					if(substrings[j].length() == 0 && j < 6)
						cout << "format error at INV " << j << endl;
					else
						sx1 << substrings[j];

					switch(j)
					{
						case 0:
							//must be INV
							break;
						case 1:
							sx1 >> invs[invi].chr;
							break;
						case 2:
							sx1 >> invs[invi].BP1_1;
							break;
						case 3:
							sx1 >> invs[invi].BP1_2;
							break;
						case 4:
							sx1 >> invs[invi].BP2_1;
							break;
						case 5:
							sx1 >> invs[invi].BP2_2;
							break;
					}
				}

				//right order
				if(invs[invi].BP1_1 > invs[invi].BP1_2)
					swap(invs[invi].BP1_1, invs[invi].BP1_2);
				if(invs[invi].BP2_1 > invs[invi].BP2_2)
					swap(invs[invi].BP2_1, invs[invi].BP2_2);
				if(invs[invi].BP1_1 > invs[invi].BP2_1)
				{
					swap(invs[invi].BP1_1, invs[invi].BP2_1);
					swap(invs[invi].BP1_2, invs[invi].BP2_2);
				}

				//at least one bp distance
				if(invs[invi].BP1_2 - invs[invi].BP1_1 < 1)
					invs[invi].BP1_2++;
				if(invs[invi].BP2_2 - invs[invi].BP2_1 < 1)
					invs[invi].BP2_2++;

				//check the distances of the BPs
				if(invs[invi].BP2_1 - invs[invi].BP1_2 < minimalDistanceofBPs)
				{
					cout << "2 breakpoints of an INV are too close: INV" << invi;
					cout  << ", bp" << invs[invi].BP1_2 << " to bp" << invs[invi].BP2_1 << endl;

					return 0;
				}

				invs[invi].results = NULL;

				invi++;
			}


			if(s1.find("TL\t") == 0)
			{
				if(substrings.size() < 9)
					cout << "format error at TL " << j << endl;

				for(j=0;j<(signed)substrings.size();j++)
				{
					sx1.str("");
					sx1.clear();

					if(substrings[j].length() == 0 && j < 9)
						cout << "format error at TL " << j << endl;
					else
						sx1 << substrings[j];


					switch(j)
					{
						case 0:
							//must be TL
							break;
						case 1:
							sx1 >> tls[tli].chrRange;
							break;
						case 2:
							sx1 >> tls[tli].from_1;
							break;
						case 3:
							sx1 >> tls[tli].from_2;
							break;
						case 4:
							sx1 >> tls[tli].to_1;
							break;
						case 5:
							sx1 >> tls[tli].to_2;
							break;
						case 6:
							sx1 >> tls[tli].chrBP;
							break;
						case 7:
							sx1 >> tls[tli].BP_1;
							break;
						case 8:
							sx1 >> tls[tli].BP_2;
							break;
					}
				}

				//right order
				if(tls[tli].from_1 > tls[tli].to_1)
				{
					swap(tls[tli].from_1, tls[tli].to_1);
					swap(tls[tli].from_2, tls[tli].to_2);
				}
				if(tls[tli].from_1 > tls[tli].from_2)
					swap(tls[tli].from_1, tls[tli].from_2);
				if(tls[tli].to_1 > tls[tli].to_2)
					swap(tls[tli].to_1, tls[tli].to_2);
				if(tls[tli].BP_1 > tls[tli].BP_2)
					swap(tls[tli].BP_1, tls[tli].BP_2);

				//at least one bp distance
				if(tls[tli].from_2 - tls[tli].from_1 < 1)
					tls[tli].from_2++;
				if(tls[tli].to_2 - tls[tli].to_1 < 1)
					tls[tli].to_2++;

				//check the distances of the BPs
				if(tls[tli].to_1 - tls[tli].from_2 < minimalDistanceofBPs)
				{
					cout << "2 breakpoints of a TL are too close: TL" << tli;
					cout  << ", bp" << tls[tli].from_2 << " to bp" << tls[tli].to_1 << endl;

					return 0;
				}
				if(tls[tli].chrBP == tls[tli].chrRange)
				{
					if(tls[tli].BP_2 < tls[tli].from_1)
					{
						if(tls[tli].from_1 - tls[tli].BP_2 < minimalDistanceofBPs)
						{
							cout << "2 breakpoints of a TL are too close: TL" << tli;
							cout  << ", bp" << tls[tli].BP_2 << " to bp" << tls[tli].from_1 << endl;

							return 0;
						}
					}
					else if(tls[tli].BP_1 > tls[tli].to_2)
					{
						if(tls[tli].BP_1 - tls[tli].to_2 < minimalDistanceofBPs)
						{
							cout << "2 breakpoints of a TL are too close: TL" << tli;
							cout  << ", bp" << tls[tli].to_2 << " to bp" << tls[tli].BP_1 << endl;

							return 0;
						}
					}
					else
					{
						cout << "At TL" << tli << " the first breakpoint is in the (maybe) inserted range\n";
						return 0;
					}
				}

				tli++;
			}
		}
	}

	file.close();
	file.clear();

	sort(invs, invs + invsn, Comparison_INV);
	sort(tls, tls + tlsn, Comparison_TL);

	return 1;
}

void SVs_class::setRange(char * chr, int from, int to, int max)
{
	int i;
	int startIndex, endIndex;

	startIndex = from - maximalDistanceFromBP;
	endIndex = to + maximalDistanceFromBP;
	if(startIndex < 0)	startIndex = 0;
	if(endIndex > max) endIndex = max - 1;

	for(i=startIndex;i<endIndex;i++)
		chr[i] = 1;
}

void SVs_class::setMask(void)
{
	int i;
	char * chr;
	int chrSize;

	for(i=0;i<invsn;i++)
	{
		switch(invs[i].chr)
		{
			case 1:
				chr = chr1mask;	chrSize = CHR1SIZE;
				break;
			case 2:
				chr = chr2mask;	chrSize = CHR2SIZE;
				break;
			case 3:
				chr = chr3mask;	chrSize = CHR3SIZE;
				break;
			case 4:
				chr = chr4mask;	chrSize = CHR4SIZE;
				break;
			case 5:
				chr = chr5mask;	chrSize = CHR5SIZE;
				break;
			default:
				continue;
		}

		setRange(chr, invs[i].BP1_1, invs[i].BP1_2, chrSize);
		setRange(chr, invs[i].BP2_1, invs[i].BP2_2, chrSize);
	}

	for(i=0;i<tlsn;i++)
	{
		switch(tls[i].chrBP)
		{
			case 1:
				chr = chr1mask;	chrSize = CHR1SIZE;
				break;
			case 2:
				chr = chr2mask;	chrSize = CHR2SIZE;
				break;
			case 3:
				chr = chr3mask;	chrSize = CHR3SIZE;
				break;
			case 4:
				chr = chr4mask;	chrSize = CHR4SIZE;
				break;
			case 5:
				chr = chr5mask;	chrSize = CHR5SIZE;
				break;
			default:
				continue;
		}

		setRange(chr, tls[i].BP_1, tls[i].BP_2, chrSize);


		switch(tls[i].chrRange)
		{
			case 1:
				chr = chr1mask;	chrSize = CHR1SIZE;
				break;
			case 2:
				chr = chr2mask;	chrSize = CHR2SIZE;
				break;
			case 3:
				chr = chr3mask;	chrSize = CHR3SIZE;
				break;
			case 4:
				chr = chr4mask;	chrSize = CHR4SIZE;
				break;
			case 5:
				chr = chr5mask;	chrSize = CHR5SIZE;
				break;
			default:
				continue;
		}

		setRange(chr, tls[i].from_1, tls[i].from_2, chrSize);
		setRange(chr, tls[i].to_1, tls[i].to_2, chrSize);
	}
}

SVs_class::CoverageCollection01 SVs_class::computeCoverageSingleReads(int chr, int from, int to, ReadPairShort * selectedReadPairs, int selectedReadPairsn)
{
	CoverageCollection01 ret;

	int i, j;
	int tempi;

	int readStart, readEnd;
	int startIndex, endIndex;


	int * coverage;
	int coveragen;
	int minCoverage;
	int sum;

	coveragen = to - from + 1;
	coverage = new int[coveragen];
	memset(coverage, 0, sizeof(int) * coveragen);

	sum = 0;

	for(i=0;i<selectedReadPairsn;i++)
	{
		if(selectedReadPairs[i].chr1 == chr)
		{
			startIndex = max(from, selectedReadPairs[i].leftmost1);
			endIndex = min(to, selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length);

			tempi = endIndex - startIndex + 1;
			if(tempi > 0)
				sum += tempi;

			readStart = selectedReadPairs[i].leftmost1 + TRIMREADS;
			readEnd = selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length - TRIMREADS;

			startIndex = max(from, readStart) - from;
			endIndex = min(to, readEnd) - from;

			for(j=startIndex;j<=endIndex;j++)
				coverage[j]++;
		}

		if(selectedReadPairs[i].chr2 == chr)
		{
			startIndex = max(from, selectedReadPairs[i].leftmost2);
			endIndex = min(to, selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length);

			tempi = endIndex - startIndex + 1;
			if(tempi > 0)
				sum += tempi;

			readStart = selectedReadPairs[i].leftmost2 + TRIMREADS;
			readEnd = selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length - TRIMREADS;

			startIndex = max(from, readStart) - from;
			endIndex = min(to, readEnd) - from;

			for(j=startIndex;j<=endIndex;j++)
				coverage[j]++;
		}

	}

	ret.meanCoverage = (double)sum / coveragen;

	minCoverage = coverage[0];

	for(i=0;i<coveragen;i++)
	{	if(minCoverage > coverage[i])
			minCoverage = coverage[i];
	}

	ret.minimalCoverage = minCoverage;

	delete [] coverage;



	return ret;
}

int SVs_class::computeCoverageReadpairs(int chr, int from, int to, ReadPairShort * selectedReadPairs, int selectedReadPairsn)
{
	int i, j;
	int tempi;

	int readStart, readEnd;
	int startIndex, endIndex;


	int * coverage;
	int coveragen;
	int minCoverage;

	coveragen = to - from + 1;
	coverage = new int[coveragen];
	memset(coverage, 0, sizeof(int) * coveragen);

	for(i=0;i<selectedReadPairsn;i++)
	{
		if(selectedReadPairs[i].chr1 == chr && selectedReadPairs[i].chr2 == chr &&
		   selectedReadPairs[i].leftmost2 - (selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length) < MAXIMALINSERT
		)
		{
			readStart = selectedReadPairs[i].leftmost1 + TRIMREADS;
			readEnd = selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length - TRIMREADS;

			startIndex = max(from, readStart) - from;
			endIndex = min(to, readEnd) - from;

			for(j=startIndex;j<=endIndex;j++)
				coverage[j]++;
		}
		else
		{
			if(selectedReadPairs[i].chr1 == chr)
			{
				readStart = selectedReadPairs[i].leftmost1 + TRIMREADS;
				readEnd = selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length - TRIMREADS;

				startIndex = max(from, readStart) - from;
				endIndex = min(to, readEnd) - from;

				for(j=startIndex;j<=endIndex;j++)
					coverage[j]++;
			}

			if(selectedReadPairs[i].chr2 == chr)
			{
				readStart = selectedReadPairs[i].leftmost2 + TRIMREADS;
				readEnd = selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length - TRIMREADS;

				startIndex = max(from, readStart) - from;
				endIndex = min(to, readEnd) - from;

				for(j=startIndex;j<=endIndex;j++)
					coverage[j]++;
			}
		}
	}

	minCoverage = coverage[0];

	for(i=0;i<coveragen;i++)
	{	if(minCoverage > coverage[i])
			minCoverage = coverage[i];
	}

	delete [] coverage;

	return minCoverage;
}

//enum RelativePos {NOWHERE, BP1_INNER, BP1_OUTER, BP2_INNER, BP2_OUTER, BP3_INNER, BP3_OUTER,};
int SVs_class::analyseINVs(ReadPairShort * selectedReadPairs, int selectedReadPairsn)
{
	int inv_i;
	int i;

	RelativePos r1relative, r2relative;

	if(inv_results != NULL)
		delete [] inv_results;

	inv_results = new INV_RESULTs[invsn];

	for(inv_i=0;inv_i<invsn;inv_i++)
	{
		//connect elements
		inv_results[inv_i].index_INV = inv_i;
		invs[inv_i].results = &inv_results[inv_i];

		for(i=0;i<selectedReadPairsn;i++)
		{
			if(selectedReadPairs[i].chr1 != invs[inv_i].chr || selectedReadPairs[i].chr2 != invs[inv_i].chr)
				continue;

			if(selectedReadPairs[i].leftmost1 >= invs[inv_i].BP1_1 && selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length <= invs[inv_i].BP1_2)
				r1relative = BP1_INNER;
			else if(selectedReadPairs[i].leftmost1 >= invs[inv_i].BP1_1 - maximalDistanceFromBP &&
					selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length <= invs[inv_i].BP1_2 + maximalDistanceFromBP)
				r1relative = BP1_OUTER;
			else if(selectedReadPairs[i].leftmost1 >= invs[inv_i].BP2_1 && selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length <= invs[inv_i].BP2_2)
				r1relative = BP2_INNER;
			else if(selectedReadPairs[i].leftmost1 >= invs[inv_i].BP2_1 - maximalDistanceFromBP &&
					selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length <= invs[inv_i].BP2_2 + maximalDistanceFromBP)
				r1relative = BP2_OUTER;
			else
				r1relative = NOWHERE;

			if(selectedReadPairs[i].leftmost2 >= invs[inv_i].BP1_1 && selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length <= invs[inv_i].BP1_2)
				r2relative = BP1_INNER;
			else if(selectedReadPairs[i].leftmost2 >= invs[inv_i].BP1_1 - maximalDistanceFromBP &&
					selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length <= invs[inv_i].BP1_2 + maximalDistanceFromBP)
				r2relative = BP1_OUTER;
			else if(selectedReadPairs[i].leftmost2 >= invs[inv_i].BP2_1 && selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length <= invs[inv_i].BP2_2)
				r2relative = BP2_INNER;
			else if(selectedReadPairs[i].leftmost2 >= invs[inv_i].BP2_1 - maximalDistanceFromBP &&
					selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length <= invs[inv_i].BP2_2 + maximalDistanceFromBP)
				r2relative = BP2_OUTER;
			else
				r2relative = NOWHERE;

			if(((r1relative == BP1_INNER && (r2relative == BP2_INNER || r2relative == BP2_OUTER)) ||
			    (r2relative == BP1_INNER && (r1relative == BP2_INNER || r1relative == BP2_OUTER)) ||
			    (r1relative == BP2_INNER && (r2relative == BP1_INNER || r2relative == BP1_OUTER)) ||
			    (r2relative == BP2_INNER && (r1relative == BP1_INNER || r1relative == BP1_OUTER))) &&
			    selectedReadPairs[i].orientation == REVERSE
			   )
			inv_results[inv_i].readsSupporting++;
		}


		//compute the coverages:
		CoverageCollection01 cc1;

		cc1 = computeCoverageSingleReads(invs[inv_i].chr, invs[inv_i].BP1_1, invs[inv_i].BP1_2, selectedReadPairs, selectedReadPairsn);
		inv_results[inv_i].meanCoverageA = cc1.meanCoverage;
		inv_results[inv_i].continuousCoverageAsingle = cc1.minimalCoverage;

		cc1 = computeCoverageSingleReads(invs[inv_i].chr, invs[inv_i].BP2_1, invs[inv_i].BP2_2, selectedReadPairs, selectedReadPairsn);
		inv_results[inv_i].meanCoverageB = cc1.meanCoverage;
		inv_results[inv_i].continuousCoverageBsingle = cc1.minimalCoverage;

		inv_results[inv_i].continuousCoverageApair = computeCoverageReadpairs(invs[inv_i].chr, invs[inv_i].BP1_1, invs[inv_i].BP1_2, selectedReadPairs, selectedReadPairsn);
		inv_results[inv_i].continuousCoverageBpair = computeCoverageReadpairs(invs[inv_i].chr, invs[inv_i].BP2_1, invs[inv_i].BP2_2, selectedReadPairs, selectedReadPairsn);
	}

	return 1;
}


ReadShort * TLstartreads;
ReadShortExtra * TLstartreadsExtra;

static bool Comparison_TLstartreadsIndex(const int &a, const int &b)
{
	if(TLstartreadsExtra[a].mask == 1 && TLstartreadsExtra[b].mask == 0)
		return 0;
	if(TLstartreadsExtra[a].mask == 0 && TLstartreadsExtra[b].mask == 1)
		return 1;

	if(TLstartreadsExtra[a].overlappingOfOthertype > TLstartreadsExtra[b].overlappingOfOthertype)
		return 0;
	if(TLstartreadsExtra[a].overlappingOfOthertype < TLstartreadsExtra[b].overlappingOfOthertype)
		return 1;

	return 0;
}



//enum RelativePos {NOWHERE, BP1_INNER, BP1_OUTER, BP2_INNER, BP2_OUTER, BP3_INNER, BP3_OUTER,};
int SVs_class::analyseTLs(ReadPairShort * selectedReadPairs, int selectedReadPairsn)
{
	int tl_i;
	int i, j;
	int start, end;

	int sum;


	//-----------------------------------------------
	//for the stacks:
	int TLstartreadsn = 0;
	int TLstartreads_i;
	int * TLstartreadsSortIndex;

	char * mappedToBP2;
	char * mappedToBP3;

	int state;
	int only2;
	int only3;


	int readLength;

	//--------------------------------------------------------------------------------------------


	RelativePos r1relative, r2relative;

	if(tl_results != NULL)
		delete [] tl_results;

	tl_results = new TL_RESULTs[tlsn];

	for(tl_i=0;tl_i<tlsn;tl_i++)
	{
		//connect elements
		tl_results[tl_i].index_TL = tl_i;
		tls[tl_i].results = &tl_results[tl_i];

		for(i=0;i<selectedReadPairsn;i++)
		{
			//the unrelated case
			if((selectedReadPairs[i].chr1 != tls[tl_i].chrBP || selectedReadPairs[i].chr2 != tls[tl_i].chrBP) &&
				(selectedReadPairs[i].chr1 != tls[tl_i].chrRange || selectedReadPairs[i].chr2 != tls[tl_i].chrRange))
			{	continue;	}

			//the position in short -----------------------------------------------------------
			r1relative = NOWHERE;
			if(tls[tl_i].chrBP == selectedReadPairs[i].chr1)
			{	if(selectedReadPairs[i].leftmost1 >= tls[tl_i].BP_1 - maximalDistanceFromBP &&
						selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length < tls[tl_i].BP_2 + maximalDistanceFromBP)
				{	r1relative = BP1;	}
			}
			if(tls[tl_i].chrRange == selectedReadPairs[i].chr1)
			{
				if(selectedReadPairs[i].leftmost1 >= tls[tl_i].from_1 &&
						selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length < tls[tl_i].from_2 + maximalDistanceFromBP)
				{	r1relative = BP2;	}

				if(selectedReadPairs[i].leftmost1 >= tls[tl_i].from_1 - maximalDistanceFromBP &&
						selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length < tls[tl_i].from_1)
				{	r1relative = BP2_OUTER;	}


				if(selectedReadPairs[i].leftmost1 >= tls[tl_i].to_1 - maximalDistanceFromBP &&
						selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length < tls[tl_i].to_2)
				{	r1relative = BP3;	}

				if(selectedReadPairs[i].leftmost1 >= tls[tl_i].to_2 &&
						selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length < tls[tl_i].to_2 + maximalDistanceFromBP)
				{	r1relative = BP3_OUTER;	}
			}


			r2relative = NOWHERE;
			if(tls[tl_i].chrBP == selectedReadPairs[i].chr2)
			{	if(selectedReadPairs[i].leftmost2 >= tls[tl_i].BP_1 - maximalDistanceFromBP &&
						selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length < tls[tl_i].BP_2 + maximalDistanceFromBP)
				{	r2relative = BP1;	}
			}
			if(tls[tl_i].chrRange == selectedReadPairs[i].chr2)
			{
				if(selectedReadPairs[i].leftmost2 >= tls[tl_i].from_1 &&
						selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length < tls[tl_i].from_2 + maximalDistanceFromBP)
				{	r2relative = BP2;	}

				if(selectedReadPairs[i].leftmost2 >= tls[tl_i].from_1 - maximalDistanceFromBP &&
						selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length < tls[tl_i].from_1)
				{	r2relative = BP2_OUTER;	}


				if(selectedReadPairs[i].leftmost2 >= tls[tl_i].to_1 - maximalDistanceFromBP &&
						selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length < tls[tl_i].to_2)
				{	r2relative = BP3;	}

				if(selectedReadPairs[i].leftmost2 >= tls[tl_i].to_2 &&
						selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length < tls[tl_i].to_2 + maximalDistanceFromBP)
				{	r2relative = BP3_OUTER;	}
			}

			//------------------------------------------------------
			//count support general

			if((r1relative == BP1 && r2relative == BP2) ||
			   (r1relative == BP2 && r2relative == BP1)
			)
			{
				tl_results[tl_i].readsSupportingFrom++;
				if(selectedReadPairs[i].orientation == FORWARD)
					tl_results[tl_i].readsSupportingForward++;
				if(selectedReadPairs[i].orientation == REVERSE)
					tl_results[tl_i].readsSupportingReverse++;
			}

			if((r1relative == BP1 && r2relative == BP3) ||
			   (r1relative == BP3 && r2relative == BP1)
			)
			{
				tl_results[tl_i].readsSupportingTo++;
				if(selectedReadPairs[i].orientation == FORWARD)
					tl_results[tl_i].readsSupportingForward++;
				if(selectedReadPairs[i].orientation == REVERSE)
					tl_results[tl_i].readsSupportingReverse++;
			}

			if( ((r1relative == BP2 || r1relative == BP2_OUTER) && (r2relative == BP3 || r2relative == BP3_OUTER)) ||
			   ((r1relative == BP3 || r1relative == BP3_OUTER) && (r2relative == BP2 || r2relative == BP2_OUTER))
			)
			{
				tl_results[tl_i].readsSupportingDeletion++;
				if(selectedReadPairs[i].orientation == FORWARD)
					tl_results[tl_i].readsSupportingForward++;
				if(selectedReadPairs[i].orientation == REVERSE)
					tl_results[tl_i].readsSupportingReverse++;
			}
		}

		//------------------------------------------------------
		//compute stacks
		tl_results[tl_i].amountOfStacks = -1;
		tl_results[tl_i].amountOfOverlappingSupportReads = -1;
		tl_results[tl_i].bestBPratioFrom = -1;
		tl_results[tl_i].bestBPratioTo = -1;
		tl_results[tl_i].bestBPC1 = -1;
		tl_results[tl_i].bestBPC2 = -1;


		if(tls[tl_i].BP_2 - tls[tl_i].BP_1 > MINIMALLENGTHFORSTACKSEARCH && tl_results[tl_i].readsSupportingFrom + tl_results[tl_i].readsSupportingTo > 0)
		{
			TLstartreadsn = tl_results[tl_i].readsSupportingFrom + tl_results[tl_i].readsSupportingTo;
			TLstartreads = new ReadShort[TLstartreadsn];

			//revisit the reads..
				TLstartreads_i = 0;
				for(i=0;i<selectedReadPairsn;i++)
				{
					//the unrelated case
					if((selectedReadPairs[i].chr1 != tls[tl_i].chrBP || selectedReadPairs[i].chr2 != tls[tl_i].chrBP) &&
						(selectedReadPairs[i].chr1 != tls[tl_i].chrRange || selectedReadPairs[i].chr2 != tls[tl_i].chrRange))
					{	continue;	}

					//the position in short -----------------------------------------------------------
					r1relative = NOWHERE;
					if(tls[tl_i].chrBP == selectedReadPairs[i].chr1)
					{	if(selectedReadPairs[i].leftmost1 >= tls[tl_i].BP_1 - maximalDistanceFromBP &&
								selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length < tls[tl_i].BP_2 + maximalDistanceFromBP)
						{	r1relative = BP1;	}
					}
					if(tls[tl_i].chrRange == selectedReadPairs[i].chr1)
					{
						if(selectedReadPairs[i].leftmost1 >= tls[tl_i].from_1 &&
								selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length < tls[tl_i].from_2 + maximalDistanceFromBP)
						{	r1relative = BP2;	}

						if(selectedReadPairs[i].leftmost1 >= tls[tl_i].from_1 - maximalDistanceFromBP &&
								selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length < tls[tl_i].from_1)
						{	r1relative = BP2_OUTER;	}


						if(selectedReadPairs[i].leftmost1 >= tls[tl_i].to_1 - maximalDistanceFromBP &&
								selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length < tls[tl_i].to_2)
						{	r1relative = BP3;	}

						if(selectedReadPairs[i].leftmost1 >= tls[tl_i].to_2 &&
								selectedReadPairs[i].leftmost1 + selectedReadPairs[i].length < tls[tl_i].to_2 + maximalDistanceFromBP)
						{	r1relative = BP3_OUTER;	}
					}


					r2relative = NOWHERE;
					if(tls[tl_i].chrBP == selectedReadPairs[i].chr2)
					{	if(selectedReadPairs[i].leftmost2 >= tls[tl_i].BP_1 - maximalDistanceFromBP &&
								selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length < tls[tl_i].BP_2 + maximalDistanceFromBP)
						{	r2relative = BP1;	}
					}
					if(tls[tl_i].chrRange == selectedReadPairs[i].chr2)
					{
						if(selectedReadPairs[i].leftmost2 >= tls[tl_i].from_1 &&
								selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length < tls[tl_i].from_2 + maximalDistanceFromBP)
						{	r2relative = BP2;	}

						if(selectedReadPairs[i].leftmost2 >= tls[tl_i].from_1 - maximalDistanceFromBP &&
								selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length < tls[tl_i].from_1)
						{	r2relative = BP2_OUTER;	}


						if(selectedReadPairs[i].leftmost2 >= tls[tl_i].to_1 - maximalDistanceFromBP &&
								selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length < tls[tl_i].to_2)
						{	r2relative = BP3;	}

						if(selectedReadPairs[i].leftmost2 >= tls[tl_i].to_2 &&
								selectedReadPairs[i].leftmost2 + selectedReadPairs[i].length < tls[tl_i].to_2 + maximalDistanceFromBP)
						{	r2relative = BP3_OUTER;	}
					}

					//------------------------------------------------------
					//find the reads again

					if((r1relative == BP1 && r2relative == BP2) ||
					   (r1relative == BP2 && r2relative == BP1)
					)
					{
						if(r1relative == BP1)
						{
							TLstartreads[TLstartreads_i].chr = selectedReadPairs[i].chr1;
							TLstartreads[TLstartreads_i].leftmost1 = selectedReadPairs[i].leftmost1;
							TLstartreads[TLstartreads_i].length = selectedReadPairs[i].length;
							TLstartreads[TLstartreads_i].direction = selectedReadPairs[i].orientation;
							TLstartreads[TLstartreads_i].pairSide = 0;

						}
						else if(r2relative == BP1)
						{
							TLstartreads[TLstartreads_i].chr = selectedReadPairs[i].chr2;
							TLstartreads[TLstartreads_i].leftmost1 = selectedReadPairs[i].leftmost2;
							TLstartreads[TLstartreads_i].length = selectedReadPairs[i].length;
							TLstartreads[TLstartreads_i].direction = selectedReadPairs[i].orientation;
							TLstartreads[TLstartreads_i].pairSide = 0;
						}
						else
						{
							TLstartreads[TLstartreads_i].chr = -1;
							TLstartreads[TLstartreads_i].leftmost1 = -1;
							TLstartreads[TLstartreads_i].length = -1;
							TLstartreads[TLstartreads_i].direction = MISSING;
							TLstartreads[TLstartreads_i].pairSide = -1;

							cout << "mark 67\n";
						}

						TLstartreads_i++;
					}

					if((r1relative == BP1 && r2relative == BP3) ||
					   (r1relative == BP3 && r2relative == BP1)
					)
					{
						if(r1relative == BP1)
						{
							TLstartreads[TLstartreads_i].chr = selectedReadPairs[i].chr1;
							TLstartreads[TLstartreads_i].leftmost1 = selectedReadPairs[i].leftmost1;
							TLstartreads[TLstartreads_i].length = selectedReadPairs[i].length;
							TLstartreads[TLstartreads_i].direction = selectedReadPairs[i].orientation;
							TLstartreads[TLstartreads_i].pairSide = 1;

						}
						else if(r2relative == BP1)
						{
							TLstartreads[TLstartreads_i].chr = selectedReadPairs[i].chr2;
							TLstartreads[TLstartreads_i].leftmost1 = selectedReadPairs[i].leftmost2;
							TLstartreads[TLstartreads_i].length = selectedReadPairs[i].length;
							TLstartreads[TLstartreads_i].direction = selectedReadPairs[i].orientation;
							TLstartreads[TLstartreads_i].pairSide = 1;
						}
						else
						{
							TLstartreads[TLstartreads_i].chr = -1;
							TLstartreads[TLstartreads_i].leftmost1 = -1;
							TLstartreads[TLstartreads_i].length = -1;
							TLstartreads[TLstartreads_i].direction = MISSING;
							TLstartreads[TLstartreads_i].pairSide = -1;

							cout << "mark 68\n";
						}

						TLstartreads_i++;
					}

				}

				sum = 0;
				for(TLstartreads_i=0;TLstartreads_i<TLstartreadsn;TLstartreads_i++)
					sum += TLstartreads[TLstartreads_i].length;

				readLength = sum / TLstartreadsn;

				sort(TLstartreads, TLstartreads + TLstartreadsn, Comparison_ReadShort);

				//----------------------------------------------------

				mappedToBP2 = new char[tls[tl_i].BP_2 - tls[tl_i].BP_1 + 2*maximalDistanceFromBP];
				mappedToBP3 = new char[tls[tl_i].BP_2 - tls[tl_i].BP_1 + 2*maximalDistanceFromBP];

				for(i=0;i<tls[tl_i].BP_2 - tls[tl_i].BP_1 + 2*maximalDistanceFromBP;i++)
				{
					mappedToBP2[i] = 0;
					mappedToBP3[i] = 0;
				}

				for(TLstartreads_i=0;TLstartreads_i<TLstartreadsn;TLstartreads_i++)
				{
					start = TLstartreads[TLstartreads_i].leftmost1 - (tls[tl_i].BP_1 - maximalDistanceFromBP);
					if(start < 0) start = 0;
					if(start >= tls[tl_i].BP_2 - tls[tl_i].BP_1 + 2*maximalDistanceFromBP) start = tls[tl_i].BP_2 + 2*maximalDistanceFromBP - tls[tl_i].BP_1 - 1;

					end = TLstartreads[TLstartreads_i].leftmost1 + TLstartreads[TLstartreads_i].length - (tls[tl_i].BP_1 - maximalDistanceFromBP);
					if(end < 0) end = 0;
					if(end >= tls[tl_i].BP_2 - tls[tl_i].BP_1 + 2*maximalDistanceFromBP) end = tls[tl_i].BP_2 - tls[tl_i].BP_1 + 2*maximalDistanceFromBP - 1;

					if(TLstartreads[TLstartreads_i].pairSide == 0)
					{
						for(j=start;j<end;j++)
							mappedToBP2[j] = 1;
					}
					else if(TLstartreads[TLstartreads_i].pairSide == 1)
					{
						for(j=start;j<end;j++)
							mappedToBP3[j] = 1;
					}
				}

				tl_results[tl_i].amountOfStacks = 0;

				state = 0;
				only2 = 0;
				only3 = 0;

				for(i=0;i<tls[tl_i].BP_2 - tls[tl_i].BP_1 + 2*maximalDistanceFromBP;i++)
				{
					if(state == 0 && mappedToBP2[i])
					{	state = 2;	tl_results[tl_i].amountOfStacks++; }
					if(state == 0 && mappedToBP3[i])
					{	state = 3;	tl_results[tl_i].amountOfStacks++;	}

					if(state == 2)
					{
						if(mappedToBP2[i])
							only3 = 0;
						else if(mappedToBP3[i])
							only3++;

						if(only3 > readLength - 2 * TRIMREADS)
						{
							state = 3;
							only2 = 0;
							tl_results[tl_i].amountOfStacks++;
						}
					}else if(state == 3)
					{
						if(mappedToBP3[i])
							only2 = 0;
						else if(mappedToBP2[i])
							only2++;

						if(only2 > readLength - 2 * TRIMREADS)
						{
							state = 2;
							only3 = 0;
							tl_results[tl_i].amountOfStacks++;
						}
					}
				}

				delete [] mappedToBP2;
				delete [] mappedToBP3;

				TLstartreadsExtra = new ReadShortExtra[TLstartreadsn];
				TLstartreadsSortIndex = new int[TLstartreadsn];

				for(i=0;i<TLstartreadsn;i++)
				{
					TLstartreadsExtra[i].overlappingOfOthertype = 0;
					TLstartreadsExtra[i].mask = 0;
					TLstartreadsSortIndex[i] = i;
				}

				for(i=0;i<TLstartreadsn;i++)
				{
					for(j=i+1;j<TLstartreadsn;j++)
					{
						if(TLstartreads[j].leftmost1 + TRIMREADS > TLstartreads[i].leftmost1 + TLstartreads[i].length - TRIMREADS)
						{	break;	}

						//overlapping reads
						if(TLstartreads[i].pairSide != TLstartreads[j].pairSide)
						{
							TLstartreadsExtra[i].overlappingOfOthertype++;
							TLstartreadsExtra[j].overlappingOfOthertype++;
						}
					}
				}

				sort(TLstartreadsSortIndex, TLstartreadsSortIndex + TLstartreadsn, Comparison_TLstartreadsIndex);


				while(TLstartreadsExtra[0].overlappingOfOthertype > 0 && TLstartreadsExtra[0].mask == 0)
				{
					for(i=0;i<TLstartreadsn;i++)
					{
						for(j=i+1;j<TLstartreadsn;j++)
						{
							if(TLstartreads[j].leftmost1 + TRIMREADS > TLstartreads[i].leftmost1 + TLstartreads[i].length - TRIMREADS)
							{	break;	}

							//overlapping reads
							if(TLstartreads[i].pairSide != TLstartreads[j].pairSide)
							{
								TLstartreadsExtra[i].overlappingOfOthertype++;
								TLstartreadsExtra[j].overlappingOfOthertype++;
							}
						}
					}

					sort(TLstartreadsSortIndex, TLstartreadsSortIndex + TLstartreadsn, Comparison_TLstartreadsIndex);

					if(TLstartreadsExtra[TLstartreadsSortIndex[0]].mask == 0 && TLstartreadsExtra[TLstartreadsSortIndex[0]].overlappingOfOthertype > 0)
					{	//remove that read from further stack calculation

						for(i=TLstartreadsSortIndex[0];i>=0;i--)
						{
							if(TLstartreads[TLstartreadsSortIndex[0]].leftmost1 + TRIMREADS > TLstartreads[i].leftmost1 + TLstartreads[i].length - TRIMREADS)
								break;

							if(TLstartreads[TLstartreadsSortIndex[0]].pairSide != TLstartreads[i].pairSide)
							{
								TLstartreadsExtra[TLstartreadsSortIndex[0]].overlappingOfOthertype--;
								TLstartreadsExtra[i].overlappingOfOthertype--;
							}
						}

						for(i=TLstartreadsSortIndex[0];i<TLstartreadsn;i++)
						{
							if(TLstartreads[TLstartreadsSortIndex[0]].leftmost1 + TLstartreads[TLstartreadsSortIndex[0]].length - TRIMREADS < TLstartreads[i].leftmost1 + TRIMREADS)
								break;

							if(TLstartreads[TLstartreadsSortIndex[0]].pairSide != TLstartreads[i].pairSide)
							{
								TLstartreadsExtra[TLstartreadsSortIndex[0]].overlappingOfOthertype--;
								TLstartreadsExtra[i].overlappingOfOthertype--;
							}
						}

						if(TLstartreadsExtra[TLstartreadsSortIndex[0]].overlappingOfOthertype != 0)
						{
							cout << "error at filtering read..";
						}

						TLstartreadsExtra[TLstartreadsSortIndex[0]].mask = 1;
					}
				}

				tl_results[tl_i].amountOfOverlappingSupportReads = 0;

				for(i=0;i<TLstartreadsn;i++)
					if(TLstartreadsExtra[i].mask)
						tl_results[tl_i].amountOfOverlappingSupportReads++;

				//------------------------------------------------------
				//only non-overlapping reads are left (others should be masked)
				/*int bestBPratioFrom;
				int bestBPratioTo;
				int bestBPC1;
				int bestBPC2;
				//*/
				//TLstartreads
				int from_forward = 0;
				int from_reverse = 0;
				int to_forward = 0;
				int to_reverse = 0;

				int stack1last = -1;
				int stack2first;
				int stack2last;

				double FRratio = 0;
				double FRratiobest = 0;

				int state = 0;

				tl_results[tl_i].bestBPratioFrom = 0;
				tl_results[tl_i].bestBPratioTo = 0;

				if(TLstartreadsn > 0)
				{
					//first read
					i=0;
					while(TLstartreadsExtra[i].mask) i++;
					if(i>=TLstartreadsn) break;

					state = TLstartreads[i].pairSide + 2;
					if(state == 2)
					{
						if(TLstartreads[i].direction == FORWARD)
							from_forward++;
						if(TLstartreads[i].direction == REVERSE)
							from_reverse++;
					}
					if(state == 3)
					{
						if(TLstartreads[i].direction == FORWARD)
							to_forward++;
						if(TLstartreads[i].direction == REVERSE)
							to_reverse++;
					}
					stack2first = TLstartreads[i].leftmost1;
					stack2last = TLstartreads[i].leftmost1 + TLstartreads[i].length;

					for(i=i+1;i<TLstartreadsn;i++)
						if(!TLstartreadsExtra[i].mask)
						{
							if(state == TLstartreads[i].pairSide + 2)
							{
								if(state == 2)
								{
									if(TLstartreads[i].direction == FORWARD)
										from_forward++;
									if(TLstartreads[i].direction == REVERSE)
										from_reverse++;
								}
								if(state == 3)
								{
									if(TLstartreads[i].direction == FORWARD)
										to_forward++;
									if(TLstartreads[i].direction == REVERSE)
										to_reverse++;
								}
								if(stack2first > TLstartreads[i].leftmost1)
									stack2first = TLstartreads[i].leftmost1;
								if(stack2last < TLstartreads[i].leftmost1 + TLstartreads[i].length)
									stack2last = TLstartreads[i].leftmost1 + TLstartreads[i].length;
							}
							else
							{
								if(stack1last != -1)
								{
									//compare
									if(from_forward + to_forward > from_reverse + to_reverse)
										FRratio = (from_forward + to_forward) / (from_forward + to_forward + from_reverse + to_reverse);
									else
										FRratio = (from_reverse + to_reverse) / (from_forward + to_forward + from_reverse + to_reverse);

									if(FRratio > FRratiobest)
									{
										if(from_forward > from_reverse)
											tl_results[tl_i].bestBPratioFrom = from_forward / (from_forward + from_reverse);
										else
											tl_results[tl_i].bestBPratioFrom = from_reverse / (from_forward + from_reverse);

										if(to_forward > to_reverse)
											tl_results[tl_i].bestBPratioTo = to_forward / (to_forward + to_reverse);
										else
											tl_results[tl_i].bestBPratioTo = to_reverse / (to_forward + to_reverse);

										tl_results[tl_i].bestBPC1 = stack1last;
										tl_results[tl_i].bestBPC2 = stack2first;

										FRratiobest = FRratio;
									}
								}

								//init
								if(state == 2)
								{
									to_forward = 0;
									to_reverse = 0;
								}
								if(state == 3)
								{
									from_forward = 0;
									from_reverse = 0;
								}

								stack1last = stack2last;

								state = TLstartreads[i].pairSide + 2;

								if(state == 2)
								{
									if(TLstartreads[i].direction == FORWARD)
										from_forward++;
									if(TLstartreads[i].direction == REVERSE)
										from_reverse++;
								}
								if(state == 3)
								{
									if(TLstartreads[i].direction == FORWARD)
										to_forward++;
									if(TLstartreads[i].direction == REVERSE)
										to_reverse++;
								}
								stack2first = TLstartreads[i].leftmost1;
								stack2last = TLstartreads[i].leftmost1 + TLstartreads[i].length;
							}
						}

				}

				//last stack
				if(stack1last != -1)
				{
					//compare
					if(from_forward + to_forward > from_reverse + to_reverse)
						FRratio = (from_forward + to_forward) / (from_forward + to_forward + from_reverse + to_reverse);
					else
						FRratio = (from_reverse + to_reverse) / (from_forward + to_forward + from_reverse + to_reverse);

					if(FRratio > FRratiobest)
					{
						if(from_forward > from_reverse)
							tl_results[tl_i].bestBPratioFrom = from_forward / (from_forward + from_reverse);
						else
							tl_results[tl_i].bestBPratioFrom = from_reverse / (from_forward + from_reverse);

						if(to_forward > to_reverse)
							tl_results[tl_i].bestBPratioTo = to_forward / (to_forward + to_reverse);
						else
							tl_results[tl_i].bestBPratioTo = to_reverse / (to_forward + to_reverse);

						tl_results[tl_i].bestBPC1 = stack1last;
						tl_results[tl_i].bestBPC2 = stack2first;

						FRratiobest = FRratio;
					}
				}


				delete [] TLstartreadsSortIndex;
				delete [] TLstartreadsExtra ;

			delete [] TLstartreads;
		}
		//------------------------------------------------------




		//compute the coverages:
		CoverageCollection01 cc1;

		cc1 = computeCoverageSingleReads(tls[tl_i].chrBP, tls[tl_i].BP_1, tls[tl_i].BP_2, selectedReadPairs, selectedReadPairsn);
		tl_results[tl_i].meanCoverageBP = cc1.meanCoverage;
		tl_results[tl_i].continuousCoverageBPsingle = cc1.minimalCoverage;

		tl_results[tl_i].continuousCoverageBPpair = computeCoverageReadpairs(tls[tl_i].chrBP, tls[tl_i].BP_1, tls[tl_i].BP_2, selectedReadPairs, selectedReadPairsn);


		cc1 = computeCoverageSingleReads(tls[tl_i].chrRange, tls[tl_i].from_1, tls[tl_i].from_2, selectedReadPairs, selectedReadPairsn);
		tl_results[tl_i].meanCoverageFROM = cc1.meanCoverage;
		tl_results[tl_i].continuousCoverageFROMsingle = cc1.minimalCoverage;

		tl_results[tl_i].continuousCoverageFROMpair = computeCoverageReadpairs(tls[tl_i].chrRange, tls[tl_i].from_1, tls[tl_i].from_2, selectedReadPairs, selectedReadPairsn);


		cc1 = computeCoverageSingleReads(tls[tl_i].chrRange, tls[tl_i].to_1, tls[tl_i].to_2, selectedReadPairs, selectedReadPairsn);
		tl_results[tl_i].meanCoverageTO = cc1.meanCoverage;
		tl_results[tl_i].continuousCoverageTOsingle = cc1.minimalCoverage;

		tl_results[tl_i].continuousCoverageTOpair = computeCoverageReadpairs(tls[tl_i].chrRange, tls[tl_i].to_1, tls[tl_i].to_2, selectedReadPairs, selectedReadPairsn);

	}

	return 1;
}


int SVs_class::saveINVresults(string filename)
{
	int i;

	if(invs == NULL && tls == NULL)
		return 0;

	ofstream out(filename.c_str());

	if(!out.is_open())
		return 0;

	//title line
	//out << "type\tnumbering\tchr\tBP1_1\tBP1_2\tBP2_1\tBP2_1\treadSupport\tcontCovSingleBP1\tcontCovPairBP1\tcontCovSingleBP2\tcontCovPairBP2\tmeanCovBP1\tmeanCovBP2\n";

	if(invs != NULL)
	{
		for(i=0;i<invsn;i++)
		{
			out << "INV" << i << "\t";
			out << invs[i].chr << "\t";
			out << invs[i].BP1_1 << "\t";
			out << invs[i].BP1_2 << "\t";
			out << invs[i].BP2_1 << "\t";
			out << invs[i].BP2_2;

			if(invs[i].results != NULL)
			{
				out << "\t";
				out << invs[i].results->readsSupporting << "\t";

				out << invs[i].results->continuousCoverageAsingle << "\t";
				out << invs[i].results->continuousCoverageApair << "\t";

				out << invs[i].results->continuousCoverageBsingle << "\t";
				out << invs[i].results->continuousCoverageBpair << "\t";

				out << invs[i].results->meanCoverageA << "\t";
				out << invs[i].results->meanCoverageB;

			}
			out << endl;
		}
	}


	out.close();
	out.clear();


	return 1;
}


int SVs_class::saveTLresults(string filename)
{
	int i;

	if(invs == NULL && tls == NULL)
		return 0;

	ofstream out(filename.c_str());

	if(!out.is_open())
		return 0;

	//title line
	//out << "type\tnumbering\tchrBP\tBP_1\tBP_2\tchrRange\tfrom_1\tfrom_2\tto_1\tto_2\t"
	//out << "readsSupportingFrom\treadsSupportingTo\treadsSupportingDeletion\treadsSupportingForward\treadsSupportingReverse\t"
	//out << "amountOfStacks\tamountOfOverlappingSupportReads\tbestStackratioFrom\tbestStackratioTo\t"
	//out << "continuousCoverageBPsingle\tcontinuousCoverageBPpair\tcontinuousCoverageFROMsingle\tcontinuousCoverageFROMpair\t"
	//out << "continuousCoverageTOsingle\tcontinuousCoverageTOpair\t"
	//out << "meanCoverageBP\tmeanCoverageFROM\tmeanCoverageTO\n"

	if(tls != NULL)
	{
		for(i=0;i<tlsn;i++)
		{
			out << "TL" << i << "\t";
			out << tls[i].chrBP << "\t";
			out << tls[i].BP_1 << "\t";
			out << tls[i].BP_2 << "\t";
			out << tls[i].chrRange << "\t";
			out << tls[i].from_1 << "\t";
			out << tls[i].from_2 << "\t";
			out << tls[i].to_1 << "\t";
			out << tls[i].to_2 << "\t";

			if(tls[i].results != NULL)
			{
				out << tl_results[i].readsSupportingFrom << "\t";
				out << tl_results[i].readsSupportingTo << "\t";
				out << tl_results[i].readsSupportingDeletion << "\t";
				out << tl_results[i].readsSupportingForward << "\t";
				out << tl_results[i].readsSupportingReverse << "\t";

				out << tl_results[i].amountOfStacks << "\t";
				out << tl_results[i].amountOfOverlappingSupportReads << "\t";
				out << tl_results[i].bestBPratioFrom << "\t";
				out << tl_results[i].bestBPratioTo << "\t";

				out << tl_results[i].continuousCoverageBPsingle << "\t";
				out << tl_results[i].continuousCoverageBPpair << "\t";
				out << tl_results[i].continuousCoverageFROMsingle << "\t";
				out << tl_results[i].continuousCoverageFROMpair << "\t";
				out << tl_results[i].continuousCoverageTOsingle << "\t";
				out << tl_results[i].continuousCoverageTOpair << "\t";

				out << tl_results[i].meanCoverageBP << "\t";
				out << tl_results[i].meanCoverageFROM << "\t";
				out << tl_results[i].meanCoverageTO << "\t";
				out << tl_results[i].bestBPC1 << "\t";
				out << tl_results[i].bestBPC2;

			}
			out << endl;
		}
	}


	out.close();
	out.clear();


	return 1;
}
