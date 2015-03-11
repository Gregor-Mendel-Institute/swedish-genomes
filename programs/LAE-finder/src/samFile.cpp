/*
 * samFile.cpp
 *
 *  Created on: Mar 7, 2012
 *      Author: alexander.platzer
 */

#include "samFile.h"
#include "SVs.h"

#include <numeric>
#include <functional>



SamFile::SamFile()
{
	singleReadCount = 0;
	selectedReadPairs = NULL;
	selectedReadPairsn = 0;
}

SamFile::~SamFile()
{
	if(file.is_open())
	{
		file.close();
		file.clear();
	}

	if(selectedReadPairs != NULL)
		delete [] selectedReadPairs;
}


int SamFile::initFile(string filename)
{
	string s;

	if(file.is_open())
	{
		file.close();
		file.clear();
	}


	file.open(filename.c_str());

	if(!file.is_open())
	{
		printf("file error\n");
		return 0;
	}

	while(!file.eof())
	{
		getline(file, s);

		substrings.clear();
		split(s, "\t", substrings);

		if(s.c_str()[0] != '@')
			break;
	}

	if(file.eof())
		return 0;

	return 1;
}

int SamFile::internNextReadPair(void)
{
	string s;
	int found = 0;

	samReadPair.firstRead.chr = 0;
	samReadPair.secondRead.chr = 0;

	if(file == NULL)
		return 0;

	while(!file.eof())
	{
		getline(file, s);

		samReadPair.name = substrings[0];

		if(s.find(samReadPair.name) == 0 && substrings.size() >= 13)
		{
			samReadPair.firstRead.chr = substrings[2].c_str()[3] - '0';
			samReadPair.firstRead.cigar = substrings[5];
			samReadPair.firstRead.flag = atoi(substrings[1].c_str());
			samReadPair.firstRead.hitTag = substrings[12];
			samReadPair.firstRead.fullLine = accumulate( substrings.begin(), substrings.end(), string("\t") );;
			samReadPair.firstRead.leftmost = atoi(substrings[3].c_str());
			samReadPair.firstRead.quality = atoi(substrings[4].c_str());
			samReadPair.firstRead.length = substrings[9].length();

			substrings.clear();
			split(s, "\t", substrings);

			samReadPair.secondRead.chr = substrings[2].c_str()[3] - '0';
			samReadPair.secondRead.cigar = substrings[5];
			samReadPair.secondRead.flag = atoi(substrings[1].c_str());
			samReadPair.secondRead.hitTag = substrings[12];
			samReadPair.secondRead.fullLine = accumulate( substrings.begin(), substrings.end(), string("\t") );;
			samReadPair.secondRead.leftmost = atoi(substrings[3].c_str());
			samReadPair.secondRead.quality = atoi(substrings[4].c_str());
			samReadPair.secondRead.length = substrings[9].length();

			if(samReadPair.firstRead.leftmost > samReadPair.secondRead.leftmost)
				swap(samReadPair.firstRead, samReadPair.secondRead);

			getline(file, s);
			substrings.clear();
			split(s, "\t", substrings);

			found = 1;
			break;
		}
		else
		{
			singleReadCount++;
			substrings.clear();
			split(s, "\t", substrings);
		}
	}

	if(file.eof() && found == 0)
		return 0;

	return 1;
}

//looking for \d+M
int SamFile::cigarFilterpattern1(string &s)
{
	int i;

	if(s.c_str()[s.length() - 1] != 'M')
		return 0;

	for(i=0;i<(signed)s.length() - 1;i++)
		if(s.c_str()[i] < '0' || s.c_str()[i] > '9')
			return 0;

	return 1;
}

int SamFile::filterReadPair(void)
{
	if(!(samReadPair.firstRead.chr >= 1 && samReadPair.firstRead.chr <= 5))
		return 0;
	if(!(samReadPair.secondRead.chr >= 1 && samReadPair.secondRead.chr <= 5))
		return 0;

	if(!(samReadPair.firstRead.quality >= minimalQ && samReadPair.secondRead.quality >= minimalQ))
		return 0;

	if(!cigarFilterpattern1(samReadPair.firstRead.cigar) || !cigarFilterpattern1(samReadPair.secondRead.cigar))
		return 0;

	//if(samReadPair.firstRead.hitTag.find("XT:A:U") == string::npos || samReadPair.secondRead.hitTag.find("XT:A:U") == string::npos)
	if(samReadPair.firstRead.fullLine.find("XT:A:U") == string::npos || samReadPair.secondRead.fullLine.find("XT:A:U") == string::npos)
		return 0;

	return 1;
}

int SamFile::toNextReadPair(void)
{
	int flag1, flag2;

	while(internNextReadPair())
	{
		if(filterReadPair())
		{
			currentReadPairShort.chr1 = samReadPair.firstRead.chr;
			currentReadPairShort.chr2 = samReadPair.secondRead.chr;

			currentReadPairShort.leftmost1 = samReadPair.firstRead.leftmost;
			currentReadPairShort.leftmost2 = samReadPair.secondRead.leftmost;

			if(samReadPair.firstRead.length != samReadPair.secondRead.length)
				cout << "unequal read lengths!\n";

			currentReadPairShort.length = samReadPair.firstRead.length;

			flag1 = (samReadPair.firstRead.flag / 16) % 2;
			flag2 = (samReadPair.secondRead.flag / 16) % 2;

			if((flag1 + flag2 + 1) % 2 == 0)
				currentReadPairShort.orientation = FORWARD;
			else
				currentReadPairShort.orientation = REVERSE;

			return 1;
		}
	}

	return 0;
}


int SamFile::selectReadPairs(string filename, SVs_class & SVs)
{
	int i, j;
	char * chrFirstReadMask = NULL;
	char * chrSecondReadMask = NULL;

	if(selectedReadPairs != NULL)
		delete [] selectedReadPairs;

	selectedReadPairs = NULL;
	selectedReadPairsn = 0;

	if(!initFile(filename))
	{
		cout << "file missing: " << filename << endl;
		return 0;
	}

	cout << "counting----------------------\n";

	j = 0;
	while(toNextReadPair())
	{
		j++;
		if(j %10000 == 0)
		{	cout << j << " readPairs\n";	cout.flush();	}

		switch(currentReadPairShort.chr1)
		{
			case 1:	chrFirstReadMask = SVs.chr1mask;
				break;
			case 2:	chrFirstReadMask = SVs.chr2mask;
				break;
			case 3:	chrFirstReadMask = SVs.chr3mask;
				break;
			case 4:	chrFirstReadMask = SVs.chr4mask;
				break;
			case 5:	chrFirstReadMask = SVs.chr5mask;
				break;
		}
		switch(currentReadPairShort.chr2)
		{
			case 1:	chrSecondReadMask = SVs.chr1mask;
				break;
			case 2:	chrSecondReadMask = SVs.chr2mask;
				break;
			case 3:	chrSecondReadMask = SVs.chr3mask;
				break;
			case 4:	chrSecondReadMask = SVs.chr4mask;
				break;
			case 5:	chrSecondReadMask = SVs.chr5mask;
				break;
		}

		if(chrFirstReadMask[currentReadPairShort.leftmost1] || chrFirstReadMask[currentReadPairShort.leftmost1 + currentReadPairShort.length] ||
		   chrSecondReadMask[currentReadPairShort.leftmost2] || chrSecondReadMask[currentReadPairShort.leftmost2 + currentReadPairShort.length])
			selectedReadPairsn++;
	}

	cout << selectedReadPairsn << "\n";

	selectedReadPairs = new ReadPairShort[selectedReadPairsn];

	initFile(filename);

	cout << "reading----------------------\n";
	i = 0;
	j = 0;
	while(toNextReadPair())
	{
		j++;
		if(j %10000 == 0)
		{	cout << j << " readPairs\n";	cout.flush();	}

		switch(currentReadPairShort.chr1)
		{
			case 1:	chrFirstReadMask = SVs.chr1mask;
				break;
			case 2:	chrFirstReadMask = SVs.chr2mask;
				break;
			case 3:	chrFirstReadMask = SVs.chr3mask;
				break;
			case 4:	chrFirstReadMask = SVs.chr4mask;
				break;
			case 5:	chrFirstReadMask = SVs.chr5mask;
				break;
		}
		switch(currentReadPairShort.chr2)
		{
			case 1:	chrSecondReadMask = SVs.chr1mask;
				break;
			case 2:	chrSecondReadMask = SVs.chr2mask;
				break;
			case 3:	chrSecondReadMask = SVs.chr3mask;
				break;
			case 4:	chrSecondReadMask = SVs.chr4mask;
				break;
			case 5:	chrSecondReadMask = SVs.chr5mask;
				break;
		}

		if(chrFirstReadMask[currentReadPairShort.leftmost1] || chrFirstReadMask[currentReadPairShort.leftmost1 + currentReadPairShort.length] ||
		   chrSecondReadMask[currentReadPairShort.leftmost2] || chrSecondReadMask[currentReadPairShort.leftmost2 + currentReadPairShort.length])
		{
			selectedReadPairs[i] = currentReadPairShort;
			i++;
		}
	}

	return 1;
}


int SamFile::saveSelectedReadPairs(string filename)
{
	int i;

	if(selectedReadPairs == NULL)
		return 0;

	ofstream out(filename.c_str());

	if(!out.is_open())
		return 0;

	out << "Chr1\tleftmost1\tChr2\tleftmost2\tlength\torientation\n";

	for(i=0;i<selectedReadPairsn;i++)
	{
		out << selectedReadPairs[i].chr1 << "\t";
		out << selectedReadPairs[i].leftmost1 << "\t";
		out << selectedReadPairs[i].chr2 << "\t";
		out << selectedReadPairs[i].leftmost2 << "\t";
		out << selectedReadPairs[i].length << "\t";
		out << selectedReadPairs[i].orientation << "\n";
	}

	out.close();
	out.clear();


	return 1;
}

int SamFile::loadSelectedReadPairs(string filename)
{
	int i, j, tempi;
	string s;
	stringstream sx1;

	if(selectedReadPairs != NULL)
		delete [] selectedReadPairs;

	ifstream in(filename.c_str());
	if(!in.is_open())
		return 0;

	getline(in, s);
	i = 0;
	while(!in.eof())
	{
		getline(in, s);
		if(s.length() > 3)
			i++;
	}
	selectedReadPairsn = i;

	selectedReadPairs = new ReadPairShort[selectedReadPairsn];

	in.close();
	in.clear();

	in.open(filename.c_str());

	getline(in, s);
	i = 0;
	while(!in.eof())
	{
		getline(in, s);

		if(s.length() > 3)
		{
			substrings.clear();
			split(s, "\t", substrings);

			for(j=0;j<(signed)substrings.size();j++)
			{
				sx1.str("");
				sx1.clear();

				sx1 << substrings[j];

				switch(j)
				{
					case 0:	sx1 >> selectedReadPairs[i].chr1;
						break;
					case 1:	sx1 >> selectedReadPairs[i].leftmost1;
						break;
					case 2:	sx1 >> selectedReadPairs[i].chr2;
						break;
					case 3:	sx1 >> selectedReadPairs[i].leftmost2;

							if(selectedReadPairs[i].leftmost2 < selectedReadPairs[i].leftmost1)
							{
								swap(selectedReadPairs[i].leftmost2,  selectedReadPairs[i].leftmost1);
								swap(selectedReadPairs[i].chr2, selectedReadPairs[i].chr1);
							}

						break;
					case 4:	sx1 >> selectedReadPairs[i].length;
						break;
					case 5:	sx1 >> tempi;
							if(tempi == 0)
								selectedReadPairs[i].orientation = FORWARD;
							else
								selectedReadPairs[i].orientation = REVERSE;
						break;
				}
			}

			i++;
		}
	}


	in.close();
	in.clear();


	return 1;
}

