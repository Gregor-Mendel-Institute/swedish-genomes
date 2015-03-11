
#include <iostream>


#include "general.h"
#include "SVs.h"
#include "samFile.h"


string delimiter = "/";

string outputDir = "output/";
//string outputDir = "c:\\Arbeit\\gmi\\AssemblyPipeline\\GivenSVs\\output\\";



SVs_class SVs;
SamFile samFile;


int main(int argc, char* argv[])
{
	int i, j;
	int tempi;

	string intermediateFilename;

	string command;
	string SVsFilename, dataFilename;
	string outputSVfilename, outputINVfilename;
	string prefix;
	string postfix = ".dat";

	if(argc != 4)
	{
		cout << "Usage:\n";
		//./LAE-finder nothing ../SVs_02.dat ../data/adal_1.9321.realigned.sort.bam.sam.lsort
		cout << "./LAE-finder <subcommand> <SV-list file> <sam file>\n";
		return 0;
	}
	else
	{
		command = argv[1];
		SVsFilename = argv[2];
		dataFilename = argv[3];
	}

	prefix = dataFilename;

	tempi = prefix.find_last_of(delimiter);
	if((unsigned)tempi != string::npos)
		prefix = prefix.erase(0, tempi + 1);

	tempi = prefix.find_last_of('.');
	if((unsigned)tempi != string::npos)
		prefix = prefix.erase(tempi);

	if(!SVs.loadFile(SVsFilename))
	{
		cout << "SV file missing or wrong\n";
		return -1;
	}


	if(command == "filter")
	{	//filter the reads for the SVs and save the result
		SVs.setMask();

		if(!samFile.selectReadPairs(dataFilename, SVs))
			return -1;

		intermediateFilename = dataFilename + ".sel";
		samFile.saveSelectedReadPairs(intermediateFilename);
	}
	else if(command == "filterAndSVcall")
	{	//filter and call
		SVs.setMask();

		if(!samFile.selectReadPairs(dataFilename, SVs))
			return -1;

		intermediateFilename = dataFilename + ".sel";
		samFile.saveSelectedReadPairs(intermediateFilename);

		//analysis:
		SVs.analyseINVs(samFile.selectedReadPairs, samFile.selectedReadPairsn);
		SVs.analyseTLs(samFile.selectedReadPairs, samFile.selectedReadPairsn);

		SVs.saveINVresults(outputDir + prefix + ".INVs.dat");
		SVs.saveTLresults(outputDir + prefix +".TLs.dat");
	}
	else if(command == "SVcall")
	{	//load the filtered reads and call
		samFile.loadSelectedReadPairs(dataFilename);

		//analysis:
		SVs.analyseINVs(samFile.selectedReadPairs, samFile.selectedReadPairsn);
		SVs.analyseTLs(samFile.selectedReadPairs, samFile.selectedReadPairsn);

		SVs.saveINVresults(outputDir + prefix + ".INVs.dat");
		SVs.saveTLresults(outputDir + prefix +".TLs.dat");
	}
	else
	{
		cout << "unknown command. Nothing done.\n";
	}


	//-------------------------------------------------------------------------
	//to generate:
	//-------------------------------------------------------------------------
	/*

	SVs.setMask();

	if(!samFile.selectReadPairs(samFilename, SVs))
		return -1;

	intermediateFilename = samFilename + ".sel";
	samFile.saveSelectedReadPairs(intermediateFilename);

	 //*/

	//-------------------------------------------------------------------------
	//to load:
/*
	samFile.loadSelectedReadPairs("c:\\Arbeit\\gmi\\AssemblyPipeline\\GivenSVs\\data\\adal_1.9321.realigned.sort.bam.sam.lsort.sel");

	//analysis:
	SVs.analyseINVs(samFile.selectedReadPairs, samFile.selectedReadPairsn);
	SVs.analyseTLs(samFile.selectedReadPairs, samFile.selectedReadPairsn);

	SVs.saveINVresults(outputDir + "adal_1.9321_INVs.dat");
	SVs.saveTLresults(outputDir + "adal_1.9321_TLs.dat");
	//*/

	return 0;
}

