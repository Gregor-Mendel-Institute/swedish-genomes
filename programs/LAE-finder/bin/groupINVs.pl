#!/usr/bin/perl


use warnings;
use strict;
local $| = 1;

use File::Find;


my $pattern = ".*.INVs.dat";
my $titleLine = "ID\tChr\tA1\tA2\tB1\tB2\t#reads supporting\t#continuousCovsingleA\t#continuousCovpairA\t#continuousCovsingleB\t#continuousCovpairB\tmeanCov A\tmeanCov B";


my $input_filename;

my $input_file;
my @input_lines;
my $input_linesn;
my @input_lines_split;
my $input_lines_splitn;


my @splitdata;
my $splitdatan;

my $joineddatas;


my $output_filename;
my $output_file;

my $i;
my $j;
my $temps;

my $id;

#my @files;
my @all_file_names;
my $all_file_namesn;

my %newFilesHash;

my %sums;
my $event;



find sub 
{
	return if -d;
	push @all_file_names, $File::Find::name;
}, '.';

$all_file_namesn = @all_file_names;

for($i=0;$i<$all_file_namesn;$i++)
{
	$temps = $all_file_names[$i];
	
	if($temps =~ /^$pattern$/)
	{
		print "at file ".$temps."\n";
		
		if($temps =~ /\/([^\/\.]+\.[^\/\.]+)\.[^\/]+$/)
		{	
			$id = $1;
		}
		else
		{
			print "wrong file namings\n";
		}
		
		if($all_file_names[$i] =~ /^.\//)
		{	$all_file_names[$i] =~	s/^.\///;	}

		$input_filename = $all_file_names[$i];
		@input_lines = ();
		open($input_file,'<'.$input_filename)  || die "Can't open".$input_filename;
		@input_lines = <$input_file>;
		$input_linesn = @input_lines;
		close($input_file);


		for($j=0;$j<$input_linesn;$j++)
		{
			if(length($input_lines[$j]) > 10)
			{
				@splitdata = split("\t", $input_lines[$j]);
				$splitdatan = @splitdata;
	
				$output_filename = $splitdata[0].".dat";
				
				if(not(exists $newFilesHash{$splitdata[0]}))
				{
					$newFilesHash{$splitdata[0]} = 1;
					open($output_file,'>'.$output_filename)  || die "Can't open".$output_filename;
					print($output_file $titleLine."\n");
					close($output_file);					
				}


				if(exists $sums{$splitdata[0]})
				{	if($splitdata[6] > 0)	{	$sums{$splitdata[0]} = 1;} else {	$sums{$splitdata[0]} = 1;	}}
				else
				{	if($splitdata[6] > 0)	{	$sums{$splitdata[0]}++;}}
				
				$splitdata[0] = $id;
				
				open($output_file,'>>'.$output_filename)  || die "Can't open".$output_filename;
				
				$joineddatas = join("\t", @splitdata);
				print($output_file $joineddatas);
				
				close($output_file);

			}
		}	
	}
}

print "\nEvents in :\n";

for $event ( keys %sums ) 
{
	if($sums{$event}) {	print $event."\n";	}
}


