#!/usr/bin/env perl
# (C) Copyright 2017 Sur Herrera Paredes
# 
# This file is part of wheelP.
# 
# wheelP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# wheelP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with wheelP.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;

# Usage:
#
#	$ percent_map.pl <infile.uc> <outfile.txt>

my $dir = shift @ARGV;
my $outfile = shift @ARGV;

if(-d $dir){
	my $subdirs_ref = get_subdirs($dir);
	my $res_ref = process_ucfiles('ref_table.uc',$subdirs_ref, $dir);
	print_res($res_ref, $outfile);
}elsif(-f $dir){
	my $subdirs_ref = [ '' ];
	my $res_ref = process_ucfiles($dir,$subdirs_ref,'');
	print_res($res_ref, $outfile);
}else{
	die "1:\nHola\n";;
}

##### subroutines
sub get_subdirs{
	my ($dir) = @_;
	opendir(DIR,$dir) or die "Can't open $dir ($!)";
	my @subdirs = grep{$_ ne '.' & $_ ne '..' &  -d "$dir/$_"} readdir DIR;
	closedir DIR;

	#print "$_\n" foreach @subdirs;

	return(\@subdirs);
}

sub process_ucfiles{
	my ($filename,$dirs_ref,$dir) = @_;

	my ($file,%Res);
	for $file (@$dirs_ref){
		my $ucfile = "$dir/$file/$filename";
		open(my $IN,$ucfile) or die "Can't open $ucfile ($!)";
		print "Processing $ucfile\n";
		uc_stats($IN,\%Res);
		close $IN;
	}

	return(\%Res);
}

sub uc_stats{
	my ($fh, $res_ref) = @_;

	my($map,$hit, $total,$read_id,$sample);
	$total = $hit = 0;
	#print "Hello\n";
	while(<$fh>){
		chomp;
		#print "Hola\n";
		next if $_ =~ /^#/;
		my (@line) = split(/\t/,$_);
		$map = $line[0];
		@line = split(/_/,$line[8]);
		$sample = $line[0];
		#print "$sample\n";

		if(exists($res_ref->{$sample})){
			$res_ref->{$sample}->{'Total'}++;
			$res_ref->{$sample}->{'Hit'}++ if $map eq 'H';
		}else{
			$res_ref->{$sample} = {'Total' => 1, 'Hit' => 0};
			$res_ref->{$sample}->{'Hit'}++ if $map eq 'H';
		}

		$hit++ if $map eq 'H';
		$total++;
		#last if $total > 100;
	}


	#my $key;
	#for $key (keys %{$res_ref}){
	#	print "\t$key\t$res_ref->{$key}->{Total}\t$res_ref->{$key}->{Hit}\n";
	#}

	print "\tProcessed $total reads...\n";
	print "\t$hit reads mapped...\n";
	print "\t". 100*$hit/$total . "% of reads mapped...\n";
}

sub print_res{
	my ($res_ref, $outfile) = @_;

	open(OUT,'>',$outfile) or die "Can't create $outfile ($!)";
	print "Printing mapping results\n";
	my ($key);
	for $key (keys %{$res_ref}){
		print OUT "$key\t$res_ref->{$key}->{Total}\t$res_ref->{$key}->{Hit}\t" . 100*$res_ref->{$key}->{Hit} / $res_ref->{$key}->{Total} . "\n";
	}
	close OUT;
}


