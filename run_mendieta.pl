#!/usr/bin/perl
use strict;
use warnings;

# define list of snapshots to process
my @snapList = (1 .. 100);

my $exec = "mendieta.x"
my $arg  = "mendieta.input"

# define (and create if needed) global output directory (which will host many sub-directories...)
my $runDir = "/mnt/sersic/marioagustin/saraswati/iden/idenfof/";
if (!-e $runDir) {system("mkdir $runDir");}

# copy code over to run directory, as a documentation ... 
if (!-e $runDir."/code") {system("mkdir $runDir"."/code");}
system("cp \*.c \*.x \*.h \*.input \*.pl Makefile $runDir"."/code/.");

# define executable file as the above copied executable (so we can link with no risk). 
my $execFile  = "$runDir"."/code/"."$exec";
my $cmdFile   = "$runDir" . "cmdFile.sh";    # instructions to run HM on all sub-boxes

open(CFILE,">$cmdFile");
# LOOP ON SNAPSHOTS
foreach (@snapList)
{
  my $i   = $_;
  my $num = sprintf("%03d",$i);

  my $dir = "$runDir"."$num";
  if(!-e $dir) {system("mkdir $dir")};

  my $nsnaps   = 32;
	my $snapdir  = "/mnt/raid/simulaciones/saraswati/outputs/";
	my $snapname = "snapshot_"."$num";
	my $foutfof  = "fof.bin";
	my $foutsub  = "sub.bin";
	my $soft     = 4.0;
	my $nfrac    = 9;

	# write input_HaloMaker.dat file
	write_input_hm($dir,$nsnaps,$snapdir,$snapname,$foutfof,$foutsub,$soft,$nfrac);

	# copy executable and stuff to running dir.
	my $cmd = "ln -s $execFile $dir"; 
	system($cmd);

	print CFILE "echo $num && cd $dir && ./$exec $arg > log.out \n";
}
close(CFILE);


sub write_input_hm {
    my $filename  = "$_[0]"."/mendieta.input";
    my $nsnaps    = $_[1];
    my $snapdir   = $_[2];
    my $snapname  = $_[3];
    my $foutfof   = $_[4];
    my $foutsub   = $_[4];
    my $softening = $_[5];
    my $nfrac     = $_[5];

    open(FILE,">$filename");

    print FILE "$nsnaps \n";
    print FILE "$snapdir \n";
    print FILE "$snapname \n";
    print FILE "$foutfof \n";
    print FILE "$foutsub \n";
    print FILE "$softening \n";
    print FILE "$nfrac \n";
    close(FILE);
}
