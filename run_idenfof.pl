#!/usr/bin/perl
use strict;
use warnings;

# define list of snapshots to process
my @snapList = (1 .. 100);

# define (and create if needed) global output directory (which will host many sub-directories...)
my $runDir = "/mnt/sersic/marioagustin/saraswati/iden/idenfof/";
if (!-e $runDir) {system("mkdir $runDir");}

# copy code over to run directory, as a documentation ... 
if (!-e $runDir."/code") {system("mkdir $runDir"."/code");}
system("cp \*.c \*.x \*.h \*.input \*.pl Makefile $runDir"."/code/.");

# define executable file as the above copied executable (so we can link with no risk). 
my $execFile  = "$runDir"."/code/idenfof.x";
my $cmdFile   = "$runDir" . "cmdFile.sh";    # instructions to run HM on all sub-boxes

open(CFILE,">$cmdFile");
# LOOP ON SNAPSHOTS
foreach (@snapList)
{
  my $i   = $_;
  my $num = "$i";
  if ($i < 100)   {$num = "0"."$num";}
  if ($i < 10)    {$num = "0"."$num";}

  my $dir = "$runDir"."$num";
  if(!-e $dir) {system("mkdir $dir")};

  my $nsnaps   = 32;
	my $snapdir  = "/mnt/raid/simulaciones/saraswati/outputs/";
	my $snapname = "snapshot_"."$num";
	my $fileout  = "idenfof.bin";
	my $soft     = 4.0;

	# write input_HaloMaker.dat file
	write_input_hm($dir,$nsnaps,$snapdir,$snapname,$fileout,$soft);

	# copy executable and stuff to running dir.
	my $cmd = "ln -s $execFile $dir"; 
	system($cmd);

	print CFILE "echo $num && cd $dir && ./idenfof.x idenfof.input  > log.out \n";
}
close(CFILE);


sub write_input_hm {
    my $filename  = "$_[0]"."/idenfof.input";
    my $nsnaps    = $_[1];
    my $snapdir   = $_[2];
    my $snapname  = $_[3];
    my $fileout   = $_[4];
    my $softening = $_[5];

    open(FILE,">$filename");

    print FILE "$nsnaps \n";
    print FILE "$snapdir \n";
    print FILE "$snapname \n";
    print FILE "$fileout \n";
    print FILE "$softening \n";
    close(FILE);
}
