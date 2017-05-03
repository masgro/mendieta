#!/usr/bin/perl
use strict;
use warnings;

# define list of snapshots to process
my @snapList = (1 .. 100);

# define (and create if needed) global output directory
# (which will host many sub-directories...)
my $runDir = "pruebarun/";
if (!-e $runDir) {system("mkdir $runDir");}

# copy code over to run directory, as a documentation
if (!-e $runDir."/code") {system("mkdir $runDir"."/code");}
system("cp \*.c \*.x \*.h \*.pl Makefile $runDir"."/code/.");

# define executable file as the above copied executable (so we can link with no risk). 
my $execFile  = "mendieta.x";
my $cmdFile   = "$runDir" . "cmdFile.sh";    # instructions to run HM on all sub-boxes

my @xcm  = (125000.0,125000.0,125000.0,125000.0,375000.0,375000.0,375000.0,375000.0);
my @ycm  = (125000.0,125000.0,375000.0,375000.0,125000.0,125000.0,375000.0,375000.0);
my @zcm  = (125000.0,375000.0,125000.0,375000.0,125000.0,375000.0,125000.0,375000.0);
my @lado = (125000.0,125000.0,125000.0,125000.0,125000.0,125000.0,125000.0,125000.0);

open(CFILE,">$cmdFile");

# LOOP ON SNAPSHOTS
foreach (@snapList)
{
  my $i   = $_;
	my $num = sprintf("%03d",$i);

  my $dir = "$runDir"."$num";
  if(!-e $dir) {system("mkdir $dir")};

	my $file = "$dir"."/run_boxes.sh";
	open(BFILE,">$file");

	for(my $count=1;$count<=8;$count++){

		my $box = sprintf "$dir/%1d", $count;
		system("mkdir $box");

    my $nsnaps   = 32;
		my $snapdir  = "/mnt/raid/simulaciones/saraswati/outputs/";
		my $snapname = "snapshot_"."$num";
		my $fileout  = "idenfof.bin";

		# write input_HaloMaker.dat file
		write_input_hm($box,$nsnaps,$snapdir,$snapname,$fileout);

		# copy executable and stuff to running dir.
		my $cmd = "ln -s "."$runDir"."/code/"."$execFile"." $box"; 
		system($cmd);

		print BFILE "date && echo $box && cd $box && nice -19 ./$execFile idenfof.input $xcm[$count-1] $ycm[$count-1] $zcm[$count-1] $lado[$count-1] > log.out \n";
	}
	close(BFILE);

	print CFILE "cd $dir && sh run_boxes.sh\n";

}
close(CFILE);

sub write_input_hm {
    my $filename  = "$_[0]"."/idenfof.input";
    my $nsnaps    = $_[1];
    my $snapdir   = $_[2];
    my $snapname  = $_[3];
    my $fileout   = $_[4];
    my $softening = "9.";

    open(FILE,">$filename");

    print FILE "$nsnaps \n";
    print FILE "$snapdir \n";
    print FILE "$snapname \n";
    print FILE "$fileout \n";
    print FILE "$softening \n";
    close(FILE);
}
