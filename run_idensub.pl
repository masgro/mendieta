#!/usr/bin/perl
use strict;
use warnings;

# define list of snapshots to process
my $nsnaps = 10;
my @snapList = (0);
for (my $k = 0; $k<$nsnaps; $k = $k + 1) {$snapList[$k] = $k + 1;}

# define (and create if needed) global output directory (which will host many sub-directories...)
#my $runDir = "/home/marioagustin/trabajo/SubhalosGoingNotts/major_mergers/RESULTS/mendieta/Dynamic/Dynamic/idensub/";
my $runDir = "/home/marioagustin/trabajo/SubhalosGoingNotts/major_mergers/RESULTS/mendieta/Static/idensub/";
if (!-e $runDir) {system("mkdir $runDir");}

# copy code over to run directory, as a documentation ... 
if (!-e $runDir."/code") {system("mkdir $runDir"."/code");}
system("cp Makefile \*.c \*.x \*.h \*.input $runDir"."/code/.");

# define executable file as the above copied executable (so we can link with no risk). 
my $execFile  = "$runDir"."code/idensub.x";
my $cmdFile   = "$runDir" . "cmdFile.sh";    # instructions to run HM on all sub-boxes

open(CFILE,">$cmdFile");
# LOOP ON SNAPSHOTS
foreach (@snapList)
  {
    my $i   = $_;
    my $num = "$i";
    #if ($i < 100)   {$num = "0"."$num";}
    #if ($i < 10)    {$num = "0"."$num";}

    my $dir = "$runDir"."$num";
    my $cmd = "mkdir $dir";
    system($cmd);

    my $nsnaps   = 1;
		#my $snapdir  = "/home/marioagustin/trabajo/SubhalosGoingNotts/major_mergers/DATA/Dynamic/";
		my $snapdir  = "/home/marioagustin/trabajo/SubhalosGoingNotts/major_mergers/DATA/Static/";
		#my $snapname = "major_merger_dynamic_"."$num";
		my $snapname = "static_merger_".$num.".gadget";
		#my $fofdir   = "/home/marioagustin/trabajo/SubhalosGoingNotts/major_mergers/RESULTS/mendieta/Dynamic/Dynamic/idenfof/"."$num"."/";
		my $fofdir   = "/home/marioagustin/trabajo/SubhalosGoingNotts/major_mergers/RESULTS/mendieta/Static/idenfof/"."$num"."/";
		my $foffile   = "idenfof.bin";
		my $fileout   = "idensub.bin";
		my $filefof   = "sub2fof.bin";
		#my $softening = "0.7";
		my $softening = "3.0";
		my $nsteps    = "9";

		# write input_HaloMaker.dat file
		write_input_hm($dir,$nsnaps,$snapdir,$snapname,$fofdir,$foffile,$fileout,$filefof,$softening,$nsteps);

		# copy executable and stuff to running dir.
		$cmd = "ln -s $execFile $dir"; 
		system($cmd);

		print CFILE "cd $dir && echo $num && ./idensub.x idensub.input  > log.out \n";
}
close(CFILE);


sub write_input_hm {
    my $filename  = "$_[0]"."/idensub.input";
    my $nsnaps    = $_[1];
    my $snapdir   = $_[2];
    my $snapname  = $_[3];
		my $fofdir    = $_[4];
		my $foffile   = $_[5];
    my $fileout   = $_[6];
		my $filefof   = $_[7];
		my $softening = $_[8];
		my $nsteps    = $_[9];

    open(FILE,">$filename");

    print FILE "$nsnaps \n";
    print FILE "$snapdir \n";
    print FILE "$snapname \n";
		print FILE "$fofdir \n";
		print FILE "$foffile \n";
    print FILE "$fileout \n";
    print FILE "$filefof \n";
		print FILE "$softening \n";
		print FILE "$nsteps \n";
    close(FILE);
}
