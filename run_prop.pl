#!/usr/bin/perl
use strict;
use warnings;

# define list of snapshots to process
my $nsnaps = 100;
my @snapList = (0);
for (my $k = 0; $k<$nsnaps; $k++) {$snapList[$k] = $k + 1;}

# define (and create if needed) global output directory (which will host many sub-directories...)
my $runDir = "/home/marioagustin/trabajo/SubhalosGoingNotts/major_mergers/RESULTS/mendieta/Dynamic/Dynamic/propiedades/";
if (!-e $runDir) {system("mkdir $runDir");}

# copy code over to run directory, as a documentation ... 
if (!-e $runDir."/code") {system("mkdir $runDir"."/code");}
system("cp \*.c \*.x \*.h \*.input $runDir"."/code/.");

# define executable file as the above copied executable (so we can link with no risk). 
my $execFile  = "$runDir"."/code/prop.x";
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
    my $cmd = "mkdir $dir";
    system($cmd);

    my $nsnaps    = 1;
		my $snapdir   =  "/home/marioagustin/trabajo/SubhalosGoingNotts/major_mergers/DATA/Dynamic/";
		my $snapname  = "major_merger_dynamic_"."$num";
		my $fofdir    = "/home/marioagustin/trabajo/SubhalosGoingNotts/major_mergers/RESULTS/mendieta/Dynamic/Dynamic/idensub/"."$num"."/";
		my $foffile   = "idensub.bin";
		my $fileout   = "propiedades.bin";
    my $filehijos = "hijos.ascii";
		my $softening = "0.7";

		# write input_HaloMaker.dat file
		write_input_hm($dir,$nsnaps,$snapdir,$snapname,$fofdir,$foffile,$fileout,$filehijos,$softening);

		# copy executable and stuff to running dir.
		$cmd = "ln -s $execFile $dir"; 
		system($cmd);

		print CFILE "cd $dir && echo $num && ./prop.x prop.input  > log.out \n";
}
close(CFILE);


sub write_input_hm {
    my $filename  = "$_[0]"."/prop.input";
    my $nsnaps    = $_[1];
    my $snapdir   = $_[2];
    my $snapname  = $_[3];
		my $fofdir    = $_[4];
		my $foffile   = $_[5];
    my $fileout   = $_[6];
    my $filehijos = $_[7];
		my $softening = $_[8];

    open(FILE,">$filename");

    print FILE "$nsnaps \n";
    print FILE "$snapdir \n";
    print FILE "$snapname \n";
		print FILE "$fofdir \n";
		print FILE "$foffile \n";
    print FILE "$fileout \n";
    print FILE "$filehijos \n";
		print FILE "$softening \n";
    close(FILE);
}
