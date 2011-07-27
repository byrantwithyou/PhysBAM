##############################
# Data conversion script
# ----------------------------
# Usage: perl convert.pl [executable] [sourceDir] [targetDir]
# Arg1: name of the converter you want to use.
# Arg2: name of the folder in which the files to be converted reside.
# Arg3: name of the folder where you want the converted files to be stored.
# 
# Notes: Put this script in the same directory as dtof and ftod and run.
# The script will create the same directory structure as the source directory
# inside the target folder.  If the target folder does not exist, it
# will be created.
##############################

use File::Find;
use File::Stat;
use Fcntl ':mode';

# get arguments

if ($#ARGV != 2) {
	print "Usage: perl convert.pl [executable] [startDir] [endDir]\n";
	exit(1);
}

if (!(-e $ARGV[0])) {
	print "Cannot execute $ARGV[0] because it does not exist in the current directory\n";
	exit(1);
}

$prog = $ARGV[0];
$startDir = $ARGV[1];
$endDir = $ARGV[2];

# check to make sure startDir is a valid directory
if (!S_ISDIR((stat($startDir))[2])) {
	print "Invalid source directory.\n";
	exit(1);
}

# make target directory if it doesn't exist

if (!S_ISDIR((stat($endDir))[2])) {
	mkdir $endDir;
}

# make directory structure for target folder
print "setting up directory structure...\n";
find({ wanted=> \&makedirs, no_chdir => 1}, $startDir);
sub makedirs() {
	$File::Find::name =~ /($startDir\/*)(.*)/;
	$name = $2;
	$mode = (stat($File::Find::name))[2];
	$isDir = S_ISDIR($mode);
	if ($isDir) {
		$newPath = "$endDir/$name";
		mkdir "$newPath";
	}
}

# convert files

find({ wanted=> \&convert, no_chdir => 1}, $startDir);
sub convert() {
	$File::Find::name =~ /($startDir\/*)(.*)/;
	$name = $2;
	
	#if it's not a directory, execute the program
	if (!S_ISDIR((stat($File::Find::name))[2])) { 
		$newPath = "$endDir/$name";
		$cmdLine = "dataconv $prog $File::Find::name $newPath";
		print "\nExecuting: $cmdLine\n";
		print `$cmdLine`;
	}
}