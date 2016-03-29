#!/usr/bin/perl -w
use strict;

# takes a (v,k,t) design from Dan Gordon's site http://www.ccrwest.org/cover/HIGH.html and forms a schema from it

my @lines;
my $maxindex=0;
my $numlines=0;
my($next_arg, $design_file, $output_directory);

sub usage {
    warn "\nUsage: $0 \n\n";
    warn "\t -d <design_file> \n";
    warn "\t -dir <output_directory> \n\n";
    exit(1);
}
if(scalar(@ARGV) == 0){
    usage();
}

# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-d"){ $design_file = shift(@ARGV); }
    elsif($next_arg eq "-dir"){ $output_directory = shift(@ARGV); }
    else { warn "\nInvalid argument: $next_arg \n"; usage(); }
}

# Verify Input
unless(defined $design_file && -e $design_file) {
    warn "\nERROR: design file $design_file does not exist \n";
    usage();
}
# output_directory
unless(defined $output_directory) {
    warn "ERROR: No output_directory specified\n";
    usage();
}
unless(-e $output_directory && -d $output_directory) {
    mkdir($output_directory, 0775);
    warn "output_directory $output_directory did not exist so it was created \n"; 
}

open( FILE, "< $design_file" ) or die "Can't open $design_file : $!";
my @entries = ();
while( <FILE> ) {
    chomp;
    my $line=$_;
    chomp $line;
    $line =~ s/^\s+//;
    if ($line !~ /\#/){
	push @lines, $line;
	@entries=split /\s+/,$line ;
	foreach my $i (@entries) {
	    if ($i>$maxindex){$maxindex=$i;}
	}
	$numlines++;
    }
}

my @values = split(/_/, $design_file);
my $tag_length = $values[1];
my $number_of_ones = $tag_length - $values[2];
my $number_of_mismatches = $values[3];
my $number_of_rows = scalar @lines;

&make_pattern($tag_length, $number_of_ones, $number_of_mismatches, $number_of_rows);

exit;

sub make_pattern {
    my($tag_length, $number_of_ones, $number_of_mismatches, $number_of_rows) = @_;
    my @entries;
    my($pattern, $line);
    my $j;
    my $outputfile = "$output_directory/schema_$tag_length" . "_$number_of_mismatches";
#    $outputfile .= "_$number_of_ones";
#    $outputfile .= "_$number_of_rows";
    unless ( open(SCHEMA, ">$outputfile") ) {
	print "Cannot open file \"$outputfile\" to write to!!\n\n";
	exit;
    }
    print SCHEMA "# $number_of_ones base index on $tag_length, $number_of_mismatches mismatches \n";
    foreach $line (@lines) {
	$pattern="";
	for($j=0; $j<$maxindex; $j++){
	    $pattern="1$pattern";
	}
	@entries=split (/\s+/,$line);
	foreach my $i (@entries) {
	    substr($pattern,$i-1,1)="0";
	}
	# Difference between Alan and Zheng's schemas
	$pattern = substr($pattern, 1, length($pattern) - 1);
	print SCHEMA"$pattern\n";
    }
    close(SCHEMA);
}
