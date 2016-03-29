#!/usr/bin/perl -w
use strict;

my $fasta_fn;

my %sbencode = (
		0=>'A',
		1=>'C',
		2=>'G',
		3=>'T',
		);

my %encode = (
    AA=>0,
    CC=>0,
    GG=>0,
    TT=>0,
    AC=>1,
    CA=>1,
    GT=>1,
    TG=>1,
    AG=>2,
    CT=>2,
    GA=>2,
    TC=>2,
    AT=>3,
    CG=>3,
    GC=>3,
    TA=>3,
	      );
my %decode = (
    A0=>'A',
    A1=>'C',
    A2=>'G',
    A3=>'T',
    C0=>'C',
    C1=>'A',
    C2=>'T',
    C3=>'G',
    G0=>'G',
    G1=>'T',
    G2=>'A',
    G3=>'C',
    T0=>'T',
    T1=>'G',
    T2=>'C',
    T3=>'A',
	      );
	

#if the input contains a number, then it is assumed to be a color string
#to convert into sequence.  Otherwise it is taken to be a sequence to turn into
#a color string

sub interactive {
    while(<STDIN>){
	chomp;
        my $line=$_;
        if($line=~/\r$/){chop $line;}

	my ($tmp,$cstr);
        my $bstr=uc $_;

	$cstr="";
        if( $line =~ /^>/){
            print "$line\n";
        }else{
	  if($bstr =~/\d+/){
	    #must be a color string
	    if($bstr =~ /^([ACTG])/){
		my $pbase=$1;
		for(my $i=1; $i<length($bstr); $i++){
		    my $tcolor=substr $bstr,$i,1;
		    $tmp=$decode{"$pbase$tcolor"};
		    $pbase=$tmp;
		    substr($cstr, length($cstr))=$tmp;
		}
		print "$cstr\n";
	    }else{
		for(my $i=0; $i<4; $i++){
		    my $pbase=$sbencode{"$i"};
		    $cstr=$pbase;
		    for(my $i=0; $i<length($bstr); $i++){
			my $tcolor=substr $bstr,$i,1;
			$tmp=$decode{"$pbase$tcolor"};
			$pbase=$tmp;
			substr($cstr, length($cstr))=$tmp;
		    }
		    print " $cstr\n";
		}
	    }
	  }else{
	    substr($cstr, length($cstr))=substr $bstr,0,1;
	    for(my $i=1; $i<length($bstr); $i++){
		$tmp=$encode{substr $bstr,$i-1,2};
		substr($cstr, length($cstr))=$tmp;
	    }
	    print " $cstr\n>";
	  }
	}
    }
}

sub xlate_fasta {
    open INFILE, "$fasta_fn";
    my ($tmp,$cstr);
    my $lcnt=0;
    my $prev_base;
    $cstr="";
    while(<INFILE>){
        my $line=$_;
        chomp $line;
	if($line=~/\r$/){chop $line;}
        if( $line =~ /^>/){
	    print "$line\n";
	    $prev_base="."; # empty
        }else{
	    if($prev_base eq "."){
		$prev_base=substr($line,0,1);
		$line=substr($line,1);
	    }
	    $tmp=substr($line,0,1);
	    $cstr=$encode{"$prev_base$tmp"};
	    $prev_base=substr($line,-1,1);
            for(my $i=1; $i<length($line); $i++){
                $tmp=$encode{substr $line,$i-1,2};
                substr($cstr, length($cstr))=$tmp;
            }
            print "$cstr\n";
        }
    }
    close(INFILE);
}


sub parse_args{
    for(my $i=0; $i<@ARGV; $i++){
        if($ARGV[$i] eq "-f"){$fasta_fn=$ARGV[++$i];}
    }
}


# MAIN starts here
 
$fasta_fn="";
parse_args();
if($fasta_fn eq ""){
    interactive();
}else{
    xlate_fasta();
}

