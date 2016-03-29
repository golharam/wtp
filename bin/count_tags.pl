
#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;
use FileHandle;


# extract read_alignment_type (unique,all,top_score,random) matched tag from a max file 
# assume the reference genome used for max file is in a multi fasta format
# use an annotation gff file containing contiguous regions of interest (exons)
# assume strand info in gff annotation file is one of '+'/'-'/'.'
# use reads aligning inside (start/end) annotated regions +/- 3bp, on the annotated strand 
# (or both strands if annotated starnd is .)
# generate gff format file with counts for each annotated region +/- 3bp
# generate coverage file (wig) for each reference sequence
# converte max to gff 
# this script represents generate_output_v0.8.pl
my ($max_file, $is_sorted, $reference_fasta_file, $reference_info_file, $tag_length, $read_alignment_type, $min_score, $countsoutputfile, $wigoutputdir, $gffreadsfile, $produce_count, $produce_wig_genome, $coverage_filter_threshold, $split_coverage_by_chromosome, $produce_gff_reads, $next_arg );

if(scalar(@ARGV) == 0){
    print "\nUsage: $0  \n\t";
    
    print "-max <max_file> \n\t";
    print "-s <is_sorted: is the max file sorted: \'yes\' or \'no\'; default no> \n\t";
    print "-r <reference_fasta_file> \n\t";
#    print "-gff <reference_info_file> \n\t";
	print "-a <annotation_gtf_file> \n\t";
    print "-t <tag_length> \n\t";
    print "-read_alignment_type <read_alignment_type: \'top_score\' \'unique\' \'all\' or \'random\'> \n\t";
    print "-min_score <min_score: only reads with aligment score at or above this value are reported> \n\t";
    print "-counts <countsoutputfile> \n\t";
    print "-wig <wigoutputdir> \n\t";
    print "-gff_reads <gffreadsfile> \n\t";
    print "-produce_count <\'yes\' or \'no\': default yes> \n\t";
    print "-produce_wig <\'yes\' or \'no\': default yes> \n\t";
    print "-coverage_filter_threshold <default 10> \n\t";
    print "-split_coverage_by_chromosome <\'yes\' or \'no\': default no> \n\t";
    print "-produce_gff_reads <\'yes\' or \'no\': default no>\n\n";
    
    exit(1);
}

# Parse the command line

$is_sorted = "no";
$read_alignment_type = "unique";
$produce_count = "yes";
$produce_wig_genome = "yes";
$coverage_filter_threshold = 10;
$split_coverage_by_chromosome = "no";
$produce_gff_reads = "no";

# residual from previous implementation; removed 01.28.2009
my $produce_wig = "no";

while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-max"){ $max_file = shift(@ARGV); }
    elsif($next_arg eq "-s"){ $is_sorted = shift(@ARGV); }
    elsif($next_arg eq "-r"){ $reference_fasta_file = shift(@ARGV); }
#    elsif($next_arg eq "-gff"){ $reference_info_file = shift(@ARGV); }
	elsif($next_arg eq "-a") { $reference_info_file = shift(@ARGV); }
    elsif($next_arg eq "-t"){ $tag_length = shift(@ARGV); }
    elsif($next_arg eq "-read_alignment_type"){ $read_alignment_type = shift(@ARGV); }
    elsif($next_arg eq "-min_score"){ $min_score = shift(@ARGV); }
    elsif($next_arg eq "-counts"){ $countsoutputfile = shift(@ARGV); }
    elsif($next_arg eq "-wig"){ $wigoutputdir = shift(@ARGV); }
    elsif($next_arg eq "-gff_reads"){ $gffreadsfile = shift(@ARGV); }
    elsif($next_arg eq "-produce_count"){ $produce_count = shift(@ARGV); }
    elsif($next_arg eq "-produce_wig_genome"){ $produce_wig_genome = shift(@ARGV); }
    elsif($next_arg eq "-coverage_filter_threshold"){ $coverage_filter_threshold = shift(@ARGV); }
    elsif($next_arg eq "-split_coverage_by_chromosome"){ $split_coverage_by_chromosome = shift(@ARGV); }
    elsif($next_arg eq "-produce_gff_reads"){ $produce_gff_reads = shift(@ARGV); }
}


&check_options($is_sorted, $read_alignment_type, $produce_count, $produce_wig_genome, $produce_gff_reads);
&check_optional_dir($produce_wig_genome, $wigoutputdir);
&check_files($max_file, $reference_fasta_file, $reference_info_file);
&check_optional_file($produce_count, $countsoutputfile);
&check_optional_file($produce_gff_reads, $gffreadsfile);
my $base_max_file_name = get_base_max_file_name($max_file);

my ($date,$time) = &time_stamp();
print "Start loading reference: time $date, $time\n";

my $window_reference_size = 50;
my %info_map = map_info_reference($reference_info_file, $window_reference_size);
my @info = get_reference_info($reference_info_file);
my @names = get_ref_sequence_names($reference_fasta_file);

if((scalar @names > 50) & ($split_coverage_by_chromosome eq "yes") & ($is_sorted eq "no")) {
    print "WARNING\:Options \n\t\-split_coverage_by_chromosome\tyes\n\t\-s\tno\nare not implemented for more than 50 sequences in the reference file\.\nPlease use sorted max file with option \-s\tyes\n\n";
    exit;
}


($date,$time) = &time_stamp();
print "Start reading max file: time $date, $time\n";

my ($i,$j,$k, $index_max_score, $temp_max_score, $info_map_temp_key, $temp_ref_seq_index, $temp_coord_start, $temp_coord_end, $temp_start, $temp_end, $temp_score, $temp_size, $temp_index_gff, $temp_strand, $wig_genome_positive_out_file, $wig_genome_negative_out_file, $wig_annotation_out_file, $number_of_keys, $temp_line, $temp_coordinate);
my (@values, @individual_values, @individual_matching_info, @temp_location, @keys, @keys_sorted, @gff_values);

my @counts;
my @reference_end;
my (%coverage_positive, %coverage_negative, %coverage);

my $count_match = 0;
my $number_of_reads = 0;
my @index_matching_to_use;
my $current_chr = 1;
my $current_number_of_bases = 0;
my $max_number_of_bases_to_keep_in_memory = 1e6;
my $temp_annotated_call = 0;

my $temp_file_name;
my (@unsort_positive_coordinate_files, @unsort_negative_coordinate_files);
my (@sort_positive_coordinate_files, @sort_negative_coordinate_files);
for($i = 0; $i < scalar @names; $i++){
    #if($split_coverage_by_chromosome eq "yes" ) {
	$temp_file_name = "$wigoutputdir\/$names[$i]\.unsorted\_positive\_matching\_locations\.txt";
	push @unsort_positive_coordinate_files, $temp_file_name;
	$temp_file_name = "$wigoutputdir\/$names[$i]\.unsorted\_negative\_matching\_locations\.txt";
	push @unsort_negative_coordinate_files, $temp_file_name;
    
	$temp_file_name = "$wigoutputdir\/$names[$i]\.sorted\_positive_matching\_locations\.txt";
	push @sort_positive_coordinate_files, $temp_file_name;
	$temp_file_name = "$wigoutputdir\/$names[$i]\.sorted\_negative_matching\_locations\.txt";
	push @sort_negative_coordinate_files, $temp_file_name;
    #}    
    $reference_end[$i] = 0;
}

# array of file handlers for unsorted coordinates files

my @unsort_positive_file_handles = ();
my @unsort_negative_file_handles = ();
if(($is_sorted eq "no")) {
    if($produce_wig_genome eq "yes" ) {
	#if($split_coverage_by_chromosome eq "yes") {
	    for($i = 0; $i < scalar @names; $i++){    
		# localize the file glob, so FILE is unique to the inner loop.
		open(my $fhp, ">", $unsort_positive_coordinate_files[$i]) or die "Can't open file \"$unsort_positive_coordinate_files[$i]\" to write to!!\n\n";

		# push the typeglobe to the end of the array
		push(@unsort_positive_file_handles, $fhp);
		#print "unsorted positive coordinates chr$i are saved in $unsort_positive_file_handles[$i], \n";

		# same for negative strand.
		open(my $fhn, ">", $unsort_negative_coordinate_files[$i]) or die "Can't open file \"$unsort_negative_coordinate_files[$i]\" to write to!!\n\n";

		# same for negative strand
		push(@unsort_negative_file_handles, $fhn);
		#print "unsorted negative coordinates chr$i are saved in $unsort_negative_file_handles[$i], \n";
	    }
	#}
    }
}

if ($is_sorted eq "yes") {

    if($produce_wig_genome eq "yes") {
	if($split_coverage_by_chromosome eq "yes") {
	    $wig_genome_positive_out_file = "$wigoutputdir\/$names[$current_chr].positive_strand.wig";
	    $wig_genome_negative_out_file = "$wigoutputdir\/$names[$current_chr].negative_strand.wig";
	}
	else {
	    $wig_genome_positive_out_file = "$wigoutputdir\/genome.positive_strand.wig";
	    $wig_genome_negative_out_file = "$wigoutputdir\/genome.negative_strand.wig";	
	}
	print "Start printing $wig_genome_positive_out_file...\n";
	unless ( open(OUT_GENOME_POSITIVE_WIG_FILE, ">$wig_genome_positive_out_file") ) {
	    print "Cannot open file \"$wig_genome_positive_out_file\" to write to!!\n\n";	
	    exit;
	}
	print OUT_GENOME_POSITIVE_WIG_FILE "browser position $names[$current_chr]:1-200000000\nbrowser hide all\nbrowser pack refGene encodeRegions\ntrack type=wiggle_0 name=\"$base_max_file_name\, $names[$current_chr] positive strand\, SOLiD coverage\" description=\"$base_max_file_name\, $names[$current_chr] positive strand\, SOLiD coverage\" visibility=full color=0,0,255 yLineMark=0 yLineOnOff=on priority=10\nvariableStep chrom=$names[$current_chr] span=1\n";

	print "Start printing $wig_genome_negative_out_file...\n";
	unless ( open(OUT_GENOME_NEGATIVE_WIG_FILE, ">$wig_genome_negative_out_file") ) {
	    print "Cannot open file \"$wig_genome_negative_out_file\" to write to!!\n\n";	
	    exit;
	}
	print OUT_GENOME_NEGATIVE_WIG_FILE "browser position $names[$current_chr]:1-200000000\nbrowser hide all\nbrowser pack refGene encodeRegions\ntrack type=wiggle_0 name=\"$base_max_file_name\, $names[$current_chr] negative strand\, SOLiD coverage\" description=\"$base_max_file_name\, $names[$current_chr] negative strand\, SOLiD coverage\" visibility=full color=0,0,255 yLineMark=0 yLineOnOff=on priority=10\nvariableStep chrom=$names[$current_chr] span=1\n";

    }
}

if($produce_gff_reads eq "yes") {

    # Output reads gff file
    unless ( open(OUT_GFF_READS_FILE, ">$gffreadsfile") ) {
	print "Cannot open file \"$gffreadsfile\" to write to!!\n\n";
	exit;
    }
}


# matching locations are in format reference seq_signed coordinate.number of errors.length of matching

open( MAX_FILE, "< $max_file" ) or die "Can't open $max_file : $!";

while(<MAX_FILE>) {
    if(/^\>/) {
	$number_of_reads++;
	chomp;
	@values = split(/\,/,$_);
	
	#$seq = <MA_FILE>;
	#chomp $seq;
    
	if(scalar @values> 1){
	    if($read_alignment_type eq "unique") {
		if(scalar @values == 2){
		    @index_matching_to_use = (1);
		}
		else {
		    @index_matching_to_use = ();		
		}
	    }
	    if($read_alignment_type eq "all") {
		@index_matching_to_use = 1..(scalar @values -1);
	    }
	    if($read_alignment_type eq "random") {
		@index_matching_to_use = (1+int(rand(scalar @values - 2)));
	    }
	    if($read_alignment_type eq "top_score") {
	    	my @index_max_score = (1);
	    	$temp_max_score = 0;
	    	
	    	for($i=1;$i < scalar @values; $i++) {
		    @individual_values = split(/\#/,$values[$i]);
		    @temp_location = split(/\./,$individual_values[0]);
		    $temp_score = $temp_location[4];
	    	    if($temp_score > $temp_max_score) {
	    		$temp_max_score = $temp_score;
	    		$index_max_score = ($i);
	    	    }
	    	    elsif($temp_score == $temp_max_score) {
	    		push @index_max_score, $i;
	    	    }
	    	}
	    		
	    	@index_matching_to_use = @index_max_score;
	    }
	}
	else {
	    @index_matching_to_use = ();
	}
	
	
	if(scalar @index_matching_to_use> 0){	
	    
	    # for each matching location
	    for($i=0;$i < scalar @index_matching_to_use; $i++) {
		
		@individual_values = split(/\#/,$values[$index_matching_to_use[$i]]);
		@temp_location = split(/\./,$individual_values[0]);
		$temp_ref_seq_index = $temp_location[0];
		
		$temp_coord_start = $temp_location[1];
		$temp_start = $temp_location[2];
		$temp_size = $temp_location[3];
		$temp_score = $temp_location[4];
		$temp_coord_end = $temp_coord_start + $temp_size - 1;

		# check if the reads alignment score is above threshold
		if($temp_score >= $min_score) {
	 	    
	 	    $temp_strand = "+";
	 	    if($temp_coord_start < 0) {
			
			##  Change produced by the new coordinate value in max file
			#$temp_coord_start = (-1)*$temp_location[1] - $temp_size + 1;
			#$temp_coord_end = (-1)*$temp_location[1];
			
			$temp_coord_start = (-1)*$temp_location[1];
			$temp_coord_end = (-1)*$temp_location[1]  + $temp_size - 1;
			$temp_strand = "-";
			
			
		    }
		
		    if($reference_end[$temp_ref_seq_index] < $temp_coord_end) {
			$reference_end[$temp_ref_seq_index] = $temp_coord_end;
		    }


		    # counts
		    if($produce_count eq "yes") {
			
			$info_map_temp_key = ($temp_coord_start - ($temp_coord_start % $window_reference_size));
			# check if read is close to an exon
			if(defined $info_map{$names[$temp_ref_seq_index]}{$info_map_temp_key}) {
			    for($j=0; $j < scalar @{$info_map{$names[$temp_ref_seq_index]}{$info_map_temp_key}}; $j++) {
				$temp_index_gff = $info_map{$names[$temp_ref_seq_index]}{$info_map_temp_key}[$j];
				
				# check strand
				if(($info[$temp_index_gff][6] eq "\.") | ($info[$temp_index_gff][6] eq $temp_strand)) {
				    
				    # check if read aligned part is inside annotated region +/- 3bp
				    if(($temp_coord_start >= ($info[$temp_index_gff][3]-4)) & ($temp_coord_end <= ($info[$temp_index_gff][4]+2))){
					$counts[$temp_index_gff]++;
				    }
				}
			    }
			}
		    
			else {
			    #$counts{$info_map_temp_key} = 1;
			}
		    }
				    

		    # coverage file(s)
		    if($produce_wig_genome eq "yes") {
		      
		      if(($is_sorted eq "yes") & ($read_alignment_type eq "all")) {	
			
			# check if a new chromosome started
			if($temp_ref_seq_index != $current_chr) {
			    
			    # print and reset the coverage
			    # Output wig file for genome rference: one file/sequence/strand in the reference fasta file
			    
			    # printing end of the previous chromosome.
			    &print_wig_file(\*OUT_GENOME_POSITIVE_WIG_FILE, $names[$current_chr], \%coverage_positive, $coverage_filter_threshold);
			    &print_wig_file(\*OUT_GENOME_NEGATIVE_WIG_FILE, $names[$current_chr], \%coverage_negative, $coverage_filter_threshold);

			    %coverage_positive = ();
			    %coverage_negative = ();
			    $current_chr = $temp_ref_seq_index;
			    $current_number_of_bases = 0;
			    
			    if($split_coverage_by_chromosome eq "yes") {
				close(OUT_GENOME_POSITIVE_WIG_FILE);
				close(OUT_GENOME_NEGATIVE_WIG_FILE);

				# open new coverage files
				$wig_genome_positive_out_file = "$wigoutputdir\/$names[$current_chr].positive_strand.wig";
				print "Start printing $wig_genome_positive_out_file...\n";
				unless ( open(OUT_GENOME_POSITIVE_WIG_FILE, ">$wig_genome_positive_out_file") ) {
				    print "Cannot open file \"$wig_genome_positive_out_file\" to write to!!\n\n";	
				    exit;
				}

				print OUT_GENOME_POSITIVE_WIG_FILE "browser position $names[$current_chr]:0-200000000\nbrowser hide all\ntrack type=wiggle_0 name=\"$base_max_file_name\, $names[$current_chr] positive strand\, SOLiD coverage\" description=\"$base_max_file_name\, $names[$current_chr] positive strand\, SOLiD coverage\" visibility=full color=0,0,255 yLineMark=0 yLineOnOff=on priority=10\nvariableStep chrom=$names[$current_chr] span=1\n";

				$wig_genome_negative_out_file = "$wigoutputdir\/$names[$current_chr].negative_strand.wig";
				print "Start printing $wig_genome_negative_out_file...\n";
				unless ( open(OUT_GENOME_NEGATIVE_WIG_FILE, ">$wig_genome_negative_out_file") ) {
				    print "Cannot open file \"$wig_genome_negative_out_file\" to write to!!\n\n";	
				    exit;
				}

				print OUT_GENOME_NEGATIVE_WIG_FILE "browser position $names[$current_chr]:0-200000000\nbrowser hide all\ntrack type=wiggle_0 name=\"$base_max_file_name\, $names[$current_chr] negative strand\, SOLiD coverage\" description=\"$base_max_file_name\, $names[$current_chr] negative strand\, SOLiD coverage\" visibility=full color=0,0,255 yLineMark=0 yLineOnOff=on priority=10\nvariableStep chrom=$names[$current_chr] span=1\n";

			    }
			    else {
				# include a new chromosome entry in the current file
    				print OUT_GENOME_POSITIVE_WIG_FILE "\nvariableStep chrom=$names[$current_chr] span=1\n";
    				print OUT_GENOME_NEGATIVE_WIG_FILE "\nvariableStep chrom=$names[$current_chr] span=1\n";
			    }
			    			    
			    &update_coverage($temp_strand, $temp_coord_start, $temp_coord_end, $names[$current_chr], \%coverage_positive, \%coverage_negative);
			    $current_number_of_bases++;

			} # close if new chromosome
			
			# coverage for the current chromosome
			else {
			
			    # check if we reach memory limit
			    if($current_number_of_bases > $max_number_of_bases_to_keep_in_memory ) {
			    
				# flush the coverage and reset it
				@keys = keys %{$coverage_positive{$names[$current_chr]}};
				$current_number_of_bases = scalar @keys;
				@keys_sorted = sort {$a <=> $b} @keys;
				for($j=0; $j < $current_number_of_bases - $tag_length; $j++){
				    
				    if($coverage_positive{$names[$current_chr]}{$keys_sorted[$j]} >= $coverage_filter_threshold) {
					$temp_coordinate = $keys_sorted[$j] + 1;
					print OUT_GENOME_POSITIVE_WIG_FILE "$temp_coordinate\t$coverage_positive{$names[$current_chr]}{$keys_sorted[$j]}\n";
				    }
				    if($coverage_negative{$names[$current_chr]}{$keys_sorted[$j]} >= $coverage_filter_threshold) {
					$temp_coordinate = $keys_sorted[$j] + 1;
					print OUT_GENOME_NEGATIVE_WIG_FILE "$temp_coordinate\t$coverage_negative{$names[$current_chr]}{$keys_sorted[$j]}\n";
				    }
				    
				    delete $coverage_positive{$names[$current_chr]}{$keys_sorted[$j]};
				    delete $coverage_negative{$names[$current_chr]}{$keys_sorted[$j]};
				}
				
				@keys = keys %{$coverage_positive{$names[$current_chr]}};
				$current_number_of_bases = scalar @keys;
			    }
			    
			    else {
				&update_coverage($temp_strand, $temp_coord_start, $temp_coord_end, $names[$current_chr], \%coverage_positive, \%coverage_negative);
				$current_number_of_bases++;
							    
			    }
			} # close current chromosome coverage update
		      } #close if is_sorted
		      
		      
		      # coverage unsorted max file: handle coverage after sorting coordinates
		      if($is_sorted eq "no") {
			#if($split_coverage_by_chromosome eq "yes") {
			  if($temp_strand eq "+") {
			      $unsort_positive_file_handles[$temp_ref_seq_index]->print("$temp_coord_start\t$temp_size\n");
			  }
			  else {
			      $unsort_negative_file_handles[$temp_ref_seq_index]->print("$temp_coord_start\t$temp_size\n");
			  }
			#}
		      }

		    
		    #close coverage file(s)
		    }
		    
		    # gff output
		    if($produce_gff_reads eq "yes") {
			@gff_values = ();
			$gff_values[0] = $names[$temp_ref_seq_index];
			$gff_values[1] = "SOLiD";
			$gff_values[2] = "read";
			$gff_values[3] = $temp_coord_start;
			$gff_values[4] = $temp_coord_end;
			$gff_values[5] = "\.";
			$gff_values[6] = $temp_strand;
			$gff_values[7] = "\.";
			# here I should put the base sequence
			$gff_values[8] = "\.";
			
			
			for($j=0;$j<scalar @gff_values - 1;$j++){
			    print OUT_GFF_READS_FILE "$gff_values[$j]\t";
			}
			print OUT_GFF_READS_FILE "$gff_values[scalar @gff_values - 1]\n";
		    }

		    $count_match++;
		
		# close if alignment score > threshold
		}
    	    # close for each matching location
    	    }
	# close if there are any matching locations
	}
    }
# close read max file
}

close(MAX_FILE);

if($produce_gff_reads eq "yes") {
    close(OUT_GFF_READS_FILE); 	
}

($date,$time) = &time_stamp();
print "Finished reading max file: time $date, $time\n";

print "\tTotal number of reads = $number_of_reads\n\tnumber of matched reads = $count_match\n";

if(($is_sorted eq "yes") & ($read_alignment_type eq "all")){
# print the last chunk of coverage

    if($produce_wig_genome eq "yes") {
	&print_wig_file(\*OUT_GENOME_POSITIVE_WIG_FILE, $names[$current_chr], \%coverage_positive, $coverage_filter_threshold);
	close(OUT_GENOME_POSITIVE_WIG_FILE);
	&print_wig_file(\*OUT_GENOME_NEGATIVE_WIG_FILE, $names[$current_chr], \%coverage_negative, $coverage_filter_threshold);
	close(OUT_GENOME_NEGATIVE_WIG_FILE);
    }
}
%coverage_positive = ();
%coverage_negative = ();

($date,$time) = &time_stamp();
print "Start generating counts file: time $date, $time\n";


#print "counts[1] = $counts[1]\n";
# counts are produced either when is_sorted=no or when is_sorted=yes AND alignment=all
if($produce_count eq "yes"){
  if((($is_sorted eq "yes") & ($read_alignment_type eq "all")) | ($is_sorted eq "no")) {

    # Output counts
    unless ( open(OUT_FILE, ">$countsoutputfile") ) {
	print "Cannot open file \"$countsoutputfile\" to write to!!\n\n";
	exit;
    }

    for($i = 1; $i < scalar @info; $i++){

	@values = @{$info[$i]};
	if(defined $counts[$i]) {
	    $values[7] = $counts[$i];
	}
	else {
	    $values[7] = 0;
	}
	
	for($j=0; $j < (scalar @values) - 1;$j++) {
	    print OUT_FILE "$values[$j]\t";
	}
	print OUT_FILE "$values[(scalar @values) - 1]\n";
    }

    close(REFERENCE_INFO_FILE);
    close(OUT_FILE);
  }
}

($date,$time) = &time_stamp();
print "Done generating counts file: time $date, $time\n";


# WIGGLE FILE if the max file was not sorted

my $file_size;
my $temp_annotated_wig_file = "";
my $temp_annotated_wig_genome_file = "";


if($is_sorted eq "no") {
    if($produce_wig_genome eq "yes" ) {
	#if($split_coverage_by_chromosome eq "yes") {
	    foreach my $fh (@unsort_positive_file_handles) {
		close $fh;
	    }
	    foreach my $fh (@unsort_negative_file_handles) {
		close $fh;
	    }
	#}
	
	($date,$time) = &time_stamp();
	print "Start generating wig file(s): time $date, $time\n";

	if($split_coverage_by_chromosome eq "no") {
	    $wig_genome_positive_out_file = "$wigoutputdir\/genome.positive_strand.wig";
	    $wig_genome_negative_out_file = "$wigoutputdir\/genome.negative_strand.wig";	
	
	    print "Start printing $wig_genome_positive_out_file...\n";
	    unless ( open(OUT_GENOME_POSITIVE_WIG_FILE, ">$wig_genome_positive_out_file") ) {
		print "Cannot open file \"$wig_genome_positive_out_file\" to write to!!\n\n";	
		exit;
	    }
	    print OUT_GENOME_POSITIVE_WIG_FILE "browser position chr1:1-200000000\nbrowser hide all\nbrowser pack refGene encodeRegions\ntrack type=wiggle_0 name=\"$base_max_file_name\, positive strand\, SOLiD coverage\" description=\"$base_max_file_name\, positive strand\, SOLiD coverage\" visibility=full color=0,0,255 yLineMark=0 yLineOnOff=on priority=10\nvariableStep chrom=chr1 span=1\n";

	    print "Start printing $wig_genome_negative_out_file...\n";
	    unless ( open(OUT_GENOME_NEGATIVE_WIG_FILE, ">$wig_genome_negative_out_file") ) {
		print "Cannot open file \"$wig_genome_negative_out_file\" to write to!!\n\n";	
		exit;
	    }
	    print OUT_GENOME_NEGATIVE_WIG_FILE "browser position chr1:1-200000000\nbrowser hide all\nbrowser pack refGene encodeRegions\ntrack type=wiggle_0 name=\"$base_max_file_name\, negative strand\, SOLiD coverage\" description=\"$base_max_file_name\, negative strand\, SOLiD coverage\" visibility=full color=0,0,255 yLineMark=0 yLineOnOff=on priority=10\nvariableStep chrom=chr1 span=1\n";
	}

	my $cmd;
	$cmd = system("rm $unsort_positive_coordinate_files[0]"); 
	$cmd = system("rm $unsort_negative_coordinate_files[0]"); 
	
	for($i=1; $i < scalar @names; $i++) {
	  
	    if($split_coverage_by_chromosome eq "yes") {
		$wig_genome_positive_out_file = "$wigoutputdir\/$names[$i].positive_strand.wig";
		$wig_genome_negative_out_file = "$wigoutputdir\/$names[$i].negative_strand.wig";

		print "Start printing $wig_genome_positive_out_file...\n";
		unless ( open(OUT_GENOME_POSITIVE_WIG_FILE, ">$wig_genome_positive_out_file") ) {
		    print "Cannot open file \"$wig_genome_positive_out_file\" to write to!!\n\n";	
		    exit;
		}
		print OUT_GENOME_POSITIVE_WIG_FILE "browser position $names[$i]:1-200000000\nbrowser hide all\nbrowser pack refGene encodeRegions\ntrack type=wiggle_0 name=\"$base_max_file_name\, $names[$i] positive strand\, SOLiD coverage\" description=\"$base_max_file_name\, $names[$i] positive strand\, SOLiD coverage\" visibility=full color=0,0,255 yLineMark=0 yLineOnOff=on priority=10\nvariableStep chrom=$names[$i] span=1\n";

		print "Start printing $wig_genome_negative_out_file...\n";
		unless ( open(OUT_GENOME_NEGATIVE_WIG_FILE, ">$wig_genome_negative_out_file") ) {
		    print "Cannot open file \"$wig_genome_negative_out_file\" to write to!!\n\n";	
		    exit;
		}
		print OUT_GENOME_NEGATIVE_WIG_FILE "browser position $names[$i]:1-200000000\nbrowser hide all\nbrowser pack refGene encodeRegions\ntrack type=wiggle_0 name=\"$base_max_file_name\, $names[$i] negative strand\, SOLiD coverage\" description=\"$base_max_file_name\, $names[$i] negative strand\, SOLiD coverage\" visibility=full color=0,0,255 yLineMark=0 yLineOnOff=on priority=10\nvariableStep chrom=$names[$i] span=1\n";
	    }
	    
	    ##  POSITIVE STRAND
	    
	    $file_size = -s $unsort_positive_coordinate_files[$i];
	    
	    if($file_size > 0){
	  
		$cmd = system("sort -n  --buffer-size=2G -T $wigoutputdir $unsort_positive_coordinate_files[$i] -o $sort_positive_coordinate_files[$i]");
		$cmd = system("rm $unsort_positive_coordinate_files[$i]"); 


		##  SORT MATCHING LOCATIONS
	    
		open( MATCHING_LOCATIONAS_SORT, "< $sort_positive_coordinate_files[$i]" ) or die "Can't open $sort_positive_coordinate_files[$i] : $!";
		$temp_line = <MATCHING_LOCATIONAS_SORT>;
		chomp $temp_line;
		@values = split(/\t/, $temp_line);
		my $reference_start = $values[0];
		close(MATCHING_LOCATIONAS_SORT);


		# WIGGLE FILE

		my $max_number_of_bases_to_keep_in_memory = 1e6;
		my (@keys, @keys_sorted);


		if($produce_wig_genome eq "yes") {

		    open( MATCHING_LOCATIONAS_SORT, "< $sort_positive_coordinate_files[$i]" ) or die "Can't open $sort_positive_coordinate_files[$i] : $!";

		    $current_number_of_bases = 0;
		    %coverage = ();
		    $coverage{$names[$i]} = ();

		    while(<MATCHING_LOCATIONAS_SORT>) {

			chomp;
			@values = split(/\t/);
			$temp_coord_start = $values[0];
			$temp_coord_end = $temp_coord_start + $values[1] - 1;

			if($current_number_of_bases > $max_number_of_bases_to_keep_in_memory ) {

			    # flush the coverage and reset it
			    @keys = keys %{$coverage{$names[$i]}};
			    @keys_sorted = sort {$a <=> $b} @keys;
			    $current_number_of_bases = scalar @keys;
			    
			    for($j=0; $j < $current_number_of_bases - $tag_length; $j++){
				if(($produce_wig_genome eq "yes") & ($coverage{$names[$i]}{$keys_sorted[$j]} >= $coverage_filter_threshold)) {
				    $temp_coordinate = $keys_sorted[$j] + 1;
				    $temp_annotated_wig_genome_file .= "$temp_coordinate\t$coverage{$names[$i]}{$keys_sorted[$j]}\n";
				}

				delete $coverage{$names[$i]}{$keys_sorted[$j]};
			    }
			    print OUT_GENOME_POSITIVE_WIG_FILE $temp_annotated_wig_genome_file;
			    $temp_annotated_wig_genome_file = "";
				
			    @keys = keys %{$coverage{$names[$i]}};
			    $current_number_of_bases = scalar @keys;
			}

			else {
			    for($j=$temp_coord_start; $j <= $temp_coord_end; $j++){
				if(defined $coverage{$names[$i]}{$j}){
				    $coverage{$names[$i]}{$j}++;
				}
				else {
				    $coverage{$names[$i]}{$j} = 1;
				    $current_number_of_bases++;
				}
			    }			    
			}

		    }

		    ##  flush remaining wig files
		    if($produce_wig_genome eq "yes") {
			&print_wig_file(\*OUT_GENOME_POSITIVE_WIG_FILE, $names[$i], \%coverage, $coverage_filter_threshold);
			if($split_coverage_by_chromosome eq "yes") {
			    close(OUT_GENOME_POSITIVE_WIG_FILE);
			}
			else {
			    if(($i + 1) < scalar @names ){
			    	print OUT_GENOME_POSITIVE_WIG_FILE "\nvariableStep chrom=$names[$i+1] span=1\n";
			    }
			}
		    }

		    %coverage = ();
		    $cmd = system("rm $sort_positive_coordinate_files[$i]");


		    print "Done positive coverage for $names[$i]...\n";
		    # close if there is some coverage of chri
		}
	    }
	    else {
		$cmd = system("rm $unsort_positive_coordinate_files[$i]"); 
	    }
	
	    ##  NEGATIVE STRAND
	    
	    $file_size = -s $unsort_negative_coordinate_files[$i];
	    
	    if($file_size > 0){
	  
		$cmd = system("sort -n  --buffer-size=2G -T $wigoutputdir $unsort_negative_coordinate_files[$i] -o $sort_negative_coordinate_files[$i]");
		$cmd = system("rm $unsort_negative_coordinate_files[$i]"); 


		##  SORT MATCHING LOCATIONS
	    
		open( MATCHING_LOCATIONAS_SORT, "< $sort_negative_coordinate_files[$i]" ) or die "Can't open $sort_negative_coordinate_files[$i] : $!";
		$temp_line = <MATCHING_LOCATIONAS_SORT>;
		chomp $temp_line;
		@values = split(/\t/, $temp_line);
		my $reference_start = $values[0];
		close(MATCHING_LOCATIONAS_SORT);


		# WIGGLE FILE

		my $max_number_of_bases_to_keep_in_memory = 1e6;
		my (@keys, @keys_sorted);


		if($produce_wig_genome eq "yes") {

		    open( MATCHING_LOCATIONAS_SORT, "< $sort_negative_coordinate_files[$i]" ) or die "Can't open $sort_negative_coordinate_files[$i] : $!";

		    $current_number_of_bases = 0;
		    %coverage = ();
		    $coverage{$names[$i]} = ();

		    while(<MATCHING_LOCATIONAS_SORT>) {

			chomp;
			@values = split(/\t/);
			$temp_coord_start = $values[0];
			$temp_coord_end = $temp_coord_start + $values[1] - 1;

			if($current_number_of_bases > $max_number_of_bases_to_keep_in_memory ) {

			    # flush the coverage and reset it
			    @keys = keys %{$coverage{$names[$i]}};
			    @keys_sorted = sort {$a <=> $b} @keys;
			    $current_number_of_bases = scalar @keys;
			    
			    for($j=0; $j < $current_number_of_bases - $tag_length; $j++){
				if(($produce_wig_genome eq "yes") & ($coverage{$names[$i]}{$keys_sorted[$j]} >= $coverage_filter_threshold)) {
				    $temp_coordinate = $keys_sorted[$j] + 1;
				    $temp_annotated_wig_genome_file .= "$temp_coordinate\t$coverage{$names[$i]}{$keys_sorted[$j]}\n";
				}

				delete $coverage{$names[$i]}{$keys_sorted[$j]};
			    }
			    print OUT_GENOME_NEGATIVE_WIG_FILE $temp_annotated_wig_genome_file;
			    $temp_annotated_wig_genome_file = "";
			    
			    @keys = keys %{$coverage{$names[$i]}};
			    $current_number_of_bases = scalar @keys;
			}

			else {
			    for($j=$temp_coord_start; $j <= $temp_coord_end; $j++){
				if(defined $coverage{$names[$i]}{$j}){
				    $coverage{$names[$i]}{$j}++;
				}
				else {
				    $coverage{$names[$i]}{$j} = 1;
				    $current_number_of_bases++;
				}
			    }			    
			}

		    }

		    ##  flush remaining wig files
		    if($produce_wig_genome eq "yes") {
			&print_wig_file(\*OUT_GENOME_NEGATIVE_WIG_FILE, $names[$i], \%coverage, $coverage_filter_threshold);
			if($split_coverage_by_chromosome eq "yes") {
			    close(OUT_GENOME_NEGATIVE_WIG_FILE);
			}
			else {
			    if(($i + 1) < scalar @names ){
				print OUT_GENOME_NEGATIVE_WIG_FILE "\nvariableStep chrom=$names[$i+1] span=1\n";
			    }
			}
		    }

		    %coverage = ();
		    $cmd = system("rm $sort_negative_coordinate_files[$i]");


		    print "Done negative coverage for $names[$i]...\n";
		    # close if there is some coverage of chri
		}
	    }
	    else {
		$cmd = system("rm $unsort_negative_coordinate_files[$i]"); 
	    }
	# close for each chromosome    
	}

	($date,$time) = &time_stamp();
	print "Done generating wig file(s): time $date, $time\n";

    }
}

exit;


sub check_files {

    my ($max_file, $reference_fasta_file, $reference_info_file) = @_;
    if(defined $max_file) {
	unless(-e $max_file) {
	    print "max file does not exist\n";
	}
    }
    else {
	print "\-max \<max_file\> needs to be defined\n";
	exit(1);
    }
    
    if(defined $reference_fasta_file) {
	unless(-e $reference_fasta_file) {
	    print "reference_fasta_file file does not exist\n";
	}
    }
    else {
	print "\-r \<reference_fasta_file\> needs to be defined\n";
	exit(1);    
    }
    
    if(defined $reference_info_file){
	unless(-e $reference_info_file) {
	    print "max file does not exist\n";
	}
    }
    else {
	print "\-a \<annotation_gtf_file\> needs to be defined\n";
	exit(1);
    }
}

sub check_optional_file {

    my($call, $filez) = @_;
    my @values;
    my $i;
    
    if($call eq "yes") {
	@values = split(/\//,$filez);
	my $temp = '';
	for (my $i = 0; $i < (scalar @values)-1; $i++) {
	    $temp .= "$values[$i]/";
	    unless(-e $temp){
		mkdir $temp or die "Can't mkdir $temp: $!\n";
            }
        }
        print "$temp directory was created\n";
    }
}

sub check_optional_dir {

    my($call, $filez) = @_;
    my @values;
    my $i;
    
    if($call eq "yes") {
	@values = split(/\//,$filez);
	my $temp = '';
	for (my $i = 0; $i < scalar @values; $i++) {
	    $temp .= "$values[$i]/";
	    unless(-e $temp){
		mkdir $temp or die "Can't mkdir $temp: $!\n";
            }
        }
        print "$temp directory was created\n";
    }
}


sub check_options {
    
    my($is_sorted, $read_alignment_type, $produce_count, $produce_wig_genome, $produce_gff_reads) = @_;
    if(($is_sorted eq "yes") & ($read_alignment_type ne "all")){
	print "\n\t\-is\_sorted $is_sorted and \-read_alignment_type $read_alignment_type IS NOT IMPLEMENTED YET !!!;\n\t USE \-is\_sorted \'yes\' \-read_alignment_type \'all\'  OR START WITH UNSORTED MAX FILE WITH \-is\_sorted \'no\'!!\n\n";
	exit;
    }
    
    if($is_sorted eq "no"){
	print "is\_sorted=$is_sorted will take some time to generate genome wiggle files...\n\n";
    }    
}

sub get_base_max_file_name {
    my ($max_file) = @_;
    my @zz = split(/\//, $max_file);
    my @temp_name = split(/\./, $zz[(scalar @zz) - 1]);
    my $name_out = $temp_name[0];
    return($name_out);
}


sub map_info_reference {
    
    my ($reference_info_file, $window_reference_size) = @_;
  
    my ($index_entry, $start, $end, $rem, $q1, $q2, $i, $temp_location);
    my %info_map;
    my @values;
    my $counter_entry = 1;
	my %archive;
	
    open( REFERENCE_INFO_FILE, "< $reference_info_file" ) or die "Can't open $reference_info_file : $!";

    while(<REFERENCE_INFO_FILE>) {
	chomp;
	if(!/^\#/ && /\texon\t/) {	
	    @values = split(/\t/,$_);
	    # values[0] = index_entry
	    # values[3] = start
	    # values[4] = end;
	    
	    $index_entry = $values[0];
	    $start = $values[3];
	    $end = $values[4];
	
		#Skip entry if you've seen it before.
		my $archive_key = join ':', $index_entry, $start, $end;
		if ($archive{$archive_key}) {
			next;
		}    
		else {
			$archive{$archive_key} = 1;
		}
	    $rem = abs($start) % $window_reference_size;
	    $q1 = (abs($start) - $rem) / $window_reference_size;

	    $rem = abs($end) % $window_reference_size;
	    $q2 = (abs($end) - $rem) / $window_reference_size;
	    
	    for($i=$q1; $i < ($q2+1); $i++) {
	        $temp_location = $i*$window_reference_size;
	        if(defined $info_map{$index_entry}{$temp_location}) {
		    $info_map{$index_entry}{$temp_location}[scalar @{$info_map{$index_entry}{$temp_location}}] = $counter_entry; 
		}
		else {
		    $info_map{$index_entry}{$temp_location}[0] = $counter_entry;
		}
	    }
	    $counter_entry++;    
	}
    }

    close(REFERENCE_INFO_FILE);

    return(%info_map);  
}

sub get_reference_info {

    my($reference_info_file) = @_;
    
    my (@info, @values);
    open( REFERENCE_INFO_FILE, "< $reference_info_file" ) or die "Can't open $reference_info_file : $!";

    @info = ["test", ".",".","0","0"];
	my %archive;
	
    while(<REFERENCE_INFO_FILE>) {
	chomp;
	if(!/^\#/ && /\texon\t/) {	
	    @values = split(/\t/,$_);
		my $archive_key = join ':', $values[0], $values[3], $values[4];
		if ($archive{$archive_key}) {
			next;
		}
		else {
			$archive{$archive_key} = 1;
		}
	    push(@info , [@values]);
	}
    }

    close(REFERENCE_INFO_FILE);

    return(@info);
}

sub get_ref_sequence_names {

    my ($reference_fasta_file) = @_;
    my @names = ("zero_name");
    my $unique_names = {};
    my @values;
    my ($seqID, $buf_name);
    
    open( REFERENCE_INFO_FILE, "< $reference_fasta_file" ) or die "Can't open $reference_fasta_file : $!";
    
    while(<REFERENCE_INFO_FILE>) {
	chomp;
	if(/^\>/) {
	    $buf_name = substr($_,1);
	    @values = split(/ /,$buf_name);
	    $seqID = $values[0];
	    if($seqID eq "") {
		$seqID = $values[1];
	    }
	    
	    push @names, $seqID;
	    if(defined $unique_names->{$seqID}) {
		print "$seqID in reference file is duplicated ..\n";
		$unique_names->{$seqID}++;
	    }
	    else {
		$unique_names->{$seqID} = 1;
	    }
	
	}
    }
    
    
    return(@names);
    
}

sub time_stamp {
  my ($d,$t);
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
        $year += 1900;
        $mon++;
        $d = sprintf("%4d-%2.2d-%2.2d",$year,$mon,$mday);
        $t = sprintf("%2.2d:%2.2d:%2.2d",$hour,$min,$sec);
        return($d,$t);
}

sub print_wig_file {
    
    my ($OUT_WIG_FILE, $chr_name, $coverage, $coverage_filter_threshold) = @_;

    my ($j, $temp_coordinate);		    
    my @keys = keys %{$$coverage{$chr_name}};
    my @keys_sort = sort {$a <=> $b} @keys;
 
			    
    for($j=0; $j<scalar @keys_sort; $j++){
	if($$coverage{$chr_name}{$keys_sort[$j]} >= $coverage_filter_threshold){
	    $temp_coordinate = $keys_sort[$j] + 1;
	    print $OUT_WIG_FILE "$temp_coordinate\t$$coverage{$chr_name}{$keys_sort[$j]}\n";
	}
	delete $$coverage{$chr_name}{$keys_sort[$j]};
    }
}

sub update_coverage {
    my ($temp_strand, $temp_coord_start, $temp_coord_end, $chr_name, $coverage_positive, $coverage_negative) = @_;
    my $j;
    
    for($j=$temp_coord_start; $j <= $temp_coord_end; $j++){
	if($temp_strand eq "+") {
	    if(defined $coverage_positive{$chr_name}{$j}){
		$$coverage_positive{$chr_name}{$j}++;
	    }
	    else {
		$$coverage_positive{$chr_name}{$j} = 1;
		$$coverage_negative{$chr_name}{$j} = 0;
	    }
	}
	else {
	    if(defined $coverage_negative{$chr_name}{$j}){
		$$coverage_negative{$chr_name}{$j}++;
	    }
	    else {
		$$coverage_negative{$chr_name}{$j} = 1;
		$$coverage_positive{$chr_name}{$j} = 0;
	    } 				    
	}
    }

}
