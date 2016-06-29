#!/usr/bin/perl -w
# select_contigs - Select a subset of contigs from a fasta input file.
#   A fasta quality file also can be processed.  Contigs can be reversed
#   and complemented, and partial contigs can be extracted.
#
#   The advantage of select_contigs over sort_contigs is that much
#   less memory is used for processing large input files, because
#   the entire file is not read into memory all at once.  This is also
#   faster for large files.  Note that the sort_contigs -o option is
#   the -n option in this program.
#
# Written by: James D. White, University of Oklahoma, Advanced Center for
#   Genome Technology
#
# Date Written: Aug 5, 2009
#
# 20100706 JDW - Add -i to ignore the list of contigs specified by -n.
# 20100204 JDW - Add very verbose mode (-V), which does what -v used
#		 to do.  -v now will only output how many contigs were
#		 written.  Fix bug when -n was used with a real file
#		 and not just '-' (STDIN).
# 20090907 JDW - Date rewritten from sort_contigs.

use strict;

our $Date_Last_Modified = "July 6, 2010";


######################## Customization Parameters ###########################

our $DEFAULT_LINE_LENGTH = 50;	# Default length of fasta output sequence lines.
our $output_line_len = $DEFAULT_LINE_LENGTH;# Length of fasta output sequence lines.
our $MIN_LINE_LENGTH = 10;	# Mimimum allowed value for $output_line_len.

our $contig_prefix = '';	# Prefix to add to output contig names
				#				(-c)
our $no_degenerate = 0;		# Input is dna.  Filter out any
				#   contigs that do not contain at
				#   least one each of A, C, G, and T.
				#				(-d)
our $empty_ok = 0;		# Empty output files with no contigs
				#   are OK?			(-e)
our $trim_ends = 0;		# Trim Ns from ends of contigs	(-g)
our $ignore_list = 0;		# Ignore contigs specified by -n (-i)
our $min_size = 1;		# Ignore input contigs < $min_size in
				#   length.			(-m)
our $MAX_QUAL = 99;		# Maximum quality value allowed in
				#   output file (-M)
our $select_file = '';		# Name of file to read for selecting
				#   contigs to be output.	(-n)
our $select_all = 1;		# Select ALL contigs in full? (not -n)
our $select_all_line = "F\t1\t0\t-";	# used with $select_all
our $Preserve_Comments = 0;	# Preserve contig header comments? (-p)
our $get_qualities = 0;		# Is a fasta quality file to be
				#   used/created? (-q)
our $QUALITY_VALUE = undef;	# Constant quality value to be used
				#   instead of reading a fasta qual
				#   file (-Q ##)
our $QUALITY_SUB = undef;	# Constant value to be subtracted from
				#   all quality values from the fasta
				#   qual file (-Q SUB##)
our $QUALITY_DIV = undef;	# All quality values from the fasta
				#   qual file are to be divided by
				#   this constant value (-q DIV#.###)
our @QUALITY_STEPS = ();	# Quality values from the fasta qual
				#   file are to be divided by stepped
				#   integer values
				#   (-Q STEPpos1[,pos2,...])
our @QUALITY_SLOPES = ();	# Quality values from the fasta qual
				#   file are to be divided by
				#   interpolated values
				#   (-Q SLOPEpos1,pos2[,pos3,...])
our $shorten_contig_name = 0;	# Is contig name to be shortened? (-s)
our $Type_Strip = '';		# Type to strip from output filename
				#   before adding ".qual"	(-t)
our $Use_Uaccno = 0;		# Use universal accession numbers
				#   (uaccno) as contig names for 454
				#   reads.			(-u)
our $Verbose_Mode = 0;		# Are some statistics to be listed?
				#				(-v/-V)
our $Extend_File = 0;		# Extend (append to) existing output
				#   file(s) (see -x)
########## Operating System specific ##########

my $directory_separator = '/';	# For Unix
#$directory_separator = '\';	# For DOS/Windows

################# End of Customization Parameters ###################


my $full_command_name = $0;
our($my_name);
if ($full_command_name =~
    m"(^|${directory_separator})([^${directory_separator}]*)$")
  {
  $my_name = $2;
  }
else
  {
  $my_name = 'select_contigs';
  }

  
################### Process flags on command line ###################

use Getopt::Std;
our($opt_c, $opt_d, $opt_e, $opt_g, $opt_h, $opt_i, $opt_l, $opt_m,
    $opt_M, $opt_n, $opt_p, $opt_q, $opt_Q, $opt_s, $opt_t, $opt_u,
    $opt_v, $opt_V, $opt_x) = ('') x 19;
display_help('') unless( getopts('c:deghil:m:M:n:pqQ:st:uvVx') );
display_more_help() if ($opt_h);		# -h	(help)

$Verbose_Mode = 1 if ($opt_v);			# -v	(verbose mode)
$Verbose_Mode = 2 if ($opt_V);			# -V	(very verbose mode)

print STDERR "\n$my_name - Last Modified: $Date_Last_Modified\n\n"
  if ($Verbose_Mode > 1);

$contig_prefix = $opt_c if ($opt_c ne '');	# -c contig_prefix
$no_degenerate = 1 if ($opt_d);			# -d
$empty_ok = 1 if ($opt_e);			# -e
$trim_ends = 1 if ($opt_g);			# -g
$ignore_list = 1 if ($opt_i ne '');		# -i
$output_line_len = $opt_l if ($opt_l ne '');	# -l output_line_length
$min_size = $opt_m if ($opt_m ne '');		# -m min_size
$MAX_QUAL = $opt_M if ($opt_M ne '');		# -M max_qual
if ($opt_n ne '')				# -n select_file
  {
  $select_file = $opt_n;
  $select_all = 0;	# we do not automatically select all contigs
  }
$Preserve_Comments = 1 if ($opt_p);		# -p
$get_qualities = 1 if ($opt_q);			# -q
my $Qmsg = '';
if ($opt_Q ne '')				# -Q
  {
  if ($opt_Q =~ /^\d{1,2}$/)		# set constant $QUALITY_VALUE
    {
    $QUALITY_VALUE = $opt_Q;
    $get_qualities = 2;			# do not read quality file
    }
  elsif ($opt_Q =~ /^S(?:UB)?(\d{1,2})$/i)
    {				# subtract $QUALITY_SUB from qualities
    $QUALITY_SUB = $1;
    }
  elsif ($opt_Q =~ /^D(?:IV)?([1-9](\.\d*)?)$/i)
    {				# divide qualities by constant values
    $QUALITY_DIV = $1;
    }
  elsif ($opt_Q =~ /^STEP([\d,]+)$/i)
    {				# divide qualities by stepped values
    @QUALITY_STEPS = map { (defined $_ && $_ ne '') ? $_ : 0 }
      split(',', $1);
    }
  elsif ($opt_Q =~ /^SLOPE([\d,]+)$/i)
    {				# divide qualities by sloped values
    @QUALITY_SLOPES = map { (defined $_ && $_ ne '') ? $_ : 0 }
      split(',', $1);
    }
  else
    {
    $Qmsg = "Invalid -Q qual_value='$opt_Q'\n";
    $get_qualities |= 4;
    }
  $get_qualities = ($get_qualities == 2) ? 2 : # adjust $get_qualities
		   ($get_qualities == 4) ? 0 : 1;
  } # end if ($opt_Q ne '')
$shorten_contig_name = 1 if ($opt_s);		# -s
$Type_Strip = $opt_t if ($opt_t ne '');		# -t type_to_remove
$Type_Strip =~ s/\./\\./g;	# make '.' mean a period
$Use_Uaccno = 1 if ($opt_u);			# -u
$Extend_File = 1 if ($opt_x);			# -x

###################### check for invalid flags ######################

my $help_msg = '';
$help_msg .= "Contig prefix='$contig_prefix' cannot contain white space\n"
  if ($contig_prefix =~ /\s/);
$help_msg .= "Invalid -c flag, contig_prefix='$contig_prefix'\n"
  if ($contig_prefix =~ /[^-\w]/);
$help_msg .= "Cannot use -i flag without -n flag\n"
  if ($ignore_list && ($select_file eq ''));
$help_msg .= "Invalid -M max_qual='$MAX_QUAL', must be in range 1-99\n"
  unless ($MAX_QUAL =~ /^\d{1,2}$/ && $MAX_QUAL > 0);
$help_msg .= $Qmsg;
$help_msg .= "Cannot use -p and -u flags together\n"
  if ($Preserve_Comments + $Use_Uaccno > 1);
$help_msg .= "Invalid value for -l flag, output_line_length='$output_line_len'\n"
  if (($output_line_len !~ /^\d+$/) ||
      ($output_line_len < $MIN_LINE_LENGTH));
$help_msg .= "Invalid value for -m flag, min_size='$min_size'\n"
  if ($min_size !~ /^\d+$/);
$min_size = 1 if ($min_size == 0);
$help_msg .= "For -n flag, missing or invalid select_file='$select_file'\n"
  if ($select_file ne '' && $select_file ne '-' && ! -f $select_file);

#################### Check command line arguments ####################

my $fasta_input_file = shift @ARGV;
$help_msg .= "Missing 'fasta_input_file' name.\n"
  if (! defined $fasta_input_file || $fasta_input_file eq '');
$help_msg .= "'fasta_input_file' name specified as '-' (STDIN), which is not\n" .
	     "    allowed when '-q' or a '-Q' modifier is specified.\n"
  if ((defined $fasta_input_file) && ($fasta_input_file eq '-') &&
      ($get_qualities == 1));

my $fasta_output_file = shift @ARGV;
$help_msg .= "Missing 'fasta_output_file' name.\n"
  if (! defined $fasta_output_file || $fasta_output_file eq '');
$help_msg .= "'fasta_output_file' name specified as '-' (STDOUT), which is not\n" .
	     "    allowed when '-q', '-Q', or '-x' is specified.\n"
  if ((defined $fasta_output_file) && ($fasta_output_file eq '-') &&
      ($get_qualities || $Extend_File));

###### if command line errors, print error(s) & usage then exit ######

display_help($help_msg) if ($help_msg ne '');


### Read and process 'select_file' and save in @SELECTIONS, if -n ####

our($Short_Input_Contigs, $Short_Contigs, $select_sequences,
    $select_input_lines, $Num_Contigs, $Degenerate_Contigs,
    $GOOD_SELECT_CONTIGS, ) = (0) x 7;
our($line);
our(%SELECTIONS, %PROCESSED, %qualities_modified, ) = ();

unless ($select_all)
  {
  our($SELECTFILE);
  if ($select_file eq '-')
    {
    $SELECTFILE = 'STDIN';
    }
  else
    {
    open($SELECTFILE, $select_file) or
      die("Can't open select_file: '$select_file'\n");
    }
  while (defined ($line = <$SELECTFILE>))
    {
    $select_input_lines++;
    chomp($line);
    $line =~ s/^\s+//;
    next if $line =~ /^#/;	# ignore comment lines
    $line =~ s/\s+#.*$//;	# remove trailing comments
    $line =~ s/\s+$//;
    next if $line eq '';
    # input is (input contig name, direction (F or C/R), begin base,
    #   length, output contig name)
    my($contig, $direction, $begin, $length, $out_contig) =
      split(' ', $line);

    if ($ignore_list)
      {
      $SELECTIONS{$contig} = 1;
      }
    else
      {
      $SELECTIONS{$contig} = [] unless (exists $SELECTIONS{$contig});
      $direction = (defined $direction && $direction =~ /^[rc]/i) ? 'R' : 'F';
      if (defined $begin)
	{
	if ($begin !~ /^\-?\d*$/i)
	  {
	  die "Invalid 'begin_base'='$begin' on select_file line $select_input_lines\n";
	  }
	}
      else
	{
	$begin = 1;
	}
      if (defined $length)
	{
	if ($length !~ /^\-?\d*$/i)
	  {
	  die "Invalid length='$length' on select_file line $select_input_lines";
	  }
	}
      else
	{
	$length = 0;
	}
      $out_contig = '-' if (! defined $out_contig);
      push @{ $SELECTIONS{$contig} },
	join("\t", $direction, $begin, $length, $out_contig);
      } # end unless ($ignore_list)
    $select_sequences++;
    } # end while (defined ($line = <$SELECTFILE>))
  close($SELECTFILE) if ($select_file ne '-');
  die "Warning: No contig names read from select_file='$select_file'." .
      "  Output file(s) not written.\n"
    unless ($empty_ok || $select_sequences);
  my $select_contigs = scalar keys %SELECTIONS;
  print "$select_sequences sequences from $select_contigs to be written according to select_file\n"
    if ($Verbose_Mode > 1);
  } # end unless ($select_all)


#################### Open input sequence file(s) #####################

my $fasta_qual_file = "${fasta_input_file}.qual";
if ($get_qualities == 1 && ! -f $fasta_qual_file)
  {
  $fasta_qual_file = $fasta_input_file;
  # Just adding .qual didn't work, so now try removing the last
  #   qualifier, e.g., xxx.fa + xxx.qual, xxx.fna + xxx.qual, or
  #   xxx.fasta + xxx.qual
  unless ($fasta_qual_file =~ s/\.f(ast|n)?a$/.qual/ &&
          -f $fasta_qual_file)
    {
    die "Can't find a fasta quality file for '$fasta_input_file'.\n";
    }
  }
if (! defined open(FASTAIN, $fasta_input_file))
  {
  die("Can't open fasta_input_file: '$fasta_input_file', $!\n");
  }
if ($get_qualities == 1)
  {
  if (! defined open(QUALIN, $fasta_qual_file))
    {
    my $err = $!;
    close(FASTAIN);
    die("Can't open input fasta quality file: '$fasta_qual_file', $err\n");
    }
  }


########################## Open output files #########################

my $append_flag = ($Extend_File) ? '>' : '';
my $output_msg = ($Extend_File) ? 'append to' : 'create';
if (!open(FASTAOUT, ">$append_flag$fasta_output_file"))
  {
  close(FASTAIN);
  close(QUALIN) if ($get_qualities);
  die("Can't $output_msg fasta_output_file: '$fasta_output_file', $!\n");
  }
if ($get_qualities)
  {
  my $qualout = $fasta_output_file;
  $qualout =~  s/$Type_Strip$// if ($Type_Strip ne '');
  $qualout .= ".qual";
  if (!open(QUALOUT, ">$append_flag$qualout"))
    {
    my $err = $!;
    close(FASTAIN);
    close(QUALIN);
    close(FASTAOUT);
    die("Can't $output_msg output fasta quality file: '$qualout', $!\n");
    exit 2;
    }
  }


################# Read and process selected contigs ##################

my($header, $contig, $comment, $sequence, ) = ('') x 4;
my($line_num, ) = (0) x 1;
my($Qline, $Qcontig, $Qheader, $Qlen, $Qcomment, $Qline_num, $Quality,
   @Qualities, );
if ($get_qualities == 1)
  {
  $Qline = <QUALIN>;
  $Qline_num = 1;
  $Qcontig = '';
  $Qlen = 0;
  }
  $Qcomment = '';
while ($line = <FASTAIN>)
  {
  chomp $line;
  $line_num++;
  if ($line =~ /^>/)
    {
    if ($contig)	# after first input line?
      {
      if ($get_qualities == 1 && length($sequence) != $Qlen)
        {
        print STDERR "Lengths of fasta sequence and quality files do not match on\n  contig='$contig'\n";
        exit 2;
        }
      process_contig($contig, $comment, $sequence, $Qcomment);
      }
    $Num_Contigs++;
    $header = $line;
    if ($header =~ m/^>(\S+)(.*)$/)
      {
      $contig = $1;
      $comment = $2;
      }
    else
      {
      print STDERR "Error: Invalid fasta_input_file format: '$fasta_input_file'\n  Fasta input line number=$line_num\n";
      exit 2;
      }
    $sequence = '';
    if ($get_qualities == 1)
      {
      if (! defined $Qline)
        {
        print STDERR "Fasta sequence file has extra contig(s) not found in quality file.\n  Contig header number $Num_Contigs\n  Sequence contig='$contig', last Quality contig='$Qcontig'\n";
        exit 2;
        }
      chomp($Qline);
      $Qheader = $Qline;
      if ($Qheader =~ m/^>(\S+)(.*)$/)
        {
        $Qcontig = $1;
        $Qcomment = $2;
        if ($contig ne $Qcontig)
          {
          print STDERR "Fasta sequence and quality files do not match on contig header number $Num_Contigs\n  Sequence contig='$contig', Quality contig='$Qcontig'\n";
          exit 2;
          }
        }
      else
        {
        print STDERR "Error: Invalid fasta_qual_input_file format: '${fasta_input_file}.qual'\n  Quality file input line number=$Qline_num,\n  Qline='$Qline'\n";
        exit 2;
        }
      $Qline = <QUALIN>;
      $Qline_num++;
      $Quality = '';
      while (defined $Qline && length($Qline) && $Qline !~ /^>/)
        {
        chomp($Qline);
        $Quality .= ' ' . $Qline;
        $Qline = <QUALIN>;
        $Qline_num++;
        }
      $Quality =~ s/^\s+//;
      $Quality =~ s/\s+$//;
      @Qualities = split(' ', $Quality);
      $Qlen = scalar @Qualities;
      } # end if ($get_qualities == 1)
					# Remove Contig name prefix?
    $contig =~ s/^\S*Contig/Contig/ if $shorten_contig_name;
    next;
    } # end if ($line =~ /^>/)

  if (!$contig)
    {
    print STDERR "Error: Invalid fasta_input_file format: '$fasta_input_file'\n";
    exit 2;
    }
  $line =~ s/\s+//g;
  $sequence .= $line;
  } # end while ($line = <FASTAIN>)
if ($contig)
  {
  process_contig($contig, $comment, $sequence, $Qcomment);
  }
else
  {
  print STDERR "Error: Empty fasta_input_file: '$fasta_input_file'\n";
  }
close(FASTAIN);
close(FASTAOUT);
if ($get_qualities)
  {
  close(QUALOUT);
  close(QUALIN);
  }


############# Print summary statistics if requested #############

if ($Verbose_Mode == 1)
  {
  print STDERR "$GOOD_SELECT_CONTIGS sequences were copied\n";
  }
elsif ($Verbose_Mode > 1)
  {
  print STDERR "$Num_Contigs contigs read\n";
  print STDERR "$Short_Input_Contigs short input contigs ignored\n"
    if $Short_Input_Contigs;
  print STDERR "$Degenerate_Contigs degenerate input contigs ignored\n"
    if $Degenerate_Contigs;
  my $UNPROCESSED_CONTIGS = 0;
  my $PROCESSED_CONTIGS = scalar keys %SELECTIONS;
  unless ($select_all)
    {
    for $contig (keys %SELECTIONS)
      {
      if (! exists $PROCESSED{$contig})
	{
	$UNPROCESSED_CONTIGS++;
	}
      } # end for $contig (keys %SEQ)
    }
  print STDERR "$GOOD_SELECT_CONTIGS sequences from $PROCESSED_CONTIGS select contigs were output\n";
  print STDERR "$UNPROCESSED_CONTIGS select contigs were not used\n"
    if $UNPROCESSED_CONTIGS;
  printf STDERR "%d short select sequences were not written\n",
    $Short_Contigs if ($Short_Contigs);
  }

exit 0;


######################################################################
# process_contig($contig, $comment, $sequence, $Qcomment) - process
#   just read contig.
######################################################################

sub process_contig
  {
  my($contig, $comment, $sequence, $Qcomment) = @_;
  
  my(@select_list);

  # If -u was specified, then look for a uaccno to become the contig
  # name.
  $contig = $1 if ($Use_Uaccno && $comment =~ /\suaccno=(\S+)/);

  # return unless this is a selected contig
  return unless ($select_all || (exists $SELECTIONS{$contig} xor $ignore_list));

  # skip this contig if contig length is < min_size
  my $seq_len = length $sequence;
  if ($seq_len < $min_size)
    {
    $Short_Input_Contigs++;
    return;
    }

  # skip this contig if -d flag and contig sequence does not contain
  # at least one each of A, C, G, and T
  if ($no_degenerate && degenerate($sequence))
    {
    $Degenerate_Contigs++;
    return;
    }
    
  if ($select_all || $ignore_list)
    {
    @select_list = ($select_all_line);
    }
  else
    {
    @select_list  = @{ $SELECTIONS{$contig} };
    }

  # contig selection info is in @{ $SELECTIONS{$contig} }
  foreach my $select_line (@select_list)
    {
    # Each select line is tab-separated values:
    #   (direction (F or R), begin base, length, output contig name)
    my($direction, $begin, $length, $out_contig) =
      split("\t", $select_line);
    my($end);

    my $rev = ($direction eq 'R') ? 1 : 0;

    # adjust begin as needed
    $begin = 1 if (($begin eq '') || ($begin eq '-'));
    if ($begin =~ /^-?\d+$/)			# integer value
      {
      $begin = $seq_len + ($begin + 1) if ($begin < 0);
      $begin = 1 if ($begin <= 0);
      }
    else # (! ($begin =~ /^-?\d+$/))		# not a valid integer
      {
      $begin = 1;
      } # end if (defined $begin && ($begin =~ /^-?\d+$/))
    if ($begin > $seq_len)	# skip contigs shorter than begin base
      {
      $Short_Contigs++;
      next;
      }

    # adjust length and compute end as necessary
    $length = 0 if (($length eq '') || ($length eq '-'));
    if ($length =~ /^-?\d+$/)			# integer value
      {
      if ($length <= 0)
	{
	$end = $seq_len + $length;
	$length = $end + 1 - $begin;
	if ($end < $begin)	# skip contigs ending before beginning
	  {
	  $Short_Contigs++;
	  next;
	  }
	}
      elsif ($begin + $length - 1 <= $seq_len)
	{
	$end = $begin + $length - 1;
	}
      else
	{
	$end = $seq_len;
	$length = $end - $begin + 1;
	}
      }
    else #  (! ($length =~ /^-?\d+$/) )		# not a valid integer
      {
      $end = $seq_len;
      $length = $end - $begin + 1;
      } # end if ((defined $length) && ($length =~ /^-?\d+$/))

    # trim Ns from ends of output contig?
    if ($trim_ends && $length >= $min_size)
      {
      my $trim_len;
      if (substr($sequence, $begin-1, $seq_len) =~ /^(N+)/i)
	{
	$trim_len = length $1;	# trim beginning of output sequence
	$begin += $trim_len;
	$length -= $trim_len;
	}
      if (($length >= $min_size) &&
	  substr($sequence, $begin-1, $seq_len) =~ /(N+)$/i)
	{
	$trim_len = length $1;	# trim end of output sequence
	$end -= $trim_len;
	$length -= $trim_len;
	}
      }

    # skip contigs shorter than min_size
    if ($length < $min_size)
      {
      $Short_Contigs++;
      next;
      }

    # begin == -1 for whole contigs
    $begin = -1 if ($begin <= 1 && $end == $seq_len);

#print "contig=$contig, rev=$rev, begin=$begin, length=$length, end=$end,\n  seq_len=$seq_len, out_contig=$out_contig\n";
    write_contig($contig, $rev, $begin, $length, $end, $seq_len,
      $out_contig, $sequence, $comment, $Qcomment);
    } # end foreach my $select_line (...)
  } # end process_contig


######################################################################
# write_contig($contig, $rev, $begin, $length, $end, $seq_len,
#   $out_contig, $sequence, $comment, $Qcomment) - Output selected
#   contig, or contig fragment.
######################################################################

sub write_contig
  {
  my($contig, $rev, $begin, $length, $end, $seq_len, $out_contig,
     $sequence, $comment, $Qcomment) = @_;
  unless ($Preserve_Comments)
    {
    ($comment, $Qcomment) = ('', '');
    }

  if ($get_qualities)
    {
    unless (exists $qualities_modified{$contig})
      {
      ########## apply quality modifiers ##########

      # constant quality value?
      @Qualities = ($QUALITY_VALUE) x length($sequence)
        if defined $QUALITY_VALUE;

      # subtract a constant value from all qualities?
      @Qualities =
        map { my $t = ($_ - $QUALITY_SUB); ($t >= 0) ? $t : 0 }
          @Qualities if defined $QUALITY_SUB;

      # divide all qualities by a constant value?
      @Qualities = map { int($_ / $QUALITY_DIV) }
        @Qualities if defined $QUALITY_DIV;

      # divide all qualities by stepped divisors?
      if (@QUALITY_STEPS)
	{
	my $step = -1;
	my $run_length = 0;
	my $divisor = 0;
	for (my $i = 0; $i < scalar @Qualities; $i++)
	  {
	  while ($step < (scalar @QUALITY_STEPS) - 1 &&
		 $run_length-- <= 0)
	    {
	    $step++;
	    $run_length = $QUALITY_STEPS[$step];
	    $divisor++;
	    }
	  $Qualities[$i] =
	    int($Qualities[$i] / $divisor);
	  } # end for (my $i = ... )
	} # end if (@QUALITY_STEPS)

      # divide all qualities by sloped divisors?
      if (@QUALITY_SLOPES)
	{
	my $step = -1;
	my $slope_length = 0;
	my $slope_step = 0;
	my $old_divisor = 1;
	my $divisor = 0;
	my $slope_divisor;
	my $num_slopes = (scalar @QUALITY_SLOPES) - 1;
	for (my $i = 0; $i < scalar @Qualities; $i++)
	  {
	  while ($step < $num_slopes && $slope_step >= $slope_length)
	    {
	    $step++;
	    $slope_length = $QUALITY_SLOPES[$step];
	    $old_divisor = $divisor || 1;
	    $divisor++;
	    $slope_step = 0;
	    }
	  if ($step >= $num_slopes && $slope_step >= $slope_length)
	    {
	    $slope_divisor = $divisor;  # used for flat divisor after all slopes
	    }
	  else
	    {		# interpolate slope from $old_divisor to $divisor
	    $slope_step++;
	    $slope_divisor = $old_divisor
	      + ($divisor - $old_divisor) * $slope_step / $slope_length;
	    }
	  $Qualities[$i] = int($Qualities[$i] / $slope_divisor);
	  } # end for (my $i = ... )
	} # end if (@QUALITY_SLOPES)

      # apply maximum allowed quality value
      @Qualities = map { ($_ > $MAX_QUAL) ? $MAX_QUAL : $_ }
        @Qualities;

      ########## end of quality modifiers ##########
      $qualities_modified{$contig} = 1;
      } # end unless (exists $qualities_modified{$contig})
    } # end if ($get_qualities)


  if ($begin > 0)	# output partial contig
    {
    if ($out_contig eq '-') # if Output Contig name not specified
      {
      $out_contig = $contig;
      if ($rev)	# reverse and complement?
        {
        $out_contig .= '.comp' if ($out_contig !~ s/\.comp$//);
        }
      $out_contig .= " (${begin}..${end})" unless ($trim_ends);
      }
    output_contig($contig, $out_contig . $comment, $rev,
		  substr($sequence, $begin - 1, $length),
		  $begin, $end, $seq_len);
    if ($get_qualities)
      {
      output_contig_qual($contig, $out_contig . $Qcomment, $rev,
			 $begin, $end,
			 @Qualities[$begin - 1 .. $end - 1]);
      } # end if ($get_qualities)
    }
  else # ($begin < 0)	# output entire input contig
    {
    if ($out_contig eq '-') # if Output Contig name not specified
      {
      $out_contig = $contig;
      if ($rev)	# reverse and complement?
        {
        $out_contig .= '.comp' if ($out_contig !~ s/\.comp$//);
        }
      }
    output_contig($contig, $out_contig . $comment, $rev, $sequence,
		  $begin, $end, $seq_len);
    if ($get_qualities)
      {
      output_contig_qual($contig, $out_contig . $Qcomment, $rev,
			 $begin, $end, @Qualities);
      }
    } # end if ($begin > 0)
  $GOOD_SELECT_CONTIGS++;
  $PROCESSED{$contig}++;
  } # end write_contig


######################################################################
# output_contig - output a contig sequence.
######################################################################

sub output_contig
  {
  my($contig, $out_contig, $rev, $sequence, $begin, $end, $contig_len) = @_;
  my($len) = length($sequence);
  my($segment, $i, $direction, $comp, $suffix, $CONTIG_POSITION);
  if ($rev)	# reverse and complement?
    {
    $sequence = reverse($sequence);
    $sequence =~ tr/acgtACGT/tgcaTGCA/;
    $direction = 'R';
    }
  else
    {
    $direction = 'F';
    }
  print FASTAOUT ">$contig_prefix$out_contig\n";
  for ($i = 0; $i < $len; $i += $output_line_len)
    {
    $segment = substr($sequence, $i, $output_line_len);
    print FASTAOUT "$segment\n";
    } # end for ($i ... )
  } # end output_contig


###########################################################################
# output_contig_qual - output contig qualities.
###########################################################################

sub output_contig_qual
  {
  my($contig, $out_contig, $rev, $begin, $end, @quals) = @_;
  my($len) = scalar @quals;
  my($i, $segment, $comp, $suffix);
  @quals = reverse(@quals) if ($rev);
  print QUALOUT ">$contig_prefix$out_contig\n";
  $segment = '';
  for ($i = 0; $i < $len; $i++)
    {
    if (length($segment) >= $output_line_len)
      {
      print QUALOUT "$segment\n";
      $segment = '';
      }
    $segment .= "$quals[$i] ";
    } # end for ($i ... )
  print QUALOUT "$segment\n" if (length($segment));
  } # end output_contig_qual


###########################################################################
# degenerate($sequence) - return 1 if $sequence does not contain at least
#   one each of A, C, G, and T;  else return 0.
###########################################################################

sub degenerate
  {
  my($sequence) = @_;
  return(($sequence !~ /A/i || $sequence !~ /C/i ||
          $sequence !~ /G/i || $sequence !~ /T/i) ? 1 : 0);
  } # end degenerate


###########################################################################
# display_help
###########################################################################

sub display_help
  {
  my($msg) = @_;
  print STDERR "\n$msg\n" if $msg;
  print STDERR <<EOF;

USAGE: $my_name [-n select_file] [-c contig_prefix] [-d] [-e]
	    [-g] [-i] [-l output_line_length] [-m min_size] [-M max_qual]
	    [-p] [-q] [-Q qual_value] [-s] [-t type_to_remove] [-u]
	    [-v] [-V] [-x] fasta_input_file fasta_output_file
              or
       $my_name -h

OPTIONS:

  -c  Contig prefix to be added to output contig names.

  -d  Input contigs are assumed to be dna.  Filter out any degenerate
      contigs that do not contain at least one each of A, C, G, and T.

  -e  Emtpy output files are OK and do not result in an error.

  -g  Remove runs of Ns from ends of contigs.  Minimum contig length
      is enforced after trimming the ends.

  -i  Ignore contigs specified in the 'select_file' (-n), and output
      the input contigs whose names are NOT listed.  May only be used
      with -n.

  -l  Specify output line length.

  -m  Specify minimum contig length to be written.

  -M  Specify a maximum quality score for all bases in the output
      quality file.

  -n  Specifies the name of a 'select_file', which contains a list of
      contig names to be output.  Each line of 'select_file' is a tab
      separated list.  The fields in the list are:  contig_name,
      direction, begin_base, and length.  If -i is specified, then
      the list is the list of contigs to be ignored, and only the
      contig_name field is used.

  -p  Preserve contig header comments.  May not be used with -u.

  -q  Also process a fasta quality file ('fasta_input_file'.qual), as
      well as the sequence file and create an output quality file
      ('fasta_output_file'.qual) in addition to the output sequence
      file.

  -Q  Specify a constant quality value to be applied to all bases in
      the output quality file or a modifier to be applied to all
      qualities from the input Fasta quality file.  If 'qual_value' is
      a simple one- or two-digit positive integer, then that value is
      used for the quality scores and the input Fasta quality file is
      not needed.  If 'qual_value' is not just a simple one- or
      two-digit integer, then it specifies a modifier to be applied to
      the values from the input Fasta quality file.

  -s  The contig name may be shortened by removing any prefix before
      the word "Contig", i.e., "gono.fasta.screen.Contig26" becomes
      "Contig26".

  -t  Specify a filetype to be removed before adding ".qual" to create
      the output quality filename.

  -u  Use universal accession numbers (uaccno) as contig names for 454
      reads.  May not be used with -p.

  -v  Verbose mode - print to STDERR the number of contigs copied.  If
      both -v and -V are specified, then -V will be used.

  -V  Verbose mode - print out some statistics to STDERR while running.
      If both -v and -V are specified, then -V will be used.

  -x  Create new or append to existing (extend) output files.

EOF

  exit 2;
  } # end display_help


###########################################################################
# display_help
###########################################################################

sub display_more_help
  {
  print STDOUT <<EOF;

$my_name - Select a subset of contigs from a fasta input file.
The output can be a set of contigs, or a single contig with contig
separators. A fasta quality file also can be processed.  If '-n' is
not specified, then all contigs are examined, but if the '-n' flag
specifies a 'select_file', then the output contig order matches the
input order, but only the contigs named in the 'select_file' are
used.  Using the '-n' flag, the selected contigs can be reversed and
complemented, and partial contigs can be extracted.  With the '-m'
flag, input and output contigs shorter than 'min_size' bases long
are discarded.

The contig names optionally may be shortened by removing everything
before the word "Contig" (-s).

Contig header comments, located after the contig name on the contig
header line, will be removed, but new comments may be added if partial
contigs are requested in 'select_file'.  If '-p' is specified, then the
original contig header comments will be added back to the contig
header (possibly after any added comments describing partial contigs).
All output lines containing sequence bases (not contig headers) are
reformatted to be 'output_line_length' bases long if possible.

If '-u' is specified, then universal accession numbers (uaccno) are
used as contig names for 454 reads instead of the rank_x_y form
present in older 454Reads files.  The uaccnos will be used in the
output contigs and in any select files used to select contigs.  The
'-p' flag may not be used with '-u'.

If '-q' is specified, then a fasta quality file named
"'fasta_input_file'.qual" also is processed, creating an additional
output fasta quality file named "'fasta_output_file'.qual".  If
'-t type_to_remove' is also specified, then the 'type_to_remove' is
removed from the end of 'fasta_output_file' before adding ".qual".


USAGE: $my_name [-n select_file] [-c contig_prefix] [-d] [-e]
	    [-g] [-i] [-l output_line_length] [-m min_size] [-M max_qual]
	    [-p] [-q] [-Q qual_value] [-s] [-t type_to_remove] [-u]
	    [-v] [-V] [-x] fasta_input_file fasta_output_file
              or
       $my_name -h           <== What you are reading

  where 'select_file' is the name of an input file containing the
            names of the contigs to be output.

        'output_line_length' is the length of the output lines
	    containing sequence data (not sequence headers).  The
	    default value is $DEFAULT_LINE_LENGTH.  'output_line_length' must be greater
	    or equal to $MIN_LINE_LENGTH.

        'max_qual' is the maximum allowed quality value for all bases
            in the output file.  Any bases with a quality value higher
            than 'maximum_quality_value' will have the quality lowered
	    to 'maximum_quality_value'.  'max_qual' must be a one- or
	    two-digit positive integer.  The default value is $MAX_QUAL.

	'min_size' is the minimum contig length to be used.  Shorter
	    contigs are discarded.  The default value is $min_size.

        'qual_value' is one of the following:

            (1) '#' - a one- or two-digit integer constant quality
	      value to be applied to all bases in the phd file.  If
	      this is specified, then an input Fasta quality file is
	      not used;  or

            (2) 'SUB#' - the word 'SUB' followed by a one- or
	      two-digit integer which is to be subtracted from all
	      scores from the input Fasta quality file;  or

            (3) 'DIV#.###' - the word 'DIV' followed by an integer or
	      real number (>= 1 and < 10), by which all of the quality
	      scores from the input Fasta quality file are to be
	      divided;  or

            (4) 'STEP#,#,...' - the word 'STEP' followed by a
	      comma-separated list of integers representing segment
	      lengths of bases.  This format divides input quality
	      scores by a divisor like the DIV format above, except
	      that the divisor is variable.  Each number gives the
	      number of base quality scores that are to be divided by
	      the same divisor.  The first number gives the number of
	      bases to be divided by 1; then the second number of base
	      qualities are to be divided by 2;  then 3;  etc.;  or

            (5) 'SLOPE#,#,...' - the word 'SLOPE' followed by a
	      comma-separated list of integers representing segment
	      lengths of bases.  This format divides input quality
	      scores by a divisor like the DIV format above, except
	      that the divisor is variable.  Each number gives the
	      number of base quality scores that are to be divided by
	      a set of divisors.  The first number gives the number of
	      base qualities to be divided by 1;  then the second
	      number of base qualities are to be divided by a linearly
	      interpolated set of divisors sloping from 1 to 2;  then
	      2 to 3;  then 3 to 4;  etc.  After the last segment, all
	      the rest of the base qualities are divided by the number
	      of segments.

            Only one of these forms may be used at a time.  If a
	    negative base quality score results, then that value is
	    replaced by zero.  When base quality scores are divided,
	    the resulting quotient is rounded down to the next lowest
	    integer.

	'type_to_remove' is a filetype to be removed before adding
	    ".qual" to create the output quality filename.  If the
	    filetype is not present, then no error or warning is
	    given.  Note:  you need to include the leading period for
	    the type to remove (e.g., '.fna').  The default type to be
	    removed is '$Type_Strip'.

        'fasta_input_file' is the name of the input sequence file in
	    Fasta format, which contains the contigs to be processed.
	    If 'fasta_input_file' is specified as '-', then standard
	    input is used.  If '-q' is specified, then
	    'fasta_input_file' may not be specified as '-'.

        'fasta_output_file' is the name of the output sequence file in
            Fasta format.  If 'fasta_output_file' is specified as '-',
            then standard output is used.  If '-q' is specified, then
	    'fasta_output_file' may not be specified as '-'.

OPTIONS:

  -c  Contig prefix to be added to output contig names.  The prefix is
      added after possible shortening by -s.  May not be used with -n.

  -d  Input contigs are assumed to be dna.  Filter out any degenerate
      contigs that do not contain at least one each of A, C, G, and T.

  -e  Emtpy output files are OK and do not result in an error.  If -e
      is not specified, and the input file is empty or an select file
      selects no contigs, then the program normally exits with an
      error.

  -g  Remove runs of Ns from ends of contigs.  Minimum contig length
      is enforced after trimming the ends.

  -i  Ignore contigs specified in the 'select_file' (-n), and output
      the input contigs whose names are NOT listed.  May only be used
      with -n.

  -l  Specify output line length.

  -m  Specify minimum contig length to be written.

  -M  Specify a maximum quality score for all bases in the output
      quality file.

  -n  Specifies the name of a 'select_file', which contains a list of
      contig names to be output.  Each line of 'select_file' is a tab
      separated list.  The fields in the list are:  contig_name,
      direction, begin_base, and length.  The fields present may be
      followed by a comment, which begins with a '#' character after a
      white space character and runs to the end of the line.  If '-i'
      is specified with -n, then the list is the list of contig names
      to be ignored, and only the contig_name field is used.  If '-n'
      is not specified, then the full sequences of all input contigs
      are used, subject to modifications by other flags.

  -p  Preserve contig header comments.  May not be used with -u.

  -q  Also process a fasta quality file, "'fasta_input_file'.qual", as
      well as the sequence file and create an output quality file,
      "'fasta_output_file'.qual", in addition to the output sequence
      file.  $my_name also allows relaxed quality file naming.
      If "'fasta_input_file'.qual" is not found, and
      'fasta_input_file' is of the form "xxx.fa", "xxx.fna", or
      "xxx.fasta", then $my_name will try to use "xxx.qual" instead.

  -Q  Specify a constant quality value to be applied to all bases in
      the output quality file or a modifier to be applied to all
      qualities from the input Fasta quality file.  If 'qual_value' is
      a simple one- or two-digit positive integer, then that value is
      used for the quality scores and the input Fasta quality file is
      not needed.  If 'qual_value' is not just a simple one- or
      two-digit integer, then it specifies a modifier to be applied to
      the values from the input Fasta quality file.  See the
      dexription of 'qual_value' above.

  -s  The contig name may be shortened by removing any prefix before
      the word "Contig", i.e., "gono.fasta.screen.Contig26" becomes
      "Contig26".

  -t  Specify a filetype to be removed before adding ".qual" to create
      the output quality filename.  See 'type_to_remove' above.  The
      '-t' flag is not used unless '-q' is also specified.

  -u  Use universal accession numbers (uaccno) as contig names for 454
      reads instead of the rank_x_y form normally present in 454Reads
      files.  May not be used with -p.

  -v  Verbose mode - print to STDERR the number of contigs copied.  If
      both -v and -V are specified, then -V will be used.

  -V  Verbose mode - print out some statistics to STDERR while running.
      If both -v and -V are specified, then -V will be used.

  -x  Create new or append to existing (extend) output files.


SELECT_FILE INPUT FORMAT:

Blank lines and comment lines (beginning with #) are allowed.  Each
non-blank, non-comment, line of 'select_file' should name one input
contig to be written to the output file.  The same input contig name
may be used more than once.  This is useful for splitting an input
contig into multiple output contigs.  Each input line may have
multiple white-space separated fields.  Leading white-space is
ignored.  Trailing fields may be omitted to use the default values.
After any fields that are present, a trailing comment may be added.

  - Input contig name.  This field is required.
  - Direction (forward or reverse & complement).  If the field begins
    with a C or an R (of either case), then that contig is reversed
    and complemented before it is output, and the Output contig name
    contains the Input contig name suffixed with '.comp' (or '.comp'
    is removed), unless Output contig name is specified explicitly.
    The default is to output the contig in the forward direction.
  - Beginning base of input contig.  The default beginning base is 1.
    A negative value indicates to start that many bases before the end
    of the sequence.  (-1 means start at the last base;  -10 means
    start at the 10th base counting from the end of the sequence,
    etc.)  A value of zero or just a minus sign '-' is treated the
    same as base 1.  A value longer than the input sequence is treated
    as a short contig and ignored.  The value for Beginning base
    refers to base positions in the original input contig before any
    base reversal and complement operation or end-trimming that may
    occur.
  - Length of sequence to be output.  A negative value indicates to
    end that many bases before the end of the sequence.  (-1 means end
    1 base before the last base;  -10 means end 10 bases before the
    last base of the sequence, etc.)  If the length is omitted or zero
    or out of range or less than Beginning base or just a minus sign
    '-', then the end of the contig is used.
  - Output contig name.  This field may be used to name the resulting
    output contig.  If this field is omitted, then an output contig
    name will be constructed from the Input contig name.
  - Trailing comment.  The comment begins with a # character and must
    be at the beginning of the line or the # must be immediately
    preceded by white space.  The # character and and all characters
    after it are ignored until the end of the line.

Assuming Contig1 is 1000 bases long, then each of the following lines
specify that all 1000 bases of Contig1 should be output in the forward
direction:

Contig1				# All of Contig1
Contig1              F		# All of Contig1 on the forward strand
Contig1 F 1 1000		# The first 1000 bases (of 1000)
    Contig1   F  1 1001		# The first 1001 bases (of only 1000)
Contig1	F	1	0	# Start at base 1 through the last base
Contig1	F	1		# Start at base 1 through the end
Contig1	f	1	-	# Start at base 1 through the end
Contig1	F 0 0 Contig1		# All of Contig1, named as Contig1
Contig1	- - - Contig1		# All of Contig1, named as Contig1

The following lines would specify as follows:

Contig9				# All of Contig9
Contig1	C			# All of Contig1 reversed and
				# complemented
Contig3	R	1	100	# First 100 bases of Contig3, then
				# reversed and complemented
Contig2	F	-100		# Last 100 bases of Contig2
Contig7	F	51	-50	# All but the first and last 50 bases
				# on the ends of Contig7
Contig5	F	-100	-50	# The 50 bases ending 50 bases before
				# the end of Contig5
Contig2	f 501 1000 Contig2:501-1000	# Second 500 bases of Contig2.
				# The output contig is named
				# "Contig2:501-1000".
Contig4	f  123 334  Fred	# Bases 123-456 of Contig4. The output
				# contig is named "Fred".
Contig6	r - - Contig6.revcomp	# All of Contig6 reversed and
				# complemented and named
				# "Contig6.revcomp"


TEXAMPLES:

For the following examples, assume 'input.fa' and 'input.fa.qual'
contain the following contigs:

Contig1 (100 bases), Contig2 (200 bases), ... Contig101 (10100 bases)

\$ echo "Contig1 F 1 1000" | $my_name -n - input.fa Contig1_1000.fa

reads the input fasta sequence file 'input.fa' and writes a new fasta
sequence file, 'Contig1_1000.fa', containing only the first 1000 bases
of contig 'Contig1, as specified from STDIN piped from the echo
command.

\$ $my_name -n select -m 2000 -v input.fa long_ones.fa

reads the input fasta sequence file 'input.fa' and outputs a new fasta
sequence file, 'long_ones.fa', containing only the selected contigs
that are at least 2000 bases long.  The output contig order is the
same as in the input fasta file, but only the requested contigs are
written.  The Verbose flag (-v) causes some statistics to be written
to Standard Output, including counts of contigs read, contigs written,
and select contigs not written to the output file, if any.

\$ $my_name -Q 20 contigs.fa contigs20.fa

reads the fasta sequence file "contigs.fa" (and does NOT use a fasta
quality file "contigs.fa.qual" or "contigs.qual") and produces new
fasta sequence and quality files named "contigs20.fa" and
"contigs20.fa.qual".  All base quality values are set to 20.

\$ $my_name -Q s20 -t .fa contigs.fa contigs_sub20.fa

reads the fasta sequence file "contigs.fa" and a fasta quality file
"contigs.fa.qual" (or "contigs.qual") and produces new fasta sequence
and quality files named "contigs_sub20.fa" and "contigs_sub20.qual".
All base quality values are reduced by 20, but any qualities which
would be negative are set to zero.  Note that the output quality file
is "*.qual", not "*.fa.qual", because of the "-t .fa" option.

\$ $my_name -Q d2.5 contigs.fa contigs_div2.5.fna

reads the fasta sequence file "contigs.fa" and the fasta quality file
"contigs.fa.qual" (or "contigs.qual") and produces new fasta sequence
and quality files named "contigs_div2.5.fna" and
"contigs_div2.5.fna.qual".  All base quality values are divided by 2.5
and rounded down to the next integer.

\$ $my_name -Q step100,50,25 contigs12.fna contigs_3step.fasta

reads the fasta sequence file "contigs12.fna" and the fasta quality
file "contigs12.fna.qual" (or "contigs12.qual") and produces new fasta
sequence and quality files named "ccontigs_3step.fasta" and
"contigs_3step.fasta.qual".  The first 100 base quality values are
divided by 1;  the next 50 are divided by 2;  the next 25 are divided
by 3;  the remaining base qualities are divided by 4. The resulting
quality scores are all rounded down to the next lowest integer.

\$ $my_name -Q slope10,5,0,4,2 contigs13.fna contigs_5slope.fasta

reads the fasta sequence file "contigs13.fna" and the fasta quality
file "contigs13.fna.qual" (or "contigs13.qual") and produces new fasta
sequence and quality files named "ccontigs_3step.fasta" and
"contigs_3step.fasta.qual".  The base quality scores are divided by:
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1.2, 1.4, 1.6, 1.8, 2.0,
  (note the jump from 2 because of the zero)
  3.25, 3.5, 3.75, 4,
  4.5, 5,
  5, 5, ...
The resulting quality scores are all rounded down to the next lowest
integer.


DATE LAST MODIFIED: $Date_Last_Modified

EOF

  exit;

  } # end display_more_help

