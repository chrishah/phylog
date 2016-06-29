#!/usr/bin/perl -w

####!/usr/local/bin/perl -w
# sort_contigs - Sort a fasta input file.  The output can be a set
#   of ordered contigs, or a single contig with contig separators. A
#   fasta quality file also can be processed.  Contigs can be reversed
#   and complemented, and partial contigs can be extracted.
#
# NOTE: If you do NOT need to change the order of the output contigs,
#   (using options -a, -b, -n, -o order_file, -S, or, -z) or join
#   output contigs as one contig (option -j join_char), then the
#   program select_contigs may be a better choice and run faster,
#   especially for large input files.
#
# Written by: James D. White, University of Oklahoma, Advanced Center for
#   Genome Technology
#
# Date Written: Jul 26, 2000
#
# 20090807 JDW - Add information about the program select_contigs as
#		 an alternative which may be faster.
# 20090625 JDW - Add -g flag for trimming Ns from ends of contigs for
#		 Genbank submissions.
# 20090518 JDW - Add missing help info for -x flag.
# 20090424 JDW - Add missing help info for -e flag.
# 20090422 JDW - Add -M and -Q flags for modifying output qualities.
# 20090414 JDW - Add -e flag to not exit with error when there are no
#		 contigs to be output.
# 20090319 JDW - Add -x flag to append to output files.
# 20090312 JDW - Add -c flag for adding a prefix to output contig names.
# 20081109 JDW - Fix help routine to exit after outputing the help.
# 20081104 JDW - Add -S flag for sorting by sequence.
# 20080908 JDW - Add -d flag to filter out degenerate contigs.
# 20080904 JDW - Allow -m flag to be used with -o flag.  Also skip
#		 contigs where begin is past end of contig or negative
#		 length value specifies an end before begin base.
# 20080731 JDW - Add -t flag to strip old file type before adding
#		 ".qual" to make output quality filename.  Add -u
#		 flag to use universal accession numbers (uaccno) as
#		 contig names for 454 reads instead of the rank_x_y
#		 form normally present in 454Reads files.  Add -p
#		 flag to preserve contig header comments.  Simplify
#		 flag processing code by using Getopt::Std.
# 20071018 JDW - Fix bogus error for missing quality file when -q
#		 was not specified.
# 20071005 JDW - Allow relaxed input quality file naming.
# 20070214 JDW - Add '-m minsize' option to exclude contigs shorter
#		 than 'minsize' long.
# 20050818 JDW - Fix bug in output qual file naming.
# 20050615 JDW - Allow output contig name on order file lines.
# 20041110 JDW - Now allows trailing comments on order file lines.

our $Date_Last_Modified = "August 7, 2009";


######################## Customization Parameters ###########################

our $DEFAULT_LINE_LENGTH = 50;	# Default length of fasta output sequence lines.
our $output_line_len = $DEFAULT_LINE_LENGTH;# Length of fasta output sequence lines.
our $MIN_LINE_LENGTH = 10;	# Mimimum allowed value for $output_line_len.

				# Output contig order?
				#   Default is to output contigs in original
				#   input order.  Order may be modified by -b.
our $alpha_sort = 0;		# Sort input contigs in alphanumeric order?
				#   (-a)  Order may be modified by -b.
our $backwards_sort = 0;	# Sort input contigs in backwards order?
				#   (May modify default or -a or -n or -o)
our $contig_prefix = '';	# Prefix to add to output contig names (-c)
our $no_degenerate = 0;		# Input is dna.  Filter out any contigs
				#   that do not contain at least one
				#   each of A, C, G, and T.  (-d)
our $empty_ok = 0;		# Empty output files with no contigs
				#   are OK?			(-e)
our $trim_ends = 0;		# Trim Ns from ends of contigs	(-g)
our $min_size = 1;		# Ignore contigs < $min_size in length. (-m)
our $MAX_QUAL = 99;		# Maximum quality value allowed in
				#   output file (-M)
our $numeric_sort = 0;		# Sort input contigs in numeric order?
				#   (-n)  Order may be modified by -b.
our $order_file = '';		# Name of file to read for contig ordering.
				#   (-o)  Order may be modified by -b.
				#   May not be used with -r.
our $size_sort = 0;		# Sort input contigs in size order from
				#   smallest to largest (-z).

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
our $reverse_all_contigs = 0;	# Reverse and complement all input contigs?
				#   (-r)  May not be used with -o.
our $shorten_contig_name = 0;	# Is contig name shortened? (-s)
our $sequence_sort = 0;		# Sort input contigs in sequence order (-S)
our $Type_Strip = '';		# Type to strip from output filename
				#   before adding ".qual" (-t)
our $Use_Uaccno = 0;		# Use universal accession numbers
				#   (uaccno) as contig names for 454
				#   reads. (-u)
our $Verbose_Mode = 0;		# Are some statistics to be listed? (see -v)
our $join_contigs = 0;		# Join the contigs into one big contig (See -j)
our $join_char = '';		#   Character to be duplicated for join string
				#   (If multi-character string is used,)
				#   (then the string will be used as is)
our $join_len = 50;			#   Length of string to join contig ends
#$join_string = "$join_char" x $join_len;   # String for joining contig ends
#@JOIN_QUALS = split(//, '0' x $join_len);  # Dummy qualities for joined string
our $Extend_File = 0;		# Extend (append to) existing output
				#   file(s) (see -x)
########## Operating System specific ##########

my $directory_separator = '/';	# For Unix
#$directory_separator = '\';	# For DOS/Windows

##################### End of Customization Parameters #######################


my $full_command_name = $0;
our($my_name);
if ($full_command_name =~ m"(^|${directory_separator})([^${directory_separator}]*)$")
  {
  $my_name = $2;
  }
else
  {
  $my_name = 'sort_contigs';
  }

  
############# Process flags on command line #############
use Getopt::Std;
our($opt_a, $opt_b, $opt_c, $opt_d, $opt_e, $opt_g, $opt_h, $opt_j,
    $opt_l, $opt_m, $opt_M, $opt_n, $opt_o, $opt_p, $opt_q, $opt_Q,
    $opt_r, $opt_s, $opt_S, $opt_t, $opt_u, $opt_v, $opt_x, $opt_z) =
  ('') x 24;
display_help('') unless( getopts('abc:deghj:l:m:M:no:pqQ:rsSt:uvxz') );
display_more_help() if ($opt_h);		# -h	(help)

$Verbose_Mode = 1 if ($opt_v);			# -v	(verbose mode)

print STDERR "\n$my_name - Last Modified: $Date_Last_Modified\n\n"
  if $Verbose_Mode;

$alpha_sort = 1 if ($opt_a);			# -a
$backwards_sort = 1 if ($opt_b);		# -b
$contig_prefix = $opt_c if ($opt_c ne '');	# -c contig_prefix
$no_degenerate = 1 if ($opt_d);			# -d
$empty_ok = 1 if ($opt_e);			# -e
$trim_ends = 1 if ($opt_g);			# -g
if ($opt_j ne '')				# -j join_char
  {
  $join_contigs = 1;
  $join_char = $opt_j;
  if (length $join_char == 1)
    {
    $join_string = "$join_char" x $join_len;  # String for joining contigs
    }
  elsif (length $join_char > 1)
    {
    $join_string = "$join_char";	# String for joining contigs
    $join_len = length $join_string;	# Adjust separator length
    }
  @JOIN_QUALS = split(//, ('0' x $join_len)); # zero quality filler
  }
$output_line_len = $opt_l if ($opt_l ne '');	# -l output_line_length
$min_size = $opt_m if ($opt_m ne '');		# -m min_size
$MAX_QUAL = $opt_M if ($opt_M ne '');		# -M max_qual
$numeric_sort = 1 if ($opt_n);			# -n
$order_file = $opt_o if ($opt_o ne '');		# -o order_file
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
  elsif ($opt_Q =~ /^S(?:UB)?(\d{1,2})$/i)	# subtract $QUALITY_SUB from qualities
    {
    $QUALITY_SUB = $1;
    }
  elsif ($opt_Q =~ /^D(?:IV)?([1-9](\.\d*)?)$/i) # divide qualities by constant values
    {
    $QUALITY_DIV = $1;
    }
  elsif ($opt_Q =~ /^STEP([\d,]+)$/i)	# divide qualities by stepped values
    {
    @QUALITY_STEPS = map { (defined $_ && $_ ne '') ? $_ : 0 }
      split(',', $1);
    }
  elsif ($opt_Q =~ /^SLOPE([\d,]+)$/i) # divide qualities by sloped values
    {
    @QUALITY_SLOPES = map { (defined $_ && $_ ne '') ? $_ : 0 }
      split(',', $1);
    }
  else
    {
    $Qmsg = "Invalid -Q qual_value='$opt_Q'\n";
    $get_qualities |= 4;
    }
  $get_qualities = ($get_qualities == 2) ? 2 :	# adjust $get_qualities
		   ($get_qualities == 4) ? 0 : 1;
  } # end if ($opt_Q ne '')
$reverse_all_contigs = 1 if ($opt_r);		# -r
$shorten_contig_name = 1 if ($opt_s);		# -s
$sequence_sort = 1 if ($opt_S);			# -S
$Type_Strip = $opt_t if ($opt_t ne '');		# -t type_to_remove
$Type_Strip =~ s/\./\\./g;	# make '.' mean a period
$Use_Uaccno = 1 if ($opt_u);			# -u
$Extend_File = 1 if ($opt_x);			# -x
$size_sort = 1 if ($opt_z);			# -z

$global_direction = $reverse_all_contigs ? 'R' : 'F';	# used if not -o

############# check for invalid flags #############
my $help_msg = '';
$help_msg = "Cannot use -a, -n, -o, -S, or -z flags with each other\n"
  if (($alpha_sort + $numeric_sort + $size_sort + $sequence_sort +
      (($order_file ne '' ? 1 : 0)) > 1)) ;
$help_msg .= "Cannot use -c with -j or -o flags\n"
  if ($contig_prefix ne '' && ($join_contigs || $order_file ne ''));
$help_msg .= "Invalid -c flag, contig_prefix='$contig_prefix'\n"
  if ($contig_prefix =~ /[^-\w]/);
$help_msg .= "Cannot use -j flag with -x flag.\n"
  if ($join_contigs && $Extend_File);
$help_msg .= "Invalid -M max_qual='$MAX_QUAL', must be in range 1-99\n"
  unless ($MAX_QUAL =~ /^\d{1,2}$/ && $MAX_QUAL > 0);
$help_msg .= $Qmsg;
$help_msg .= "Cannot use -p and -u flags together\n"
  if ($Preserve_Comments + $Use_Uaccno > 1);
$help_msg .= "Cannot use -o flag with -r\n"
  if (($order_file ne '') && $reverse_all_contigs);
$help_msg .= "Invalid value for -l flag, output_line_length='$output_line_len'\n"
  if (($output_line_len !~ /^\d+$/) ||
      ($output_line_len < $MIN_LINE_LENGTH));
$help_msg .= "Invalid value for -m flag, min_size='$min_size'\n"
  if ($min_size !~ /^\d+$/);
$min_size = 1 if ($min_size == 0);
$help_msg .= "For -o flag, missing or invalid order_file='$order_file'\n"
  if ($order_file ne '' && $order_file ne '-' && ! -f $order_file);

############# Check command line arguments #############
$fasta_input_file = shift @ARGV;
$help_msg .= "Missing 'fasta_input_file' name.\n"
  if (! defined $fasta_input_file || $fasta_input_file eq '');
$help_msg .= "'fasta_input_file' name specified as '-' (STDIN), which is not\n" .
	     "    allowed when '-q' or a '-Q' modifier is specified.\n"
  if ((defined $fasta_input_file) && ($fasta_input_file eq '-') &&
      ($get_qualities == 1));

$fasta_output_file = shift @ARGV;
$help_msg .= "Missing 'fasta_output_file' name.\n"
  if (! defined $fasta_output_file || $fasta_output_file eq '');
$help_msg .= "'fasta_output_file' name specified as '-' (STDOUT), which is not\n" .
	     "    allowed when '-j', '-q', '-Q', or '-x' is specified.\n"
  if ((defined $fasta_output_file) && ($fasta_output_file eq '-') &&
      ($get_qualities || $join_contigs || $Extend_File));

########### if errors, print error(s) and usage, then exit ###########
display_help($help_msg) if ($help_msg ne '');


######## Read and process 'order_file' and save in @ORDERING, if -o ########
@ORDERING = ();
if ($order_file ne '')
  {
  open(ORDERIN, $order_file) || die("Can't open order_file: '$order_file'\n");
  $order_lines = 0;
  while (defined ($line = <ORDERIN>))
    {
    chomp($line);
    $line =~ s/^\s+//;
    next if $line =~ /^#/;	# ignore comment lines
    $line =~ s/\s+#.*$//;	# remove trailing comments
    $line =~ s/\s+$//;
    next if $line eq '';
    push @ORDERING, $line;
    $order_lines++;
    } # end while (<ORDERIN>)
  close(ORDERIN);
  print "$order_lines order lines read from order_file\n" if $Verbose_Mode;
  }

############# Read contigs and save them in hashes %SEQ and %QUAL #############
########## Also save input contig order in @ORDERING if no 'order_file' #######
$fasta_qual_file = "${fasta_input_file}.qual";
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

%SEQ = ();
%QUAL = ();
%COMMENTS = ();
%Short_Contigs = ();
$Num_Contigs = 0;
$Short_Contigs = 0;
$Degenerate_Contigs = 0;
$header = '';
$contig = '';
$comment = '';
$sequence = '';
$line_num = 0;
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
      save_contig($contig, $comment, $sequence, $Qcomment);
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
      my @Qualities = split(' ', $Quality);
      $Qualref = \@Qualities;
      $Qlen = scalar @Qualities;
      } # end if ($get_qualities == 1)
					# Remove Contig name prefix?
    $contig =~ s/^.*Contig/Contig/ if $shorten_contig_name;
    next;
    } # end if ($line =~ /^>/)

  if (!$contig)
    {
    print STDERR "Error: Invalid fasta_input_file format: '$fasta_input_file'\n";
    exit 2;
    }
  $line =~ s/\s//g;
  $sequence .= $line;
  } # end while ($line = <FASTAIN>)
if ($contig)
  {
  save_contig($contig, $comment, $sequence, $Qcomment);
  }
else
  {
  print STDERR "Error: Empty fasta_input_file: '$fasta_input_file'\n";
  }
close(FASTAIN);
close(QUALIN) if ($get_qualities == 1);
if ($Verbose_Mode)
  {
  print STDERR "$Num_Contigs contigs read\n";
  print STDERR "$Short_Contigs short input contigs ignored\n"
    if $Short_Contigs;
  print STDERR "$Degenerate_Contigs degenerate input contigs ignored\n"
    if $Degenerate_Contigs;
  }

####################### Check ordering values #######################
# Input is in @ORDERING: Each line is white-space separated values:
#   (contig name, direction (F or C/R), begin base, length).  All but
#   contig name may be omitted.
@ORDERING1 = ();
$Short_Contigs = 0;
$BAD_ORDER_CONTIGS = 0;
for $contig_line (@ORDERING)
  {
# input is (input contig name, direction (F or C/R), begin base, length,
#   output contig name)
  my($contig, $direction, $begin, $length, $out_contig) =
    split(' ', $contig_line);
  $out_contig = '-' if $join_contigs || ! defined $out_contig;
  if (exists $SEQ{$contig})	# requested contig exists
    {
    my($end, $rev);
    my $seq_len = length $SEQ{$contig};

    $rev = (defined $direction && $direction =~ /^[rc]/i) ? 1 : 0;

    # adjust begin as needed
    $begin = 1 if ((! defined $begin) || ($begin eq '-'));
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
    $length = 0 if ((! defined $length) || ($length eq '-'));
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
      if (substr($SEQ{$contig}, $begin-1, $seq_len) =~ /^(N+)/i)
	{
	$trim_len = length $1;	# trim beginning of output sequence
	$begin += $trim_len;
	$length -= $trim_len;
	}
      if (($length >= $min_size) &&
	  substr($SEQ{$contig}, $begin-1, $seq_len) =~ /(N+)$/i)
	{
	$trim_len = length $1;	# trim end of output sequence
	$end -= $trim_len;
	$length -= $trim_len;
	}
      }

    # save it if long enough
    if ($length >= $min_size)
      {
      $begin = -1 if ($begin <= 1 && $end == $seq_len);
      $line = join("\t", $contig, $rev, $begin, $length, $end,
        $seq_len, $out_contig);
      push @ORDERING1, $line;
      }
    else
      {
      $Short_Contigs++;
      }
    }
  elsif (exists $Short_Contigs{$contig})
    {			# requested contig is discarded short contig
    $Short_Contigs++;
    }
  else # (! exists $SEQ{$contig})
    {				# requested contig does not exist
    print STDERR "Invalid contig name ($contig) read from 'order_file'\n";
    $BAD_ORDER_CONTIGS++;
    } # end if (exists $SEQ{$contig})
  } # end for $contig_line (@ORDERING)
printf STDERR "%d short%s contigs not written\n",
  $Short_Contigs, ($no_degenerate ? ' or degenerate' : '')
  if ($Short_Contigs && $Verbose_Mode);
# Output is in @ORDERING1.  Each line is tab-separated values: (Input
#   contig name, Reverse? (0=F, 1=R), Begin base, Output sequence length,
#   End base, Input sequence length, Output Contig Name).  All fields
#   are defined.  Only valid contig names remain.  Reverse? = 0/1.
#   Begin base is valid, or Begin base = -1 to mean entire input contig.
#   Other fields are forced to valid values.  Output contig name is set
#   to '-' if it is missing or if the '-j' flag is specified.  If Output
#   contig name is '-', then Output contig name is to be contructed from
#   Input contig name.

####################### Change order of contigs? #######################
if ($backwards_sort)
  {
  $backwards_phrase = 'reverse-';
  }
else
  {
  $backwards_phrase = '';
  }
if ($alpha_sort)
  {
  @ORDERING = sort @ORDERING1;
  print STDERR "Contigs names ${backwards_phrase}sorted alphanumerically\n" if ($Verbose_Mode);
  }
elsif ($numeric_sort)
  {
  @ORDERING = sort by_contignumber @ORDERING1;
  print STDERR "Contigs names ${backwards_phrase}sorted numerically\n" if ($Verbose_Mode);
  }
elsif ($order_file ne '')
  {
  @ORDERING = @ORDERING1;
  print STDERR "Contigs will be output in ${backwards_phrase}order of 'order_file'\n" if ($Verbose_Mode);
  }
elsif ($size_sort)
  {
  @ORDERING = sort by_size @ORDERING1;
  print STDERR "Contigs ${backwards_phrase}sorted by size\n" if ($Verbose_Mode);
  }
elsif ($sequence_sort)
  {
  @ORDERING = sort by_sequence @ORDERING1;
  print STDERR "Contigs ${backwards_phrase}sorted by sequence\n" if ($Verbose_Mode);
  }
else # default is to copy original input order
  {
  @ORDERING = @ORDERING1;
  print STDERR "Contigs will be output in ${backwards_phrase}order from 'fasta_input_file'\n" if ($Verbose_Mode);
  }
if ($backwards_sort)
  {
  @ORDERING = reverse @ORDERING;
  }
# Output is in @ORDERING.  Each line is tab-separated values: (Contig
#   name, Reverse? (0=F, 1=R), Begin base, Output sequence length,
#   End base, Input sequence length).  All fields are defined.  Only
#   valid contig names remain.  Reverse? = 0/1.  Begin base is valid,
#   or Begin base = -1 to mean entire input contig.  Other fields are
#   forced to valid values.  List is now sorted.
unless ($empty_ok || scalar @ORDERING)
  {
  print "ERROR:  No valid contigs to be output.  Program terminating.\n";
  exit 2;
  }

######################## Open output files #######################
my $append_flag = ($Extend_File) ? '>' : '';
my $output_msg = ($Extend_File) ? 'append to' : 'create';
if (!open(FASTAOUT, ">$append_flag$fasta_output_file"))
  {
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
    close(FASTAOUT);
    die("Can't $output_msg output fasta quality file: '$qualout', $!\n");
    exit 2;
    }
  }
if ($join_contigs)
  {
  if (!open(TOC, ">${fasta_output_file}.toc"))
    {
    my $err = $!;
    close(FASTAOUT);
    close(QUALOUT);
    die("Can't create output table-of-contents file: '${fasta_output_file}.toc', $!\n");
    exit 2;
    }
  print TOC "#Joined contigs for input file: '$fasta_input_file'\n";
  print TOC "#Contig_Name	Direction	Output_Offset	Output_Length	Begin_Base	End_Base	Input_Contig_Length\n";
  }

################## Output contigs in selected order ##################
%PROCESSED = ();
$GOOD_ORDER_CONTIGS = 0;
$FIRST_CONTIG = 1;
%qualities_modified = ();
for $contig_line (@ORDERING)
  {
  my($contig, $rev, $begin, $length, $end, $seq_len, $out_contig) =
    split("\t", $contig_line);
  my $sequence = $SEQ{$contig};
  if ($Preserve_Comments)
    {
    ($comment, $Qcomment) = @{ $COMMENTS{$contig} };
    }
  else
    {
    ($comment, $Qcomment) = ('', '');
    }

  if ($get_qualities)
    {
    unless (exists $qualities_modified{$contig})
      {
      ########## apply quality modifiers ##########

      # constant quality value?
      @{ $QUAL{$contig} } = ($QUALITY_VALUE) x length($sequence)
        if defined $QUALITY_VALUE;

      # subtract a constant value from all qualities?
      @{ $QUAL{$contig} } =
        map { my $t = ($_ - $QUALITY_SUB); ($t >= 0) ? $t : 0 }
          @{ $QUAL{$contig} } if defined $QUALITY_SUB;

      # divide all qualities by a constant value?
      @{ $QUAL{$contig} } = map { int($_ / $QUALITY_DIV) }
        @{ $QUAL{$contig} } if defined $QUALITY_DIV;

      # divide all qualities by stepped divisors?
      if (@QUALITY_STEPS)
	{
	my $step = -1;
	my $run_length = 0;
	my $divisor = 0;
	for (my $i = 0; $i < scalar @{ $QUAL{$contig} }; $i++)
	  {
	  while ($step < (scalar @QUALITY_STEPS) - 1 &&
		 $run_length-- <= 0)
	    {
	    $step++;
	    $run_length = $QUALITY_STEPS[$step];
	    $divisor++;
	    }
	  ${ $QUAL{$contig} }[$i] =
	    int(${ $QUAL{$contig} }[$i] / $divisor);
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
	for (my $i = 0; $i < scalar @{ $QUAL{$contig} }; $i++)
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
	  ${ $QUAL{$contig} }[$i] = int(${ $QUAL{$contig} }[$i] / $slope_divisor);
	  } # end for (my $i = ... )
	} # end if (@QUALITY_SLOPES)

      # apply maximum allowed quality value
      @{ $QUAL{$contig} } = map { ($_ > $MAX_QUAL) ? $MAX_QUAL : $_ }
        @{ $QUAL{$contig} };

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
			 @{ $QUAL{$contig} }[$begin - 1 .. $end - 1]);
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
			 $begin, $end, @{ $QUAL{$contig} });
      }
    } # end if ($begin > 0)
  $GOOD_ORDER_CONTIGS++;
  $PROCESSED{$contig}++;
  $FIRST_CONTIG = 0;
  } # end for $contig_line (@ORDERING)
close(FASTAOUT);
close(QUALOUT) if ($get_qualities);


############# Print summary statistics if requested #############
if ($Verbose_Mode)
  {
  my $UNPROCESSED_CONTIGS = 0;
  for $contig (keys %SEQ)
    {
    if (! exists $PROCESSED{$contig})
      {
      $UNPROCESSED_CONTIGS++;
      }
    } # end for $contig (keys %SEQ)
  print STDERR "$GOOD_ORDER_CONTIGS contigs output\n";
  print STDERR "$BAD_ORDER_CONTIGS invalid contig names from 'order_file' were skipped\n" if $BAD_ORDER_CONTIGS;
  print STDERR "$UNPROCESSED_CONTIGS input contigs were not used\n" if $UNPROCESSED_CONTIGS;
  } # end if ($Verbose_Mode)

exit 0;


###########################################################################
# save_contig($contig, $comment, $sequence, $Qcomment) - save just read
#   contig into hashes %SEQ and %QUAL.  Save contig name in list @ORDERING
#   unless -o was specified.
###########################################################################

sub save_contig
  {
  my($contig, $comment, $sequence, $Qcomment) = @_;

  # If -u was specified, then look for a uaccno to become the contig
  # name.
  $contig = $1 if ($Use_Uaccno && $comment =~ /\suaccno=(\S+)/);

  # return and do not save if contig length is < min_size
  if (length($sequence) < $min_size)
    {
    $Short_Contigs++;
    $Short_Contigs{$contig} = 1;
    return;
    }

  # return and do not save if -d flag and contig sequence does not
  # contain at least one each of A, C, G, and T
  if ($no_degenerate && degenerate($sequence))
    {
    $Degenerate_Contigs++;
    $Short_Contigs{$contig} = 1;
    return;
    }

  # if -o was not specified, then save list of contigs in input order
  push @ORDERING, "$contig\t$global_direction" if ($order_file eq '');

  # now save the contig info
  $SEQ{$contig} = $sequence;
  $QUAL{$contig} = $Qualref if ($get_qualities == 1);
  $Qcomment = $comment if ($get_qualities == 2);
  $COMMENTS{$contig} = [$comment, $Qcomment] if $Preserve_Comments;
  } # end save_contig


###########################################################################
# by_contignumber - sort contigs by numerically, not alphanumerically
###########################################################################

sub by_contignumber
  {
  my($contiga) = split("\t", $a);
  my($contigb) = split("\t", $b);
  my($numa, $numb);
  my($num) = 1;
  if ($contiga =~ /(\d+)\D*$/)	# get contig number from $a
    {
    $numa = $1;
    }
  else
    {
    $num = 0;
    }
  if ($contigb =~ /(\d+)\D*$/)	# get contig number from $b
    {
    $numb = $1;
    }
  else
    {
    $num = 0;
    }
			# if both are numeric, compare contig numbers
  return $numa <=> $numb if $num && $numa <=> $numb;
  return $contiga cmp $contigb;	# if same numbers, compare names
  } # end by_contignumber


###########################################################################
# by_size - sort contigs by size of input contig
###########################################################################

sub by_size
  {
  my($contiga, $reva, $begina, $lena, $enda, $seq_lena) = split("\t", $a);
  my($contigb, $revb, $beginb, $lenb, $endb, $seq_lenb) = split("\t", $b);
					# first compare contig sizes
  return $seq_lena <=> $seq_lenb if $seq_lena <=> $seq_lenb;
  return $contiga cmp $contigb;	# if same size, compare names
  } # end by_size


###########################################################################
# by_sequence - sort contigs by contig sequence
###########################################################################

sub by_sequence
  {
  my($contiga, ) = split("\t", $a);
  my($contigb, ) = split("\t", $b);
  my $seqa = $SEQ{$contiga};
  my $seqb = $SEQ{$contigb};
  # first compare contig sequences
  # if same size, compare names
  return $seqa cmp $seqb or $contiga cmp $contigb;
  } # end by_sequence


###########################################################################
# output_contig - output a contig sequence.
###########################################################################

sub output_contig
  {
  my($contig, $out_contig, $rev, $sequence, $begin, $end, $contig_len) = @_;
  my($len) = length($sequence);
  my($segment, $i, $direction, $comp, $suffix);
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
  if ($join_contigs)
    {
    $begin = 1 if ($begin <= 0);
    if ($FIRST_CONTIG)
      {
      print FASTAOUT ">$fasta_input_file\n";
      $CONTIG_POSITION = 1;
      }
    else
      {
      print FASTAOUT "$join_string\n";
      $CONTIG_POSITION += $join_len;
      }
    print TOC "$contig\t$direction\t$CONTIG_POSITION\t$len\t$begin\t$end\t$contig_len\n";
    $CONTIG_POSITION += $len;
    }
  else # (! $join_contigs)
    {
    print FASTAOUT ">$contig_prefix$out_contig\n";
    }
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
  if ($join_contigs)
    {
    if ($FIRST_CONTIG)
      {
      print QUALOUT ">$fasta_input_file\n";
      }
    else
      {
      printf QUALOUT "%s\n", join(' ', @JOIN_QUALS);
      }
    }
  else # (! $join_contigs)
    {
    print QUALOUT ">$contig_prefix$out_contig\n";
    }
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

USAGE: $my_name [-o order_file / -a / -n / -S / -z] [-b]
	    [-c contig_prefix] [-d] [-e] [-g] [-j join_char]
	    [-l output_line_length] [-m min_size] [-M max_qual] [-p]
	    [-q] [-Q qual_value] [-r] [-s] [-t type_to_remove] [-u]
	    [-v] [-x] fasta_input_file fasta_output_file
              or
       $my_name -h

OPTIONS:

  -a  Sort input contigs in alphanumeric order by contig name.

  -b  Output the list of contigs in backwards order.

  -c  Contig prefix to be added to output contig names.

  -d  Input contigs are assumed to be dna.  Filter out any degenerate
      contigs that do not contain at least one each of A, C, G, and T.

  -e  Emtpy output files are OK and do not result in an error.

  -g  Remove runs of Ns from ends of contigs.  Minimum contig length
      is enforced after trimming the ends.

  -j  Join contigs into a single contig.  In addition, a table of contents
      file ('fasta_output_file'.toc) is created.  May not be used with
      -c or -x.

  -l  Specify output line length.

  -m  Specify minimum contig length to be used.  May not be specified with
      -o.

  -M  Specify a maximum quality score for all bases in the output
      quality file.

  -n  Sort input contigs in numeric order by contig name.

  -o  Specifies the name of an 'order_file', which contains a list of contig
      names to be output.  Each line of 'order_file' is a tab separated list.
      The fields in the list are:  contig_name, direction, begin_base, and
      length.  May not be used with -a, -m, -n, -r, or -z.

  -p  Preserve contig header comments.  May not be used with -u.

  -q  Also process a fasta quality file ('fasta_input_file'.qual), as well
      as the sequence file and create an output quality file
      ('fasta_output_file'.qual) in addition to the output sequence file.

  -Q  Specify a constant quality value to be applied to all bases in
      the output quality file or a modifier to be applied to all
      qualities from the input Fasta quality file.  If 'qual_value' is
      a simple one- or two-digit positive integer, then that value is
      used for the quality scores and the input Fasta quality file is
      not needed.  If 'qual_value' is not just a simple one- or
      two-digit integer, then it specifies a modifier to be applied to
      the values from the input Fasta quality file.

  -r  Each output contig is to be reversed and complemented.  May not be
      used with -o.

  -s  The contig name may be shortened by removing any prefix before the
      word "Contig", i.e., "gono.fasta.screen.Contig26" becomes "Contig26"
      (or "Contig26r" and "Contig26f" if the two ends are output as separate
      files.)

  -S  Sort input contigs in sequence order.

  -t  Specify a filetype to be removed before adding ".qual" to create
      the output quality filename.

  -u  Use universal accession numbers (uaccno) as contig names for 454
      reads.  May not be used with -p.

  -v  Verbose mode - print out some statistics while running.

  -x  Create new or append to existing (extend) output files.  May not
      be used with -j.

  -z  Sort input contigs by contig size.

EOF

  exit 2;
  } # end display_help


###########################################################################
# display_help
###########################################################################

sub display_more_help
  {
  print STDOUT <<EOF;

$my_name - Sort contigs in an input fasta sequence file, possibly with
an accompanying fasta quality file.  Contigs are left unsorted, if a sort
order is not specified.

NOTE: If you do NOT need to change the order of the output contigs,
(using options -a, -b, -n, -o order_file, -S, or, -z) or join output
contigs as one contig (option -j join_char), then the program
select_contigs may be a better choice and run faster, especially for
large input files.

If none of -a, -b, -n, -o, or -z is specified, the output contig order
matches the input order.  If -a or -n is specified, then the contigs are
sorted by contig name either alphanumerically (-a) or numerically (-n).
If -z is specified, then the contigs are sorted by contig size.  If -o
is specified, then the contigs are output in the order specified by the
lines in 'order_file', and only the contigs named in 'order_file' are
output.  If -o is not specified, then all input contigs are output.  If
-b is specified, then the order of contigs to be output is backwards from
what it would otherwise be;  this applies whether -a, -n, -o, -z, or none
of these is specified.  Contigs shorter than 'min_size' bases long are
discarded.

If -j is specified, then all of the contigs are joined with sequences of
'join_char's and output as a single contig, and a table-of-contents file
named 'fasta_output_file'.toc is also created.  If -j is not requested,
the contigs are output as separate contigs.

The contig names optionally may be shortened by removing everything before
the word "Contig" (-s).

Contig header comments, located after the contig name on the contig header
line, will be removed, but new comments may be added if partial contigs
are requested in 'order_file'.  If '-p' is specified, then the original
contig header comments will be added back to the contig header (possibly
after any added comments describing partial contigs).  All output lines
containing sequence bases (not contig headers) are reformatted to be
'output_line_length' bases long, with the exception that contig separator
sequences (specified by -j) are output on separate output lines.

If '-u' is specified, then universal accession numbers (uaccno) are
used as contig names for 454 reads instead of the rank_x_y form
normally present in 454Reads files.  The uaccnos will be used in the
output contigs and in any order files used to select contigs.  The '-p'
flag may not be used with '-u'.

If '-q' is specified, then a fasta quality file named
"'fasta_input_file'.qual" also is processed, creating an additional
output fasta quality file named "'fasta_output_file'.qual".  If
'-t type_to_remove' is also specified, then the 'type_to_remove' is
removed from the end of 'fasta_output_file' before adding ".qual".

All output contigs may be reversed and complemented if requested (-r), or
on an individual basis as specified on each line in the 'order_file' (-o).
Both -r and -o may not be specified together.


USAGE: $my_name [-o order_file / -a / -n / -S / -z] [-b]
	    [-c contig_prefix] [-d] [-e] [-g] [-j join_char]
	    [-l output_line_length] [-m min_size] [-M max_qual] [-p]
	    [-q] [-Q qual_value] [-r] [-s] [-t type_to_remove] [-u]
	    [-v] [-x] fasta_input_file fasta_output_file
              or
       $my_name -h           <== What you are reading

  where 'order_file' is the name of an input file containing the output
            ordering for the contigs.

        'join_char' - The output contigs are joined together and output as
            a single contig, named the same as the input sequence file,
            'fasta_input_file'.  If a single character is specified for
            'join_char', then each pair of input contigs will be separated
            by $join_len 'join_char's.  If a multi-character string is
            specified, then that string will be used as the separator
            without duplication.  If '-j' is not specified, the contigs
            are output as separate contigs.

        'output_line_length' is the length of the output lines containing
            sequence data (not sequence headers).  The default value is $DEFAULT_LINE_LENGTH.
            'output_line_length' must be greater or equal to $MIN_LINE_LENGTH.

        'max_qual' is the maximum allowed quality value for all bases
            in the output file.  Any bases with a quality value higher
            than 'maximum_quality_value' will have the quality lowered to
            'maximum_quality_value'.  'max_qual' must be a one- or
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

        'fasta_input_file' is the name of the input sequence file in Fasta
            format, which contains the contigs to be processed.  If
            'fasta_input_file' is specified as '-', then standard input
            is used.  If '-q' is specified, then 'fasta_input_file' may
            not be specified as '-'.

        'fasta_output_file' is the name of the output sequence file in
            Fasta format.  If 'fasta_output_file' is specified as '-',
            then standard output is used.  If '-j' or '-q' is specified,
            then 'fasta_output_file' may not be specified as '-'.

OPTIONS:

  -a  Sort input contigs in alphanumeric order by contig name.  May be
      modified by -b.  May not be used with -n, -o, -S, or -z.

  -b  Output the list of contigs in backwards order.

  -c  Contig prefix to be added to output contig names.  The prefix is
      added after possible shortening by -s.  May not be used with -j
      or -o.

  -d  Input contigs are assumed to be dna.  Filter out any degenerate
      contigs that do not contain at least one each of A, C, G, and T.

  -e  Emtpy output files are OK and do not result in an error.  If -e
      is not specified, and the input file is empty or an order file
      selects no contigs, then the program normally exits with an
      error.

  -j  Join contigs into a single contig, named the same as the
      'fasta_input_file'.  In addition, a table of contents file,
      'fasta_output_file'.toc, is created.  If both -j and -o are specified,
      then the output_contig_name fields in the order file are not used.
      May not be used with -c or -x.

  -l  Specify output line length.

  -m  Specify minimum contig length to be used.

  -M  Specify a maximum quality score for all bases in the output
      quality file.

  -n  Sort input contigs in numeric order by contig name.  May be modified
      by -b.  The program searches each contig name and looks for the last
      string of numeric characters within the contig name.  May not be
      used with -a, -o, -S, or -z.

  -o  Specifies the name of an 'order_file', which contains a list of
      contig names to be output.  Each line of 'order_file' is a
      white-space separated list of one or more fields.  The fields in
      the list are:  input_contig_name, direction, begin_base, length,
      and output_contig_name.  The fields present may be followed by a
      comment, which begins with a '#' character after a white space
      character and runs to the end of the line.  May be modified by
      -b.  May not be used with -a, -c, -n, -r, -S, or -z.  If both -j
      and -o are specified, then the output_contig_name fields in the
      order file are not used.

  -p  Preserve contig header comments.  May not be used with -u.

  -q  Also process a fasta quality file, "'fasta_input_file'.qual", as well
      as the sequence file and create an output quality file,
      "'fasta_output_file'.qual", in addition to the output sequence file.
      $my_name also allows relaxed quality file naming.  If
      "'fasta_input_file'.qual" is not found, and 'fasta_input_file' is of
      the form "xxx.fa", "xxx.fna", or "xxx.fasta", then $my_name
      will try to use "xxx.qual" instead.

  -Q  Specify a constant quality value to be applied to all bases in
      the output quality file or a modifier to be applied to all
      qualities from the input Fasta quality file.  If 'qual_value' is
      a simple one- or two-digit positive integer, then that value is
      used for the quality scores and the input Fasta quality file is
      not needed.  If 'qual_value' is not just a simple one- or
      two-digit integer, then it specifies a modifier to be applied to
      the values from the input Fasta quality file.  See the
      dexription of 'qual_value' above.

  -r  Each output contig is to be reversed and complemented.  May not be
      used with -o.

  -s  The contig name may be shortened by removing any prefix before the
      word "Contig", i.e., "gono.fasta.screen.Contig26" becomes "Contig26".

  -S  Sort input contigs in sequence order.  May be modified by -b.
      May not be used with -a, -n, -o, or -z.

  -t  Specify a filetype to be removed before adding ".qual" to create
      the output quality filename.  See 'type_to_remove' above.  The
      '-t' flag is not used unless '-q' is also specified.

  -u  Use universal accession numbers (uaccno) as contig names for 454
      reads instead of the rank_x_y form normally present in 454Reads
      files.  May not be used with -p.

  -v  Verbose mode - print out some statistics while running.

  -x  Create new or append to existing (extend) output files.  May not
      be used with -j.

  -z  Sort input contigs by contig size.  May be modified by -b.  May not
      be used with -a, -n, -o, or -S.


ORDER_FILE INPUT FORMAT:

Blank lines and comment lines (beginning with #) are allowed.  Each
non-blank, non-comment, line of 'order_file' should name one input
contig to be written to the output file.  The same input contig name
may be used more than once.  This is useful for splitting an input
contig into multiple output contigs.  Each input line may have multiple
white-space separated fields.  Leading white-space is ignored.
Trailing fields may be omitted to use the default values.  After any
fields that are present, a trailing comment may be added.

  - Input contig name.  This field is required.
  - Direction (forward or reverse & complement).  If the field begins
    with a C or an R (of either case), then that contig is reversed and
    complemented before it is output, and the Output contig name
    contains the Input contig name suffixed with '.comp' (or '.comp' is
    removed), unless Output contig name is specified explicitly.  The
    default is to output the contig in the forward direction.
  - Beginning base of input contig.  The default beginning base is 1.  A
    negative value indicates to start that many bases before the end of
    the sequence.  (-1 means start at the last base;  -10 means start at
    the 10th base counting from the end of the sequence, etc.)  A value
    of zero or just a minus sign '-' is treated the same as base 1.  A
    value longer than the input sequence is treated as a short contig
    and ignored.  The value for Beginning base refers to base positions
    in the original input contig before any base reversal and complement
    operation that may occur.
  - Length of sequence to be output.  A negative value indicates to end
    that many bases before the end of the sequence.  (-1 means end 1
    base before the last base;  -10 means end 10 bases before the last
    base of the sequence, etc.)  If the length is omitted or zero or
    out of range or less than Beginning base or just a minus sign '-',
    then the end of the contig is used.
  - Output contig name.  This field may be used to name the resulting
    output contig.  If -j is specified, then this field is not used.
    If this field is omitted, then an output contig name will be
    constructed from the Input contig name.
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
Contig1	C			# All of Contig1 reversed and complemented
Contig3	R	1	100	# First 100 bases of Contig3, then reversed
				# and complemented
Contig2	F	-100		# Last 100 bases of Contig2
Contig7	F	51	-50	# All but the first and last 50 bases on
				# the ends of Contig7
Contig5	F	-100	-50	# The 50 bases ending 50 bases before the
				# end of Contig5
Contig2	f 501 1000 Contig2:501-1000	# Second 500 bases of Contig2.  The
				# output contig is named "Contig2:501-1000"
Contig4	f  123 334  Fred	# Bases 123-456 of Contig4.  The output
				# contig is named "Fred"
Contig6	r - - Contig6.revcomp	# All of Contig6 reversed and complemented
				# and named "Contig6.revcomp"


TABLE-OF-CONTENTS FILE OUTPUT FORMAT

The Table-of-Contents file is produced only when -j is specified and is
named "'fasta_output_file'.toc".  Two header lines (comments beginning
with #) begin the file.  Each subsequent line of the table of contents
file contains the following tab-separated fields:

  - Contig name,
  - Direction (F or R),
  - Beginning base offset in the output contig (starting at 1),
  - Output sequence length,
  - Beginning base from input contig,
  - Ending base from input contig, and
  - Total input contig length.

If the order file or the -r flag specifies to reverse and complement an
input contig, then Direction is 'R' and the Contig name has '.comp'
added (or removed) as a suffix.  "Beginning base from input contig" and
"Ending base from input contig" refer to base positions before the
reverse and complement, if Direction is 'R'.


EXAMPLES:

For the following examples, assume 'input.fa' and 'input.fa.qual'
contain the following contigs:

Contig1 (100 bases), Contig2 (200 bases), ... Contig101 (10100 bases)

\$ cat input.fa | $my_name -n - - > numeric.fa
\$ $my_name -n alpha.fa numeric.fa

Either of the above lines reads the input fasta sequence file 'input.fa',
sorts the contigs by contig number, and writes a new fasta sequence file,
'numeric.fa'.  The output contig order is: Contig1, Contig2, ... Contig101.

\$ $my_name -z -b input.fa big_first

sorts the contigs by decreasing contig size and writes a new fasta
sequence file, 'big_first.fa'.  The output contig order is: Contig101,
Contig100, ... Contig1.

\$ $my_name -r input.fa revcomp.fa

reads the input fasta sequence file 'input.fa', reverses and complements
all input contigs, and writes a new fasta sequence file, 'revcomp.fa'.
The output contig order is unchanged from the input.

\$ $my_name -q -a input.fa alpha.fa

reads the input fasta sequence file 'input.fa' and the input fasta
quality file 'input.fa.qual', sorts the contigs in alphanumeric order,
and writes a new fasta sequence file, 'alpha.fa' and a new fasta
quality file, 'alpha.fa.qual'.  The output contig order is: Contig1,
Contig10, Contig100, Contig101, Contig11, Contig12, ... Contig19,
Contig2, Contig20, Contig21, ... Contig99.

\$ $my_name -v -o order input.fa ordered.fa

reads the input fasta sequence file 'input.fa' and writes a new fasta
sequence file, 'ordered.fa', containing the contigs or parts of contigs,
as specified in the order_file 'order'.  The output contig order is as
specified in the file 'order', and only the requested contigs are written.
The Verbose flag (-v) causes some statistics to be written to Standard
Output, including counts of contigs read, contigs written, and contigs
not written to the output file, if any.

\$ $my_name -j X -n -q input.fa one_big_one.fa

reads the input fasta sequence file 'input.fa' and the input fasta
quality file 'input.fa.qual', sorts the contigs in numeric order,
and writes all of the contigs as one big contig in a new fasta sequence
file, 'one_big_one.fa' and a new fasta quality file, 'one_big_one.fa.qual'.
The output contig order is: Contig1, Contig2, ... Contig101, but all of
these contigs are contained in one long contig named the same as the input
file, 'input.fa'.  Separating consecutive input contigs is a line of $join_len
Xs in the sequence file ($join_len zero qualities in the quality file).
A table of contents file, 'one_big_one.fa.toc', is also created listing
all output contigs and their sizes and positions.

\$ $my_name -r input.fa long_ones.fa -m 2000

reads the input fasta sequence file 'input.fa' and outputs a new fasta
sequence file, 'long_ones.fa', containing only the contigs at least 2000
bases long.  The output contig order is unchanged from the input.

$my_name -Q 20 contigs.fa contigs20.fa

reads the fasta sequence file "contigs.fa" (and does NOT use a fasta
quality file "contigs.fa.qual" or "contigs.qual") and produces new
fasta sequence and quality files named "contigs20.fa" and
"contigs20.fa.qual".  All base quality values are set to 20.

$my_name -Q s20 -t .fa contigs.fa contigs_sub20.fa

reads the fasta sequence file "contigs.fa" and a fasta quality file
"contigs.fa.qual" (or "contigs.qual") and produces new fasta sequence
and quality files named "contigs_sub20.fa" and "contigs_sub20.qual".
All base quality values are reduced by 20, but any qualities which
would be negative are set to zero.  Note that the output quality file
is "*.qual", not "*.fa.qual", because of the "-t .fa" option.

$my_name -Q d2.5 contigs.fa contigs_div2.5.fna

reads the fasta sequence file "contigs.fa" and the fasta quality file
"contigs.fa.qual" (or "contigs.qual") and produces new fasta sequence
and quality files named "contigs_div2.5.fna" and
"contigs_div2.5.fna.qual".  All base quality values are divided by 2.5
and rounded down to the next integer.

$my_name -Q step100,50,25 contigs12.fna contigs_3step.fasta

reads the fasta sequence file "contigs12.fna" and the fasta quality
file "contigs12.fna.qual" (or "contigs12.qual") and produces new fasta
sequence and quality files named "ccontigs_3step.fasta" and
"contigs_3step.fasta.qual".  The first 100 base quality values are
divided by 1;  the next 50 are divided by 2;  the next 25 are divided
by 3;  the remaining base qualities are divided by 4. The resulting
quality scores are all rounded down to the next lowest integer.

$my_name -Q slope10,5,0,4,2 contigs13.fna contigs_5slope.fasta

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
The resulting quality scores are all rounded down to the next lowest integer.


DATE LAST MODIFIED: $Date_Last_Modified

EOF

  exit;

  } # end display_more_help

