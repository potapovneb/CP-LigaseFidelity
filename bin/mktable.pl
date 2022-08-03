#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Text::ParseWords;

my @bases = qw/A C G T/;

### command-line options
my $o_size   = 4;
my $o_prefix = "table";

GetOptions(
    "size=f"   => \$o_size,
    "prefix=s" => \$o_prefix,
    );

### command-line usage
if( @ARGV == 0 )
{
    print "usage: $0 overhangs.csv\n";
    print "\n";
    print "options\n";
    print "  --size\toverhang size (default '$o_size')\n";
    print "  --prefix\tprefix for output files (default '$o_prefix')\n";
    print "\n";
    exit;
}

### command-line arguments
my @files = @ARGV;

my $data = load_overhangs(\@files);

### generate all k-mers: AAA, AAT, AAG, AAC, ATA, ATT, ATG, ...
my @array = ();
my @LIST = ();
deep(\@array,0,["A","C","G","T"],$o_size,\@LIST);

if ($o_size == 3 || $o_size == 4)
{
    open(T,">","${o_prefix}-mismatch-by-position.csv");
    mismatch_by_position(\*T,$data);
    close(T);
}

open(T,">","${o_prefix}-ligation_fidelity.csv");
ligation_fidelity(\*T,$data);
close(T);

open(T,">","${o_prefix}-overhang_matrix.csv");
overhang_matrix(\*T,$data);
close(T);

################################################################################
#                             SUBROUTINES                                      #
################################################################################

sub load_overhangs {
    my ($files) = @_;
    
    my %data = ();
    my $total = 0;

    for my $file ( @$files )
    {
    my ($loc_data,$loc_head) = read_csv($file);

    for my $entry ( @$loc_data )
    {
        my $o1 = $$entry{"O1"};
        my $o2 = $$entry{"O2"};
        my $n  = $$entry{"Count"};

        ($o1,$o2) = sort($o1,$o2);

        ### double counts (=> symmetric)
        $data{$o1}{$o2} += $n;
        $data{$o2}{$o1} += $n;
        $total += 2 * $n;
    }
    }

    return \%data;
}

sub deep {
    my ($array,$level,$dict,$depth,$result) = @_;
    
    for( my $i = 0; $i < @$dict; $i++ )
    {
    $$array[$level] = $$dict[$i];
    
    if( $level < $depth-1 )
    {
        deep($array,$level+1,$dict,$depth,$result);
    }
    else
    {
        push @$result, join("",@$array);
    }
    }
}

sub mismatch_by_position {
    my ($fh,$data) = @_;

    my %out = ();
    my %in = ();

    for my $b1 ( @bases )
    {
    for my $b2 ( @bases )
    {
        $out{$b1}{$b2} = 0;
        $in{$b1}{$b2} = 0;
    }
    }

    for my $o1 ( keys %$data )
    {
    for my $o2 ( keys %{$$data{$o1}} )
    {
        my @top = split(//,$o1);
        my @bot = split(//,scalar reverse($o2));

        if( $o_size == 3 )
        {
        ####   0  1  2
        ### p-N1 N2 N3   : top
        ###   N3 N2 N1-p : bot
        
        ### N3 / N1
        $out{ $bot[0] }{ $top[0] } += $$data{$o1}{$o2};

        ### N2 / N2
        $in{  $top[1] }{ $bot[1] } += $$data{$o1}{$o2};
        }
        else
        {
        ####   0  1  2  3
        ### p-N1 N2 N3 N4   : top
        ###   N4 N3 N2 N1-p : bot
        
        ### N4 / N1
        $out{ $bot[0] }{ $top[0] } += $$data{$o1}{$o2};

        ### N3 / N2
        $in{  $bot[1] }{ $top[1] } += $$data{$o1}{$o2};
        }
    }
    }

    print_mismatch_by_position($fh,\%out,"edge");
    print_mismatch_by_position($fh,\%in,"middle");
}

sub print_mismatch_by_position {
    my ($fh,$mat,$text) = @_;

    my $total = 0;

    for my $b1 ( @bases )
    {
    for my $b2 ( @bases )
    {
        $total += $$mat{$b1}{$b2};
    }
    }

    print $fh ",\n";
    print $fh "Count ($text)\n";
    print $fh join(",","",@bases), "\n";

    for my $b1 ( @bases )
    {
    my @values = ( $b1 );
    
    for my $b2 ( @bases )
    {
        push @values, $$mat{$b1}{$b2};
    }

    print $fh join(",",@values), "\n";
    }

    print $fh ",\n";
    print $fh "Fraction ($text)\n";
    print $fh join(",","",@bases), "\n";

    for my $b1 ( @bases )
    {
    my @values = ( $b1 );
    
    for my $b2 ( @bases )
    {
        push @values, $$mat{$b1}{$b2} / $total;
    }

    print $fh join(",",@values), "\n";
    }
}

sub ligation_fidelity {
    my ($fh,$data) = @_;

    my %mm = ();
    my %stat = ();

    for my $o ( @LIST )
    {
    $stat{$o} = {
        "Total" => 0,
        "Correct" => 0,
        "Mismatch" => 0
    };

    for( my $i = 1; $i <= 4; $i++ )
    {
        $mm{$o}{$i}{"Total"} = 0;
    }
    }

    my $total = 0;  # total number of overhangs

    ### check all overhang pairs
    for my $o1 ( keys %$data )
    {
    for my $o2 ( keys %{$$data{$o1}} )
    {
        if( $o1 eq scalar reverse(complement($o2)) )
        {
        ### WC ligations
        $stat{$o1}{"Correct"} += $$data{$o1}{$o2};
        }
        else
        {
        ### mismtach ligations
        $stat{$o1}{"Mismatch"} += $$data{$o1}{$o2};

        ### number of mismatches
        my ($mm,$seq) = number_of_mismatches_2($o1,$o2);

        ### record number of overhangs with 1, 2, 3, 4 mismatches
        $mm{$o1}{$mm}{"Total"} += $$data{$o1}{$o2};

        ### keep record of mismatch overhangs and their counts
        $mm{$o1}{$mm}{"Data"}{$seq} += $$data{$o1}{$o2};
        }

        ### total number of ligation for a given overhang
        $stat{$o1}{"Total"} += $$data{$o1}{$o2};

        $total += $$data{$o1}{$o2};
    }
    }

    print $fh join(",", qw/Overhang GC Total Correct Mismatch Fidelity fTotal fCorrect fMismatch fFidelity MM1 MM2 MM3 MM4 MM1_Overhangs MM2_Overhangs MM3_Overhangs MM4_Overhangs MO1 MO2 MO3 MO4 MO5 MO6 MO7 MO8 MO9 M10 Mismatch_Overhangs/ ), "\n";

    for my $o ( @LIST )
    {
    ### --------------------------------------------------------------------
    ### compute fractions
    $stat{$o}{"fTotal"}    = $stat{$o}{"Total"}    / $total;
    $stat{$o}{"fCorrect"}  = $stat{$o}{"Correct"}  / $total;
    $stat{$o}{"fMismatch"} = $stat{$o}{"Mismatch"} / $total;

    ### --------------------------------------------------------------------
    ### compute fidelity
    if( $stat{$o}{"Total"} > 0 )
    {
        $stat{$o}{"Fidelity"}  = $stat{$o}{"Correct"}  / $stat{$o}{"Total"};
        $stat{$o}{"fFidelity"} = $stat{$o}{"fCorrect"} / $stat{$o}{"fTotal"};
    }
    else
    {
        $stat{$o}{"Fidelity"}  = 0;
        $stat{$o}{"fFidelity"} = 0;
    }

    ### --------------------------------------------------------------------
    ### mismatch overhangs
    my %sets = ();

    for( my $i = 1; $i <= 4; $i++ )
    {
        ### get list of all mismatch overhangs
        my @set = sort { $mm{$o}{$i}{"Data"}{$b} <=> $mm{$o}{$i}{"Data"}{$a}} keys %{$mm{$o}{$i}{"Data"}};

        ### keep at most top three
        if( scalar @set > 3 )
        {
        @set = @set[0..2];
        }

        ### format output string
        my $str = join(";", map { sprintf( "%s=%i", $_, $mm{$o}{$i}{"Data"}{$_} ) } @set );

        ### store output string
        $sets{$i} = $str;
    }

    ### --------------------------------------------------------------------
    ## sort overhangs by ligation count
    my @sorted = sort { $$data{$o}{$b} <=> $$data{$o}{$a} } keys %{$$data{$o}};

    ## ligation counts (top 10 mismatch overhangs)
    my @counts = (0) x 10;
    
    my $i = 0;
    my $sum = 0;

    for my $o2 ( @sorted )
    {
        ## non WC overhnags only
        if( $o ne scalar reverse(complement($o2)) )
        {
        $counts[$i] = $$data{$o}{$o2};
        $sum += $$data{$o}{$o2};
        $i++;

        last if( $i >= 9 );
        }
    }

    $counts[9] = $stat{$o}{"Mismatch"} - $sum;

    ### --------------------------------------------------------------------
    if( @sorted > 3 )
    {
        @sorted = @sorted[0..3];
    }

    my $mismatch_overhangs = join("; ", map { ($o ne scalar reverse(complement($_)) ? sprintf("%s (%.0f%%)", $_, $$data{$o}{$_} / $stat{$o}{"Mismatch"} * 100) : ()) }  @sorted);
    
    ### --------------------------------------------------------------------
    ### ptint out values
    print $fh join(",",
               $o,
               gc_content($o) * 100,
               $stat{$o}{"Total"},
               $stat{$o}{"Correct"},
               $stat{$o}{"Mismatch"},
               $stat{$o}{"Fidelity"},
               $stat{$o}{"fTotal"},
               $stat{$o}{"fCorrect"},
               $stat{$o}{"fMismatch"},
               $stat{$o}{"fFidelity"},
               $mm{$o}{"1"}{"Total"},
               $mm{$o}{"2"}{"Total"},
               $mm{$o}{"3"}{"Total"},
               $mm{$o}{"4"}{"Total"},
               $sets{1},
               $sets{2},
               $sets{3},
               $sets{4},
               @counts,
               $mismatch_overhangs
        ), "\n";
    }
}

sub gc_content {
    my ($seq) = @_;

    my $num = 0;

    for( my $i = 0; $i < length($seq); $i++ )
    {
    my $b = substr($seq,$i,1);

    if( uc($b) eq "G" || uc($b) eq "C" )
    {
        $num++;
    }
    }

    return $num/length($seq);
}

sub number_of_mismatches {
    my ($o1,$o2) = @_;

    my $mm = 0;
    my $str = "";

    my $o2_rc = scalar reverse(complement($o2));
    my $o2_r  = scalar reverse($o2);

    for( my $i = 0; $i < length($o1); $i++ )
    {
    if( substr($o1,$i,1) ne substr($o2_rc,$i,1) )
    {
        $mm++;
        $str .= lc(substr($o2_r,$i,1));
    }
    else
    {
        $str .= uc(substr($o2_r,$i,1));
    }
    }

    return ($mm,$str);
}

sub number_of_mismatches_2 {
    my ($o1,$o2) = @_;

    my $mm = 0;
    my $str = "";

    my $o2_rc = scalar reverse(complement($o2));
    my $o2_r  = $o2;

    for( my $i = 0; $i < length($o1); $i++ )
    {
    if( substr($o1,$i,1) ne substr($o2_rc,$i,1) )
    {
        $mm++;
        $str .= substr($o2_r,$i,1);
    }
    else
    {
        $str .= substr($o2_r,$i,1);
    }
    }

    return ($mm,$str);
}

sub overhang_matrix {
    my ($fh,$data) = @_;

    my @LIST2 = ();

    for my $o ( @LIST )
    {
    push @LIST2, scalar reverse(complement($o));
    }

    print $fh join(",","Overhang",@LIST), "\n";

    for my $o1 ( @LIST2 )
    {
    my @values = map { exists $$data{$o1}{$_} ? $$data{$o1}{$_} : 0 } @LIST;
        
    print $fh join(",", $o1, @values), "\n";
    }
}

sub complement {
    my ($seq) = @_;

    my %bases = (
    "A" => "T",
    "C" => "G",
    "G" => "C",
    "T" => "A",
    "a" => "t",
    "c" => "g",
    "g" => "c",
    "t" => "a",
    "-" => "-",
    );

    my $complement = join("",map { $bases{$_} } split(//,$seq));

    return $complement;
}

sub read_csv {
    my ($file) = @_;

    my @head = ();
    my @array = ();

    open(CSV,$file) || die "Can't open '$file'";

    while( my $line = <CSV> )
    {
    chomp($line);

    my @tokens = quotewords(',',0,$line);
    
    if( @head == 0 )
    {
        @head = @tokens;
    }
    else
    {
        my %entry = ();

        for( my $i = 0; $i < @tokens; $i++ )
        {
        $entry{$head[$i]} = $tokens[$i];
        }

        push @array, \%entry;
    }
    }

    close(CSV);

    return (\@array,\@head);
}
