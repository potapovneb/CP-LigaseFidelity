#!/usr/bin/perl -w

use strict;
use Text::ParseWords;

### command-line usage
if( @ARGV == 0 )
{
    print "usage: $0 barcodes.csv\n\n";
    exit;
}

### command-line arguments
my $file = shift @ARGV;

### read all barcodes and their counts
my ($data, $head) = read_csv($file);

### total number of bases and per position
my %bctot = ();
my %bcpos = ();

### barcode length
my $bclen = length($$data[0]{"Barcode"});

for my $entry ( @$data )
{
    my $barcode = $$entry{"Barcode"};
    my $count   = $$entry{"Count"};

    for( my $k = 0; $k < $bclen; $k++ )
    {
	my $base = substr($barcode, $k, 1);

	### ignore deletions (rare)
	next if( $base eq "-" );

	### number of bases at each position
	$bcpos{$base}{$k+1} += $count;
	$bcpos{"*"}{$k+1}   += $count;

	### total number of bases
	$bctot{$base} += $count;
	$bctot{"*"}   += $count;
    }
}

### print table headers
print join( ",",
	    "Base",
	    (map { "c$_" } 1..$bclen, "ANY"),
	    (map { "f$_" } 1..$bclen, "ANY") ), "\n";

### print table values
for my $base ( "A", "C", "G", "T", "*" )
{
    ### init
    my @data = ($base);
    
    ### print base counts (per position)
    for my $pos ( sort { $a <=> $b } keys %{$bcpos{$base}} )
    {
	push @data, $bcpos{$base}{$pos};
    }

    ### print base counts (total)
    push @data, $bctot{$base};
    
    ### print base fractions (per position)
    for my $pos ( sort { $a <=> $b } keys %{$bcpos{$base}} )
    {
	push @data, $bcpos{$base}{$pos} / $bcpos{"*"}{$pos};
    }
    
    ### print base fractions (total)
    push @data, $bctot{$base} / $bctot{"*"};
    
    ### actually print data
    print join(",",@data), "\n";
}

################################################################################
#                             SUBROUTINES                                      #
################################################################################

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
