use strict ; 
use warnings ; 

my %transcript ; 
my $min_daf = 0.1 ; 
my $max_daf = 0.9 ;

while (<STDIN>) { 	
	
	chomp $_ ; 
	my @split = split ( /\t/, $_ ) ; 
	
	my @transcripts = split ( /\,/, $split[$#split] ) ; 
	pop @transcripts ; 

	foreach my $trs ( @transcripts ) { 
		if ( !exists( $transcript{$trs} ) ) { 
			push @{ $transcript{$trs} }, (0,0,0,0) ; 
		}
		foreach ( 2..3 ) { 
			${ $transcript{$trs} }[$_-2] += $split[$_] ; 
		}
		foreach ( 4..5 ) { 
			if ( $split[7] ne "NA" && $split[7] >= $min_daf && $split[7] <= $max_daf ) { 
				${ $transcript{$trs} }[$_-2] += $split[$_] ;
			}
		}
	}
}

foreach my $trs ( sort keys %transcript ) { 
	print $trs ; 
	foreach ( 0..3 ) { 
		print "\t", ${ $transcript{$trs} }[$_] ;
	}
	print "\n" ;
}
