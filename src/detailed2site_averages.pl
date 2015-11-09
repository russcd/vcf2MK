use strict ; 
use warnings ; 

my %sites ; 
my %freqs ; 

while (<STDIN>) { 	
	
	chomp $_ ; 
	my @split = split ( /\t/, $_ ) ; 

	if ( !exists( $sites{$split[0]}{$split[1]} ) ) { 
		push @{ $sites{$split[0]}{$split[1]} }, (0,0,0,0,0) ; 
	}

	foreach ( 2..5 ) { 
		${ $sites{$split[0]}{$split[1]} }[$_-2] += $split[$_] ;
	}	
	${ $sites{$split[0]}{$split[1]} }[4] ++ ;
	$freqs{$split[0]}{$split[1]} = $split[7] ; 
}

foreach my $contig ( sort keys %sites ) {
	foreach my $position ( sort { $a<=>$b } keys %{ $sites{$contig} } ) {  
		print $contig, "\t", $position ; 	
		foreach ( 0..3 ) { 
			print "\t", ${ $sites{$contig}{$position} }[$_]/${ $sites{$contig}{$position} }[4] ;
		}
		print "\t", ${ $sites{$contig}{$position} }[4], "\t", $freqs{$contig}{$position}, "\n" ;
	}
}
