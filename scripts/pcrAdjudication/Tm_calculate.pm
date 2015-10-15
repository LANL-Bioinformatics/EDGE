# Po-E (Paul) Li
# B-11 LANL
# 09/10/2008
#
# Tm_calculate.pm - A module is used to calculate Tm
#
# 20131006
#   - supported Tm calculation of internal loop and bugles
# 20130930
#   - bugfix: Fixed bug of Tm calculation using Santalucia's method
# 2012
#   - support single mismatc and dangling ends
#
# USAGE:
# 
# Tm(<OLIGO> [,<TARGET>]) is the main function to claculate Tm. Use a dash and a base (5' -N and/or -N 3')
# at the flanking sequence to represent the dangling bases which N can be A,T,C or G. Puting "-" inside of
# oligo/target also means a bugle at its opposite strand in a duplex.
# 
# Exampels:
# 1) Exact match of DNA oligo complex
#
#     oligo: TATTCCTCGCCTC
#            |||||||||||||
#    target: TATTCCTCGCCTC
#
#    Tm( "TATTCCTCGCCTC" );
#
# 2) DNA oligo complex with mismatches
#
#     oligo: 5' TATTCCTCGCCTC
#               ||| |||| ||||
#    target:    TATTCCTCTCCTC 3'
#
#    Tm( "TATTCCTCGCCTC", "TATTCCTCTCCTC" );
#
# 3) DNA oligo binds to a dangling target:
#
#     oligo: 5'          TATTCCTCGCCTC
#                        |||||||||||||
#    target:    ATCGGTAGCTATTCCTCGCCTCCGTAGCT 3'
#
#    Tm( "TATTCCTCGCCTC", "-CTATTCCTCGCCTCC-");
#
# 4) DNA oligo partially bind on target
#
#            5'      GTT
#     oligo:            TATTCCTCGCCTC
#                       |||||||||||||
#    target:   ATCGGTAGCTATTCCTCGCCTCCGTAGCT 3'
#
#    Tm( "-TTATTCCTCGCCTC", "-CTATTCCTCGCCTCC-");
#
# 5) DNA oligo binds to a dangling target:
#
#     oligo: 5'          TATTCCTGTAGCT
#                        |||||||||||||
#    target:    ATCGGTAGCTATTCCTGTAGCTTGAGAGGA 3'
#                             GA
#                            T  G
#                             AT
#
#    Tm( "TATTCC------TGTAGCT", "-CTATTCCGTATGATGTAGCTT-");
#

package Tm_calculate;

use strict;
use Exporter;

our ( @ISA, @EXPORT, $VERSION );
$VERSION = 0.9;
@ISA = qw( Exporter );
@EXPORT  = qw( Tm_parameters expandDegen gc com_rev equivalent_na_conc Tm );

# Parameters
my $T = 310.15; # =273.15+37
#my $T = 298.15; # =273.15+25
my $r = 1.987;  # molar gas constant;

my %param_H;
my %param_S;
my %param_loop_terminal_H;
my %param_loop_terminal_S;
my %param_loop_S;
my %param_bulge_S;

# default settings
my ( $oligo_conc, $salt_conc, $mg_conc, $dntp_conc, $bound, $dplex_path, $salt_correct_formula, $debug ) = ( 5e-9, 0.08, 0.004, 0.0002, 0.5, 'dplex', 1, 0 );

sub Tm_parameters {
	my ( $class, %args ) = @_;
	
	my ( $param_H, $param_S) = &buildThermoParameters();
	%param_H = %$param_H;
	%param_S = %$param_S;
	
	foreach my $argument ( keys %args ) {
		if ( $argument =~ /^-oligo/i ) {
			$oligo_conc = $args{$argument};
		}
		elsif ( $argument =~ /^-salt/i ) {
			$salt_conc = $args{$argument};
		}
		elsif ( $argument =~ /^-mg/i ) {
			$mg_conc = $args{$argument};
		}
		elsif ( $argument =~ /^-bound/i ) {
			$bound = $args{$argument};
		}
		elsif ( $argument =~ /^-dntp/i ) {
			$dntp_conc = $args{$argument};
		}
		elsif ( $argument =~ /^-dplex/i ) {
			$dplex_path = $args{$argument};
		}
		elsif ( $argument =~ /^-sc/i ) {
			# 1: Owczarzy 2004
			# 2: SantaLucia 1998
			$salt_correct_formula = $args{$argument};
		}
		elsif ( $argument =~ /^-debug/i ) {
			$debug = 1 if $args{$argument};
		}
	}
	
	return 1;
}

sub Tm_approx {

	# Approximate Tm calculated function used on sequence length > 60 (DNA/DNA)
	my ( $seq ) = @_;
	$salt_conc = 0.05 unless ($salt_conc);

	# Formula from MELTING approximate mode
	my $tm = 81.5
			+ 16.6 * log( $salt_conc / ( 1.0 + 0.7 * $salt_conc ) ) / log(10)
			+ 0.41 * gc($seq) - 500.0 / length($seq);

	# Could change to Primer3 formula w/o salt correction
	# my $tm= 81.5+ ( 16.6* ( log($salt_conc)/log(10) ) ) + ( 0.41 * gc($seq) ) - ( 600 / length($seq) );

	return $tm;
}

sub _dG {
	my ($dH,$dS,$tm) = @_;
	$tm = $T unless $tm;
	return $dH - $tm*$dS/1000;
}

sub _K {
	my ($dG,$T) = @_;
	my $r = 1.987;
	my $e = 2.718281828;
	return $e**($dG/-$r*$T);
}

sub _dH_dS {
	my ( $sequence, $target, $type ) = @_;
	
	$type = "" unless ( defined $type );
	
	$sequence = uc $sequence;      # change sequence to upper case
	$sequence =~ s/\+(\w)/\l$1/g;  # change bases go with "+" to lower case as LNA
	
	# default target is reverse-complimentary sequence
	unless ( defined $target ){
		$target = uc &com_rev( $sequence, 'c' );
		$target =~ s/[\+\*]//g;
	}
	else{
		$target = uc &com_rev( $target, 'c' ); 
	}
	
	my @snn = split( //, $sequence );
	my @tnn = split( //, $target );
	my $dH = 0;
	my $dS = 0;
	my @dH_debug;
	my @dS_debug;
	my $default_enthalpy = 0.0;
	my $default_entropy	 = 0.0;

    my ($sbg, $tbg, $mm) = (0, 0, 0);

	for ( my $i=0; defined $snn[$i+1]; $i++ ) {
		my $cur = $snn[$i].$tnn[$i];
		my $nex = $snn[$i+1].$tnn[$i+1];
		
        # dH/dS of loop terminal 
        my $ildh_val = defined $param_loop_terminal_H{$cur}{$nex} ? $param_loop_terminal_H{$cur}{$nex} : $default_enthalpy;
        my $ilds_val = defined $param_loop_terminal_S{$cur}{$nex} ? $param_loop_terminal_S{$cur}{$nex} : $default_entropy;
      
        # dH/dS of internal loop
        if( $mm > 1 && $snn[$i] eq &com_rev($tnn[$i],'c') ){
            $dS += $param_loop_S{$mm};
            $mm = 0;
        }

        $mm++ if( $cur ne "-" && $nex ne "-" && $snn[$i] ne &com_rev($tnn[$i],'c') );

        # deal with bulge
        if( $sbg > 1 && $snn[$i] ne "-" ){
            $dS += $param_bulge_S{$sbg};
            $sbg = 0;
        }
        if( $tbg > 1 && $tnn[$i] ne "-" ){
            $dS += $param_bulge_S{$tbg};
            $tbg = 0;
        }

        $sbg ++ if $snn[$i] eq "-";
        $tbg ++ if $tnn[$i] eq "-";

        $dH += defined $param_H{$cur}{$nex} ? $param_H{$cur}{$nex} : $ildh_val;
		$dS	+= defined $param_S{$cur}{$nex} ? $param_S{$cur}{$nex} : $ilds_val;
 
		push @dS_debug, $param_S{$cur}{$nex};
	}
	
	return ($dH,$dS);
}
sub equivalent_na_conc {
	# Convert concentration of divalent cations to concentration of monovalent cations (Ahsen et al., 2001)
	# [Monovalent cations] = [Monovalent cations] + 120*(Ã([divalent cations] - [dNTP]))
	# Note that the concentrations are in mMol/L.
	# if [dntp conc] > [Mg++], notice from Primer3 manul " According to the formula concentration of desoxynucleotide triphosphate [dNTP]
	#                                                      must be smaller than concentration of divalent cations. The concentration of dNTPs
	#                                                      is included to the formula beacause of some magnesium is bound by the dNTP."
	
	my $salt_conc = shift;
    if( $mg_conc > $dntp_conc ){
        my $salt_mm = $salt_conc * 1000;
        my $mg_mm   = $mg_conc * 1000;
        my $dntp_mm = $dntp_conc * 1000;
        $salt_mm += 120 * sqrt( $mg_mm - $dntp_mm );

        #$salt_conc += 120 * sqrt( $mg_conc - $dntp_conc );

        $salt_conc = $salt_mm/1000;
    }

    return $salt_conc;
}

sub Tm_main {

	# Input sequence, salt concentration, oligo concentration and oligo concentration correaction factor as parameters
	# return Tm using nearest-neighbor parameters if sequence length <= 60 or
	# return an approximate Tm if sequence length > 60
	my ( $sequence, $target ) = @_;

	my ($length,$gc,$seq,$tar) = _prep_seq($sequence, $target);
	my $tm;
	
	# Use approximate mode if sequence length > 60
	return Tm_approx($seq) if ( $length > 60 );

	# retrieve dH, dS
	my ( $dH, $dS ) = _dH_dS( $seq, $tar );

	_debug( "dH", $dH ) if $debug;
	_debug( "dS(orig)", $dS ) if $debug;
	_debug( "GC%", $gc ) if $debug;
	_debug( "Length", $length ) if $debug;
	
	if ( $salt_correct_formula =~ /2/ ) {
		my $salt_correction = 0;
		
		# Convert concentration of divalent cations to concentration of monovalent cations (Ahsen et al., 2001)
		my $salt_conc_cor = equivalent_na_conc($salt_conc);
	
		# Correction for the concentration of salt:
		# SantaLucia et al. 1996 i.e. 12.5*log[Na+]
        #$salt_correction = 12.5 * (log($salt_conc)/log(10));
		
		# SantaLucia 1998
		$dS += 0.368 * ( $length - 1 ) * log($salt_conc_cor);
		
		# For selfcomplementary oligonucleotide duplexes
		$tm = $dH * 1000 / ( $dS + $r * log($oligo_conc/4) ) + $salt_correction;
		
	}
	else{
		# Tricky.... 
		# To get the same result as IDT web tool, Tm(1M) need to minus 1 before salt correction.
		
		$tm = $dH * 1000 / ( $dS + $r * log( (1-$bound)/$bound*$oligo_conc ) );
		
		#
		#debug
		#
		my $salt_conc_debug = equivalent_na_conc($salt_conc);
		$dS += 0.368 * ( $length - 1 ) * log($salt_conc_debug);
		_debug( "dS(w/ SantaLucia 1998)", $dS ) if $debug;

		# Deoxynucleoside triphosphates (dNTPs) correction
		my $dG = _dG( $dH, $dS );
		my $Ka = 3*(10**4);

		_bound( $dH, $dS, $tm );	

		_debug( "Na+", $salt_conc ) if $debug;
		_debug( "Na+ (w/ Ahsen Mg2+ correction)", equivalent_na_conc($salt_conc) ) if $debug;
		_debug( "Mg2+", $mg_conc ) if $debug;

		my $D = ($Ka*$dntp_conc - $Ka*$mg_conc + 1)**2 + 4*$Ka*$mg_conc;
		my $mg_conc_c = ( sqrt($D) - ( $Ka*$dntp_conc - $Ka*$mg_conc + 1 ) )/(2*$Ka);
		
		_debug( "Tm(1M,orig K)", $tm ) if $debug;
		_debug( "Tm(1M,orig C)", $tm-273.15 ) if $debug;
		_debug( "dG", $dG ) if $debug;
		_debug( "D", $D ) if $debug;
		_debug( "Mg2+(corrected)", $mg_conc_c ) if $debug;

		if ( sqrt($mg_conc_c)/$salt_conc < 0.22 || $salt_conc == 0 ) {
			# IDTDNA OligoAnalyzer employs the improved quadratic TM salt correction function (Owczarzy,R. et al., Biochemistry, 43, 3537)
			# 1/Tm[Na+] = 1/Tm[1M Na+] + [( 4.29 * f(CG) - 3.95 ) * ln[Na+] + ( 0.94 * ln^2[Na+] )]*1e-5
			# where f(GC) is the fraction of GC base pairs
						
			my $salt_correction = ( ( 4.29 * $gc - 3.95 )*log($salt_conc) + 0.94*(log($salt_conc)**2) )*1e-5;
			$tm = 1 / ( 1/$tm + $salt_correction );
			
			_debug( "Owczarzy salt correction", "R < 0.22 OR [NA+] = 0" ) if $debug;
		}
		else{
			# IDTDNA Tm correction for divalent (Mg++) ions (Owczarzy, R. et al., Biochemistry, 47, 5336),
			my ($a, $d, $g) = (3.92, 1.42, 8.31);
			
			# R value is in the range from 0.22 to 6.0, using a, d, and g coefficients that varies with Na+ concentration
			if( sqrt($mg_conc_c)/$salt_conc <= 6 ){
				$a = 3.92 * ( 0.843 - 0.352 * sqrt($salt_conc) * log($salt_conc) );
				$d = 1.42 * ( 1.279 - 4.03 * 1e-3 * log($salt_conc) - 8.03 * 1e-3 * log($salt_conc)**2 );
				$g = 8.31 * ( 0.486 - 0.258 * log($salt_conc) + 5.25 * 1e-3 * log($salt_conc)**3 );		
			}
			
			my $salt_correction = ( $a
									- 0.911*log($mg_conc_c) 
									+ $gc*( 6.26 + $d*log($mg_conc_c) )
									+ ( 1/2/( $length - 1 ) )*( -48.2 + 52.5*log($mg_conc_c) + $g*(log($mg_conc_c)**2) )
								)*1e-5;
	
			$tm = 1 / ( 1/$tm + $salt_correction );
			
			_debug( "Owczarzy salt correction", "R > 0.22" ) if $debug;
			_debug( "a", $a ) if $debug;
			_debug( "d", $d ) if $debug;
			_debug( "g", $g ) if $debug;
		}
	}
	
	# Convert degree K to degree C
	$tm =  $tm - 273.15;
	
	return $tm;
}

sub _prep_seq {
    my ($seq,$tar) = @_;
    my ($length,$gc);

    #oligo
    if( $seq =~ /^-/ ){
        $seq =~ s/^-//;
    }
    else{
        $seq = "E$seq";
    }

    if( $seq =~ /-$/ ){
        $seq =~ s/-$//;
    }
    else{
        $seq = "${seq}E";
    }
    
    #target
    if( $tar =~ /^-/ ){
        $tar =~ s/^-//;
    }
    else{
        $tar = "E$tar";
    }

    if( $tar =~ /-$/ ){
        $tar =~ s/-$//;
    }
    else{
        $tar = "${tar}E";
    }

    #calc gc%
    my ($real_seq) = $seq =~ /^\w(.+)\w$/;
    $real_seq =~ s/-//g;
    $length = length($real_seq);
    $gc = gc($real_seq);

    if ( length($seq) != length($tar) ){
        print STDERR "WARNING: The length of oligo/target doesn't match. Ignored target sequence.\n";
        $tar = $seq;
    }

    return ($length, $gc, $seq, $tar);
}

sub _bound {
	my ( $dH, $dS, $tm ) = @_;
	
	my $dG = _dG( $dH, $dS, $tm );
	my $K = exp( -$dG/($r*$tm) );
	$K = exp( -12100/($r*335.24) );
	_debug( "K", $K ) if $debug;
}

sub _debug {
	my ( $text, $value ) = @_;
	print STDERR "$text = $value\n";
}

sub Tm1 {
	my ( $sequence ) = @_;
	
	my $salt_mg_conc = $salt_conc + ( 120 * ( ( $mg_conc > $dntp_conc ) ? sqrt( ( $mg_conc - $dntp_conc ) * 1000 ) : 0 ) ) / 1000;
	
	my $mo = $sequence;
	$mo = "1 $mo";
	my @dplex_result = qx(echo "$mo" | $dplex_path --salt $salt_mg_conc --strand $oligo_conc --tm 90 --no-opt) or die "ERROR: dplex error.";
	
	foreach my $line (@dplex_result){
		if ($line =~ /Min Tm = ([^ ]+) C/){
			my $tm = $1;
			return $tm;
		}
	}
}

sub Hair_Tm {
	my ( $sequence ) = @_;
	
	my $salt_mg_conc = $salt_conc + ( 120 * ( ( $mg_conc > $dntp_conc ) ? sqrt( ( $mg_conc - $dntp_conc ) * 1000 ) : 0 ) ) / 1000;
	
	my $mo = $sequence;
	$mo = "1 $mo";
	my @dplex_result = qx(echo "$mo" | $dplex_path --salt $salt_mg_conc --strand $oligo_conc --tm 90 --no-opt) or die "ERROR: dplex error.";
	
	foreach my $line (@dplex_result){
		if ($line =~ /Max Hairpin Tm = ([^ ]+) C/){
			my $tm = $1;
			return $tm;
		}
	}
}

sub Tm {
	# Return a list (MAX, AVG, MIN) Tms for degenerated sepquence or a Tm value for regular sequence.
	my ( $seq, $target ) = @_;
	
	$seq =~ s/\s//g if defined $seq;
	$target =~ s/\s//g if defined $target;
    $target = $seq unless defined $target;

    my @seqPool = &expandDegen($seq);    # expand input sequence w or w/o degenerate base(s)
	my @tarPool = &expandDegen($target);
	my @tmPool;

	# Calculating Tm for all non-degenerated seuquencces
	foreach my $exp_seq (@seqPool) {
		foreach my $exp_tar (@tarPool) {
			push @tmPool, Tm_main($exp_seq, $exp_tar);
		}
	}

	@tmPool = sort(@tmPool);
	my $avg = ( eval join "+", @tmPool ) / ( $#tmPool + 1 );
	#return ( $tmPool[0], $avg, pop(@tmPool) );
	
    $avg = 0 if $avg < 0;
    return sprintf "%.2f",$avg;
}

sub expandDegen
{
	# Original sequence as a parameter and return an array
	my $seq = shift;
	my @seqPool;

	# Convert IUB codes to bases:
	if ( $seq =~ /[RYKMSWHBVDN]/ ) {
		$seq =~ s/R/[AG]/g;
		$seq =~ s/Y/[CT]/g;
		$seq =~ s/K/[GT]/g;
		$seq =~ s/M/[AC]/g;
		$seq =~ s/S/[CG]/g;
		$seq =~ s/W/[AT]/g;
		$seq =~ s/H/[ACT]/g;
		$seq =~ s/B/[CGT]/g;
		$seq =~ s/V/[ACG]/g;
		$seq =~ s/D/[AGT]/g;
		$seq =~ s/N/[ACGT]/g;
	}
	else { push( @seqPool, $seq ); }    # no degenerate bases found

	if ( $seq =~ /^(.*)\[\w+\]/ ) {
		push @seqPool, $1;
		$seq =~ s/^$1//;
	}

	while ( $seq =~ /^\[(\w+)\]([^\[]*)/ ) {
		my @seqPoolTemp;
		my @degBase = split( '', $1 );
		my $seqFrag = $2;

		foreach my $tempSeq (@seqPool) {
			foreach my $base (@degBase) {    #degenerate bases
				my $str = $tempSeq . $base . $seqFrag;
				push @seqPoolTemp, $str;
			}
		}
		$seq =~ s/^\[\w+\][^\[]*//;
		@seqPool = @seqPoolTemp;
	}
	return @seqPool;
}

sub com_rev {
	my ( $seq, $arg ) = @_;
	if ( $arg =~ /c/i ) {    #complement
		$seq =~ tr/ACGTURYSWKMBDHVNacgturyswkmbdhvn/TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn/;
	}

	if ( $arg =~ /r/i ) {    #reserve
		$seq = reverse($seq);
	}

	return $seq;
}

sub gc {
	# calculate the gc content of the chunk of sequence passed as a parameter
	my $seq = shift;
	return ( $seq =~ tr/gGcC// ) / length($seq);
}

sub buildThermoParameters {
	my %param_H;
	my %param_S;
	
	# Build a hash with the thermodynamic values
	# Nearest-neighbor parameters all97a (Allawi et al 1997), J/mol-l

	#_dH
	my $AT_AT = -7.9;
	my $AT_CG = -8.4;
	my $AT_GC = -7.8;
	my $AT_TA = -7.2;
	my $CG_AT = -8.5;
	my $CG_CG = -8.0;
	my $CG_GC = -10.6;
	my $GC_AT = -8.2;
	my $GC_CG = -9.8;
	my $TA_AT = -7.2;
	
	$param_H{AT}{AT} = $param_H{TA}{TA} = $AT_AT; 	# AA/TT
	$param_H{AT}{CG} = $param_H{GC}{TA} = $AT_CG; 	# AC/TG
	$param_H{AT}{GC} = $param_H{CG}{TA} = $AT_GC; 	# AG/TC
	$param_H{AT}{TA}                    = $AT_TA;   # AT/TA
	$param_H{CG}{AT} = $param_H{TA}{GC} = $CG_AT; 	# CA/GT
	$param_H{CG}{CG} = $param_H{GC}{GC} = $CG_CG; 	# CC/GG
	$param_H{CG}{GC}                    = $CG_GC;	# CG/GC
	$param_H{GC}{AT} = $param_H{TA}{CG} = $GC_AT;	# GA/CT
	$param_H{GC}{CG}                    = $GC_CG;	# GC/CG
	$param_H{TA}{AT}                    = $TA_AT;	# TA/AT
	
	$param_S{AT}{AT} = $param_S{TA}{TA} = -22.2;	# AA/TT
	$param_S{AT}{CG} = $param_S{GC}{TA} = -22.4;	# AC/TG
	$param_S{AT}{GC} = $param_S{CG}{TA} = -21.0;	# AG/TC
	$param_S{AT}{TA}                    = -20.4;	# AT/TA
	$param_S{CG}{AT} = $param_S{TA}{GC} = -22.7;	# CA/GT
	$param_S{CG}{CG} = $param_S{GC}{GC} = -19.9;	# CC/GG
	$param_S{CG}{GC}                    = -27.2;	# CG/GC
	$param_S{GC}{AT} = $param_S{TA}{CG} = -22.2;	# GA/CT
	$param_S{GC}{CG}                    = -24.4;	# GC/CG
	$param_S{TA}{AT}                    = -21.3;	# TA/AT
	
	my $AE_AT =  0.2;
	my $AE_CG = -6.3;
	my $AE_GC = -3.7;
	my $AE_TA = -2.9;
	my $CE_AT =  0.6;
	my $CE_CG = -4.4;
	my $CE_GC = -4.0;
	my $CE_TA = -4.1;
	my $GE_AT = -1.1;
	my $GE_CG = -5.1;
	my $GE_GC = -3.9;
	my $GE_TA = -4.2;
	my $TE_AT = -6.9;
	my $TE_CG = -4.0;
	my $TE_GC = -4.9;
	my $TE_TA = -0.2;
	
	$param_H{AE}{AT} = $param_H{TA}{EA} = $AE_AT;	# AA/ET
	$param_H{AE}{CG} = $param_H{GC}{EA} = $AE_CG;	# AC/EG
	$param_H{AE}{GC} = $param_H{CG}{EA} = $AE_GC;	# AG/EC
	$param_H{AE}{TA} = $param_H{AT}{EA} = $AE_TA;	# AT/EA
	$param_H{CE}{AT} = $param_H{TA}{EC} = $CE_AT;	# CA/ET
	$param_H{CE}{CG} = $param_H{GC}{EC} = $CE_CG;	# CC/EG
	$param_H{CE}{GC} = $param_H{CG}{EC} = $CE_GC;	# CG/EC
	$param_H{CE}{TA} = $param_H{AT}{EC} = $CE_TA;	# CT/EA
	$param_H{GE}{AT} = $param_H{TA}{EG} = $GE_AT;	# GA/ET
	$param_H{GE}{CG} = $param_H{GC}{EG} = $GE_CG;	# GC/EG
	$param_H{GE}{GC} = $param_H{CG}{EG} = $GE_GC;	# GG/EC
	$param_H{GE}{TA} = $param_H{AT}{EG} = $GE_TA;	# GT/EA
	$param_H{TE}{AT} = $param_H{TA}{ET} = $TE_AT;	# TA/ET
	$param_H{TE}{CG} = $param_H{GC}{ET} = $TE_CG;	# TC/EG
	$param_H{TE}{GC} = $param_H{CG}{ET} = $TE_GC;	# TG/EC
	$param_H{TE}{TA} = $param_H{AT}{ET} = $TE_TA;	# TT/EA
	
	$param_S{AE}{AT} = $param_S{TA}{EA} = _ENTROPY(-0.51, $AE_AT);	# AA/ET
	$param_S{AE}{CG} = $param_S{GC}{EA} = _ENTROPY(-0.96, $AE_CG);	# AC/EG
	$param_S{AE}{GC} = $param_S{CG}{EA} = _ENTROPY(-0.58, $AE_GC);	# AG/EC
	$param_S{AE}{TA} = $param_S{AT}{EA} = _ENTROPY(-0.50, $AE_TA);	# AT/EA
	$param_S{CE}{AT} = $param_S{TA}{EC} = _ENTROPY(-0.42, $CE_AT);	# CA/ET
	$param_S{CE}{CG} = $param_S{GC}{EC} = _ENTROPY(-0.52, $CE_CG);	# CC/EG
	$param_S{CE}{GC} = $param_S{CG}{EC} = _ENTROPY(-0.34, $CE_GC);	# CG/EC
	$param_S{CE}{TA} = $param_S{AT}{EC} = _ENTROPY(-0.02, $CE_TA);	# CT/EA
	$param_S{GE}{AT} = $param_S{TA}{EG} = _ENTROPY(-0.62, $GE_AT);	# GA/ET
	$param_S{GE}{CG} = $param_S{GC}{EG} = _ENTROPY(-0.72, $GE_CG);	# GC/EG
	$param_S{GE}{GC} = $param_S{CG}{EG} = _ENTROPY(-0.56, $GE_GC);	# GG/EC
	$param_S{GE}{TA} = $param_S{AT}{EG} = _ENTROPY( 0.48, $GE_TA);	# GT/EA
	$param_S{TE}{AT} = $param_S{TA}{ET} = _ENTROPY(-0.71, $TE_AT);	# TA/ET
	$param_S{TE}{CG} = $param_S{GC}{ET} = _ENTROPY(-0.58, $TE_CG);	# TC/EG
	$param_S{TE}{GC} = $param_S{CG}{ET} = _ENTROPY(-0.61, $TE_GC);	# TG/EC
	$param_S{TE}{TA} = $param_S{AT}{ET} = _ENTROPY(-0.10, $TE_TA);	# TT/EA
	
	my $EA_AT = -0.7;
	my $EA_CG = -2.1;
	my $EA_GC = -5.9;
	my $EA_TA = -0.5;
	my $EC_AT =  4.4;
	my $EC_CG = -0.2;
	my $EC_GC = -2.6;
	my $EC_TA =  4.7;
	my $EG_AT = -1.6;
	my $EG_CG = -3.9;
	my $EG_GC = -3.2;
	my $EG_TA = -4.1;
	my $ET_AT =  2.9;
	my $ET_CG = -4.4;
	my $ET_GC = -5.2;
	my $ET_TA = -3.8;
	
	$param_H{EA}{AT} = $param_H{TA}{AE} = $EA_AT;	# EA/AT
	$param_H{EA}{CG} = $param_H{GC}{AE} = $EA_CG;	# EC/AG
	$param_H{EA}{GC} = $param_H{CG}{AE} = $EA_GC;	# EG/AC
	$param_H{EA}{TA} = $param_H{AT}{AE} = $EA_TA;	# ET/AA
	$param_H{EC}{AT} = $param_H{TA}{CE} = $EC_AT;	# EA/CT
	$param_H{EC}{CG} = $param_H{GC}{CE} = $EC_CG;	# EC/CG
	$param_H{EC}{GC} = $param_H{CG}{CE} = $EC_GC;	# EG/CC
	$param_H{EC}{TA} = $param_H{AT}{CE} = $EC_TA;	# ET/CA
	$param_H{EG}{AT} = $param_H{TA}{GE} = $EG_AT;	# EA/GT
	$param_H{EG}{CG} = $param_H{GC}{GE} = $EG_CG;	# EC/GG
	$param_H{EG}{GC} = $param_H{CG}{GE} = $EG_GC;	# EG/GC
	$param_H{EG}{TA} = $param_H{AT}{GE} = $EG_TA;	# ET/GA
	$param_H{ET}{AT} = $param_H{TA}{TE} = $ET_AT;	# EA/TT
	$param_H{ET}{CG} = $param_H{GC}{TE} = $ET_CG;	# EC/TG
	$param_H{ET}{GC} = $param_H{CG}{TE} = $ET_GC;	# EG/TC
	$param_H{ET}{TA} = $param_H{AT}{TE} = $ET_TA;	# ET/TA
	
	$param_S{EA}{AT} = $param_S{TA}{AE} = _ENTROPY(-0.48, $EA_AT);	# EA/AT
	$param_S{EA}{CG} = $param_S{GC}{AE} = _ENTROPY(-0.92, $EA_CG);	# EC/AG
	$param_S{EA}{GC} = $param_S{CG}{AE} = _ENTROPY(-0.82, $EA_GC);	# EG/AC
	$param_S{EA}{TA} = $param_S{AT}{AE} = _ENTROPY(-0.12, $EA_TA);	# ET/AA
	$param_S{EC}{AT} = $param_S{TA}{CE} = _ENTROPY(-0.19, $EC_AT);	# EA/CT
	$param_S{EC}{CG} = $param_S{GC}{CE} = _ENTROPY(-0.23, $EC_CG);	# EC/CG
	$param_S{EC}{GC} = $param_S{CG}{CE} = _ENTROPY(-0.31, $EC_GC);	# EG/CC
	$param_S{EC}{TA} = $param_S{AT}{CE} = _ENTROPY( 0.28, $EC_TA);	# ET/CA
	$param_S{EG}{AT} = $param_S{TA}{GE} = _ENTROPY(-0.50, $EG_AT);	# EA/GT
	$param_S{EG}{CG} = $param_S{GC}{GE} = _ENTROPY(-0.44, $EG_CG);	# EC/GG
	$param_S{EG}{GC} = $param_S{CG}{GE} = _ENTROPY(-0.01, $EG_GC);	# EG/GC
	$param_S{EG}{TA} = $param_S{AT}{GE} = _ENTROPY(-0.01, $EG_TA);	# ET/GA
	$param_S{ET}{AT} = $param_S{TA}{TE} = _ENTROPY(-0.29, $ET_AT);	# EA/TT
	$param_S{ET}{CG} = $param_S{GC}{TE} = _ENTROPY(-0.35, $ET_CG);	# EC/TG
	$param_S{ET}{GC} = $param_S{CG}{TE} = _ENTROPY(-0.52, $ET_GC);	# EG/TC
	$param_S{ET}{TA} = $param_S{AT}{TE} = _ENTROPY( 0.13, $ET_TA);	# ET/TA
	
	my $AT_AG = -0.6;
	my $AT_GA = -0.7;
	my $CG_AG = -0.7;
	my $CG_GA = -4.0;
	my $GC_AG = -0.6;
	my $GC_GA =  0.5;
	my $TA_AG =  0.7;
	my $TA_GA =  3.0;

	$param_H{AE}{AT} = $param_H{TA}{EA} = $AE_AT;	# AA/ET
	$param_H{AE}{CG} = $param_H{GC}{EA} = $AE_CG;	# AC/EG
	$param_H{AE}{GC} = $param_H{CG}{EA} = $AE_GC;	# AG/EC
	$param_H{AE}{TA} = $param_H{AT}{EA} = $AE_TA;	# AT/EA
	$param_H{CE}{AT} = $param_H{TA}{EC} = $CE_AT;	# CA/ET
	$param_H{CE}{CG} = $param_H{GC}{EC} = $CE_CG;	# CC/EG
	$param_H{CE}{GC} = $param_H{CG}{EC} = $CE_GC;	# CG/EC
	$param_H{CE}{TA} = $param_H{AT}{EC} = $CE_TA;	# CT/EA
	$param_H{GE}{AT} = $param_H{TA}{EG} = $GE_AT;	# GA/ET
	$param_H{GE}{CG} = $param_H{GC}{EG} = $GE_CG;	# GC/EG
	$param_H{GE}{GC} = $param_H{CG}{EG} = $GE_GC;	# GG/EC
	$param_H{GE}{TA} = $param_H{AT}{EG} = $GE_TA;	# GT/EA
	$param_H{TE}{AT} = $param_H{TA}{ET} = $TE_AT;	# TA/ET
	$param_H{TE}{CG} = $param_H{GC}{ET} = $TE_CG;	# TC/EG
	$param_H{TE}{GC} = $param_H{CG}{ET} = $TE_GC;	# TG/EC
	$param_H{TE}{TA} = $param_H{AT}{ET} = $TE_TA;	# TT/EA
	
	$param_S{AE}{AT} = $param_S{TA}{EA} = _ENTROPY(-0.51, $AE_AT);	# AA/ET
	$param_S{AE}{CG} = $param_S{GC}{EA} = _ENTROPY(-0.96, $AE_CG);	# AC/EG
	$param_S{AE}{GC} = $param_S{CG}{EA} = _ENTROPY(-0.58, $AE_GC);	# AG/EC
	$param_S{AE}{TA} = $param_S{AT}{EA} = _ENTROPY(-0.50, $AE_TA);	# AT/EA
	$param_S{CE}{AT} = $param_S{TA}{EC} = _ENTROPY(-0.42, $CE_AT);	# CA/ET
	$param_S{CE}{CG} = $param_S{GC}{EC} = _ENTROPY(-0.52, $CE_CG);	# CC/EG
	$param_S{CE}{GC} = $param_S{CG}{EC} = _ENTROPY(-0.34, $CE_GC);	# CG/EC
	$param_S{CE}{TA} = $param_S{AT}{EC} = _ENTROPY(-0.02, $CE_TA);	# CT/EA
	$param_S{GE}{AT} = $param_S{TA}{EG} = _ENTROPY(-0.62, $GE_AT);	# GA/ET
	$param_S{GE}{CG} = $param_S{GC}{EG} = _ENTROPY(-0.72, $GE_CG);	# GC/EG
	$param_S{GE}{GC} = $param_S{CG}{EG} = _ENTROPY(-0.56, $GE_GC);	# GG/EC
	$param_S{GE}{TA} = $param_S{AT}{EG} = _ENTROPY( 0.48, $GE_TA);	# GT/EA
	$param_S{TE}{AT} = $param_S{TA}{ET} = _ENTROPY(-0.71, $TE_AT);	# TA/ET
	$param_S{TE}{CG} = $param_S{GC}{ET} = _ENTROPY(-0.58, $TE_CG);	# TC/EG
	$param_S{TE}{GC} = $param_S{CG}{ET} = _ENTROPY(-0.61, $TE_GC);	# TG/EC
	$param_S{TE}{TA} = $param_S{AT}{ET} = _ENTROPY(-0.10, $TE_TA);	# TT/EA
	
	# self-deplex dangling
	$param_H{EE}{AT} = $param_H{EE}{TA} = $param_H{AT}{EE} = $param_H{TA}{EE} = 2.3;	# EA/ET
	$param_H{EE}{CG} = $param_H{EE}{GC} = $param_H{CG}{EE} = $param_H{GC}{EE} = 0.1;	# EC/EG
	
	$param_S{EE}{AT} = $param_S{EE}{TA} = $param_S{AT}{EE} = $param_S{TA}{EE} = 4.1;	# EA/ET
	$param_S{EE}{CG} = $param_S{EE}{GC} = $param_S{CG}{EE} = $param_S{GC}{EE} = -2.8;	# EC/EG
	
	# Single G-A mismatch
	# Allawi et. al. Biochemistry 1998, 37, 2170-2179
	$param_H{AT}{AG} = $param_H{GA}{TA} = $AT_AG;	# Aa/Tg
	$param_H{AT}{GA} = $param_H{AG}{TA} = $AT_GA;	# Ag/Ta
	$param_H{CG}{AG} = $param_H{GA}{GC} = $CG_AG;	# Ca/Gg
	$param_H{CG}{GA} = $param_H{AG}{GC} = $CG_GA;	# Cg/Ga
	$param_H{GC}{AG} = $param_H{GA}{CG} = $GC_AG;	# Ga/Cg
	$param_H{GC}{GA} = $param_H{AG}{CG} = $GC_GA;	# Gg/Ca
	$param_H{TA}{AG} = $param_H{GA}{AT} = $TA_AG;	# Ta/Ag
	$param_H{TA}{GA} = $param_H{AG}{AT} = $TA_GA;	# Tg/Aa
	
	$param_S{AT}{AG} = $param_S{GA}{TA} = _ENTROPY( 0.14, $AT_AG);	# Aa/Tg
	$param_S{AT}{GA} = $param_S{AG}{TA} = _ENTROPY( 0.02, $AT_GA);	# Ag/Ta
	$param_S{CG}{AG} = $param_S{GA}{GC} = _ENTROPY( 0.03, $CG_AG);	# Ca/Gg
	$param_S{CG}{GA} = $param_S{AG}{GC} = _ENTROPY( 0.11, $CG_GA);	# Cg/Ga
	$param_S{GC}{AG} = $param_S{GA}{CG} = _ENTROPY(-0.25, $GC_AG);	# Ga/Cg
	$param_S{GC}{GA} = $param_S{AG}{CG} = _ENTROPY(-0.52, $GC_GA);	# Gg/Ca
	$param_S{TA}{AG} = $param_S{GA}{AT} = _ENTROPY( 0.42, $TA_AG);	# Ta/Ag
	$param_S{TA}{GA} = $param_S{AG}{AT} = _ENTROPY( 0.74, $TA_GA);	# Tg/Aa
	
	my $AT_CT =  0.7;
	my $AT_TC = -1.2;
	my $CG_CT = -0.8;
	my $CG_TC = -1.5;
	my $GC_CT =  2.3;
	my $GC_TC =  5.2;
	my $TA_CT =  1.2;
	my $TA_TC =  1.0;
	
	# Single C-T mismatch
	# Allawi et. al. Nucleic Acids Research, 1998, 26, 11, 2694-2701
	$param_H{AT}{CT} = $param_H{TC}{TA} = $AT_CT;	# Ac/Tt
	$param_H{AT}{TC} = $param_H{CT}{TA} = $AT_TC;	# At/Tc
	$param_H{CG}{CT} = $param_H{TC}{GC} = $CG_CT;	# Cc/Gt
	$param_H{CG}{TC} = $param_H{CT}{GC} = $CG_TC;	# Ct/Gc
	$param_H{GC}{CT} = $param_H{TC}{CG} = $GC_CT;	# Gc/Ct
	$param_H{GC}{TC} = $param_H{CT}{CG} = $GC_TC;	# Gt/Cc
	$param_H{TA}{CT} = $param_H{TC}{AT} = $TA_CT;	# Tc/At
	$param_H{TA}{TC} = $param_H{CT}{AT} = $TA_TC;	# Tt/Ac
	
	$param_S{AT}{CT} = $param_S{TC}{TA} = _ENTROPY( 0.64, $AT_CT);	# Ac/Tt
	$param_S{AT}{TC} = $param_S{CT}{TA} = _ENTROPY( 0.73, $AT_TC);	# At/Tc
	$param_S{CG}{CT} = $param_S{TC}{GC} = _ENTROPY( 0.62, $CG_CT);	# Cc/Gt
	$param_S{CG}{TC} = $param_S{CT}{GC} = _ENTROPY( 0.40, $CG_TC);	# Ct/Gc
	$param_S{GC}{CT} = $param_S{TC}{CG} = _ENTROPY( 0.62, $GC_CT);	# Gc/Ct
	$param_S{GC}{TC} = $param_S{CT}{CG} = _ENTROPY( 0.98, $GC_TC);	# Gt/Cc
	$param_S{TA}{CT} = $param_S{TC}{AT} = _ENTROPY( 0.97, $TA_CT);	# Tc/At
	$param_S{TA}{TC} = $param_S{CT}{AT} = _ENTROPY( 0.75, $TA_TC);	# Tt/Ac
	
	my $AT_AC = 2.3;
	my $AT_CA = 5.3;
	my $CG_AC = 1.9;
	my $CG_CA = 0.6;
	my $GC_AC = 5.2;
	my $GC_CA =-0.7;
	my $TA_AC = 3.4;
	my $TA_CA = 7.6;
	
	# Single A-C mismatch
	# Allawi et. al. Biochemistry 1998, 37, 9435-9444
	$param_H{AT}{AC} = $param_H{CA}{TA} = $AT_AC; 	# Aa/Tc
	$param_H{AT}{CA} = $param_H{AC}{TA} = $AT_CA;	# Ac/Ta
	$param_H{CG}{AC} = $param_H{CA}{GC} = $CG_AC;	# Ca/Gc
	$param_H{CG}{CA} = $param_H{AC}{GC} = $CG_CA;	# Cc/Ga
	$param_H{GC}{AC} = $param_H{CA}{CG} = $GC_AC;	# Ga/Cc
	$param_H{GC}{CA} = $param_H{AC}{CG} = $GC_CA;	# Gc/Ca
	$param_H{TA}{AC} = $param_H{CA}{AT} = $TA_AC;	# Ta/Ac
	$param_H{TA}{CA} = $param_H{AC}{AT} = $TA_CA;	# Tc/Aa
	
	$param_S{AT}{AC} = $param_S{CA}{TA} = _ENTROPY( 0.88, $AT_AC);	# Aa/Tc
	$param_S{AT}{CA} = $param_S{AC}{TA} = _ENTROPY( 0.77, $AT_CA);	# Ac/Ta
	$param_S{CG}{AC} = $param_S{CA}{GC} = _ENTROPY( 0.75, $CG_AC);	# Ca/Gc
	$param_S{CG}{CA} = $param_S{AC}{GC} = _ENTROPY( 0.79, $CG_CA);	# Cc/Ga
	$param_S{GC}{AC} = $param_S{CA}{CG} = _ENTROPY( 0.81, $GC_AC);	# Ga/Cc
	$param_S{GC}{CA} = $param_S{AC}{CG} = _ENTROPY( 0.47, $GC_CA);	# Gc/Ca
	$param_S{TA}{AC} = $param_S{CA}{AT} = _ENTROPY( 0.92, $TA_AC);	# Ta/Ac
	$param_S{TA}{CA} = $param_S{AC}{AT} = _ENTROPY( 1.33, $TA_CA);	# Tc/Aa
	
	my $AT_GT =  1.0;
	my $AT_TG = -2.5;
	my $CG_GT = -4.1;
	my $CG_TG = -2.8;
	my $GC_GT =  3.3;
	my $GT_GT =  5.8;
	my $GC_TG = -4.4;
	my $GT_TG =  4.1;
	my $TA_GT = -0.1;
	my $TG_GT = -1.4;
	my $TA_TG = -1.3;
	
	# Single G-T mismatch
	# Allawi et. al. Biochemistry 1997, 36, 10581-10594
	$param_H{AT}{GT} = $param_H{TG}{TA} = $AT_GT;	# Ag/Tt
	$param_H{AT}{TG} = $param_H{GT}{TA} = $AT_TG;	# At/Tg
	$param_H{CG}{GT} = $param_H{TG}{GC} = $CG_GT;	# Cg/Gt
	$param_H{CG}{TG} = $param_H{GT}{GC} = $CG_TG;	# Ct/Gg
	$param_H{GC}{GT} = $param_H{TG}{CG} = $GC_GT;	# Gg/Ct
	$param_H{GT}{GT} = $param_H{TG}{TG} = $GT_GT;	# gg/tt <-- double mismatch!
	$param_H{GC}{TG} = $param_H{GT}{CG} = $GC_TG;	# Gt/Cg
	$param_H{GT}{TG}                    = $GT_TG;	# gt/tg <-- double mismatch!
	$param_H{TA}{GT} = $param_H{TG}{AT} = $TA_GT;	# Tg/At
	$param_H{TG}{GT}                    = $TG_GT;	# tg/gt <-- double mismatch!
	$param_H{TA}{TG} = $param_H{GT}{AT} = $TA_TG;	# Tt/Ag
	
	$param_S{AT}{GT} = $param_S{TG}{TA} = _ENTROPY( 0.71, $AT_GT);	# Ag/Tt
	$param_S{AT}{TG} = $param_S{GT}{TA} = _ENTROPY( 0.07, $AT_TG);	# At/Tg
	$param_S{CG}{GT} = $param_S{TG}{GC} = _ENTROPY(-0.47, $CG_GT);	# Cg/Gt
	$param_S{CG}{TG} = $param_S{GT}{GC} = _ENTROPY(-0.32, $CG_TG);	# Ct/Gg
	$param_S{GC}{GT} = $param_S{TG}{CG} = _ENTROPY( 0.08, $GC_GT);	# Gg/Ct
	$param_S{GT}{GT} = $param_S{TG}{TG} = _ENTROPY( 0.74, $GT_GT);	# gg/tt <-- double mismatch!
	$param_S{GC}{TG} = $param_S{GT}{CG} = _ENTROPY(-0.59, $GC_TG);	# Gt/Cg
	$param_S{GT}{TG}                    = _ENTROPY( 1.15, $GT_TG);	# gt/tg <-- double mismatch!
	$param_S{TA}{GT} = $param_S{TG}{AT} = _ENTROPY( 0.43, $TA_GT);	# Tg/At
	$param_S{TG}{GT}                    = _ENTROPY( 0.52, $TG_GT);	# tg/gt <-- double mismatch!
	$param_S{TA}{TG} = $param_S{GT}{AT} = _ENTROPY( 0.34, $TA_TG);	# Tt/Ag
	
	my $AT_AA =  1.2;
	my $CG_AA = -0.9;
	my $GC_AA = -2.9;
	my $TA_AA =  4.7;
	
	# Single A-A, C-C, G-G or T-T mismatch
	# Peyret et. al. Biochemistry 1999, 38, 3468-3477
	$param_H{AT}{AA} = $param_H{AA}{TA} = $AT_AA;	# Aa/Ta
	$param_H{CG}{AA} = $param_H{AA}{GC} = $CG_AA;	# Ca/Ga
	$param_H{GC}{AA} = $param_H{AA}{CG} = $GC_AA;	# Ga/Ca
	$param_H{TA}{AA} = $param_H{AA}{AT} = $TA_AA;	# Ta/Aa
	
	$param_S{AT}{AA} = $param_S{AA}{TA} = _ENTROPY( 0.61, $AT_AA);	# Aa/Ta
	$param_S{CG}{AA} = $param_S{AA}{GC} = _ENTROPY( 0.43, $CG_AA);	# Ca/Ga
	$param_S{GC}{AA} = $param_S{AA}{CG} = _ENTROPY( 0.17, $GC_AA);	# Ga/Ca
	$param_S{TA}{AA} = $param_S{AA}{AT} = _ENTROPY( 0.69, $TA_AA);	# Ta/Aa
	
	my $AT_CC =  0.0;
	my $CG_CC = -1.5;
	my $GC_CC =  3.6;
	my $TA_CC =  6.1;
	
	$param_H{AT}{CC} = $param_H{CC}{TA} = $AT_CC;	# Ac/Tc
	$param_H{CG}{CC} = $param_H{CC}{GC} = $CG_CC;	# Cc/Gc
	$param_H{GC}{CC} = $param_H{CC}{CG} = $GC_CC;	# Gc/Cc
	$param_H{TA}{CC} = $param_H{CC}{AT} = $TA_CC;	# Tc/Ac
	
	$param_S{AT}{CC} = $param_S{CC}{TA} = _ENTROPY(1.33, $AT_CC);	# Ac/Tc
	$param_S{CG}{CC} = $param_S{CC}{GC} = _ENTROPY( 0.70, $CG_CC);	# Cc/Gc
	$param_S{GC}{CC} = $param_S{CC}{CG} = _ENTROPY( 0.79, $GC_CC);	# Gc/Cc
	$param_S{TA}{CC} = $param_S{CC}{AT} = _ENTROPY(1.05, $TA_CC);	# Tc/Ac
	
	my $AT_GG = -3.1;
	my $CG_GG = -4.9;
	my $GC_GG = -6.0;
	my $TA_GG =  1.6;
	
	$param_H{AT}{GG} = $param_H{GG}{TA} = $AT_GG;	# Ag/Tg
	$param_H{CG}{GG} = $param_H{GG}{GC} = $CG_GG;	# Cg/Gg
	$param_H{GC}{GG} = $param_H{GG}{CG} = $GC_GG;	# Gg/Cg
	$param_H{TA}{GG} = $param_H{GG}{AT} = $TA_GG;	# Tg/Ag
	
	$param_S{AT}{GG} = $param_S{GG}{TA} = _ENTROPY(-0.13, $AT_GG);	# Ag/Tg
	$param_S{CG}{GG} = $param_S{GG}{GC} = _ENTROPY(-0.11, $CG_GG);	# Cg/Gg
	$param_S{GC}{GG} = $param_S{GG}{CG} = _ENTROPY(-1.11, $GC_GG);	# Gg/Cg
	$param_S{TA}{GG} = $param_S{GG}{AT} = _ENTROPY( 0.44, $TA_GG);	# Tg/Ag
	
	my $AT_TT = -2.7;
	my $CG_TT = -5.0;
	my $GC_TT = -2.2;
	my $TA_TT =  0.2;
	
	$param_H{AT}{TT} = $param_H{TT}{TA} = $AT_TT;	# At/Tt
	$param_H{CG}{TT} = $param_H{TT}{GC} = $CG_TT;	# Ct/Gt
	$param_H{GC}{TT} = $param_H{TT}{CG} = $GC_TT;	# Gt/Ct
	$param_H{TA}{TT} = $param_H{TT}{AT} = $TA_TT;	# Tt/At
		
	$param_S{AT}{TT} = $param_S{TT}{TA} = _ENTROPY( 0.69, $AT_TT);	# At/Tt
	$param_S{CG}{TT} = $param_S{TT}{GC} = _ENTROPY(-0.12, $CG_TT);	# Ct/Gt
	$param_S{GC}{TT} = $param_S{TT}{CG} = _ENTROPY( 0.45, $GC_TT);	# Gt/Ct
	$param_S{TA}{TT} = $param_S{TT}{AT} = _ENTROPY( 0.68, $TA_TT);	# Tt/At
	
	# LNA modifications
	$param_H{aT}{AT} = -6.682;
	$param_H{aT}{TA} = -6.123;
	$param_H{aT}{GC} = -8.634;
	$param_H{aT}{CG} = -7.21;
	$param_H{tA}{AT} = -5.982;
	$param_H{tA}{TA} = -6.823;
	$param_H{tA}{GC} = -9.334;
	$param_H{tA}{CG} = -7.01;
	$param_H{gC}{AT} = -6.982;
	$param_H{gC}{TA} = -7.323;
	$param_H{gC}{GC} = -8.834;
	$param_H{gC}{CG} = -8.61;
	$param_H{cG}{AT} = -7.282;
	$param_H{cG}{TA} = -6.723;
	$param_H{cG}{GC} = -11.434;
	$param_H{cG}{CG} = -6.81;
	$param_H{AT}{aT} = -6.775;
	$param_H{AT}{tA} = -6.075;
	$param_H{AT}{gC} = -6.675;
	$param_H{AT}{cG} = -7.275;
	$param_H{TA}{aT} = -5.527;
	$param_H{TA}{tA} = -6.227;
	$param_H{TA}{gC} = -6.827;
	$param_H{TA}{cG} = -6.527;
	$param_H{GC}{aT} = -8.715;
	$param_H{GC}{tA} = -8.915;
	$param_H{GC}{gC} = -8.515;
	$param_H{GC}{cG} = -10.315;
	$param_H{CG}{aT} = -8.132;
	$param_H{CG}{tA} = -7.432;
	$param_H{CG}{gC} = -10.232;
	$param_H{CG}{cG} = -7.632;
	
	$param_S{aT}{AT} = -17.474;
	$param_S{aT}{TA} = -16.149;
	$param_S{aT}{GC} = -21.851;
	$param_S{aT}{CG} = -17.776;
	$param_S{tA}{AT} = -16.574;
	$param_S{tA}{TA} = -17.949;
	$param_S{tA}{GC} = -23.551;
	$param_S{tA}{CG} = -17.576;
	$param_S{gC}{AT} = -17.474;
	$param_S{gC}{TA} = -18.149;
	$param_S{gC}{GC} = -20.751;
	$param_S{gC}{CG} = -19.776;
	$param_S{cG}{AT} = -17.974;
	$param_S{cG}{TA} = -16.749;
	$param_S{cG}{GC} = -28.051;
	$param_S{cG}{CG} = -15.276;
	$param_S{AT}{aT} = -17.28;
	$param_S{AT}{tA} = -15.48;
	$param_S{AT}{gC} = -16.08;
	$param_S{AT}{cG} = -17.48;
	$param_S{TA}{aT} = -15.384;
	$param_S{TA}{tA} = -16.284;
	$param_S{TA}{gC} = -16.784;
	$param_S{TA}{cG} = -16.284;
	$param_S{GC}{aT} = -22.072;
	$param_S{GC}{tA} = -22.272;
	$param_S{GC}{gC} = -19.772;
	$param_S{GC}{cG} = -24.272;
	$param_S{CG}{aT} = -20.914;
	$param_S{CG}{tA} = -19.214;
	$param_S{CG}{gC} = -25.414;
	$param_S{CG}{cG} = -18.114;

    #internal loop
    $param_loop_terminal_H{AT}{AA} = -3.20;
    $param_loop_terminal_S{AT}{AA} = -0.00806061583104949;
    $param_loop_terminal_H{AA}{TA} = -3.20;
    $param_loop_terminal_S{AA}{TA} = -0.00806061583104949;
    
    $param_loop_terminal_H{AT}{AC} = -0.90;
    $param_loop_terminal_S{AT}{AC} = -0.00193454779945188;
    $param_loop_terminal_H{CA}{TA} = -0.90;
    $param_loop_terminal_S{CA}{TA} = -0.00193454779945188;
    
    $param_loop_terminal_H{AT}{AG} = -2.30;
    $param_loop_terminal_S{AT}{AG} = -0.00580364339835563;
    $param_loop_terminal_H{GA}{TA} = -2.30;
    $param_loop_terminal_S{GA}{TA} = -0.00580364339835563;
    
    $param_loop_terminal_H{AT}{AT} = 0.00;
    $param_loop_terminal_S{AT}{AT} = 0.00;
    $param_loop_terminal_H{TA}{TA} = 0.00;
    $param_loop_terminal_S{TA}{TA} = 0.00;
    
    $param_loop_terminal_H{AT}{CA} = -2.20;
    $param_loop_terminal_S{AT}{CA} = -0.00515879413187168;
    $param_loop_terminal_H{AC}{TA} = -2.20;
    $param_loop_terminal_S{AC}{TA} = -0.00515879413187168;
    
    $param_loop_terminal_H{AT}{CC} = -0.50;
    $param_loop_terminal_S{AT}{CC} = -0.000967273899725939;
    $param_loop_terminal_H{CC}{TA} = -0.50;
    $param_loop_terminal_S{CC}{TA} = -0.000967273899725939;
    
    $param_loop_terminal_H{AT}{CG} = 0.00;
    $param_loop_terminal_S{AT}{CG} = 0.00;
    $param_loop_terminal_H{GC}{TA} = 0.00;
    $param_loop_terminal_S{GC}{TA} = 0.00;
    
    $param_loop_terminal_H{AT}{CT} = -1.20;
    $param_loop_terminal_S{AT}{CT} = -0.00290182169917782;
    $param_loop_terminal_H{TC}{TA} = -1.20;
    $param_loop_terminal_S{TC}{TA} = -0.00290182169917782;
    
    $param_loop_terminal_H{AT}{GA} = -2.70;
    $param_loop_terminal_S{AT}{GA} = -0.00677091729808157;
    $param_loop_terminal_H{AG}{TA} = -2.70;
    $param_loop_terminal_S{AG}{TA} = -0.00677091729808157;
    
    $param_loop_terminal_H{AT}{GC} = 0.00;
    $param_loop_terminal_S{AT}{GC} = 0.00;
    $param_loop_terminal_H{CG}{TA} = 0.00;
    $param_loop_terminal_S{CG}{TA} = 0.00;
    
    $param_loop_terminal_H{AT}{GG} = -1.30;
    $param_loop_terminal_S{AT}{GG} = -0.00290182169917782;
    $param_loop_terminal_H{GG}{TA} = -1.30;
    $param_loop_terminal_S{GG}{TA} = -0.00290182169917782;
    
    $param_loop_terminal_H{AT}{GT} = -2.90;
    $param_loop_terminal_S{AT}{GT} = -0.00773819119780751;
    $param_loop_terminal_H{TG}{TA} = -2.90;
    $param_loop_terminal_S{TG}{TA} = -0.00773819119780751;
    
    $param_loop_terminal_H{AT}{TA} = 0.00;
    $param_loop_terminal_S{AT}{TA} = 0.00;
    $param_loop_terminal_H{AT}{TA} = 0.00;
    $param_loop_terminal_S{AT}{TA} = 0.00;
    
    $param_loop_terminal_H{AT}{TC} = -2.80;
    $param_loop_terminal_S{AT}{TC} = -0.00806061583104949;
    $param_loop_terminal_H{CT}{TA} = -2.80;
    $param_loop_terminal_S{CT}{TA} = -0.00806061583104949;
    
    $param_loop_terminal_H{AT}{TG} = -3.50;
    $param_loop_terminal_S{AT}{TG} = -0.00967273899725939;
    $param_loop_terminal_H{GT}{TA} = -3.50;
    $param_loop_terminal_S{GT}{TA} = -0.00967273899725939;
    
    $param_loop_terminal_H{AT}{TT} = -2.40;
    $param_loop_terminal_S{AT}{TT} = -0.00644849266483959;
    $param_loop_terminal_H{TT}{TA} = -2.40;
    $param_loop_terminal_S{TT}{TA} = -0.00644849266483959;
    
    $param_loop_terminal_H{CG}{AA} = -2.80;
    $param_loop_terminal_S{CG}{AA} = -0.00580364339835563;
    $param_loop_terminal_H{AA}{GC} = -2.80;
    $param_loop_terminal_S{AA}{GC} = -0.00580364339835563;
    
    $param_loop_terminal_H{CG}{AC} = -2.00;
    $param_loop_terminal_S{CG}{AC} = -0.00386909559890376;
    $param_loop_terminal_H{CA}{GC} = -2.00;
    $param_loop_terminal_S{CA}{GC} = -0.00386909559890376;
    
    $param_loop_terminal_H{CG}{AG} = -3.00;
    $param_loop_terminal_S{CG}{AG} = -0.00677091729808157;
    $param_loop_terminal_H{GA}{GC} = -3.00;
    $param_loop_terminal_S{GA}{GC} = -0.00677091729808157;
    
    $param_loop_terminal_H{CG}{AT} = 0.00;
    $param_loop_terminal_S{CG}{AT} = 0.00;
    $param_loop_terminal_H{TA}{GC} = 0.00;
    $param_loop_terminal_S{TA}{GC} = 0.00;
    
    $param_loop_terminal_H{CG}{CA} = -2.40;
    $param_loop_terminal_S{CG}{CA} = -0.00515879413187167;
    $param_loop_terminal_H{AC}{GC} = -2.40;
    $param_loop_terminal_S{AC}{GC} = -0.00515879413187167;
    
    $param_loop_terminal_H{CG}{CC} = -1.40;
    $param_loop_terminal_S{CG}{CC} = -0.00290182169917782;
    $param_loop_terminal_H{CC}{GC} = -1.40;
    $param_loop_terminal_S{CC}{GC} = -0.00290182169917782;
    
    $param_loop_terminal_H{CG}{CG} = 0.00;
    $param_loop_terminal_S{CG}{CG} = 0.00;
    $param_loop_terminal_H{GC}{GC} = 0.00;
    $param_loop_terminal_S{GC}{GC} = 0.00;
    
    $param_loop_terminal_H{CG}{CT} = -2.40;
    $param_loop_terminal_S{CG}{CT} = -0.00548121876511366;
    $param_loop_terminal_H{TC}{GC} = -2.40;
    $param_loop_terminal_S{TC}{GC} = -0.00548121876511366;
    
    $param_loop_terminal_H{CG}{GA} = -5.10;
    $param_loop_terminal_S{CG}{GA} = -0.0132194099629212;
    $param_loop_terminal_H{AG}{GC} = -5.10;
    $param_loop_terminal_S{AG}{GC} = -0.0132194099629212;
    
    $param_loop_terminal_H{CG}{GC} = 0.00;
    $param_loop_terminal_S{CG}{GC} = 0.00;
    $param_loop_terminal_H{CG}{GC} = 0.00;
    $param_loop_terminal_S{CG}{GC} = 0.00;
    
    $param_loop_terminal_H{CG}{GG} = -2.90;
    $param_loop_terminal_S{CG}{GG} = -0.00644849266483959;
    $param_loop_terminal_H{GG}{GC} = -2.90;
    $param_loop_terminal_S{GG}{GC} = -0.00644849266483959;
    
    $param_loop_terminal_H{CG}{GT} = -2.90;
    $param_loop_terminal_S{CG}{GT} = -0.00612606803159761;
    $param_loop_terminal_H{TG}{GC} = -2.90;
    $param_loop_terminal_S{TG}{GC} = -0.00612606803159761;
    
    $param_loop_terminal_H{CG}{TA} = 0.00;
    $param_loop_terminal_S{CG}{TA} = 0.00;
    $param_loop_terminal_H{AT}{GC} = 0.00;
    $param_loop_terminal_S{AT}{GC} = 0.00;
    
    $param_loop_terminal_H{CG}{TC} = -3.10;
    $param_loop_terminal_S{CG}{TC} = -0.00806061583104949;
    $param_loop_terminal_H{CT}{GC} = -3.10;
    $param_loop_terminal_S{CT}{GC} = -0.00806061583104949;
    
    $param_loop_terminal_H{CG}{TG} = -5.90;
    $param_loop_terminal_S{CG}{TG} = -0.016121231662099;
    $param_loop_terminal_H{GT}{GC} = -5.90;
    $param_loop_terminal_S{GT}{GC} = -0.016121231662099;
    
    $param_loop_terminal_H{CG}{TT} = -5.30;
    $param_loop_terminal_S{CG}{TT} = -0.0141866838626471;
    $param_loop_terminal_H{TT}{GC} = -5.30;
    $param_loop_terminal_S{TT}{GC} = -0.0141866838626471;
    
    $param_loop_terminal_H{GC}{AA} = -6.00;
    $param_loop_terminal_S{GC}{AA} = -0.016121231662099;
    $param_loop_terminal_H{AA}{CG} = -6.00;
    $param_loop_terminal_S{AA}{CG} = -0.016121231662099;
    
    $param_loop_terminal_H{GC}{AC} = -4.00;
    $param_loop_terminal_S{GC}{AC} = -0.0106400128969853;
    $param_loop_terminal_H{CA}{CG} = -4.00;
    $param_loop_terminal_S{CA}{CG} = -0.0106400128969853;
    
    $param_loop_terminal_H{GC}{AG} = -3.40;
    $param_loop_terminal_S{GC}{AG} = -0.00838304046429147;
    $param_loop_terminal_H{GA}{CG} = -3.40;
    
    $param_loop_terminal_H{GC}{AT} = 0.00;
    $param_loop_terminal_S{GC}{AT} = 0.00;
    $param_loop_terminal_H{TA}{CG} = 0.00;
    $param_loop_terminal_S{TA}{CG} = 0.00;
    
    $param_loop_terminal_H{GT}{AA} = 0.00;
    $param_loop_terminal_S{GT}{AA} = 0.0016121231662099;
    $param_loop_terminal_H{AA}{TG} = 0.00;
    $param_loop_terminal_S{AA}{TG} = 0.0016121231662099;
    
    $param_loop_terminal_H{GT}{AC} = 0.00;
    $param_loop_terminal_S{GT}{AC} = 0.000644849266483959;
    $param_loop_terminal_H{CA}{TG} = 0.00;
    $param_loop_terminal_S{CA}{TG} = 0.000644849266483959;
    
    $param_loop_terminal_H{GT}{AG} = 0.00;
    $param_loop_terminal_S{GT}{AG} = 0.0016121231662099;
    $param_loop_terminal_H{GA}{TG} = 0.00;
    $param_loop_terminal_S{GA}{TG} = 0.0016121231662099;
    
    $param_loop_terminal_H{GT}{AT} = -3.70;
    $param_loop_terminal_S{GT}{AT} = -0.00999516363050137;
    $param_loop_terminal_H{TA}{TG} = -3.70;
    $param_loop_terminal_S{TA}{TG} = -0.00999516363050137;
    
    $param_loop_terminal_H{GC}{CA} = -2.80;
    $param_loop_terminal_S{GC}{CA} = -0.00580364339835563;
    $param_loop_terminal_H{AC}{CG} = -2.80;
    $param_loop_terminal_S{AC}{CG} = -0.00580364339835563;
    
    $param_loop_terminal_H{GC}{CC} = -2.50;
    $param_loop_terminal_S{GC}{CC} = -0.00612606803159761;
    $param_loop_terminal_H{CC}{CG} = -2.50;
    $param_loop_terminal_S{CC}{CG} = -0.00612606803159761;
    
    $param_loop_terminal_H{GC}{CG} = 0.00;
    $param_loop_terminal_S{GC}{CG} = 0.00;
    $param_loop_terminal_H{GC}{CG} = 0.00;
    $param_loop_terminal_S{GC}{CG} = 0.00;
    
    $param_loop_terminal_H{GC}{CT} = -3.30;
    $param_loop_terminal_S{GC}{CT} = -0.00838304046429147;
    $param_loop_terminal_H{TC}{CG} = -3.30;
    $param_loop_terminal_S{TC}{CG} = -0.00838304046429147;
    
    $param_loop_terminal_H{GT}{CA} = 0.00;
    $param_loop_terminal_S{GT}{CA} = 0.000644849266483959;
    $param_loop_terminal_H{AC}{TG} = 0.00;
    $param_loop_terminal_S{AC}{TG} = 0.000644849266483959;
    
    $param_loop_terminal_H{GT}{CC} = 0.00;
    $param_loop_terminal_S{GT}{CC} = 0.000644849266483959;
    $param_loop_terminal_H{CC}{TG} = 0.00;
    $param_loop_terminal_S{CC}{TG} = 0.000644849266483959;
    
    $param_loop_terminal_H{GT}{CG} = -4.50;
    $param_loop_terminal_S{GT}{CG} = -0.0116072867967113;
    $param_loop_terminal_H{GC}{TG} = -4.50;
    $param_loop_terminal_S{GC}{TG} = -0.0116072867967113;
    
    $param_loop_terminal_H{GT}{CT} = 0.00;
    $param_loop_terminal_S{GT}{CT} = 0.000644849266483959;
    $param_loop_terminal_H{TC}{TG} = 0.00;
    $param_loop_terminal_S{TC}{TG} = 0.000644849266483959;
    
    $param_loop_terminal_H{GC}{GA} = -4.00;
    $param_loop_terminal_S{GC}{GA} = -0.00967273899725939;
    $param_loop_terminal_H{AG}{CG} = -4.00;
    $param_loop_terminal_S{AG}{CG} = -0.00967273899725939;
    
    $param_loop_terminal_H{GC}{GC} = 0.00;
    $param_loop_terminal_S{GC}{GC} = 0.00;
    $param_loop_terminal_H{CG}{CG} = 0.00;
    $param_loop_terminal_S{CG}{CG} = 0.00;
    
    $param_loop_terminal_H{GC}{GG} = -5.30;
    $param_loop_terminal_S{GC}{GG} = -0.0138642592294051;
    $param_loop_terminal_H{GG}{CG} = -5.30;
    $param_loop_terminal_S{GG}{CG} = -0.0138642592294051;
    
    $param_loop_terminal_H{GC}{GT} = -3.70;
    $param_loop_terminal_S{GC}{GT} = -0.00935031436401741;
    $param_loop_terminal_H{TG}{CG} = -3.70;
    $param_loop_terminal_S{TG}{CG} = -0.00935031436401741;
    
    $param_loop_terminal_H{GT}{GA} = 0.00;
    $param_loop_terminal_S{GT}{GA} = 0.0016121231662099;
    $param_loop_terminal_H{AG}{TG} = 0.00;
    $param_loop_terminal_S{AG}{TG} = 0.0016121231662099;
    
    $param_loop_terminal_H{GT}{GC} = -5.90;
    $param_loop_terminal_S{GT}{GC} = -0.016121231662099;
    $param_loop_terminal_H{CG}{TG} = -5.90;
    $param_loop_terminal_S{CG}{TG} = -0.016121231662099;
    
    $param_loop_terminal_H{GT}{GG} = 0.00;
    $param_loop_terminal_S{GT}{GG} = 0.0016121231662099;
    $param_loop_terminal_H{GG}{TG} = 0.00;
    $param_loop_terminal_S{GG}{TG} = 0.0016121231662099;
    
    $param_loop_terminal_H{GT}{GT} = -2.00;
    $param_loop_terminal_S{GT}{GT} = -0.0048363694986297;
    $param_loop_terminal_H{TG}{TG} = -2.00;
    $param_loop_terminal_S{TG}{TG} = -0.0048363694986297;
    
    $param_loop_terminal_H{GC}{TA} = 0.00;
    $param_loop_terminal_S{GC}{TA} = 0.00;
    $param_loop_terminal_H{AT}{CG} = 0.00;
    $param_loop_terminal_S{AT}{CG} = 0.00;
    
    $param_loop_terminal_H{GC}{TC} = -2.50;
    $param_loop_terminal_S{GC}{TC} = -0.00612606803159761;
    $param_loop_terminal_H{CT}{CG} = -2.50;
    $param_loop_terminal_S{CT}{CG} = -0.00612606803159761;
    
    $param_loop_terminal_H{GC}{TG} = -4.50;
    $param_loop_terminal_S{GC}{TG} = -0.0116072867967113;
    $param_loop_terminal_H{GT}{CG} = -4.50;
    $param_loop_terminal_S{GT}{CG} = -0.0116072867967113;
    
    $param_loop_terminal_H{GC}{TT} = -6.10;
    $param_loop_terminal_S{GC}{TT} = -0.0167660809285829;
    $param_loop_terminal_H{TT}{CG} = -6.10;
    $param_loop_terminal_S{TT}{CG} = -0.0167660809285829;
    
    $param_loop_terminal_H{GT}{TA} = -3.50;
    $param_loop_terminal_S{GT}{TA} = -0.00967273899725939;
    $param_loop_terminal_H{AT}{TG} = -3.50;
    $param_loop_terminal_S{AT}{TG} = -0.00967273899725939;
    
    $param_loop_terminal_H{GT}{TC} = 0.00;
    $param_loop_terminal_S{GT}{TC} = 0.000644849266483959;
    $param_loop_terminal_H{CT}{TG} = 0.00;
    $param_loop_terminal_S{CT}{TG} = 0.000644849266483959;
    
    $param_loop_terminal_H{GT}{TG} = -2.00;
    $param_loop_terminal_S{GT}{TG} = -0.0048363694986297;
    $param_loop_terminal_H{GT}{TG} = -2.00;
    $param_loop_terminal_S{GT}{TG} = -0.0048363694986297;
    
    $param_loop_terminal_H{GT}{TT} = 0.00;
    $param_loop_terminal_S{GT}{TT} = 0.000644849266483959;
    $param_loop_terminal_H{TT}{TG} = 0.00;
    $param_loop_terminal_S{TT}{TG} = 0.000644849266483959;
    
    $param_loop_terminal_H{TA}{AA} = -2.70;
    $param_loop_terminal_S{TA}{AA} = -0.00677091729808157;
    $param_loop_terminal_H{AA}{AT} = -2.70;
    $param_loop_terminal_S{AA}{AT} = -0.00677091729808157;
    
    $param_loop_terminal_H{TA}{AC} = -2.60;
    $param_loop_terminal_S{TA}{AC} = -0.00709334193132355;
    $param_loop_terminal_H{CA}{AT} = -2.60;
    $param_loop_terminal_S{CA}{AT} = -0.00709334193132355;
    
    $param_loop_terminal_H{TA}{AG} = -2.40;
    $param_loop_terminal_S{TA}{AG} = -0.00612606803159761;
    $param_loop_terminal_H{GA}{AT} = -2.40;
    $param_loop_terminal_S{GA}{AT} = -0.00612606803159761;
    
    $param_loop_terminal_H{TA}{AT} = 0.00;
    $param_loop_terminal_S{TA}{AT} = 0.00;
    $param_loop_terminal_H{TA}{AT} = 0.00;
    $param_loop_terminal_S{TA}{AT} = 0.00;
    
    $param_loop_terminal_H{TG}{AA} = 0.00;
    $param_loop_terminal_S{TG}{AA} = 0.0016121231662099;
    $param_loop_terminal_H{AA}{GT} = 0.00;
    $param_loop_terminal_S{AA}{GT} = 0.0016121231662099;
    
    $param_loop_terminal_H{TG}{AC} = 0.00;
    $param_loop_terminal_S{TG}{AC} = 0.000644849266483959;
    $param_loop_terminal_H{CA}{GT} = 0.00;
    $param_loop_terminal_S{CA}{GT} = 0.000644849266483959;
    
    $param_loop_terminal_H{TG}{AG} = 0.00;
    $param_loop_terminal_S{TG}{AG} = 0.0016121231662099;
    $param_loop_terminal_H{GA}{GT} = 0.00;
    $param_loop_terminal_S{GA}{GT} = 0.0016121231662099;
    
    $param_loop_terminal_H{TG}{AT} = -2.30;
    $param_loop_terminal_S{TG}{AT} = -0.00580364339835563;
    $param_loop_terminal_H{TA}{GT} = -2.30;
    $param_loop_terminal_S{TA}{GT} = -0.00580364339835563;
    
    $param_loop_terminal_H{TA}{CA} = -2.60;
    $param_loop_terminal_S{TA}{CA} = -0.00677091729808157;
    $param_loop_terminal_H{AC}{AT} = -2.60;
    $param_loop_terminal_S{AC}{AT} = -0.00677091729808157;
    
    $param_loop_terminal_H{TA}{CC} = -0.50;
    $param_loop_terminal_S{TA}{CC} = -0.000967273899725939;
    $param_loop_terminal_H{CC}{AT} = -0.50;
    $param_loop_terminal_S{CC}{AT} = -0.000967273899725939;
    
    $param_loop_terminal_H{TA}{CG} = 0.00;
    $param_loop_terminal_S{TA}{CG} = 0.00;
    $param_loop_terminal_H{GC}{AT} = 0.00;
    $param_loop_terminal_S{GC}{AT} = 0.00;
    
    $param_loop_terminal_H{TA}{CT} = -2.70;
    $param_loop_terminal_S{TA}{CT} = -0.00709334193132355;
    $param_loop_terminal_H{TC}{AT} = -2.70;
    $param_loop_terminal_S{TC}{AT} = -0.00709334193132355;
    
    $param_loop_terminal_H{TG}{CA} = 0.00;
    $param_loop_terminal_S{TG}{CA} = 0.000644849266483959;
    $param_loop_terminal_H{AC}{GT} = 0.00;
    $param_loop_terminal_S{AC}{GT} = 0.000644849266483959;
    
    $param_loop_terminal_H{TG}{CC} = 0.00;
    $param_loop_terminal_S{TG}{CC} = 0.000644849266483959;
    $param_loop_terminal_H{CC}{GT} = 0.00;
    $param_loop_terminal_S{CC}{GT} = 0.000644849266483959;
    
    $param_loop_terminal_H{TG}{CG} = -3.70;
    $param_loop_terminal_S{TG}{CG} = -0.00935031436401741;
    $param_loop_terminal_H{GC}{GT} = -3.70;
    $param_loop_terminal_S{GC}{GT} = -0.00935031436401741;
    
    $param_loop_terminal_H{TG}{CT} = 0.00;
    $param_loop_terminal_S{TG}{CT} = 0.000644849266483959;
    $param_loop_terminal_H{TC}{GT} = 0.00;
    $param_loop_terminal_S{TC}{GT} = 0.000644849266483959;
    
    $param_loop_terminal_H{TA}{GA} = -1.90;
    $param_loop_terminal_S{TA}{GA} = -0.00419152023214574;
    $param_loop_terminal_H{AG}{AT} = -1.90;
    $param_loop_terminal_S{AG}{AT} = -0.00419152023214574;
    
    $param_loop_terminal_H{TA}{GC} = 0.00;
    $param_loop_terminal_S{TA}{GC} = 0.00;
    $param_loop_terminal_H{CG}{AT} = 0.00;
    $param_loop_terminal_S{CG}{AT} = 0.00;
    
    $param_loop_terminal_H{TA}{GG} = -1.50;
    $param_loop_terminal_S{TA}{GG} = -0.00354667096566178;
    $param_loop_terminal_H{GG}{AT} = -1.50;
    $param_loop_terminal_S{GG}{AT} = -0.00354667096566178;
    
    $param_loop_terminal_H{TA}{GT} = -2.30;
    $param_loop_terminal_S{TA}{GT} = -0.00580364339835563;
    $param_loop_terminal_H{TG}{AT} = -2.30;
    $param_loop_terminal_S{TG}{AT} = -0.00580364339835563;
    
    $param_loop_terminal_H{TG}{GA} = 0.00;
    $param_loop_terminal_S{TG}{GA} = 0.0016121231662099;
    $param_loop_terminal_H{AG}{GT} = 0.00;
    $param_loop_terminal_S{AG}{GT} = 0.0016121231662099;
    
    $param_loop_terminal_H{TG}{GC} = -2.90;
    $param_loop_terminal_S{TG}{GC} = -0.00612606803159761;
    $param_loop_terminal_H{CG}{GT} = -2.90;
    $param_loop_terminal_S{CG}{GT} = -0.00612606803159761;
    
    $param_loop_terminal_H{TG}{GG} = 0.00;
    $param_loop_terminal_S{TG}{GG} = 0.0016121231662099;
    $param_loop_terminal_H{GG}{GT} = 0.00;
    $param_loop_terminal_S{GG}{GT} = 0.0016121231662099;
    
    $param_loop_terminal_H{TG}{GT} = -2.00;
    $param_loop_terminal_S{TG}{GT} = -0.0048363694986297;
    $param_loop_terminal_H{TG}{GT} = -2.00;
    $param_loop_terminal_S{TG}{GT} = -0.0048363694986297;
    
    $param_loop_terminal_H{TA}{TA} = 0.00;
    $param_loop_terminal_S{TA}{TA} = 0.00;
    $param_loop_terminal_H{AT}{AT} = 0.00;
    $param_loop_terminal_S{AT}{AT} = 0.00;
    
    $param_loop_terminal_H{TA}{TC} = -1.40;
    $param_loop_terminal_S{TA}{TC} = -0.00354667096566178;
    $param_loop_terminal_H{CT}{AT} = -1.40;
    $param_loop_terminal_S{CT}{AT} = -0.00354667096566178;
    
    $param_loop_terminal_H{TA}{TG} = -3.70;
    $param_loop_terminal_S{TA}{TG} = -0.00999516363050137;
    $param_loop_terminal_H{GT}{AT} = -3.70;
    $param_loop_terminal_S{GT}{AT} = -0.00999516363050137;
    
    $param_loop_terminal_H{TA}{TT} = -2.30;
    $param_loop_terminal_S{TA}{TT} = -0.00644849266483959;
    $param_loop_terminal_H{TT}{AT} = -2.30;
    $param_loop_terminal_S{TT}{AT} = -0.00644849266483959;
    
    $param_loop_terminal_H{TG}{TA} = -2.90;
    $param_loop_terminal_S{TG}{TA} = -0.00773819119780751;
    $param_loop_terminal_H{AT}{GT} = -2.90;
    $param_loop_terminal_S{AT}{GT} = -0.00773819119780751;
    
    $param_loop_terminal_H{TG}{TC} = 0.00;
    $param_loop_terminal_S{TG}{TC} = 0.000644849266483959;
    $param_loop_terminal_H{CT}{GT} = 0.00;
    $param_loop_terminal_S{CT}{GT} = 0.000644849266483959;
    
    $param_loop_terminal_H{TG}{TG} = -2.00;
    $param_loop_terminal_S{TG}{TG} = -0.0048363694986297;
    $param_loop_terminal_H{GT}{GT} = -2.00;
    $param_loop_terminal_S{GT}{GT} = -0.0048363694986297;
    
    $param_loop_terminal_H{TG}{TT} = 0.00;
    $param_loop_terminal_S{TG}{TT} = 0.000644849266483959;
    $param_loop_terminal_H{TT}{GT} = 0.00;
    $param_loop_terminal_S{TT}{GT} = 0.000644849266483959;

    #############################################################

    my $NC_R                 = 1.9872e-3;  # Kcal/(Mol . K) -- McQuarrie Stat. Mech.
    my $MAX_SEQUENCE_LENGTH  = 256;
    my $MAX_LOOP_LENGTH      = $MAX_SEQUENCE_LENGTH/2;
    my $MAX_BULGE_LENGTH     = $MAX_SEQUENCE_LENGTH/2;
    my $MAX_HAIRPIN_LENGTH   = $MAX_SEQUENCE_LENGTH/2;

    ##/////////////////////////////////////////////////////
    ## From Table 4 in the above reference:
    ##/////////////////////////////////////////////////////
    ## Internal Loop
    $param_loop_S{0} = 0.0;
    $param_loop_S{1} = 0.0;
    $param_loop_S{2} = 0.0;
    $param_loop_S{3} = _ENTROPY(3.2, 0.0);
    $param_loop_S{4} = _ENTROPY(3.6, 0.0);
    $param_loop_S{5} = _ENTROPY(4.0, 0.0);
    $param_loop_S{6} = _ENTROPY(4.4, 0.0);
    $param_loop_S{7} = _ENTROPY(4.6, 0.0);
    $param_loop_S{8} = _ENTROPY(4.8, 0.0);
    $param_loop_S{9} = _ENTROPY(4.9, 0.0);
    $param_loop_S{10} = _ENTROPY(4.9, 0.0);
    
    $param_loop_S{12} = _ENTROPY(5.2, 0.0);
    
    $param_loop_S{14} = _ENTROPY(5.4, 0.0);
    
    $param_loop_S{16} = _ENTROPY(5.6, 0.0);
    
    $param_loop_S{18} = _ENTROPY(5.8, 0.0);
    
    $param_loop_S{20} = _ENTROPY(5.9, 0.0);
    
    $param_loop_S{25} = _ENTROPY(6.3, 0.0);
    
    $param_loop_S{30} = _ENTROPY(6.6, 0.0);
    
    # Jacobson-Stockmayer _ENTROPY extrapolation
    # dS(n) = dS(x) + 2.44*R*ln(n/x)
    
    my $extrapolation_x = 30;
    
    my $dS_loop_x = _ENTROPY(6.6, 0.0);
    
    $param_loop_S{11} = $dS_loop_x + 2.44*$NC_R*(-1.00);
    $param_loop_S{13} = $dS_loop_x + 2.44*$NC_R*(-8.36e-1);
    $param_loop_S{15} = $dS_loop_x + 2.44*$NC_R*(-6.93e-1);
    $param_loop_S{17} = $dS_loop_x + 2.44*$NC_R*(-5.68e-1);
    $param_loop_S{19} = $dS_loop_x + 2.44*$NC_R*(-4.57e-1);
    $param_loop_S{21} = $dS_loop_x + 2.44*$NC_R*(-3.57e-1);
    $param_loop_S{22} = $dS_loop_x + 2.44*$NC_R*(-3.10e-1);
    $param_loop_S{23} = $dS_loop_x + 2.44*$NC_R*(-2.66e-1);
    $param_loop_S{24} = $dS_loop_x + 2.44*$NC_R*(-2.23e-1);
    $param_loop_S{26} = $dS_loop_x + 2.44*$NC_R*(-1.43e-1);
    $param_loop_S{27} = $dS_loop_x + 2.44*$NC_R*(-1.05e-1);
    $param_loop_S{28} = $dS_loop_x + 2.44*$NC_R*(-6.90e-2);
    $param_loop_S{29} = $dS_loop_x + 2.44*$NC_R*(-3.39e-2);

    for( my $i=31; $i<$MAX_LOOP_LENGTH; $i++ ){
        $param_loop_S{$i} = $dS_loop_x + 2.44*$NC_R*( log($i/$extrapolation_x) );
    }

    #############################################################
    ## Bulge
    $param_bulge_S{0}  = 0.0;
    $param_bulge_S{1}  = _ENTROPY(4.0, 0.0);
    $param_bulge_S{2}  = _ENTROPY(2.9, 0.0);
    $param_bulge_S{3}  = _ENTROPY(3.1, 0.0);
    $param_bulge_S{4}  = _ENTROPY(3.2, 0.0);
    $param_bulge_S{5}  = _ENTROPY(3.3, 0.0);
    $param_bulge_S{6}  = _ENTROPY(3.5, 0.0);
    $param_bulge_S{7}  = _ENTROPY(3.7, 0.0);
    $param_bulge_S{8}  = _ENTROPY(3.9, 0.0);
    $param_bulge_S{9}  = _ENTROPY(4.1, 0.0);
    $param_bulge_S{10} = _ENTROPY(4.3, 0.0);
    
    $param_bulge_S{12} = _ENTROPY(4.5, 0.0);
    
    $param_bulge_S{14} = _ENTROPY(4.8, 0.0);
    
    $param_bulge_S{16} = _ENTROPY(5.0, 0.0);
    
    $param_bulge_S{18} = _ENTROPY(5.2, 0.0);
    
    $param_bulge_S{20} = _ENTROPY(5.3, 0.0);
    
    $param_bulge_S{25} = _ENTROPY(5.6, 0.0);
    
    $param_bulge_S{30} = _ENTROPY(5.9, 0.0);
    
    ## Jacobson-Stockmayer _ENTROPY extrapolation
    ## dS(n) = dS(x) + 2.44*R*ln(n/x)
    ##
    $extrapolation_x = 30;
    my $dS_bulge_x      = _ENTROPY(5.9, 0.0);
    
    $param_bulge_S{11} = $dS_bulge_x + 2.44*$NC_R*(-1.00);
    $param_bulge_S{13} = $dS_bulge_x + 2.44*$NC_R*(-8.36e-1);
    $param_bulge_S{15} = $dS_bulge_x + 2.44*$NC_R*(-6.93e-1);
    $param_bulge_S{17} = $dS_bulge_x + 2.44*$NC_R*(-5.68e-1);
    $param_bulge_S{19} = $dS_bulge_x + 2.44*$NC_R*(-4.57e-1);
    $param_bulge_S{21} = $dS_bulge_x + 2.44*$NC_R*(-3.57e-1);
    $param_bulge_S{22} = $dS_bulge_x + 2.44*$NC_R*(-3.10e-1);
    $param_bulge_S{23} = $dS_bulge_x + 2.44*$NC_R*(-2.66e-1);
    $param_bulge_S{24} = $dS_bulge_x + 2.44*$NC_R*(-2.23e-1);
    $param_bulge_S{26} = $dS_bulge_x + 2.44*$NC_R*(-1.43e-1);
    $param_bulge_S{27} = $dS_bulge_x + 2.44*$NC_R*(-1.05e-1);
    $param_bulge_S{28} = $dS_bulge_x + 2.44*$NC_R*(-6.90e-2);
    $param_bulge_S{29} = $dS_bulge_x + 2.44*$NC_R*(-3.39e-2);
    
    for( my $i = 31; $i < $MAX_BULGE_LENGTH; $i++){
        $param_bulge_S{i} = $dS_bulge_x + 2.44*$NC_R*( log( $i/$extrapolation_x) );
    }

    return ( \%param_H, \%param_S );
}

sub _ENTROPY {
	my ( $_dG, $_dH) = @_;
	return -($_dG-$_dH)*1000/(273.15+37);
}

sub _dimerCheck {
	my ( $sequence, $target ) = @_;
	
	$sequence = uc $sequence;      # change sequence to upper case
	$sequence = uc $sequence;      # change sequence to upper case
}

sub _dimerCheck_func {
	my ( $sequence, $target ) = @_;
}

1;
