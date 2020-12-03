        my $mapping_file=$ARGV[0];;
    my %hash;
        my %unique;
        my $sample_num=0;
        my $ntc_sample_num=0;
        open (my $fh, $mapping_file) or die "Cannot read $mapping_file $!\n";
        my @header;
        while(<$fh>){
                chomp;
                @header = split /\t/,$_ if (/SampleID/);
                next if (/^#/);
                next if (/^\n/);
                my @array = split /\t/,$_;
		my $ntc=0;
                for my $i (1..$#array){
                        if ($header[$i] !~ /SampleID|Barcode|Linker|Description/){
                                $unique{$header[$i]}->{$array[$i]}++ if ( $array[$i] ne 'NTC' && $array[$i] !~ /no template/i );
				$ntc=1 if ( $array[$i] eq 'NTC' || $array[$i] =~ /no template/i );
                        }
                        $hash{$array[0]}->{$header[$i]}=$array[$i];
                }
                $sample_num++;
                $ntc_sample_num++ if (  $ntc );
        }
        close $fh;
        my @category_for_analysis;
        foreach my $feature (keys %unique){
                my $unique_feature_num = scalar (keys %{$unique{$feature}});
                push @category_for_analysis, $feature if ($unique_feature_num > 1 && $unique_feature_num < ($sample_num - $ntc_sample_num));
		print $feature,"\t",$unique_feature_num,"\t",$ntc_sample_num,"\t",$sample_num,"\n";
        }
	print join(",",@category_for_analysis),"\n";

