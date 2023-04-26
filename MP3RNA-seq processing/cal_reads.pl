open I, "$ARGV[0]" or die;;
while(<I>)
	{
	chomp;
	@ss=split /\t/,$_;
	$hash{$ss[15]}++;
	}
open K, "$ARGV[1]";
while(<K>)
	{
	chomp;
	@kk=split;
	if(exists $hash{$kk[3]})
		{
		print "$kk[3]\t$hash{$kk[3]}\n";
		}
	else
		{
		print "$kk[3]\t0\n";
		}

	}
