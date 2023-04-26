open K, "$ARGV[0]" or die;
while(<K>)
        {
        chomp;
        @ss=split;
        $header=$ss[0];
	$header=~s/@//g;
        $file2=<K>;
        $file3=<K>;
        $file4=<K>;
        chomp($file2);
        chomp($file3);
        chomp($file4);
        $file2=~s/(^.{6})//g;
        $hash{$header}=$1;
	}
open I, "samtools view -h $ARGV[1]|" or die;
while(<I>)
	{
	chomp;
	@ss=split;
	$ID=$hash{$ss[0]};
	if(/^@/)
		{
		print "$_\n";
		}
	else
		{
		if(exists $hash_ID{"$ID\t$ss[2]\t$ss[3]"})
			{
			next;
			}
		else
			{
			$hash_ID{"$ID\t$ss[2]\t$ss[3]"}++;
			print "$_\n";
			}
		}
	}
