open I, "$ARGV[0]" or die;
while(<I>)
	{
	chomp;
	@ss=split;
	$file2=<I>;
	$file3=<I>;
	$file4=<I>;
	chomp($file2);
	chomp($file3);
	chomp($file4);
	if($file2=~/[ATCG](A{5,})/)
		{
		@mm=split /$1/,$file2;
		next if(length($mm[0])<10);
$len=length($mm[0]);
		$file4=~/(^.{$len})/;
		print "$_\n$mm[0]\n$file3\n$1\n";
		}
	else
		{
		print "$_\n$file2\n$file3\n$file4\n";
		}
	}
