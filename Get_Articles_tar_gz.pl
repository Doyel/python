use LWP::Simple;
use warnings;
use Cwd;

# Download protein records corresponding to a list of PMC numbers.

$num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: Get_Articles_tar_gz.pl <File containing PMC ids> <xmlFTPfiles_dir>\n";
    exit;
}

$PMID_filename = $ARGV[0];
#print "$PMID_filename\n";

open($fh, '<:encoding(UTF-8)', $PMID_filename)
or die "Could not open file '$PMID_filename' $!";
$pwd = cwd();
$xmlFTP_dir = $ARGV[1];
chdir($xmlFTP_dir); 
while ($row = <$fh>) 
{
	chomp $row;

	#assemble the oa URL
	$base = 'http://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi';
	$url = $base . "?id=PMC$row";

	#post the oa URL
	$docsums = get($url);
    die "Couldn't get it!" unless defined $docsums;
    #$code = getprint($url);
    #print "$code\n";

    print "Processing: PMC$row\n";
	$filename = "File_Details_PMC" . $row . ".xml";
	open($fh1, '>:encoding(UTF-8)', $filename) or die "Could not open file '$filename'";
	print $fh1 $docsums;	
    close $fh1
}
chdir($pwd);


