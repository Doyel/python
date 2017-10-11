use LWP::Simple;
use File::Spec;
use warnings;

# Download protein records corresponding to a list of PMC numbers.

$db = 'pubmed'; # Use 'pubmed' for getting bibliographic information for PubMed articles. Use 'pmc' for getting bibliographic information for PubMed Central articles
#$id_list = "";

$num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: Get_Biblio.pl <File containing PMC ids> <Dir for writing eSummary>\n";
    exit;
}

$PMID_filename = $ARGV[0];
#print "$PMID_filename\n";

$eSum_dir = $ARGV[1];

open($fh, '<:encoding(UTF-8)', $PMID_filename)
or die "Could not open file '$PMID_filename' $!";
while ($row = <$fh>) {
chomp $row;
#$id_list = $id_list . ',' . $row;
print "Processing: PMID$row\n";
$filename = "eSummary_Results_$row.xml";

#assemble the esummary URL
$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "esummary.fcgi?db=$db&id=$row";

#post the esummary URL
$docsums = get($url);
die "Couldn't get it!" unless defined $docsums;
#$code = getprint($url);
#print "$code\n";
$output_path = File::Spec->catfile($eSum_dir, $filename);
open($ofh, '>:encoding(UTF-8)', $output_path) or die "Could not open file '$output_path'";
print $ofh $docsums;
close $ofh;
}



