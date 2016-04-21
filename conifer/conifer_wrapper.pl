#!/usr/bin/env perl

# Execute CoNIFER plotcalls and
# returns a HTML page with links to PNG plots

use strict;
use warnings;
use Getopt::Long;
use File::Basename; 
use File::Path qw(make_path remove_tree);

my $command;
my $dir=$ENV{'CONIFER_PATH'};

our ($multiple, $input, $regions, $sample, $window, $html_file, $html_folder, $verbose);

GetOptions('multiple'=>\$multiple, 'input=s'=>\$input, 'regions=s'=>\$regions,
    'sample:s'=>\$sample, 'window:i'=>\$window, 'verbose'=>\$verbose,
    'html_file=s'=>\$html_file, 'html_folder=s'=>\$html_folder);

make_path($html_folder);

# Build command
if ($multiple){
    # Reformat file with regions as required by CoNIFER plotcalls
    system("awk '{print \$5,\$1,\$2,\$3,\$4}' OFS=\"\t\" $regions > regions_sorted");
    
    $command = "python ".$dir."/conifer.py plotcalls --input $input --calls regions_sorted --window $window --outputdir $html_folder 2>&1";
}else{
    my $sample_command = ($sample eq "") ? "" : "--sample $sample";
    my $plot_name = $regions;
    $plot_name =~ s/[:-]/_/g;
    $command = "python ".$dir."/conifer.py plot --input $input --region $regions $sample_command --output $html_folder/$plot_name.png 2>&1";
}

# Run CoNIFER
system($command);
$verbose and print $command,"\n";

# Write HTML file
open(HTML, ">$html_file");
print HTML "<html><head><title>CoNIFER: Copy Number Analysis for Targeted Resequencing</title></head><body><h3>CoNIFER Output Files:</h3><p><ul>\n";
opendir(DIR, $html_folder);

my @FILES= grep { /png$/ }  readdir(DIR); 
closedir(DIR);
foreach my $file (@FILES) {
    print HTML "<li><a href=$file>$file</a><img src=\"$file\" height=\"50\" width=\"100\"></li>\n";
}
print HTML "</ul></p></body></html>\n";
close(HTML);
