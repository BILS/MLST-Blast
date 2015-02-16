#!/usr/bin/perl 
#===============================================================================
=pod

=head1 FILE

mlst-blast.pl

=head1 USAGE

mlst-blast.pl -q query.fas [options]

=head1 DESCRIPTION

Perform (nucleotide) Blast against allele-type database, and report ST-type(s).

=head1 OPTIONS

=over 8

=item B<-q, --query>=F<query.fas>

Use (DNA) fasta file F<query.fas> as query.

=item B<-d, --db>=I<blastdb>

Name of blast database. Defaults to I<blastdb>.

=item B<-p, --profiles>=F<PROFILES.txt>

Name of special ST-profiles file. Needed if reporting ST-types.

=item B<-c, --createdb>=F<blastdb.fas>

Create new blast data base using fasta file F<blastdb.fas>.

=item B<-l, --length>=I<n>

Min length of matching sequence alignment (e.g. 200).

=item B<-s, --similarity>=I<n>

Min matching percentage (e.g. 95) in blast.

=item B<-r, --report>

Parse blast results and report ST-types (default). Use B<--noreport> to do blast only.

=item B<-k, --keepblastfile>

Print full blast report to file.

=item B<-h, --help>

Print help message.

=item B<-v, --verbose>

Be verbose (or not, B<--noverbose>).

=back

=head1 EXAMPLES

=over 8

=item mlst-blast.pl -q query.fas

=item mlst-blast.pl -c mlst-data/BLASTDB/blastdb.fas

=item mlst-blast.pl -q query.fas -d myfolder/mydb -p myfolder/myprofiles.txt

=back

=head1 REQUIREMENTS

Programs I<blastn> and I<makeblastdb> (L<ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST>) needs to be installed in the PATH.

Blast data base, preferrably created from F<mlst-data/BLASTDB/blastdb.fas> (see below).

Special file F<PROFILES.txt> with ST-type definitions. Preferrably taken from F<mlst-data/PROFILES/PROFILES.txt>


=head1 BUGS

---

=head1 NOTES

Will take nucleotide sequences in a fasta formatted file and 
perform blastn against default blast database I<blastdb>.
The usual blast database files (ending in .nhr, .nin, .nsq)
needs to be present.

The script saves the blast report to F<queryfile.blast.out.tab>, 
and, optionally, a more detailed blast-report to
F<queryfile.blast.out.raw>.

ST-types for queries are given in file queryfile.mlst-blast.out.

The script uses flat files to keep information, and those files 
can easily be created using the companion script
F<download-mlst-data.pl>. If supplied manually, they need to
conform to a strict format.

=over 8

=item Format of F<query.fas>

In order to allow iterating over any number 
of query strains, the queries needs to be named systematically,
starting with a common strain identifier followed by underscore,
then any unique identifier. For example: ">B250_001",
">B250_002". The script will issue a warning if several species
are found for one strain.

=item Format of F<blastdb.fas>

In order to allow several species in the 
blast database, a special format needs to be used for the fasta
file to be used as blast database. The sequence identifier
needs to start with the species, followed by a pipe symbol, then 
gene name, followed by strain number. For example,
">achromobacter|nusA_1", ">achromobacter|nusA_2", etc.

=item Format of F<PROFILES.txt>

In order to do successful parsing and 
genotyping of the blast results, the special file F<PROFILES.txt>
needs to be present. The path to this file can be set using B<-p>.
The format of this file needs to follow a "standard" ST-profiles 
file, except that species names needs to be added to each entry.
The file is tab separated with ST type in the first column,
followed by locus name and locus type. For example (content below on three lines):

ecoli|ST_1 ecoli|adk_4 ecoli|fumC_2 ecoli|gyrB_2 ecoli|icd_4 ecoli|mdh_4 ecoli|purA_4 ecoli|recA_4

ecoli|ST_2 ecoli|adk_5 ecoli|fumC_3 ecoli|gyrB_2 ecoli|icd_6 ecoli|mdh_5 ecoli|purA_5 ecoli|recA_4

ecoli|ST_3 ...

=back

=head1 AUTHOR

Johan Nylander (JN), Johan.Nylander@bils.se

=head1 COMPANY

BILS/NRM

=head1 VERSION

1.0

=head1 CREATED

12/02/2014 02:13:45 PM

=head1 REVISION

---

=head1 TODO

Handle when No hits: report ST type for the loci that have hits, but give a warning about No hits for certain seqs? Adjust min_percent. How stringent? Perhaps use info from all-against-all on db. Adjust min_length. Set to length of query seq?

=head1 DOWNLOAD

https://github.com/nylander/MLST-Blast

=cut

#===============================================================================

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;

exec('perldoc', $0) unless (@ARGV);

## Global variables
my $data_dir          = 'mlst-data';                    # Default data directory
my $db_dir            = $data_dir . '/BLASTDB/';        # Default blastdb directory
my $dbname            = 'blastdb';                      # Default blast database
my $dbname_def        = $db_dir . $dbname;              # Default blast database with path
my $min_percent       = 95;                             # Default similarity cut-off for blast
my $min_length        = 200;                            # Default cut-off for for blast 
my $blastoutfile      = q{};
my $blastoutfile_suf  = '.blast.out.tab';               # Default suffix for blast output file
my $blastrawfile      = q{};
my $blastrawfile_suf  = '.blast.out.raw';               # Default suffix for blast raw report file
my $profiles          = q{};
my $profiles_dir      = $data_dir . '/PROFILES/';       # Default profiles directory
my $profiles_file     = 'PROFILES.txt';                 # Default special ST-profiles file
my $profiles_file_def = $profiles_dir . $profiles_file; # Default special ST-profiles file with path
my $outfile           = q{};
my $outfile_suf       = '.mlst-blast.out';              # Default suffix for ST-typing out file
my $outformat         = 7;                              # Tabular with comment lines
my $keepblastfile     = 0;
my $report            = 1;
my $VERBOSE           = 1;
my $createdb          = 0;
my $query             = q{};
my $db                = q{};
my $length            = q{};
my $similarity        = q{};
my %profiles_hash     = ();
my %HoH               = ();
my $blastn            = find_prog('blastn');
my $makeblastdb       = find_prog('makeblastdb');

## Get options
my $r = GetOptions(
    "createdb=s"     => \$createdb,
    "db=s"           => \$db,
    "keepblastfile!" => \$keepblastfile,
    "length=i"       => \$length,
    "outfile=s"      => \$outfile,
    "profiles=s"     => \$profiles,
    "query=s"        => \$query,
    "report!"        => \$report,
    "similarity=i"   => \$similarity,
    "verbose!"       => \$VERBOSE,
    "help"           => sub { exec("perldoc", $0); exit(0); },
);
if (!$query) {
    die "Error! Need a query (fasta) file.\n" unless ($createdb);
}
$blastrawfile = $query . $blastrawfile_suf;
if (!$blastoutfile) {
    $blastoutfile = $query . $blastoutfile_suf;
}
if (!$outfile) {
    $outfile = $query . $outfile_suf;
}
if ($length) {
    $min_length = $length;
}
if ($similarity) {
    $min_percent = $similarity;
}

## Set up BLAST db
if ($createdb) { # create data base with fasta file as argument
    my $file = $createdb;
     print Dumper($file);warn "\nfile (hit return to continue)\n" and getc(); 
    if (-e $file) {
        my $name = fileparse($file, qr/\.[^.]*$/);
        my $dir = dirname($file);
        my $out = $dir . '/' . $name;
        print STDERR "Creating blast database from file $file\n" if ($VERBOSE);
        system("$makeblastdb -dbtype nucl -in $file -title $name -out $out"); # Need to put the path to destination in -out 
        $dbname = $out;
    }
    else {
        die "Error! Could not find file $file: $!\n";
    }
}
elsif ($db) { # specify name of data base, with or without path
    $dbname = $db;
}
else { # If no data base provided, check cwd for blastdb.nsq, then mlst-data/BLASTDB/blastdb.nsq
    my $nsq = $dbname . '.nsq';
    if (-e $nsq) {
        # Just checking cwd for nsq file
    }
    else {
        $dbname = $dbname_def;
    }
}

## If still no query file, stop here
if (! $query) {
    print STDERR "\nEnd of script.\n" if ($VERBOSE);
    exit;
}

## Check for nsq file. $dbname needs a path
my $nsq = $dbname . '.nsq'; # Case sensitive
die "Error! cannot find file $nsq associated with $dbname.\nDo you need to create blast database?\n" unless (-f $nsq);
print STDERR "Using blast database $dbname\n" if ($VERBOSE);

## Perform blastn
print STDERR "Performing blastn using $query as query file\n" if ($VERBOSE);
system("$blastn -db $dbname -query $query -outfmt $outformat -max_target_seqs 1 -out $blastrawfile");

## Parse blast output.
die "Error! Cannot parse blast result.\n" unless (-s $blastrawfile);
open my $BLAST, "<", $blastrawfile or die "Error! Could not open file $blastrawfile: $!\n"; 
open my $BLASTOUTFILE, ">", $blastoutfile or die "Error! Could not write to $blastoutfile: $!\n";
print $BLASTOUTFILE "Query\tHit\tLength\tPercent_id\n";
{
    local $/ = "# Query: ";
    foreach my $chunk (<$BLAST>) {
        chomp($chunk);
        $chunk =~ s/\n/#/g;
        my ($query, $query_label);
        my ($species, $allele_allele_type);
        my ($allele, $allele_type);
        if ($chunk =~ /^(.+?)##/) { # JN: Watch this on exotic fasta header formats
            $query = $1;
            my (@splits) = split /_/, $query; # B250_purA
            $query_label = $splits[0];
        }
        my $hits_found = '';
        if ($chunk =~ /(\d+) hits found/) {
            $hits_found = $1;
        }
        if ($hits_found) {
            if ($chunk =~ /#($query.+)##/) {
                my $result = $1;
                my ($query_id, $subject_id, $identity, $alignment_length, $mismatches, $gap_opens, $q_start, $q_end, $s_start, $s_end, $evalue, $bit_score) = split/\t/, $result;
                ($species, $allele_allele_type) = split /\|/, $subject_id; # 'ecoli|fumC_40', or 'abaumannii|Oxf_cpn60_39'
                if ($allele_allele_type =~ /^(\S+?)_(\d+)$/) {
                    $allele = $1;
                    $allele_type = $2;
                }
                else {
                    die "Error! Could not parse allele and allele type from blast file.\n";
                }

                ## Filter hit on length and percent
                if (($alignment_length > $min_length) && ($identity >= $min_percent)) {
                    print $BLASTOUTFILE $query_id, "\t", $subject_id, "\t", $alignment_length, "\t", $identity, "\n";
                    $HoH{$query_label}->{$species}->{$allele}->{$allele_type} = $subject_id;
                }
                else {
                    print $BLASTOUTFILE $query_id, "\t", "No hits", "\t", "\t", "\n";
                }
            }
        }
        elsif ($hits_found eq 0) {
            print $BLASTOUTFILE $query, "\t", "No hits", "\t", "\t", "\n";
        }
    }
}

## If report, parse the special file PROFILES.txt and create a reduced representation based on the blast results
if ($report) {
    if ($profiles) {
        $profiles_file = $profiles;
    }
    elsif (-e $profiles_file) {
        # Just checking if default profile file exists in cwd
    }
    else {
        $profiles_file = $profiles_file_def;
    }
    open my $PF, "<", $profiles_file or die "Error! Could not open profile file $profiles_file: $!.\n";
    while (<$PF>) {
        my $line = $_;
        chomp($line);
        my ($st_line, @allele_strings) = split /\t/, $line;
        my ($sp, $st_id) = split /\|/, $st_line; # achromobacter|ST_1
        foreach my $query (keys %HoH) {
            if (exists($HoH{$query}->{$sp})) { # If species was hit for query
                my @alleles_array = ();
                foreach my $allele_string (@allele_strings) {
                    my $allele;
                    if ($allele_string =~ /\|(\S+?)_\d+$/) {
                        $allele = $1;
                    }
                    else {
                        die "Error! Could not parse allele and allele type from blast file.\n";
                    }
                    if(exists($HoH{$query}{$sp}{$allele})) { # If allele was hit for query 
                        push @alleles_array, $allele_string; 
                    }
                }
                ## Build key
                @alleles_array = sort(@alleles_array);
                my $key = join('', @alleles_array);
                ## The key might point to several ST-types due to the fact that we might not have all genes
                push @{$profiles_hash{$query}->{$key}}, $st_line;
            }
            else {
                next;
            }
        }
    }
    close($PF);
}

## Parse hits and Report
open my $OUTFILE, ">", $outfile or die "Error! Could not open outfile $outfile for writing:$!\n";
foreach my $query (sort keys %HoH) { # B250
    print STDERR "Parsing blast report for query: $query\n" if($VERBOSE);
    my @spp = keys %{$HoH{$query}}; # acronym(s)
    if (scalar(@spp) > 1) {
        print STDERR "# Warning: More than one species found (@spp) for query $query!\n" if ($VERBOSE);
    }
    else {
        ## Rebuild hit
        my @key_array = ();
        foreach my $acronym (sort keys %{$HoH{$query}}) {
            foreach my $locus (sort keys %{ $HoH{$query}->{$acronym} }) {
                my ($allele_type) = keys % {$HoH{$query}->{$acronym}->{$locus}};
                my $hit = $acronym . '|' . $locus . '_' . $allele_type;
                push(@key_array, $hit);
            }
        }
        my $key = join('', @key_array);

        print $OUTFILE "# Query,Species,ST-type,";
        ## See if key is found in profiles, and report
        if (exists($profiles_hash{$query}{$key})) {
            ## Check if more than one STtype
            my $nst = scalar(@{$profiles_hash{$query}{$key}});
            if ($nst > 1) {
                print STDERR "# Note: More than one ST-type for query $query!\n" if ($VERBOSE);
                print $OUTFILE "$1," while ($key =~ /\|(\S+?)_\d+/g); # allele
                foreach my $res (@{$profiles_hash{$query}{$key}}) {
                    print $OUTFILE "\n";
                    print $OUTFILE "$query,";
                    my ($sp, $string) = split /\|/, $res;
                    my ($rest, $st) = split /_/, $string; 
                    print $OUTFILE $sp, ",", $st, ",";
                    print $OUTFILE "$1," while ($key =~ /\|\S+?_(\d+)/g); # allele_type
                }
                print $OUTFILE "\n";
            }
            else {
                print $OUTFILE "$1," while ($key =~ /\|(\S+?)_\d+/g);
                print $OUTFILE "\n";
                print $OUTFILE "$query,";
                my $res = shift(@{$profiles_hash{$query}{$key}});
                my ($sp, $string) = split /\|/, $res;
                my ($rest, $st) = split /_/, $string; 
                print $OUTFILE $sp, ",", $st, ",";
                print $OUTFILE "$1," while ($key =~ /\|\S+?_(\d+)/g);
                print $OUTFILE "\n";
            }
        }
        else {
            ## Report this to out or not?
            print STDERR "# Warning. Key does not exist = $key\n" if ($VERBOSE);
            print $OUTFILE "$1," while ($key =~ /\|(\S+?)_\d+/g);
            print $OUTFILE "\n";
            print $OUTFILE "$query,";
            if ($key =~ /^(\S+?)\|/) {
                my $species = $1;
                print $OUTFILE "$species,";
            }
            print $OUTFILE "NA,"; # ST-type
            print $OUTFILE "$1," while ($key =~ /\|\S+?_(\d+)/g);
            print $OUTFILE "\n";
        }
    }
}

## Cleanup
close($BLASTOUTFILE);
if (-f $blastoutfile) {
    print STDERR "Printed summary blast results to file $blastoutfile\n" if ($VERBOSE);
}
else {
    die "Error! Could not write to $blastoutfile: $!\n";
}
if ($keepblastfile) {
    if (-f $blastrawfile) {
        print STDERR "Printed detailed blast results to file $blastrawfile\n" if ($VERBOSE);
    }
    else {
        die "Error! Could not save raw blast output $blastrawfile: $!.\n";
    }
}
else {
    unlink($blastrawfile);
}
close($OUTFILE);
if (-f $outfile) {
    print STDERR "Printed results of ST-typing to file $outfile\n" if ($VERBOSE);
}
else {
    die "Error! Could not write to $outfile: $!\n";
}

print STDERR "\nEnd of script.\n" if ($VERBOSE);

## Subroutines
sub find_prog {
    my $PROG = shift;
    my $prog = q{};
    if ( -x $PROG ) {
        $prog = $PROG;
    }
    else {
        FIND_PROG:
        foreach ( split( /:/, $ENV{PATH} ) ) {
            if (-x "$_/$PROG") {
                $prog = "$_/$PROG";
                last FIND_PROG;
            }
        }
    }
    if ($prog eq '') {
        die "Can not find executable file \'$PROG\' in the PATH.";
    }
    else {
        return $prog;
    }
} 

__END__
# Query,Species,ST-type,adk,fumC,gyrB,icd,mdh,purA,recA,
Apa,ecoli,131,53,40,47,13,36,28,29,
