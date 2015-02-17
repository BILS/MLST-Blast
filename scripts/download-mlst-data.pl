#!/usr/bin/perl 
#===============================================================================
=pod

=head1 FILE

download-mlst-data.pl

=head1 USAGE

=over 8

=item download-mlst-data.pl [options]

=back

=head1 DESCRIPTION

Download the whole, or parts of, the data base at L<http://pubmlst.org>.

The script will create a new folder in the current working directory named 
I<mlst-data>. The data folder will also contain a I<README.txt> file with 
detailed information of the folder content.

Importantly, the I<mlst-data> folder also contains a fasta file 
(I<mlst-data/BLASTDB/blastdb.fas>) formatted to be used as a blast data base, 
and a file (I<mlst-data/PROFILES/PROFILES.txt>) with ST-profile definitions.

=head1 OPTIONS

=over 8

=item B<-u, --url>=L<http://pubmlst.org/data/dbases.xml>

I<URL> to XML file with information for file paths.
               
=item B<-d, --datafolder>=I<name>

Use I<name> as data folder instead of the default (I<mlst-data>).

=item B<-s, --species>=I<name,name>

Specify which species (I<name>, one or more) to download as a comma separated list.

=item B<-q, --query>

Query URL about available species and genes. Prints a I<REPORT.txt> file.

=item B<-f, --fixfasta>, B<-nofixfasta>

Replace dashes to underscores in downloaded fasta files.

=item B<-v, --verbose>, B<-noverbose>

Print some progress info.

=item B<-c, --createblast>, B<-nocreateblast>

Create Fasta file for BLAST.

=item B<-h, --help>

Print help message.

=back

=head1 EXAMPLES

=over

=item download-mlst-data.pl

=item download-mlst-data.pl -q

=item download-mlst-data.pl -s efaecalis,efaecium -d Enterococcus

=back

=head1 REQUIREMENTS

Uses perldoc and Perl modules XML::Simple, LWP::Simple

=head1 BUGS

---

=head1 NOTES

Access to the xml file at pubmlst.org (L<http://pubmlst.org/data/dbases.xml>) is crucial for functionality.

=head1 AUTHOR

Johan Nylander (JN), johan.nylander@bils.se

=head1 COMPANY

BILS/NRM

=head1 VERSION

2.0

=head1 CREATED

11/28/2014 12:50:20 PM

=head1 REVISION

---

=head1 TODO

Subset genes. Dynamically set cut off for BLAST in order to allow new allele types.

=head1 DOWNLOAD

https://github.com/nylander/MLST-Blast

=cut

#===============================================================================

use strict;
use warnings;
use Cwd 'abs_path';
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use LWP::Simple;
use List::Util qw(first);
use XML::Simple;

## Global parameters
my $url              = "http://pubmlst.org/data/dbases.xml"; # make sure this is correct
my $createblast      = 1;
my $profiles         = 1;
my $data_folder      = q{};
my @genes            = ();
my @species          = ();
my %sp_hash          = ();
my %get_hash         = ();
#my %gene_hash        = ();
my $data_folder_abs  = q{};
my $data_folder_base = 'mlst-data';
my $date             = Date();
my $help             = q{};
my $VERBOSE          = 1;
my $fixfasta         = 1;
my $query            = 0;
my $readme           = 'README.txt';
my $species_file     = 'SPECIES.txt';
my $blastdb_folder   = 'BLASTDB';
my $blastdb_file     = 'blastdb.fas';
my $profiles_folder  = 'PROFILES';
my $profiles_file    = 'PROFILES.txt';
my $query_file       = 'REPORT.txt';
my $report_all       = 1;

## Get options
my $r = GetOptions(
                   "createblast!" => \$createblast,
                   "datafolder=s" => \$data_folder,
                   "fixfasta!"    => \$fixfasta,
                   #"genes=s"      => \@genes,
                   "species=s"    => \@species,
                   "profiles!"    => \$profiles,
                   "query"        => \$query,
                   "url=s"        => \$url,
                   "verbose!"     => \$VERBOSE,
                   "help"         => sub { exec("perldoc", $0); exit(0); },
                  );

## Get XML file from url, keep content in memory
my $xmlfile = get($url);
die "Error! Couldn't get content of $url" unless defined $xmlfile;
my $xml = new XML::Simple;
my $data = $xml->XMLin($xmlfile);

## Get all species and gene data
foreach my $sp (@{$data->{'species'}}) {
    my $species_name = $sp->{'content'};
    chomp($species_name);
    $species_name =~ s/\// or /g; # species name might have '/'
    my $url_profile = $sp->{'mlst'}->{'database'}->{'profiles'}->{'url'};
    chomp($url_profile);
    my $acronym = basename($url_profile, ".txt");
    $sp_hash{$acronym}->{'species_name'} = $species_name;
    foreach my $x (@{$sp->{'mlst'}->{'database'}->{'loci'}->{'locus'}}) {
        my $gene_name =  $x->{'content'};
        chomp($gene_name);
        $gene_name = lc($gene_name);
        #push @{$gene_hash{$gene_name}}, $acronym;
        push @{$sp_hash{$acronym}->{'genes'}}, $gene_name;
    }
}

## If using subsets
if (@species) {
    my @s = split(/,/, join(',', @species));
    @species = sort(@s);
    foreach my $sp (@species) {
        if(exists($sp_hash{$sp})) {
            $get_hash{$sp} = 1;
        }
        else {
            die "Error! Species acronym \"$sp\" was not found.\n";
        }
    }
}
#if (@genes) {
#    my @g = split(/,/,join(',',@genes));
#    @genes = sort(@g);
#    ## Check against %gene_hash
#}

## Query and Report?
if ($query) {
    ## Create REPORT file in cwd
    open RF, ">", $query_file or die "Error! Could not open file $query_file: $!\n";
    if (@species) {
        foreach my $sp (@species) {
            if (exists($sp_hash{$sp})) {
                print RF "$sp,$sp_hash{$sp}->{'species_name'}";
                foreach my $g (@{$sp_hash{$sp}->{'genes'}}) {
                        print RF ",$g";
                    }
                print RF "\n";
            }
            else {
                die "Error! Species acronym \"$sp\" is not found.\n"; 
            }
        }
    }
    else {
        foreach my $sp (sort keys %sp_hash) {
            print RF "$sp,$sp_hash{$sp}->{'species_name'}";
                foreach my $g (sort @{$sp_hash{$sp}->{'genes'}}) {
                        print RF ",$g";
                    }
                print RF "\n";
        }
    }
    close(RF);
    if (-e $query_file) {
        print STDERR "Created $query_file with a list of available species and loci.\n" if ($VERBOSE);
    }
    else {
        die "Error! Could not create $query_file: $!\n";
    }
    exit;
}
else {
    ## Create data folder
    if (!$data_folder) {
        $data_folder = $data_folder_base;
    }
    if (-e $data_folder) {
        die "Covardly refuses to overwrite existing data folder \"$data_folder\". Quitting.\n";
    }
    else {
        mkdir($data_folder);
        if (-e $data_folder) {
            print STDERR "Created folder $data_folder.\n" if $VERBOSE;
        }
        else {
            die "Failed creating folder $data_folder: $!\n";
        }
    }
    $data_folder_abs = abs_path($data_folder);

    ## Create README file in data folder
    my $rf = $data_folder_abs . '/' . $readme;
    open README, ">", $rf or die "Error! Could not open file $rf for writing: $!\n";
    print README "Data files downloaded $date from $url.\n";
    print README "\nFasta files and profiles are saved in folder $data_folder.\n";
    print README "\nSpecies names and acronyms are given in $data_folder/$species_file.\n";

    ## Create BLASTDB folder and file
    if ($createblast) {
        my $bfolder = $data_folder_abs . '/' . $blastdb_folder;
        mkdir($bfolder);
        die "Failed creating folder $bfolder: $!\n" unless (-d $bfolder);
        my $bfile = $bfolder . '/' . $blastdb_file;
        open BF, ">", $bfile or die "Error! Could not create file $bfile: $!\n"; 
        print README "\nFasta file $blastdb_file for blast was created in $data_folder/$blastdb_folder.\n";
    }

    ## Create PROFILES folder and file
    if ($profiles) {
        my $pfolder = $data_folder_abs . '/' . $profiles_folder;
        mkdir($pfolder);
        die "Failed creating folder $pfolder: $!\n" unless (-d $pfolder);
        my $pfile = $pfolder . '/' . $profiles_file;
        open PF, ">", $pfile or die "Error! Could not create file $pfile: $!\n"; 
        print README "\nA $profiles_file with all ST-profiles was saved in $data_folder/$profiles_folder.\n";
    }

    ## Create SPECIES file in data folder
    my $sf = $data_folder_abs . '/' . $species_file;
    open SPECIES, ">", $sf or die "Error! Could not open file $rf for writing: $!\n";
}

## Get data for species
print STDERR "Retrieving species data...\n" if $VERBOSE;
my $dashes = 0;
foreach my $sp (@{$data->{'species'}}) {
    ## Get species name and acronym
    my $species_name = $sp->{'content'};
    chomp($species_name);
    $species_name =~ s/\// or /g; # species name might have '/'
    my $url_profile = $sp->{'mlst'}->{'database'}->{'profiles'}->{'url'};
    chomp($url_profile);
    my $acronym = basename($url_profile, ".txt");
    ## Check if acronym is among the species
    if (@species && !exists($get_hash{$acronym})) {
        next;
    }
    else {
        ## Get profile count to make sure to download or not (error on website 12/01/2014 02:25:36 PM: C. fetus no data, but fixed now.)
        my $nprofiles = $sp->{'mlst'}->{'database'}->{'profiles'}->{'count'};
        if ($nprofiles < 1 ) {
            print STDERR "Skipping species $species_name\n" and next;
        }

        ## Create species dir
        my $species_dir = $data_folder_abs . '/' . $acronym . '/';
        mkdir($species_dir);
        die "Error! Could not create folder $species_dir: $!\n" unless (-d $species_dir);
        print STDERR "\n# Species: $species_name [$acronym]\n" if $VERBOSE;

        ## Profile file
        my $profile_file = $species_dir . '/'. basename($url_profile);
        my $rc = getstore($url_profile, $profile_file); ## TODO: need to retrieve data without downloading if report?
        if (is_error($rc)) {
            die "Error! Could not save file $profile_file: $!\n";
        }
        print STDERR "  Profile: $url_profile.\n" if $VERBOSE;
        open my $P, "<", $profile_file or die "Error! Could not open file $profile_file: $!\n";
        my @pheader = ();
        my $clonal_index;
        my $search = 'clonal_complex';
        while (<$P>) {
            my $line = $_;
            chomp($line);
            if (/^ST/) {
                @pheader = split /\t/, $line;
                $clonal_index = first { $pheader[$_] eq $search } 0 .. $#pheader; # Some have clonal_complex, others don't, so skip 
            }
            else {
                my @pparts = split /\t/, $line;
                my $length = scalar(@pparts);
                if($clonal_index) {
                    $length = $clonal_index;
                }
                for (my $i = 0; $i < $length; $i = $i + 1) {
                    $pparts[$i] = $acronym . '|' . $pheader[$i] . '_' . $pparts[$i];
                    print PF $pparts[$i];
                    if ($i < $length-1) {
                        print PF "\t";
                    }
                    else {
                        print PF "\n";
                    }
                }
            }
        }
        close($P);

        ## Fasta files
        foreach my $x (@{$sp->{'mlst'}->{'database'}->{'loci'}->{'locus'}}) {
            my $gene_name =  $x->{'content'};
            chomp($gene_name);
            my $url_gene_fasta =  $x->{'url'};
            chomp($url_gene_fasta);
            my $gene_file = $species_dir . '/'. basename($url_gene_fasta);
            my $file_content = get($url_gene_fasta);
            my $dash = 0;
            if ($file_content =~ /(\S+)-(\d+)$/m) {
                $dash = 1;
                $dashes = 1;
            }
            if ($fixfasta) {
                $file_content =~ s/(\S+)-(\d+)$/$1_$2/mg; # Replace all '-' with '_'
                $file_content =~ s/>[\w|\|]+\n\s*\n//g;   # Remove empty sequences
                $file_content =~ s/>[\w|\|]+\n>/>/g;      # -"-
                $file_content =~ s/>[\w|\|]+\n*\Z//;      # -"-
            }
            if ($createblast) {
                $file_content =~ s/>(\S+)/>$acronym\|$1/mg;
                print BF $file_content;
            }
            open my $GF, ">", $gene_file or die "Error! Could not save file $gene_file: $!\n";
            print $GF $file_content;
            close($GF);
            if ($dash && $fixfasta) {
                print STDERR " *$gene_name: $url_gene_fasta\n" if $VERBOSE;
            }
            else {
                print STDERR "  $gene_name: $url_gene_fasta\n" if $VERBOSE;
            }
        }

        ## Print to SPECIES file
        print SPECIES $acronym, "\t", $species_name, "\n";
    }
}

## End of script
print STDERR "\nEnd of script.\n" if $VERBOSE;
if ($dashes && $fixfasta) {
    print STDERR "*Note that dashes in the fasta headers have been replaced with underscores for files marked above.\n"; 
    print README "\nNote that dashes in the fasta headers have been replaced with underscores for some files.\n"; 
}
close(SPECIES);
close(README);
if($createblast) {
    close(BF);
}
if ($profiles) {
    close(PF);
}

## Get date, e.g. "11-03-2014"
sub Date {
    @_ = localtime( shift || time );
    return ( sprintf( "%02d-%02d-%04d", $_[4] + 1, $_[3], $_[5] + 1900 ) );
}

__END__

