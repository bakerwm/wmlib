#!/usr/bin/env perl

# retrieve the intron sequences of input genes
# input: gene_name, (CSNK2A1)

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::SeqDumper;

sub geneid2slice {
    # return Bio::EnsEMBL::Gene
    # eg: CSNK2A1
    my $registry = shift(@_);
    my $gene_name = shift(@_);
    my $gene_adaptor = $registry->get_adaptor('Human', 'Core', 'Gene');
    my $gene = $gene_adaptor->fetch_by_display_label($gene_name);
    return $gene;
}

sub get_features {
    # input geneAdaptor, or transcript Adaptor
    # return feautre
    my $in = shift(@_);
    my $group = shift(@_);
    my @hits = ();
    if ($group eq "transcript") {
        push @hits, @{$in->get_all_Transcripts()}
    } elsif ($group eq "intron") {
        push @hits, @{$in->get_all_Introns()};
    } elsif ($group eq "exon") {
        push @hits, @{$in->get_all_Exons()};
    } else {
        die("unknown group: " + $group);
    }
    return @hits;
}

sub parse_id{
    my $fn = shift(@_);
    open my $FN, "<$fn" or die "cannot open file: $fn\n";
    my @ids = ();
    while(<$FN>) {
        chomp;
        next if(/(^\s*$)|(^\#)/);
        my $id = (split /\t/)[0];
        push @ids, $id;
    }
    close $FN;
    return @ids;
}

if (scalar(@ARGV) != 2) {
    die("Usage: perl $0 <id.list> <out.fa>\n");
}

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    '-host' => 'ensembldb.ensembl.org',
    '-user' => 'anonymous'
);

# my $slice_adaptor = $registry->get_adaptor('Human', 'Core', 'Slice');
# my $slice = $slice_adaptor->fetch_by_gene_stable_id('ENSG00000101266');
# my @transcripts = $slice->get_all_Transcripts();

my $dumper = Bio::EnsEMBL::Utils::SeqDumper->new();

my $flag_na = 1;
my $id_file = shift(@ARGV);
my $fa_out = shift(@ARGV);
my @gene_names = parse_id($id_file);
my @introns = ();
foreach my $gene_name (@gene_names) {
    my $gene = geneid2slice($registry, $gene_name);
    if (defined $gene) {
        push @introns, get_features($gene, "intron")
    } else {
        print $flag_na, ". ", $gene_name, " - unknown gene_name\n";
        $flag_na ++;
    }
}

foreach my $intron (@introns) {
    $dumper->dump($intron->feature_Slice(), 'FASTA', $fa_out);
}


## EOF



# sub gene2tx {
#     # return reference to a list of Bio::EnsEMBL::Transcripts
#     my $gene = shift(@_);
#     my @tx_ids = ();    
#     push @tx_ids, @{$gene->get_all_Transcripts()};
#     return @tx_ids;
# }

# sub gene2intron {
#     # return listref to Bio::EnsEMBL::Intron objects
#     my $gene = shift(@_);
#     my @introns = ();
#     push @introns, @{$gene->get_all_Introns()};
#     return @introns;
# }

# sub tx2intron {
#     # return listref to Bio::EnsEMBL::Intron objects
#     my $tx = shift(@_);
#     my @introns = ();
#     push @introns, @{$tx->get_all_Introns()};
#     return @introns;
# }