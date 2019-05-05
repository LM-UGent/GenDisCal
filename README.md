# GenDisCal

GenDisCal is a k-mer based Genome Distance Calculator. It works based on a two-step process.
First a signature is computed for a genome, then signatures for all input files are compared
using a distance method.

## Installation

To install GenDisCal on a GNU-Linux system, simply download the source code (with git clone
or from here, as a compressed file) and run `make`. This should produce an obj directory and
a bin directory. The runnable application will be located in the bin directory. GenDisCal
can be compiled on windows, but no makefile is provided for this purpose at the moment.

## Usage

GenDisCal treats unnamed arguments as input files and outputs to stdout by default. The output
consists of comma-separated values (.csv).

    GenDisCal file1.fasta file2.fasta

Alternatively, if the number of files is large, they may be listed in a file. Such files are
specified using the `--filelist` option. Only one such file is permitted, but additional
single files may still be added as unnamed arguments.

    GenDisCal --filelist manyfiles.list extrafile.fasta

It is possible to specify the output file using the `--output` option. "stderr" and "stdout"
are interpreted as the standard output streams of the same name.

    GenDisCal file1.fasta file2.fasta --output output.csv

The output format is a list of comparisons, but this may be changed by adding one the
`--distancematrix` or `--histogram` option (mutually exclusive).

    GenDisCal --filelist manyfiles.list --distancematrix
    GenDisCal --filelist manyfiles.list --histogram
    GenDisCal --filelist manyfiles.list --histogram <bin width>

You may want to limit the output to a specific range of values. The `--below` and `--above`
can be used for this purpose

    GenDisCal --filelist manyfiles.list --above 0.20 --below 0.40

It is possible to perform one-against-all comparisons using the `--search` option. If this
option is absent, GenDisCal will perform all-against-all comparions, and if it is present but
does not have an argument, the last genome to have been added is compared to the rest (note:
genomes specified by the `--filelist` argument are added first, in the specified order)

    GenDisCal --filelist reference_genomes.list -s query.fasta

By default, distance computations are done with the PaSiT4 method. This can be changed using
the `--basis`, `--method`, and `--preset` options.

    GenDisCal file1.fasta file2.fasta --basis freq 4 --method manhattan
    GenDisCal file1.fasta file2.fasta --preset PaSiT6

Additional help is available through

    GenDisCal --help
    GenDisCal --help <command>
