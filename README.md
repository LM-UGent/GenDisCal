# GenDisCal

GenDisCal is a k-mer based Genome Distance Calculator. It works based on a two-step process.
First a signature is computed for a genome, then signatures for all input files are compared
using a distance method.

## Installation

To install GenDisCal on a GNU-Linux system, simply download the source code (with git clone
or from here, as a compressed file) and run `make`. This should produce an obj directory and
a bin directory. The runnable application will be located in the bin directory. GenDisCal
can also be compiled on windows by opening and building "GenDisCal.vcxproj" with Microsoft
Visual Studio, though we recommend using the pre-built exectuable instead.


Alternatively, you may use the "release" libraries, which also countains the graphical user
interface (the ".jar" file). In order to use the GUI, please download the binary that is
suitable for your operating system as well as the ".jar" file. You will need to link the
binary before the GUI can produce results. This can be done through `File > Set GenDisCal Path`.


## Usage (GUI)

Before you can use the GUI to run GenDisCal, you will need to indicate the location of the
executable (GenDisCal_vX_X, GenDisCal_vX_X.exe, or GenDisCal_st_vX_X.exe) from the `File`menu.

You can add genomes in the fasta format using the `Add Files` button. The full path to each
file added in this manner will be displayed in the box below. The `Import File List` and 
`Export File List` simplify this process by reading/writing the path of all the files from/to
another file. It is also possible to files containing signatures, however please make sure that
you then select the approppriate `Basis`. (see next paragraph)

The GUI offers a set of preset parameters, which can be selected with the appropriate button in
the top-centre of the window. Alternatively, it is possible to manually select the signature
type (`Basis`) and distance calculation methods (`Method`) below.

GenDisCal can be run in of two modes: `All Against All` and `Query Against All`. The first
mode will compare all the genomes/signatures against each other, whereas the second mode
will compare one genome/signature (specified using the text box below) against all genomes/signatures
added using `Add Files` and/or `Import File List`.

If the names of the added files are not explicit, it possible to specify their taxonomic lineage
or add an alias for each file. These are comma-separated files with 2 columns - the first containing
the name of the file to be compared (or just the GCF number in case of genomes downloaded from RefSeq),
and the second containing either the lineage (for the taxonomy file) or text used to more easily
identify a genome (for the alias file). There are three possible formats for taxonomy: `PCOF_S`, `GTDB`,
and `NCBI`. Further documentation is in the `--help` of the commandline version.

The GUI also enables to filter the added files (list is particularly useful when these are used as a
reference across multiple analyses). First, only `Include`'d names (meaning either file names, 
aliases, or part of a genome's taxonomic lineage) are kept from the input list, then `Exclude`'d names
are removed from that list.

The final step before running the program is to select an output file and the output format.
The output formats `Comparison List`, `Distance Matrix`, and `Histogram` , should be somewhat
self-explanatory. If `Signature Database` or `Readable Signatures` is selected, then the comparison
of signatures will not be perfomed, and a new file containing the signatures of all input files 
(Added files AND Query, if that `Query against All` is selected) will be generated. The file
generated when using `Signature Database` will then be usable as any regular input file, whereas
the file generated using `Reable Signatures` is intended for informative purposes.


## Usage (command line)

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
