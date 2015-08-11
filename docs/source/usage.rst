================
Usage convention
================


------
Naming
------

For a list of tools, type::

  jsa

This will output the version of Japsa installed, the java version was used to
compile the package, together with a list of tools and a brief description for
each tool.


As you probably noticed from looking at the list of tools, every tool's name
starts with *jsa.*, followed by group (*e.g.*, hts) and the specific tool
function. In a shell that allows auto-completion, one can hit the tab key to see
the list of tools.


-------------
General usage
-------------

For the usage of a tool, type the tool name followed by --help. E.g.,::

   jsa.seq.sort --help

which will print out::

   Sort sequences based on their lengths

   Usage: jsa.seq.sort [options]
   Options:
     --input=s       Name of the input file, - for standard input
                     (REQUIRED)
     --output=s      Name of the output file, - for standard output
                     (REQUIRED)
     --alphabet=s    Alphabet of the input file. Options: DNA (DNA=DNA16), DNA4
                     (ACGT), DNA5(ACGTN), DNA16 and Protein
                     (default='DNA')
     --number        Add the order number to the beginning of contig name
                     (default='false')
     --reverse       Reverse sort order
                     (default='false')
     --sortKey=s     Sort key
                     (default='length')
     --help          Display this usage and exit
                     (default='false')

To specify an option, one can a single dash (-) or a double dash (--). One can
even shorten the option to a non-ambiguous prefix. For example::
   
    jsa.seq.sort -i=input.fas -o output.fasta --a DNA --sortKey=length --reverse -n true
     
In the above example, the option *input* is shorten to -i, *output* to -a. One can
have an equal side (-i=input.fas) or without (-o output.fas). For boolean
option, by presence of the option specifies a true value (e.g., --reverse), or
one can change it with true/false.

If there are more than one options sharing the same prefix, for example *input*
and *infile*, one has to use longer prefixes to identify the two (--inp and --inf).









