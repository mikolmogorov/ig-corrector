# DESCRIPTION

Llama-fixer is a tool for immunoglobulines NGS data processing.
It provides both duplicates removal and error correction.


# DATA REQUEREMENTS

Reads should be obtained by 454 amplicon sequencing protocol.
They are expected to fully cover amplicons (immunoglobuline regions).
Also, piplne relies on CDR3 region extracion, so some markers for this
extracion should be provided (see below).

# USAGE

blah blah blah

# INSTALL

make

# DEPENDENCIES

* muscle (should be in PATH directory)
* Biopython
* Cython
* NumPy

# USAGE

./llama-fixer.py [-c config_file] -o output_dir [-q] reads_file

# THIRD PARTY

* align python module: https://github.com/FredrikAppelros/align