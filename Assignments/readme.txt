Submit "A4_Lastname_Name" .py or .ipynb in your branch or repo on GitHub.
Make sure all authors are mentioned in the comments.
    Make sure Biopython is available on your machine
        Check via "import Bio"
    Visit the "Biopython Tutorial and Cookbook"
    Choose 1 gene from the "Assignment 1 - Research and usage of biological DBs online"
    Programmatically, get nucleotide sequences for 5 different species that have this gene and:
        Save the basic info: accession number, title, organism, and length of the sequence
        Determine the GC percentage
        Perform a sequence alignment
        Create a phylogenetic tree with visualisation according to the scores
         of the sequence alignment via NJ (Neighbor-Joining) algorithm
        Using: Entrez, Bio.SeqIO, Bio.Align, Bio.Phylo
The code must automatically generate:
    A table with: accession number and title, organism, length of the sequence, GC percentage
    Phylogenetic tree visualisation with accession numbers and distances mentioned on the graph
In case of ValueError: Sequences must all be the same length see
https://stackoverflow.com/questions/32833230/biopython-alignio-valueerror-says-strings-must-be-same-length
