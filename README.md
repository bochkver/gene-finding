# gene-finding
Implementation of a hidden Markov model for gene finding and apply it to E. coli data.

We will use hidden Markov models to find genes in the sequence.
Also we will use the Viterbi algorithm to find parts of a sequence, which are likely a gene. As known, each gene starts with a start codon and terminates with a stop codon. The number of nucleotides in between must be divisible by 3. This fact can be used to create a hidden Markov model that searches for start codon, then accepts several internal codons terminated by a stop codon.
