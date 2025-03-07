# EDGE Phylogenetics - Data and Methods

Phylogenetics data and code for the 2019 publication in Journal of Ecology:  https://doi.org/10.1111/1365-2745.13252

1) Gather all species and genes.. make sure that all species share at least one gene (see SGS_gene_availability_EDGE.xlsx). You may be required to search other names (outdated, synonyms - although GenBank is generally good with this).

Genes used:
- internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2
- ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit (rbcL) gene
- tRNA-Leu (trnL) gene
- trnL-trnF intergenic spacer
- maturase K (matK) gene
- photosystem II protein D1 (psbA) gene
- NADH dehydrogenase subunit F (ndhF) gene

2) Format genes.

- Make sure .fasta is correctly formatted and species names are the same
- Remove any species missing sequences
- Align genes using MAFFT : https://www.ebi.ac.uk/Tools/msa/. Can play around with MAFFT options but not really necessary. If there is one alignment that's particularly gap-y, especially in protein coding gene, check if the reverse complement is better: http://www.cellbiol.com/scripts/complement/dna_sequence_reverse_complement.php
- Save output as a new .fasta
- Use MEGA https://www.megasoftware.net/ to visualize your aligned sequences. Check for gaps of 3 in coding sequences, also use translate portion to make sure things make sense
- Use GBlocks to trim extra stuff off ends: http://molevol.cmima.csic.es/castresana/Gblocks_server.html. Do not trim on noncoding genes (ITSs, TrnL-TrnF spacer)
- Fill in any missing species with length of trimmed gene. Use dashes (easy to do in excel with =rept("-",701) for example)

3) combine all genes, use phyutility : https://code.google.com/archive/p/phyutility/

E.g. :

`java -jar phyutility.jar -concat -in ITS_sequences_SGS_aligned.fasta rbcL_sequences_SGS_aligned_trimmed.fasta TrnL_sequences_SGS_aligned_trimmed.fasta matK_sequences_SGS_aligned_trimmed.fasta TrnL-TrnFspacer_sequences_SGS_aligned.fasta psbA_sequences_SGS_aligned_trimmed.fasta ndhF_sequences_SGS_aligned_trimmed.fasta -out SGS_allgenes_nexus.fasta`

outputs a nexus file, may need to reformat for other applications.

4) Format as .fasta using custom python script (nexus.to.fasta.py). Make sure to edit directories in script as needed.

5) Use CIPRES and RAxML to develop best final tree. https://www.phylo.org/portal2/task!display.action?id=1680048

Once you get output from Cipres, look at the ‘best tree’ file (e.g., RAxML_bestTree.result) add this to the top of the file:

#NEXUS

begin trees;
trees tree1 = 

And this to the bottom:

end;

6) Process remainder using R (see script)