# Mutate-protein
In this project, I want to check what is the probability of a truncated variant of the **dystrophin** protein arising from a **spontaneous mutation** in this gene.

Only **single nucleotide substitutions** in the exon region of the gene were studied here. I took the gene sequence and data on exon regions from the NCBI website.

You can start the experiment from the [dystrophin_mutation.py](https://github.com/AlenaSt97/Mutate-protein/blob/main/dystrophin_mutation.py) program, while all the functions are loaded from the [gene_classes.py](https://github.com/AlenaSt97/Mutate-protein/blob/main/gene_classes.py) file. First, I take a raw gene sequence, then I cut and splice exons from that sequence. After that, I translate the mRNA sequence into the protein sequence.

I **introduced mutations in two ways**: the first method assumes that the **probabilities of replacing one nucleotide with any other are the same**. In the second method, I take into account the fact that the **probabilities of transition from one nucleotide to another are not equal**, I took data on these probabilities from a textbook on biophysics.

As expected, the probabilities of the appearance of a truncated variant of the dystrophin protein were not equal: in the **first variant, 385 variants of the truncated protein were obtained per 10,000 mutation events** in the random position of the exon. In the **second case, 135 variants of the truncated protein per 10,000 events**.
