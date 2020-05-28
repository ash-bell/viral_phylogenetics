# Viral Phylogeny
---
Without a set of common genes shared between all viruses, it becomes difficult to deduce phylogeny of viruses. Other method of gene-sharing maps exist to combat this problem and are preferred. However, if you would like a classic tree, you can use this instead. This method uses a series of custom genes that the user selects to base a phylogenetic tree on. 

## Before you start:
---
*   [HMMER](http://hmmer.org/) installed `conda install -c bioconda hmmer` 
*   [MUSCLE](https://drive5.com/muscle/) installed `conda install -c bioconda muscle`
*   [raxml-ng](https://github.com/amkozlov/raxml-ng) installed if you want to make a phylogenetic tree with it `conda install -c bioconda raxml-ng`

## Files:
*   all_SAR11_phages.fna --> file of all my phage genomes

```
$ cat all_SAR11_phages.fna
>VIRSorter_AAAM_phage_000000000001
ATGAGTACAGATAGATTAGAAGTAAAACCAGATGACAATAACAATGAAACACTTGAACAGTCAGCAGAAAAATTAAAAGA
>VIRSorter_AAAO_prophage_000000000001
ATGAGTACAGATAGATTAGAAGTAAAACCAGATGACAATAACAATGAAACACTTGAACAGTCAGCAGAAAAATTAAAAGA
>VIRSorter_AACP_phage_000000000001
ATGAGTACAGATAGATTAGAAGTAAAACCAGATGACAATAACAATGAAACACTTGAACAGTCAGCAGAAAAATTAAAAGA
```

*   Gene called your organisms and BLAST them against a database to get an annotated gene function for each. Then select genes that are specific to your taxa of interest and find them within your annotations. For example, for a viral phlogenetics I chose capsid, DNA polymerase, exonuclease, endonuclease, helicase, head-tail-connector and RNA polymerase.

## Tutorial
---
1. Create a file for each catagory of genes you have chosen.

```
$ cat capsid_genes.fna
>VIRSorter_AACP_phage_000000000014_capsid
ATGAGTACAGATAGATTAGAAGTAAAACCAGATGACAATAACAATGAAACACTTGAACAGTCAGCAGAAAAATTAAAAGAACAAGGTGTTGATATTAATAAAGACCTTAGCGTCAATCCAAATGGAGAAGGAATTGAAGTTAGAGAACCGAAAGTTGAGACAGAAAGTACAGAAAAAAGACCAGAGTGGTTGCCAGAAAAATTTCAAAACGCAGAAGAACTGGCTAAAGCGTATGGGACTTTGGAAAAAGAATTTTCAGGTAGGACTAAAGAAGAAGTTAAACCTACTGAAGAAGTAAAAGCTGATGAAGCTCCACAAACAGGTTTAGATAAATACTATGAGGAGTTTGCTGATAAAGGAGAACTTGCAGAAAAGAGTTATTCAGAATTAGCTAAATTAGGTTTAGATAAAAACTTAGTTGATACTTATATTGAAGGACAAAAATTAGTTTCAGAAACAAACACTAAAGCAATTCAAGATATTGCAGGTGGTAAAGAAGAATATACTGAACTAGTTGAATGGGCAGGTAAGAACTTATCCGAAGCAGAAACTAAAGTCTTTAATGACATGGTTGATGGTGGAGATATTGAGACAGCTAAATTTGCTGTTCAAGGTCTTATGGCTAAATCAGGTGCAAATCCAAAACAACCTTCTTTATATGAGGGTACTAGCGATACAGTTTCTAAAGATGCTTTTGCTAGTGTTGCACAAGTTACAGAAGCAATGAATGACCCAAGATATGACAGCGACCCTGCATATAGACAACTAGTAGAAGACAAAATTGGGAGAAGCACTGTTCTT
>VIRSorter_AALP_phage_000000000004_capsid
ATGAGTACAGATAGAGTAGAAATAAAACCAGACGCAATAGAAAAATCTTTAGAACAGTCTGCAGAAGATTTAAAAACTAATGATGGAGTTGATGTCAGCAAGGATGTTGCAGTTAATAAAAGTGGCGAAGGAGCAAGAATTAACGAATTGTCACAAGAAGACTTACAACAAAGCACAGAAGATAAACCTGACTGGTTGCCTTCGAAGTTTAAAAATGCTGAAGAACTTGCAAAAGCATATGGTGAACTTGAAAAATCTTTTTCTTCTAGAAAACAAGAAGAAGCACCAAAAGAAACTTCAATACCAGAAGTTAGAAAACCTACAGAAGGTCAAGAAGCATTAGGTAAATTTTATGATGAGTATGCAGCTAATAGTTCTTTATCAGATAAATCTTATGAAGAATTAGCAACTAAACATGGTTTATCTAAAGAACTTGTAGATGGTTATATTGAAGGTCAAAAAGCTATTGGTGATAATCAAACTAAAGCCATACACGATTTAACTGGTGGCGGTGAGAAATATACCGAACTAATGGATTGGGCAGGTAAAAACTTATCTGAACAAGAACAACAAGCATACAATAGTATGGTTGACAGTGGGAATGTTGATGCTGCAAAACTTGCAGTTCAAGGTCTCATGTCTAAAGCAGGTGTTAATTATAATCCTAAACAACCAGAATTATTTGAAGGTGGAGACCAAATACCTAATGATAGTTTTGAAAGTGTATCTCAGGTTACTGAAGCAATGAATGACCCAAGATATGCAAAAGACCCTGCGTATAGAAAAAAAGTTACTGATAAGATTGCACGAAGTTCAGTAATC
[...]
```

2. Run `muscle` on the genes to generate a multiple sequence alignment.

```
$ muscle -in capsid_genes.fna -out capsid_genes.msa
```
It should look like this:
```
$ cat capsid_genes.msa
>VIRSorter_AALP_phage_000000000009_capsid
ATGTCAAACGCAAACACATCAAGACTGGGTCTAGTAAATGATACTGGCACATCTTATGAT
GCTTTGTTCCTAAAAGTATTCTCAGGGGAAGTATTAGCAAGTTTTGGTCGAGAGAACAAA
ATGCTTGATATGACTACTGTTAGAACTATAGGTTCAGGCAAAAGTGCGACTTTTCCTGTT
ACAGGAACAATCGCATCAAGTTACCACACAGCAGGAAATGAAATTCTTGGGACTGCGGTA
AATCACAACGAAAGAATAATCAATATAGACGATATGTTGATTTCACATGCGTTCATAGCA
GAAATAGATGAA-------------TTAAAAAATCATTT----------CGATGCTAGAA
GC---ATTTATAGTAATGAAATGGGCAGAGCATTAAGTAACAAGGTCGACCAACATTTAG
TTCAGTTAATGGTACTAGCTAGTCAAGCAAGTGCTACTATAACTGGCGGTAATGGTGGTT
ATGAAATTACTGATGCAGACGCAAAAACAAATGCAGAAAGTTTAATCGCTTCAGTATTTG
ACGCAAATCAGAAACTAGACGAAAATGATGT------------ACCAACATCCGACAGAT
ATGTTGTCGTTACACCAGATGTCTACTACCAATTATGTAATGTAGAC-AAATTAATT---
------TCTAGAGATTTCTCTAAAAATAATGGTGATTTCGGTTCAGGCACAGTTGTTTCA
ATCGCAGGAATACCAGTAGTTAAATCAAACTCTGCTAAACTTGC---TTATGATGACAAC
TCTGGAGCAATTTCTGGCACAAATAATACATACAACGTGGATGCACAGCACATAGCGGCT
ACTGTATTTCACAAAAGTGCAGTTGGAACAGTGAAACTAAAAGATTTG------------
--------GTTTTAGAAAACACTTATGACCCAAGA-AGACTAGGTCACTTGCTTACAGCA
AGACTAGCGCTTGGACATAATATTCTGAGACCAGAAGCAGCAGTTTCAATCAAACAAGCA
>VIRSorter_AAJR_prophage_000000000025_capsid
ATG---GTT---GATAAAATAGAGATTAATCAAGATGAAAATAATATT---TCACTTGA-
--------------------------GGAACAATCAAAC--------------AAACAAG
A--------------A-----ACAAACTCACAGACTACAC--------C-----------
--AGAAGC-------TCAAACTA---------AAG-------------AGACTTCTAGT-
-------------GAGAAACCTTCATGGCTTCCAGATAAATTTTCTA-----ATGCAGAA
GAATTAGCTAAAGCCTATAGTGAGCTTGAGAAGAAGTTTTCTTCACCTTCAATAGAAGAA
[...]
```

4. Build your Hidden Markov Model (HMM) for searching for your genes

```
$ hmmbuild capsid_genes.hmm capsid_genes.msa
```
It should look like this:
```
$ cat capsid_genes.hmm
HMMER3/f [3.3 | Nov 2019]
NAME  capsid
LENG  836
MAXL  1191
ALPH  DNA
RF    no
MM    no
CONS  yes
CS    no
MAP   yes
DATE  Thu May 28 10:29:29 2020
NSEQ  4
EFFN  4.000000
CKSUM 1368186433
STATS LOCAL MSV      -11.9181  0.69604
STATS LOCAL VITERBI  -13.4641  0.69604
STATS LOCAL FORWARD   -5.7680  0.69604
HMM          A        C        G        T   
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   1.07870  1.73134  1.49329  1.35371
          1.38629  1.38629  1.38629  1.38629
          0.03279  4.12713  4.12713  1.46634  0.26236  0.00000        *
      1   0.09846  3.58835  3.37856  3.44116      1 A - - -
          1.38629  1.38629  1.38629  1.38629
          0.03279  4.12713  4.12713  1.46634  0.26236  1.09861  0.40547
[...]
```

4b. Optionally you can check through your gene calls and see and genes that may have a capsid that is simular to the HMM you have created. It should look like this:

```
$ nhmmer capsid_genes.hmm all_gene_calls.fna

# nhmmer :: search a DNA model, alignment, or sequence against a DNA database
# HMMER 3.3 (Nov 2019); http://hmmer.org/
# Copyright (C) 2019 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query file:                      capsid.hmm
# target sequence database:        gene_calls.fna
# number of worker threads:        2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       capsid  [M=836]
Scores for complete hits:
    E-value  score  bias  Sequence                                                   start    end  Description
    ------- ------ -----  --------                                                   -----  -----  -----------
   1.1e-208  694.5  78.7  VIRSorter_AALP_phage_000000000004                              1    821 
   1.2e-201  671.2  81.7  VIRSorter_AACP_phage_000000000014                              1    800 
   1.9e-163  544.8  84.2  VIRSorter_AAJR_prophage_000000000025                           2    764 
   3.3e-139  464.6  51.0  VIRSorter_AALP_phage_000000000009                              7    928 
   9.6e-101  337.4  71.7  VIRSorter_SCGC_AAA795_D22_000000000027_phage_000000000044      2    751 
    6.5e-75  251.9  86.2  HTVC011P_000000000025                                          7    789 
    3.2e-74  249.7  86.2  Eyrgjafa_VP5_000000000030                                      1    779 
    3.2e-74  249.7  86.2  Gjalp_VP5_000000000048                                         1    779 
    7.4e-70  235.2  60.6  HTVC025P_000000000026                                        151    769 
    3.7e-44  150.2  48.3  HTVC109P_000000000032                                         74    743 
    1.5e-41  141.6  59.5  Eistla_VP4_000000000040                                       73    732 

[...]
```
4c. Or create a multi HMM model that looks for capsids, head-tail conntectors and integrase. It should look like this:

```
$ cat capsid_genes.hmm head_tail_genes.hmm integrase_genes.hmm > viral_genes.hmms
$ hmmpress viral_genes.hmms

Working...    done.
Pressed and indexed 3 HMMs (3 names).
Models pressed into binary file:   viral_genes.hmms.h3m
SSI index for binary model file:   viral_genes.hmms.h3i
Profiles (MSV part) pressed into:  viral_genes.hmms.h3f
Profiles (remainder) pressed into: viral_genes.hmms.h3p

$ nhmmscan viral_genes.hmms all_genes_calls.fna

# nhmmscan :: search DNA sequence(s) against a DNA profile database
# HMMER 3.3 (Nov 2019); http://hmmer.org/
# Copyright (C) 2019 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query sequence file:             gene_calls.fna
# target HMM database:             viral_genes.hmms
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       Aegir_VP1_000000000001  [L=2286]
Scores for complete hit:
    E-value  score  bias  Model     start    end  Description
    ------- ------ -----  --------  -----  -----  -----------

   [No hits detected that satisfy reporting thresholds]


Annotation for each hit  (and alignments):

   [No targets detected that satisfy reporting thresholds]


Internal pipeline statistics summary:
-------------------------------------
Query sequence(s):                         1  (4572 residues searched)

[...]
```

5. However, what we want is an MSA of all the genes against of the genomes. This is done using `hmmalign` which extracts the region of the genome that matches the HMM model. Important: The HMM must be run against all of the genomes one by one, not the other way around. I.e capsid_genes.hmm needs to be run against phage_1 phage_2 before head_tail_connector.hmm is run against phage 1 and phage_2. Otherwise you will have and MSA with an unequal length as it needs to include insertion/delation events in all tested species.

    For every hmm model I have created in my working directory, run `hmmalign` and only include the location of the genome that my HMM matches `(--trim)` and output this in the Pfam MSA format with the name being `capsid_genes_global.msa`.

    This will take about 2-3 minutes per HMM

```
$ for i in *hmm; do hmmalign --trim --outformat Pfam -o $(basename $i .hmm)_global.msa $i all_SAR11_phages.fna; done
```

Then for every `global.msa` I have created, trim it so I only include the useful information and put it into a `trimmed` file

```
$ for i in *global*; do egrep -v "#|//|^$" $i | sed -E -e 's/[[:blank:]]+/\n/g' > $(basename $i _global.msa)_trm.msa; done
```

Then for every phage I have in my file of all phage genomes, remove the MSA of each gene for each phage, concatinate it into one long MSA and append the viral name to the top, calling it `phage_1_new.msa`. Then concat the `_new.msa` files into one file.

```
$ for i in $(grep ">" all_SAR11_phages.fna | sed 's/>//g'); do MSA=$(grep $i -A 1 -h *_trm.msa | sed 's/^--$//g; /^$/d' | awk 'NR%2==0' ORS=''); echo ">"$i$'\n'$MSA > ${i}_new.msa; done
$ cat *_new.msa > viral.msa
```

```
$ cat viral.msa

>phage_1
--------------------...-..---------.....---G.........A...........TGA....AA.TG.........AAT.........A.........AACCA................ACTGA..........................
............................................................ATCAT......C..AaaaGA.........GGAA......GATAA......AACATTT.....................................
...................GAAAATGA.AG...CTaaG......AA.....C.......T..AATT......A.......
>phage_2
..........AGG.AACT....AAA.........T.GAG.......T.......ATAA-.......-.....A..TG....A....T..GAtatact...............................................................
..........atctagA-G.AGC...........................TGA.A............A.C...AAcgtgctcta..........ttaagtctttACGA.A.......A....AGA.AACT....AG..A..tA.A.TG.....
-....-.--..............--..GA.....T............T.........AG..........-T...AAC.T.
>phage_3
T--CA....................GCGTT.....A......A.....AGAA...Tg...............gT...A.......A...........CA........AA........................A.............G...........
T.................G..........T........A.........A...-.......-.-.......A...T..GC..................TT.....T.......................AA...........................A..
.TGtcatccC.......ATA.....A.....................A....AAAC--......-............-...

[...]
```

6. Check the MSA file is formatted correctly and remove duplicate phages if you wish using RAxML. Follow the onscreen suggestions and if you're happy with your MSA file (Or new one as I have used) you can start making your phylogenetic tree. 

```
$ raxml-ng --parse --msa viral.msa --model GTR+G --prefix T1
$ raxml-ng --parse --msa T1.raxml.reduced.phy --model GTR+G --prefix T2
```

6b. Other alternatives include `muscle` which can trim and refine your MSA but is VERY computationally intensive. Or you can use `trimal`. You could also use IQ-TREE to make your phylogenetic tree. 

```
$ muscle -in viral.msa -out viral.afa -refine
$ trimal -in viral.msa -out viral.afa -automated1

$ iqtree -s viral.afa -bb 1000 -m MFP 
```
