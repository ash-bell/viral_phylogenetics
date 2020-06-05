# Viral Phylogeny
---
Without a set of common genes shared between all viruses, it becomes difficult to deduce phylogeny of viruses. Other method of gene-sharing maps exist to combat this problem and are preferred. However, if you would like a classic tree, you can use this instead. This method uses a series of custom genes that the user selects to base a phylogenetic tree on. 

## Before you start:
---
*   [HMMER](http://hmmer.org/) installed `conda install -c bioconda hmmer` 
*   [MUSCLE](https://drive5.com/muscle/) installed `conda install -c bioconda muscle`
*   [trimal](http://trimal.cgenomics.org/) installed `conda install -c bioconda trimal`
*   [raxml-ng](https://github.com/amkozlov/raxml-ng) installed if you want to make a phylogenetic tree with it `conda install -c bioconda raxml-ng`

## Files:
*   all_gene_calls.faa --> file of all my phage genomes gene calls as amimoacids
*   I like to genecall them with glimmer3 and use prodigal to convert the nuclotides to amino acids

```
$ long-orfs -n -t 1.15 REPLACE_ME run3.longorfs
$ extract -t REPLACE_ME run3.longorfs > run3.train
$ build-icm -r run3.icm < run3.train
$ glimmer3 -o50 -g110 -t30 REPLACE_ME run3.icm run3.run1
$ tail -n +2 run3.run1.predict > run3.coords
$ upstream-coords.awk 25 0 run3.coords | extract REPLACE_ME - > run3.upstream
$ elph run3.upstream LEN=6 | get-motif-counts.awk > run3.motif
$ startuse="$(start-codon-distrib -3 REPLACE_ME run3.coords)"
$ glimmer3 -o50 -g110 -t30 -b run3.motif -P $startuse REPLACE_ME run3.icm FINAL
$ extract -t REPLACE_ME FINAL.predict > $(basename REPLACE_ME .fasta)_Glimmer.fnn
$ anvi-script-reformat-fasta -o $(basename REPLACE_ME .fasta)_rfmt.fna --simplify-names --prefix $(basename REPLACE_ME .fasta) $(basename REPLACE_ME .fasta)_Glimmer.fnn
$ prodigal -i $(basename REPLACE_ME .fasta)_rfmt.fna -a $(basename REPLACE_ME .fasta)_Glimmer.faa -q -p meta
$ rm run3* FINAL*
```

This should result in a file that looks like this
```
$ cat all_SAR11_phages.fna

>Aegir_VP1_000000000001_1 # 1 # 2286 # 1 # ID=1_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.287
MLTETLSKLETLGPVGIKIRSIMNERYDFLNKESNGPLKSGIVRMIDKHKIDYNMMIQLSYSIVATGTTEGQNLIQLAIAIGNRIRNYYKLPTKAELDLRLGVFVLNAYALNNMLIIKLVNDFKSFNKAKTVYKVYMGYNRSDMRKLMKEFNEVHDPFKPLFTKAPEWEFGTVLNPNGESIKLIKQSKIDTIARINKHNTPIVLNATNKKQAVAFYVNPEVYEIYKWALQTNQHCFEHNSVDTIAKERKEAKKAEALQVLKATELYVGKKFYQQYTCDSRGRFYPLSAYLNELNSDNAKGMLSFFEGKPLGNNGKAQLFHHIANMWGEDKLSHADRVKFVEDNYYEFVTYGSKPKDARGWMQAAEPIQFLSAVIELAKLDKHFVANGTVEDFISHTICYRDGSNNGLQWLFSLVKDNKNGHLVNIKPTTDNKPGDMYNHVAVSVKDIMHNKAKEEDNLSLDYYNLYFKSIEKIRNRWRIAELNNDKNVENKKRLIKWYQKRYRTELKLTDIIYWDKAKFTIKEWRKIVKRNVMTYGYSATKQGMGEQIIQDTRDIDNVYLSNKQHSAARALGSLVYTTIETEFPEVAAAMKLFKDNCAAYMKKHNKQYSHNTLISNFPFTQHYVKYKSARVKLTDGLYIMNNDKSFKWVNRVDFVIKTELPILNIGKAKAAISPNSIHNLDSLHLMLVIDECDFDIVSAHDSYGAHACNVNMMQKVIRSQFKRIIDANPLQHNLNETGNLVPMIKQGQLDSSEILKSEFAFA
>Aegir_VP1_000000000002_1 # 1 # 144 # 1 # ID=2_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.222
MFDKYVYPILEKIGEGIEKIYWFCSNNKQEVTWFSFGFIVCALINLII
```

*   Gene called your organisms and BLAST them against a database to get an annotated gene function for each. Then select genes that are specific to your taxa of interest and find them within your annotations. For example, for a viral phlogenetics I chose capsid, DNA polymerase, exonuclease, endonuclease, helicase, head-tail-connector and RNA polymerase. I like to upload my files to [eggnog-mapper](http://eggnog-mapper.embl.de/).

## Tutorial
---
1. Create a file for each catagory of genes you have chosen.

```
$ capsid.faa

>VIRSorter_AAJR_prophage_000000000025_1 # 1 # 765 # 1 # ID=25_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.318
MVDKIEINQDENNISLEEQSNKQETNSQTTPEAQTKETSSEKPSWLPDKFSNAEELAKAYSELEKKFSSPSIEENKSYENKKTDLTIKQSEEAEQGKLDKFYNEYAETGKLTDTSYNELNKLGLDRQVVDGYIDGQTALSEQKATSIMSTVGGKEQYSEMISWASKNLAPEEVQAFNHTIDSGSLEQAQLAIAGVQAKYNQNNAEPNLFSGNKTNSNLGYRSVGEMLKDINDPRYSTDSAFRADVEEKVKLSNAI
>VIRSorter_AACP_phage_000000000014_1 # 1 # 801 # 1 # ID=14_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.353
MSTDRLEVKPDDNNNETLEQSAEKLKEQGVDINKDLSVNPNGEGIEVREPKVETESTEKRPEWLPEKFQNAEELAKAYGTLEKEFSGRTKEEVKPTEEVKADEAPQTGLDKYYEEFADKGELAEKSYSELAKLGLDKNLVDTYIEGQKLVSETNTKAIQDIAGGKEEYTELVEWAGKNLSEAETKVFNDMVDGGDIETAKFAVQGLMAKSGANPKQPSLYEGTSDTVSKDAFASVAQVTEAMNDPRYDSDPAYRQLVEDKIGRSTVL
>VIRSorter_AALP_phage_000000000004_1 # 1 # 822 # 1 # ID=4_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.347
MSTDRVEIKPDAIEKSLEQSAEDLKTNDGVDVSKDVAVNKSGEGARINELSQEDLQQSTEDKPDWLPSKFKNAEELAKAYGELEKSFSSRKQEEAPKETSIPEVRKPTEGQEALGKFYDEYAANSSLSDKSYEELATKHGLSKELVDGYIEGQKAIGDNQTKAIHDLTGGGEKYTELMDWAGKNLSEQEQQAYNSMVDSGNVDAAKLAVQGLMSKAGVNYNPKQPELFEGGDQIPNDSFESVSQVTEAMNDPRYAKDPAYRKKVTDKIARSSVI

[...]
```

2. Run `muscle` and `trimal` on the genes to generate a multiple sequence alignment and refine it

```
$ muscle -in capsid.faa -out capsid.msa
$ trimal -in capsid.msa -out capsid.afa -automated1 -keepseqs

```
It should look something like this:
```
$ cat capsid.afa

>VIRSorter_AALP_phage_000000000009_1 # 1 # 948 # 1 # ID=9_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.373
MSNANTSRLGLVND-TGTSYDALFLKVFSGE--------VLASFGRENKMLDMTTV--RT
IGSGKSA----TFPVTGTIASSYHTAGNEILGTAVNHNE----RIINIDDMLISHAFIAE
IDELKNHFDARSIYSNEMGRALSNK--VDQHLVQLMVLASQASATITGGNGGYEITDADA
KTNAESLIASVFDANQKLDENDVPTSDRYVVVTPDVYYQLCN--VDKLISRDFSKNNGDF
GSGTVVSIAGIPVVKSNSAKLAYDDNSGAISGTNNTYNVDAQHIAATVFH--KSAVGTVK
(py37) [agb214@login01 amino_acids]$ head capsid.msa -n 20
>VIRSorter_AALP_phage_000000000009_1 # 1 # 948 # 1 # ID=9_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.373
MSNANTSRLGLVND-TGTSYDALFLKVFSGE--------VLASFGRENKMLDMTTV--RT
IGSGKSA----TFPVTGTIASSYHTAGNEILGTAVNHNE----RIINIDDMLISHAFIAE
IDELKNHFDARSIYSNEMGRALSNK--VDQHLVQLMVLASQASATITGGNGGYEITDADA
KTNAESLIASVFDANQKLDENDVPTSDRYVVVTPDVYYQLCN--VDKLISRDFSKNNGDF
GSGTVVSIAGIPVVKSNSAKLAYDDNSGAISGTNNTYNVDAQHIAATVFH--KSAVGTVK
LKDL--VLENTYDPR-----RLGHLLTARLALGHNILRPEAAVSIKQA

[...]
```

4. Build your Hidden Markov Model (HMM) for searching for simular genes

```
$ hmmbuild capsid.hmm capsid.afa
```
It should look like this:
```
$ cat capsid.hmm

HMMER3/f [3.3 | Nov 2019]
NAME  capsid
LENG  221
ALPH  amino
RF    no
MM    no
CONS  yes
CS    no
MAP   yes
DATE  Mon Jun  1 11:22:15 2020
NSEQ  5
EFFN  0.893555
CKSUM 1592742671
STATS LOCAL MSV      -10.4788  0.70418
STATS LOCAL VITERBI  -11.2667  0.70418
STATS LOCAL FORWARD   -4.8939  0.70418
HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y   
            m->m     m->i     m->d     i->m     i->i     d->m     d->d
  COMPO   2.46880  4.69472  2.73572  2.49653  3.47874  2.83545  3.82868  2.96866  2.56685  2.56273  3.82548  2.90495  3.54610  3.04849  3.03048  2.55907  2.83349  2.75135  4.74168  3.42452
          2.68618  4.42225  2.77519  2.73123  3.46354  2.40513  3.72494  3.29354  2.67741  2.69355  4.24690  2.90347  2.73739  3.18146  2.89801  2.37887  2.77519  2.98518  4.58477  3.61503

[...]
```

5. Now we want an MSA of all the genes against of the genomes. This is done using `hmmalign` which extracts the region of the genome that matches the HMM model. 
**Important:** The HMM must be run against all of the genomes one by one, not the other way around. I.e capsid.hmm needs to be run against phage_1 phage_2 before head_tail_connector.hmm is run against phage 1 and phage_2. Otherwise you will have and MSA with an unequal length as it needs to include insertion/delation events in all tested species.

    For every hmm model I have created in my working directory, run `hmmalign` and only include the location of the genome that my HMM matches `(--trim)` and output this in the Pfam MSA format with the name being `capsid_genes_global.msa`.

    This will take about 30 seconds per HMM

```
$ for i in *.hmm; do echo $i; hmmalign --trim --outformat Pfam -o $(basename $i .hmm)_global.afa $i SAR11_phage.faa; done
```

Then for every `global.msa` I have created, trim it so I only include the useful information and put it into a `trimmed` file

```
$ for i in *global*; do egrep -v "#|//|^$" $i | sed -E -e 's/[[:blank:]]+/\n/g; s/^/>/g' > $(basename $i _global.msa)_trm.msa; done
```

Now this provides the last oppertunity to trim MSAs before they are concatinated into one file. So we refined them with `muscle` and `trimal`

```
$ for i in *_trm.msa; do muscle -in $i -out $(basename $i .msa).afa; trimal -in $(basename $i .msa).afa -out $(basename $i .msa)_final.msa -automated1 -keepseqs; done
```

Then for every phage I have in my file of all phage genomes, remove the MSA of each gene for each phage, concatinate it into one long MSA and append the viral name to the top, then concat the files into one file.

```
$ for i in *_trm_final.msa; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $i > $(basename $i .msa).afa; done # convert multi-contig fasta into one contig

$ for i in $(grep ">" SAR11_phage.faa); do MSA=$(grep $i -A 1 -h *_trm_final.afa | sed 's/^--$//g; /^$/d' | awk 'NR%2==0' ORS=''); echo $i$'\n'$MSA >> viral.msa; done
```

```
$ cat viral.msa

>Bylgja_VP3
KTFENESRPEWLPEKFKSPEDMAKAYGELKEETNKAEKAVENAGLNMSSLQEEYNEGGQLKYVFDLESNGLYTIHCIVLKDIETNKIIQDVNEALKLLSEAELIIGHNIIKYDIPVLKKLSVKNFKNFLYLAWKHLSLPEPTEIQYDLADFLQEPNKRIVIEAFRGVGKSWITSAFVCHQTLHAEETQGKYLVASEIEDLKEKLGADKVIVALTDKNNFRKDVLPTYKDNPLRKFLIDEYKDTAENRYESLSEIKEHYLDRGRECSELTIP--TLIPEQHQTQSSDFQSVGSRGVNNLASELVQGPLIKRKINKETVTKFNYQTGGKNV--QIANYYDKLVAQKLRY-PDKFQWIGDSKLLSGVIWDSMLSSARLGMSFLQNCAKVLSKNGHAIRWTNPVGFPVIQDYPEFKSMRVKTRMTELQAINTMLSFIGESPVSSITGNIGTDAKNILDETSMSVQSQGWFFN-RNVSTTRDTSN
>Eistla_VP4
ETQSTQSKPEGLPEKFNSVEDLAKSYAELEPKTETAEKAVADAGLDMSSLQQEYSEKGELKYIFDIETDGLMKIHCLVLKDVDSNQILSSVEEALDKLSKADVIIGHNIIKFDIPVIKKLRIKNFKNFLYLCWKHLNLPNPTPIQYDIADYLQSPEKRLVIEAFRGVGKSWITSAFVCHQ-----------------------LEADDYVVALTDSKNFRKDVLPTYKNNALREYVIKKHQQSAKQRYETLKQHREHFLDRAQECSELTIP--SLIPPDGFHSSTDFQSVGARGVNNLASDFISGALSKRNIDFNTAQKFNYQTGGRPC--QIANYYDKLVAQKLRY-PDKFQWLGDAKLLAAVIWDSILKSARIGMDYLQTIARIVAKEQLPVHWVTPVGFPVYQSYPEMKSKRVKAMLTELQAINIMLSVIGEAPVNSITGTTSVDAKNILDETSMSIQSQGWHFNTHYTSLALDQDN
>Eyrgjafa_VP5
DKSTDDVRPSWLPEKFKSAEDMAKAYSEL-EKKQSQMRADAEAGEGMDKFYTEYQDKGELKIILDLETNGFLVIHCIVCKDIETHQVYKNLNDCLELLNKVEAIIGHNVLGFDLPVLKRCKINDFRNFLYLTWKHLRLPEPTPIQYDIADYLANGSTRCIISAFRGVGKSWITASYILWRTLHSDLGKGKTLLQQAINYYKEKTKSKEVILAFSDSKNFRKEFDKTYKSHPLRKWAEQNYKSLIESQYTKMEVDREQYLERARELAKLTIP--HLYPPKGANEATEYQSVGSRGVTNLASNMIDGALPSRKINSDTCKKFNYQTGGQPV--HIANYYDKIVAQKLRF-QDKFLWLGDVDLLYGLVFQQKLVSQMINDDIKQMTMGKAGNHQTGLKIICQCSLIILKNVIDGISNEITINKSELEAVNTILSTIGESPLNTLSGSLPVDAKNVLSEVSREVQSQGWHFNTHKVTLSRDTDN

[...]
```

6. Check the MSA file is formatted correctly and remove duplicate phages if you wish using RAxML. Follow the onscreen suggestions and if you're happy with your MSA file (Or new one as I have used) you can start making your phylogenetic tree. 

```
$ raxml-ng --check --msa viral.msa --model LG+G8+F --prefix T1 --data-type AA
$ raxml-ng --parse --msa viral.msa --model LG+G8+F --prefix T2 --data-type AA
```

6b. Other alternatives include `muscle` which can trim and refine your MSA or you can use `trimal`. You could also use IQ-TREE to make your phylogenetic tree. 

```
$ muscle -in viral.msa -out viral.afa
$ trimal -in viral.msa -out viral.afa -automated1 -keepseqs # or not keep seqs

$ iqtree -s viral.afa -bb 1000 -m MFP 
```
7. Lets make the tree with raxml-ng. Here are some options:
The tree you should use should be under {prefix}.raxml.support

```
# Just make me a standard tree with bootstrappings
$ raxml-ng --all --msa viral.msa --model GTR+G --prefix T3 --threads 1

# make a quick tree to check its rough structure
$ raxml-ng --search1 --msa viral.msa --model GTR+G --prefix T4 --threads 1 

# make a tree with 50 random starting+most pasimonious trees, with 2000 rounds of bootstrapping 
$ raxml-ng --all --msa viral.msa --model GTR+G --prefix T5 --threads 1 --tree pars{50},rand{50} --bs-trees 2000
```
After the tree has been compiled, I would recomend [iTOL](https://itol.embl.de/) to view and edit the tree.
---

7b. Because we are using a multi gene tree, genes can mutate and different rates #biology. Each gene may have mutations specific to them, therefore if we use the same model for all genes, this may be incorrect. To get around this, we partition up an MSA by genes, and run a model on each, combining the results after. **IMPORTANT:** This does not affect the structure of the tree drastically (if at all) so is not recomended, but I have included this just in case. To do this, you need to make a file like below. You may need to do some re-jiggling of the numbers if your MSA file has been trimmed. **IMPORTANT:** If you are making a tree of organisms that are very closely related, you need the dev verion of [raxml-ng](https://github.com/amkozlov/raxml-ng/wiki/Installation#building-development-branch) 

```
$ for i in *_trm_final.afa; do echo $i $(grep -v ">" $i | awk 'NR%2==0' | head -n 1 | wc -c); done

capsid_trm_final.afa 181
dna_pol_trm_final.afa 488
endonuclease_trm_final.afa 458
exonuclease_trm_final.afa 192
head_tail_trm_final.afa 433
helicase_trm_final.afa 408
rna_pol_trm_final.afa 176
tube_trm_final.afa 135

$ cat viral.part

LG+G8+F, capsid=1-180
LG+G8+F, dna_pol=181-667
LG+G8+F, endonuclease=668-1124
LG+G8+F, exonuclease=1125-1315
LG+G8+F, head_tail=1316-1747
LG+G8+F, helicase=1748-2154
LG+G8+F, rna_pol=2155-2329
LG+G8+F, tube=2330-2463
```

```
$ raxml-ng --all --msa viral.afa --model viral.part --prefix PARTITION_MODEL --threads 16 --tree pars{50},rand{50} --bs-trees 1000
```
If you want to learn more about the options for RAxML, I highly recomend their [tutorial](https://github.com/amkozlov/raxml-ng/wiki/Tutorial)

**NOTE:** If any of the raxml-ng commands fail, rerun the command, they have [checkpoints](https://github.com/amkozlov/raxml-ng/wiki/Advanced-Tutorial#checkpointing) and so your progess won't be lost. 
