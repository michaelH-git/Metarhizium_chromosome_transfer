# Metarhizium Chromosome Transfer

This is the repository of sripts that were used to determine the horizontal transfer of an accessory chromosome within _Metarhizium robertsii_ and between _M. robertsii_ and _M. guizhouense_. 
The following software tools are required:

# Dependencies

•	[all2vcf]( https://github.com/MatteoSchiavinato/all2vcf)  

•	[bcftools]( https://samtools.github.io/bcftools/bcftools.html)	(version: 1.14) 

•	[bedtools](https://github.com/arq5x/bedtools2)	(version: 2.25.0)

•	[BioKIT]( https://jlsteenwyk.com/BioKIT/)	(version: 0.1.3)

•	[bowtie2]( https://github.com/BenLangmead/bowtie2)	(version: 2.4.4)

•	[braker2]( https://github.com/Gaius-Augustus/BRAKER)	(version: 2.1.6)

•	[busco]( https://busco.ezlab.org/)	(version: 5.5.0)

•	[bwa-mem2]( https://github.com/bwa-mem2/bwa-mem2)	(version: 2.2.1)

•	[Canu]( https://canu.readthedocs.io/en/latest/)	(version: 2.1.1)

•	[fastqc]( https://github.com/s-andrews/FastQC)	(version: 0.12.1)

•	[Flye]( https://github.com/fenderglass/Flye)	(version: 2.8.3)

•	[FMLRC2]( https://github.com/HudsonAlpha/fmlrc2)	(version: 0.1.4)

•	[gffread]( https://github.com/gpertea/gffread)	(version: 0.12.7)

•	[medaka]( https://github.com/nanoporetech/medaka)	(version: 1.4.3)

•	[minimap2]( https://github.com/lh3/minimap2)	(version: 2.24-r1122)

•	[mosdepth]( https://github.com/brentp/mosdepth)	(version: 0.3.4)

•	[mummer]( https://github.com/mummer4/mummer)	(version: 4.0.0rc1)

•	[nanofilt]( https://github.com/wdecoster/nanofilt)	(version: 2.3.0)

•	[picard]( https://github.com/wdecoster/nanofilt)	(version: 2.24.0)

•	[pilon]( https://github.com/broadinstitute/pilon/wiki)	(version: 1.24)

•	[racon]( https://github.com/broadinstitute/pilon/wiki)	(version: 1.4.20)

•	[REPET3]( https://urgi.versailles.inra.fr/Tools/REPET)	(version: 3.0)

•	[rust]( https://www.rust-lang.org/learn)	(version: 1.74.1)

•	[samtools]( https://www.htslib.org/)	(version: 1.3.1)

•	[Sniffels2](https://github.com/fritzsedlazeck/Sniffles)	(version: 2.2)

•	[tapestry]( https://github.com/johnomics/tapestry)	(version: 1.0.0)

•	[trimmomatic]( http://www.usadellab.org/cms/?page=trimmomatic)	(version: 0.39)

•	[whatshap](https://github.com/whatshap/whatshap)	(version: 1.6)


# Data

Sequencing reads have been deposited in the Sequence Read Archive and are available under the BioProject PRJNA1017668. The Nanopore-based assemblies and Gene annotations were deposited at NCBI under the BioProjects PRJNA1015426, PRJNA1015429, PRJNA1015431. 
The IDs of the previously published sequencing reads and assemblies are given below:


| Species	Strain |	Location	| BioSample |	Assembly | Accession | Sequencing run |
|----------------|------------|-----------|----------|-----------|----------------|
| M.robertsii	| ARSEF 23	| USA	| SAMN20219665	| 	| SRR15182435 |
| M.robertsii	| ARSEF 2575	| USA	| SAMN20219666	| 	| SRR15182434 |
| M.robertsii	| ESALQ 1426	| Brazil	| SAMN20243868	| 	| SRR15182431 |
| M.robertsii	| ESALQ 1635	| Brazil	| SAMN20243869	| 	| SRR15182430 |
| M.robertsii	| ESALQ 5168	| Brazil	| SAMN20219667	| 	| SRR15182429 |
| M.robertsii	| S1-CTAB-1	| China	| SAMN21033591	| 	| SRR15665281 |
| M.anisopliae	| 15R	| Korea	| SAMN26650411	| 	| SRR18320042 |
| M.anisopliae	| ARSEF-549	| Brazil	| SAMN03268434	| GCA_000814975.1	|  |
| M.anisopliae	| BRIP-53284	| 	| SAMN02981522	| GCA_000426985.1	|  |
| M.anisopliae	| BRIP-53293	| Australia	| SAMN02981521	| GCA_000426965.1	|  |
| M.anisopliae	| CQMa421-2	| China	| SAMN15446573	| 	| SRR12224770 |
| M.anisopliae	| E6	| Brazil	| SAMN02840975	| GCA_000739145.1	|  |
| M.anisopliae	| ESALQ 1076	| Brazil	| SAMN20243862	| 	| SRR15182425 |
| M.anisopliae	| ESALQ 1116	| Brazil	| SAMN20219669	| 	| SRR15182427 |
| M.anisopliae	| ESALQ 1175	| Brazil	| SAMN20243863	| 	| SRR15182424 |
| M.anisopliae	| ESALQ 1604	| 	| SAMN20243864	| 	| SRR15182423 |
| M.anisopliae	| ESALQ 1641	| Brazil	| SAMN20219670	| 	| SRR15182426 |
| M.anisopliae	| ESALQ 43	| Brazil	| SAMN20219668	| 	| SRR15182428 |
| M.anisopliae	| JEF-290	| South Korea	| SAMN11321365	| GCA_013305495.1	|  |
| M.anisopliae	| TNAU-MA-GDU	| India	| SAMN26879068	| GCA_023212845.1	|  |
| M.brunneum	| ARSEF-3297	| Mexico	| SAMN03268432	| GCF_000814965.1	|  |
| M.brunneum	| ARSEF 4556	| USA	| SAMN14166897	| 	| SRR11149706 |
| M.brunneum	| ESALQ 5022	| Brazil	| SAMN20243865	| 	| SRR15182422 |
| M.brunneum	| ESALQ 5181	| Brazil	| SAMN20243866	| 	| SRR15182433 |
| M.brunneum	| ESALQ 5286	| Brazil	| SAMN20243867	| 	| SRR15182432 |
| M. guizhouense	| ARSEF-977	| France	| SAMN03268433	| GCA_000814955.1	|  |
| M. majus	| ARSEF-297	| Samoa	| SAMN03268436	| GCF_000814945.1	|  |
| M acridum	| ARSEF-324	| Australia	| SAMN18235592	| GCA_019434415.1	|  |
| M. acridum	| CQMa-102	| 	| SAMN02981259	| GCF_000187405.1	|  |
| M. album	| ARSEF-1941	| Philipines	| SAMN03268435	| GCF_000804445.1	|  |



