#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It corrects the consecutives lines, which represent the same read, such that they look like read which has been split.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2015 Daniel Nicorici

This file is part of FusionCatcher.

FusionCatcher is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FusionCatcher is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FusionCatcher (see file 'COPYING.txt').  If not, see
<http://www.gnu.org/licenses/>.

By default, FusionCatcher is running BLAT aligner
<http://users.soe.ucsc.edu/~kent/src/> but it offers also the option to disable
all its scripts which make use of BLAT aligner if you choose explicitly to do so.
BLAT's license does not allow to be used for commercial activities. If BLAT
license does not allow to be used in your case then you may still use
FusionCatcher by forcing not use the BLAT aligner by specifying the option
'--skip-blat'. Fore more information regarding BLAT please see its license.

Please, note that FusionCatcher does not require BLAT in order to find
candidate fusion genes!

This file is not running/executing/using BLAT.

"""



# info SAM format
"""
http://samtools.github.io/hts-specs/SAMv1.pdf


Col   Field   Type     Regexp/Range              Brief description
--------------------------------------------------------------------------------
1     QNAME   String   [!-?A-~]{1,255}           Query template NAME
2     FLAG    Int      [0,2^16-1]                bitwise FLAG
3     RNAME   String   \*|[!-()+-<>-~][!-~]*     Reference sequence NAME
4     POS     Int      [0,2^31-1]                1-based leftmost mapping POSition
5     MAPQ    Int      [0,2^8-1]                 MAPping Quality
6     CIGAR   String   \*|([0-9]+[MIDNSHPX=])+   CIGAR string
7     RNEXT   String   \*|=|[!-()+-<>-~][!-~]*   Ref. name of the mate/next read
8     PNEXT   Int      [0,2^31-1]                Position of the mate/next read
9     TLEN    Int      [-2^31+1,2^31-1]          observed Template LENgth
10    SEQ     String   \*|[A-Za-z=.]+            segment SEQuence
11    QUAL    String   [!-~]+                    ASCII of Phred-scaled base QUALity+33


CIGAR: CIGAR string. The CIGAR operations are given in the following table (set `*' if unavailable):

Op     BAM    Description
--------------------------------------------------------------------------------
M      0      alignment match (can be a sequence match or mismatch)
I      1      insertion to the reference
D      2      deletion from the reference
N      3      skipped region from the reference
S      4      soft clipping (clipped sequences present in SEQ)
H      5      hard clipping (clipped sequences NOT present in SEQ)
P      6      padding (silent deletion from padded reference)
=      7      sequence match
X      8      sequence mismatch

Notes:
------
* H can only be present as the first and/or last operation.
* S may only have H operations between them and the ends of the CIGAR string.
* For mRNA-to-genome alignment, an N operation represents an intron. For other
  types of alignments, the interpretation of N is not defined.
* Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.

"""

"""
EXAMPLE conversion SAM to PSL:

SAM:
----
@SQ	SN:ENSG00000075624|ENSG00000122566|36634	LN:48237

QNAME                                               POS                     RNEXT       SEQ
        FLAG                                                MAPQ                PNEXT                                                                                                           QUAL
            RNAME                                               CIGAR               TLEN                                                                                                                                                                                                                TAGS
000CG	16	ENSG00000075624|ENSG00000122566|36634	41506	3	62M151N39M	*	0	0	TGCAGAAATACCATACCATCAATGGTCATAATGCAGAAGTAAGAAAGGCTTTGTCTAGACAAGAAATTTCGGACCAGGACCAGGAAGTAACTTTAGAGGAG	cbcdeddbbbddddceeeeeeggggghiiiiiiiiihiiiiiiiiihhgihhiiihiiihiihiiiiiihggeihhhgghiiiiiiiigggggeeeeebbb	NH:i:2	HI:i:1	AS:i:91	nM:i:0
000CG	272	ENSG00000122566|ENSG00000075624|11603	4872	3	62M151N39M	*	0	0	TGCAGAAATACCATACCATCAATGGTCATAATGCAGAAGTAAGAAAGGCTTTGTCTAGACAAGAAATTTCGGACCAGGACCAGGAAGTAACTTTAGAGGAG	cbcdeddbbbddddceeeeeeggggghiiiiiiiiihiiiiiiiiihhgihhiiihiiihiihiiiiiihggeihhhgghiiiiiiiigggggeeeeebbb	NH:i:2	HI:i:2	AS:i:91	nM:i:0
002MX	16	ENSG00000066468|ENSG00000075624|120125	155339	0	81M	*	0	0	TCCCCCATGCCATCCTGCTTCTGGATCTGCCTGGCCGGGACCTGACTAACTACCTCATGAAGAACCACACCGAGCGCGGCT	LWQJNQTKKGKKGVOYLGGGLLGMGGGUGG\TFQT\\MMHHRUNOHHHZMHMT[WYUOHMMHGHPGGWLZ[OGOGGZ\ZRZ	NH:i:6	HI:i:1	AS:i:67	nM:i:6
002MX	272	ENSG00000075624|ENSG00000066468|36634	35214	0	81M	*	0	0	TCCCCCATGCCATCCTGCTTCTGGATCTGCCTGGCCGGGACCTGACTAACTACCTCATGAAGAACCACACCGAGCGCGGCT	LWQJNQTKKGKKGVOYLGGGLLGMGGGUGG\TFQT\\MMHHRUNOHHHZMHMT[WYUOHMMHGHPGGWLZ[OGOGGZ\ZRZ	NH:i:6	HI:i:2	AS:i:67	nM:i:6
002MX	272	ENSG00000075624|ENSG00000122566|36634	35214	0	81M	*	0	0	TCCCCCATGCCATCCTGCTTCTGGATCTGCCTGGCCGGGACCTGACTAACTACCTCATGAAGAACCACACCGAGCGCGGCT	LWQJNQTKKGKKGVOYLGGGLLGMGGGUGG\TFQT\\MMHHRUNOHHHZMHMT[WYUOHMMHGHPGGWLZ[OGOGGZ\ZRZ	NH:i:6	HI:i:3	AS:i:67	nM:i:6
002MX	272	ENSG00000122566|ENSG00000075624|11603	46817	0	81M	*	0	0	TCCCCCATGCCATCCTGCTTCTGGATCTGCCTGGCCGGGACCTGACTAACTACCTCATGAAGAACCACACCGAGCGCGGCT	LWQJNQTKKGKKGVOYLGGGLLGMGGGUGG\TFQT\\MMHHRUNOHHHZMHMT[WYUOHMMHGHPGGWLZ[OGOGGZ\ZRZ	NH:i:6	HI:i:4	AS:i:67	nM:i:6
002MX	272	ENSG00000075624|ENSG00000182472|36634	35214	0	81M	*	0	0	TCCCCCATGCCATCCTGCTTCTGGATCTGCCTGGCCGGGACCTGACTAACTACCTCATGAAGAACCACACCGAGCGCGGCT	LWQJNQTKKGKKGVOYLGGGLLGMGGGUGG\TFQT\\MMHHRUNOHHHZMHMT[WYUOHMMHGHPGGWLZ[OGOGGZ\ZRZ	NH:i:6	HI:i:5	AS:i:67	nM:i:6
002MX	272	ENSG00000182472|ENSG00000075624|39718	74932	0	81M	*	0	0	TCCCCCATGCCATCCTGCTTCTGGATCTGCCTGGCCGGGACCTGACTAACTACCTCATGAAGAACCACACCGAGCGCGGCT	LWQJNQTKKGKKGVOYLGGGLLGMGGGUGG\TFQT\\MMHHRUNOHHHZMHMT[WYUOHMMHGHPGGWLZ[OGOGGZ\ZRZ	NH:i:6	HI:i:6	AS:i:67	nM:i:6
003IV	16	ENSG00000196924|ENSG00000189143|26115	17130	3	1S59M	*	0	0	GGGTGACGAGATCCCCTTCACCCCGCACCGCGCGCGTGCCGTGCCCACCGGGGATGTCAG	GJWJGGGGGJEVVOZWJGOEOFFKFF\FFTFFGRWGGGPG_HG^^___LYGHMQHGGGGL	NH:i:2	HI:i:1	AS:i:48	nM:i:5
003IV	272	ENSG00000189143|ENSG00000196924|33152	50282	3	1S59M	*	0	0	GGGTGACGAGATCCCCTTCACCCCGCACCGCGCGCGTGCCGTGCCCACCGGGGATGTCAG	GJWJGGGGGJEVVOZWJGOEOFFKFF\FFTFFGRWGGGPG_HG^^___LYGHMQHGGGGL	NH:i:2	HI:i:2	AS:i:48	nM:i:5
0046S	0	ENSG00000075624|ENSG00000182472|36634	36247	0	32M1I68M	*	0	0	CGAGGACTTTGATTGCACATTGTTGTTTTTTTTAATAGTCATTCCAAATATGAGATGCATTGTTACAGGAAGTCCCTTGCCATCCTAAAAGCCACCCCACT	bbbeeeeegggggiiihiiiiihiigiiiiiiifhhfhhiiiiiihihihiihiiiiiiiiiggggggeeeeeddbddcccccccccccccccccacccc_	NH:i:6	HI:i:1	AS:i:92	nM:i:1
0046S	256	ENSG00000075624|ENSG00000066468|36634	36247	0	32M1I68M	*	0	0	CGAGGACTTTGATTGCACATTGTTGTTTTTTTTAATAGTCATTCCAAATATGAGATGCATTGTTACAGGAAGTCCCTTGCCATCCTAAAAGCCACCCCACT	bbbeeeeegggggiiihiiiiihiigiiiiiiifhhfhhiiiiiihihihiihiiiiiiiiiggggggeeeeeddbddcccccccccccccccccacccc_	NH:i:6	HI:i:2	AS:i:92	nM:i:1
0046S	256	ENSG00000075624|ENSG00000122566|36634	36247	0	32M1I68M	*	0	0	CGAGGACTTTGATTGCACATTGTTGTTTTTTTTAATAGTCATTCCAAATATGAGATGCATTGTTACAGGAAGTCCCTTGCCATCCTAAAAGCCACCCCACT	bbbeeeeegggggiiihiiiiihiigiiiiiiifhhfhhiiiiiihihihiihiiiiiiiiiggggggeeeeeddbddcccccccccccccccccacccc_	NH:i:6	HI:i:3	AS:i:92	nM:i:1
0046S	256	ENSG00000182472|ENSG00000075624|39718	75965	0	32M1I68M	*	0	0	CGAGGACTTTGATTGCACATTGTTGTTTTTTTTAATAGTCATTCCAAATATGAGATGCATTGTTACAGGAAGTCCCTTGCCATCCTAAAAGCCACCCCACT	bbbeeeeegggggiiihiiiiihiigiiiiiiifhhfhhiiiiiihihihiihiiiiiiiiiggggggeeeeeddbddcccccccccccccccccacccc_	NH:i:6	HI:i:4	AS:i:92	nM:i:1
0046S	256	ENSG00000122566|ENSG00000075624|11603	47850	0	32M1I68M	*	0	0	CGAGGACTTTGATTGCACATTGTTGTTTTTTTTAATAGTCATTCCAAATATGAGATGCATTGTTACAGGAAGTCCCTTGCCATCCTAAAAGCCACCCCACT	bbbeeeeegggggiiihiiiiihiigiiiiiiifhhfhhiiiiiihihihiihiiiiiiiiiggggggeeeeeddbddcccccccccccccccccacccc_	NH:i:6	HI:i:5	AS:i:92	nM:i:1
0046S	256	ENSG00000066468|ENSG00000075624|120125	156372	0	32M1I68M	*	0	0	CGAGGACTTTGATTGCACATTGTTGTTTTTTTTAATAGTCATTCCAAATATGAGATGCATTGTTACAGGAAGTCCCTTGCCATCCTAAAAGCCACCCCACT	bbbeeeeegggggiiihiiiiihiigiiiiiiifhhfhhiiiiiihihihiihiiiiiiiiiggggggeeeeeddbddcccccccccccccccccacccc_	NH:i:6	HI:i:6	AS:i:92	nM:i:1
0061C	16	ENSG00000066468|ENSG00000137309|120125	126030	3	19M672N55M1311N27M	*	0	0	GGGTGCTGCCAAGACCCGGAAAACCACCACAACTCCAGGAAGGAAACCAAGGGGCAGACCCAAAAAAACTGGAGAAGGAGGAAGAGGAGGGCATCTCGCAG	ca`^W_]``bcbXOaaab][[[F_^Xa^]bZZTG`Z_bbabdaeedd]b_cgdcge\Thhhifihgecfhgfihffhhhiiihgfeafaggcgeeeee___	NH:i:2	HI:i:1	AS:i:90	nM:i:4
0061C	272	ENSG00000137309|ENSG00000066468|9359	5905	3	19M672N55M1311N27M	*	0	0	GGGTGCTGCCAAGACCCGGAAAACCACCACAACTCCAGGAAGGAAACCAAGGGGCAGACCCAAAAAAACTGGAGAAGGAGGAAGAGGAGGGCATCTCGCAG	ca`^W_]``bcbXOaaab][[[F_^Xa^]bZZTG`Z_bbabdaeedd]b_cgdcge\Thhhifihgecfhgfihffhhhiiihgfeafaggcgeeeee___	NH:i:2	HI:i:2	AS:i:90	nM:i:4

PSL:
----

matches					tNumInsert					qEnd
	misMatches				tBaseInsert					tName
		repMatches				strand															tSize						blockSizes
			nCount					qName																tStart						qStarts
				qNumInsert					qSize																tEnd						tStarts
					qBaseInsert					qStart																	blockCount
252	0	0	0	0	0	0	151	-	000CG	101	0	101	ENSG00000075624|ENSG00000122566|36634	48237	41505	41757	2	62,39,	0,62,	41505,41718,
252	0	0	0	0	0	0	151	-	000CG	101	0	101	ENSG00000122566|ENSG00000075624|11603	48237	4871	5123	2	62,39,	0,62,	4871,5084,
81	0	0	0	0	0	0	0	-	002MX	81	0	81	ENSG00000066468|ENSG00000075624|120125	156759	155338	155419	1	81,	0,	155338,
81	0	0	0	0	0	0	0	-	002MX	81	0	81	ENSG00000075624|ENSG00000066468|36634	156759	35213	35294	1	81,	0,	35213,
81	0	0	0	0	0	0	0	-	002MX	81	0	81	ENSG00000075624|ENSG00000122566|36634	48237	35213	35294	1	81,	0,	35213,
81	0	0	0	0	0	0	0	-	002MX	81	0	81	ENSG00000122566|ENSG00000075624|11603	48237	46816	46897	1	81,	0,	46816,
81	0	0	0	0	0	0	0	-	002MX	81	0	81	ENSG00000075624|ENSG00000182472|36634	76352	35213	35294	1	81,	0,	35213,
81	0	0	0	0	0	0	0	-	002MX	81	0	81	ENSG00000182472|ENSG00000075624|39718	76352	74931	75012	1	81,	0,	74931,
59	0	0	0	0	0	0	0	-	003IV	60	1	60	ENSG00000196924|ENSG00000189143|26115	59267	17129	17188	1	59,	1,	17129,
59	0	0	0	0	0	0	0	-	003IV	60	1	60	ENSG00000189143|ENSG00000196924|33152	59267	50281	50340	1	59,	1,	50281,
100	0	0	0	0	1	0	0	+	0046S	101	0	101	ENSG00000075624|ENSG00000182472|36634	76352	36246	36346	2	32,68,	0,33,	36246,36278,
100	0	0	0	0	1	0	0	+	0046S	101	0	101	ENSG00000075624|ENSG00000066468|36634	156759	36246	36346	2	32,68,	0,33,	36246,36278,
100	0	0	0	0	1	0	0	+	0046S	101	0	101	ENSG00000075624|ENSG00000122566|36634	48237	36246	36346	2	32,68,	0,33,	36246,36278,
100	0	0	0	0	1	0	0	+	0046S	101	0	101	ENSG00000182472|ENSG00000075624|39718	76352	75964	76064	2	32,68,	0,33,	75964,75996,
100	0	0	0	0	1	0	0	+	0046S	101	0	101	ENSG00000122566|ENSG00000075624|11603	48237	47849	47949	2	32,68,	0,33,	47849,47881,
100	0	0	0	0	1	0	0	+	0046S	101	0	101	ENSG00000066468|ENSG00000075624|120125	156759	156371	156471	2	32,68,	0,33,	156371,156403,
2084	0	0	0	0	0	0	1983	-	0061C	101	0	101	ENSG00000066468|ENSG00000137309|120125	129484	126029	128113	3	19,55,27,	0,19,74,	126029,126720,128086,
2084	0	0	0	0	0	0	1983	-	0061C	101	0	101	ENSG00000137309|ENSG00000066468|9359	129484	5904	7988	3	19,55,27,	0,19,74,	5904,6595,7961,


Some random example of a valid PSL file:
-----------------------------------------
matches					tNumInsert					qEnd
	misMatches				tBaseInsert					tName
		repMatches				strand															tSize						blockSizes
			nCount					qName																tStart						qStarts
				qNumInsert					qSize																tEnd						tStarts
					qBaseInsert					qStart																	blockCount
101	0	0	0	0	0	1	151	-	000CG/1	101	0	101	ENSG00000075624|ENSG00000122566|36634	48237	41505	41757	2	62,39,	0,62,	41505,41718,
101	0	0	0	0	0	1	151	-	000CG/1	101	0	101	ENSG00000122566|ENSG00000075624|11603	48237	4871	5123	2	62,39,	0,62,	4871,5084,
101	0	0	0	0	0	1	42416	-	00ABP/1	101	0	101	ENSG00000227143|ENSG00000066468|7267	127392	7704	50221	2	60,41,	0,60,	7704,50180,
101	0	0	0	0	0	1	42416	-	00ABP/1	101	0	101	ENSG00000184009|ENSG00000066468|13877	134002	14314	56831	2	60,41,	0,60,	14314,56790,
101	0	0	0	0	0	1	42416	-	00ABP/1	101	0	101	ENSG00000168038|ENSG00000066468|715832	835957	716269	758786	2	60,41,	0,60,	716269,758745,
101	0	0	0	0	0	1	42416	-	00ABP/1	101	0	101	ENSG00000168036|ENSG00000066468|65260	185385	65697	108214	2	60,41,	0,60,	65697,108173,
101	0	0	0	0	0	1	42416	-	00ABP/1	101	0	101	ENSG00000157593|ENSG00000066468|3795	123920	4232	46749	2	60,41,	0,60,	4232,46708,
101	0	0	0	0	0	1	42416	-	00ABP/1	101	0	101	ENSG00000137309|ENSG00000066468|9359	129484	9796	52313	2	60,41,	0,60,	9796,52272,
101	0	0	0	0	0	1	42416	-	00ABP/1	101	0	101	ENSG00000136238|ENSG00000066468|29455	149580	29892	72409	2	60,41,	0,60,	29892,72368,
101	0	0	0	0	0	1	42416	-	00ABP/1	101	0	101	ENSG00000120451|ENSG00000066468|41074	161199	41511	84028	2	60,41,	0,60,	41511,83987,


Other tags in the SAM files
------------------------------
AS:Alignment score generated by aligner. For example, if you use bowtie2, the score "can be negative. Can be greater than 0 in --local mode (but not in --end-to-end mode)."
NM: number of mismatches (Edit distance to the reference), including mismatches in xxI in CIGAR line, and mismatches in xxM of CIGAR line (which is in MD string). For example,
   if CIGAR=22M3I4M and MD=25T0, then NM=4; if CIGAR=6M1D23M, MD=1T4^T1C21, then NM=3.
CC: chromosome name of next alignment, '=' if on the same chr.
CP: start position of next alignment
HI: hit index (increasing from from 0 to NH-1). The best hit does not have to be the first one. See example below.

NM:i:1 indicates there's a mismatch in base-space, as does the MD:Z:4A45 which indicates that the mismatch occurs at the 5th position

$grep HWUSI-EAS1533_0026_FC:6:71:14711:9709#0 ~/scratch/mouse_adult_wt_smallRNAseq_76nt_strand_ox/accepted_hits.sam
HWUSI-EAS1533_0026_FC:6:71:14711:9709#0 256 chr1 18909796 0 20M2D9M * 0 0 TAGCCTCTGTCAGCACTCCTGAGTTCAGA B;8:98@@;@=+?=8BBD=3=B=CBDDD= AS:i:-16 XN:i:0 XM:i:1 XO:i:1 XG:i:2 NM:i:3 MD:Z:20^GG7A1 YT:Z:UU NH:i:10 CC:Z:chr10 CP:i:73888540 HI:i:0


Other examples:
---------------
Reference: ...CTTCTATTATCCTT...  M       =/X         MD
     Read:    CTTCTATTATCCTT     14M     14=         14       // example 1
     Read:    CTTATATTATCCTT     14M     3=1X10=     3C10     // example 2
     Read:    CTTATATTGGCCTT     14M     3=1X4=2X4=  3C4AT4
     Read:    CTTCTATTGGCCTT     14M     8=2X4=      8AT4
     Read:    TTTATATTATCCTG     14M     1X12=1X     0C12T0


UCSC PSL format

By ZHENGUO ZHANG on June 21, 2012 4:54 PM | 0 Comments | 0 TrackBacks
The PSL format is one important format used by the UCSC genome database. Although
the UCSC has provided the document for this format
http://genome.ucsc.edu/FAQ/FAQformat#format2,  I can put some more summaries
here to make use of the data in the UCSC database efficiently.

1. PSL format uses the 0-based right-half open coordinates, it means that all
   the starts (qStart, qStarts, tStart, tStarts) need plus 1 to get the normal
   coordinates (in the genome browser).

2. all the blockSizes, qStarts, tStarts follow the order of query sequence if
   the aligned region of query is from its plus strand or the order of reversed
   query sequence if aligned region is on query minus strand.

3. The qStart, qEnd, tStart, tEnd are all determined in the point of view of its
   plus strand regardless the aligned regions are on the minus or plus strand.

4. If the aligned region for target sequence is from minus strand, then tStarts
   is based on the reversed strand of the target, that is, the first base
   position starts from the last base of the plus strand. The coordinate can be
   converted with the following formula: pos_on_minus = tSize - pos_on_plus.

5. The point in 4 is also applied to the query sequence if the aligned region
   is on the minus strand of query.

References:
1. http://genome.ucsc.edu/FAQ/FAQformat#format2
2. https://lists.soe.ucsc.edu/pipermail/genome/2005-October/008731.html

"""

"""

Example of SAM input:
=====================

59f/1__00b	179	ENSG00000182944|ENSG00000151702|32517	19096	255	33M	=	153165	134116	AGCCAACAGAGCAGCAGCTACGGGCAGCAGAGT	C@<A@@@;..((>;.);;:-@EHE=:CD=3EC;	XA:i:0	MD:Z:33	NM:i:0
59f/1__00a	115	ENSG00000182944|ENSG00000151702|32517	153165	255	47M	=	19096	-134116	CCTCCCCTTGGAGGGGCACAAACGATCAGTAAGAATACAGAGCAACG	F@8<?)?:CC???9CGCAB;BHCF@AEHDAEA?E<<<CC4?ADD@?8	XA:i:0	MD:Z:47	NM:i:0

0K8/2__01a	67	ENSG00000182944|ENSG00000123268|32517	19066	255	65M	=	83092	64062	GGATCCTACAGCCAAGCTCCAAGTCAATATAGCCAACAGAGCAGCAGCTACGGGCAGCAGACTGC	??;BAB?AFHCFHBCE7AA:CHGIGIEHIGGFGHBFGGGCFFAFEHEC@HHD@AAGBHHB@>CCE	XA:i:2	MD:Z:61G2A0	NM:i:2
0K8/2__01b	131	ENSG00000182944|ENSG00000123268|32517	83092	255	36M	=	19066	-64062	ATCAGGAGATATGCAAACATATCAGATCCGAACTAC	HEEEB@CCC;A@AC@ACBB>CDCEDC:CC?BBBBCC	XA:i:0	MD:Z:36	NM:i:0

0K8/2__00a	67	ENSG00000182944|ENSG00000123268|32517	19066	255	61M	=	83088	64062	GGATCCTACAGCCAAGCTCCAAGTCAATATAGCCAACAGAGCAGCAGCTACGGGCAGCAGA	??;BAB?AFHCFHBCE7AA:CHGIGIEHIGGFGHBFGGGCFFAFEHEC@HHD@AAGBHHB@	XA:i:0	MD:Z:61	NM:i:0
0K8/2__00b	131	ENSG00000182944|ENSG00000123268|32517	83088	255	40M	=	19066	-64062	CTGCATCAGGAGATATGCAAACATATCAGATCCGAACTAC	>CCEHEEEB@CCC;A@AC@ACBB>CDCEDC:CC?BBBBCC	XA:i:0	MD:Z:40	NM:i:0


179 =
    read paired
    read mapped in proper pair
    read reverse strand
    mate reverse strand
    second in pair

115 =
    read paired
    read mapped in proper pair
    read reverse strand
    mate reverse strand
    first in pair

67 =
    read paired
    read mapped in proper pair
    first in pair

131 =
    read paired
    read mapped in proper pair
    second in pair

Example of SAM output:
=====================

59f/1__00b	179	ENSG00000182944|ENSG00000151702|32517	19096	255	47S33M	=	153165	134116	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaAGCCAACAGAGCAGCAGCTACGGGCAGCAGAGT	C@<A@@@;..((>;.);;:-@EHE=:CD=3EC;	XA:i:0	MD:Z:33	NM:i:0
59f/1__00a	115	ENSG00000182944|ENSG00000151702|32517	153165	255	47M33S	=	19096	-134116	CCTCCCCTTGGAGGGGCACAAACGATCAGTAAGAATACAGAGCAACGaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	F@8<?)?:CC???9CGCAB;BHCF@AEHDAEA?E<<<CC4?ADD@?8	XA:i:0	MD:Z:47	NM:i:0

0K8/2__01a	67	ENSG00000182944|ENSG00000123268|32517	19066	255	65M36S	=	83092	64062	GGATCCTACAGCCAAGCTCCAAGTCAATATAGCCAACAGAGCAGCAGCTACGGGCAGCAGACTGCaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	??;BAB?AFHCFHBCE7AA:CHGIGIEHIGGFGHBFGGGCFFAFEHEC@HHD@AAGBHHB@>CCE	XA:i:2	MD:Z:61G2A0	NM:i:2
0K8/2__01b	131	ENSG00000182944|ENSG00000123268|32517	83092	255	65S36M	=	19066	-64062	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaATCAGGAGATATGCAAACATATCAGATCCGAACTAC	HEEEB@CCC;A@AC@ACBB>CDCEDC:CC?BBBBCC	XA:i:0	MD:Z:36	NM:i:0


"""

import os
import sys
import optparse
import gc


cigar_set = (['M','I','D','N','S','H','P','=','X'])
# SAM columns
sam_QNAME = 0
sam_FLAG = 1
sam_RNAME = 2
sam_POS = 3
sam_MAPQ = 4
sam_CIGAR = 5
sam_RNEXT = 6
sam_PNEXT = 7
sam_TLEN = 8
sam_SEQ = 9
sam_QUAL = 10
sam_TAG = 11


dna = "a"*10000

############################
def merge_sam(file_in, file_ou):
    # It merges consecutive lines which represent one read

    fin = None
    if file_in == '-':
        fin = sys.stdin
    else:
        fin = open(file_in,'r')

    fou = None
    if file_ou == '-':
        fou = sys.stdout
    else:
        fou = open(file_ou,'w')

    data = []
    size = 10**6
    last = ['_']
    limit = 10**5
    while True:
        lines = fin.readlines(size)
        if not lines:
            break

        for line in lines:
            if line.startswith('@'):
                data.append(line)
            else:
                temp = line.split("\t")
                if last[sam_QNAME][:-1] != temp[sam_QNAME][:-1] or last[sam_QNAME][-1:] == temp[sam_QNAME][-1:]:
                    last = temp
                else:
                    sa = "-" if int(last[sam_FLAG]) & 0x10 else "+"
                    sb = "-" if int(temp[sam_FLAG]) & 0x10 else "+"
                    if sa != sb:
                        last = temp
                    else:
                        if temp[sam_QNAME][-1:] == "a":
                            (last,temp) = (temp,last)

                        la = len(last[sam_SEQ])
                        lb = len(temp[sam_SEQ])

                        if sa == "+":
                            last[sam_CIGAR] = last[sam_CIGAR] + str(lb)+"H"
                            temp[sam_CIGAR] = str(la)+"H" + temp[sam_CIGAR]
                        else:
                            last[sam_CIGAR] = str(lb)+"H" + last[sam_CIGAR]
                            temp[sam_CIGAR] = temp[sam_CIGAR] + str(la)+"H"

                        last[sam_QNAME] = last[sam_QNAME][:-1]
                        temp[sam_QNAME] = temp[sam_QNAME][:-1]

                        #last[sam_QNAME] = last[sam_QNAME].partition("__")[0]
                        #temp[sam_QNAME] = temp[sam_QNAME].partition("__")[0]


                        data.append("\t".join(last))
                        data.append("\t".join(temp))

                        if len(data) > limit:
                            fou.writelines(data)
                            data = []
                        last = ['_']
    if data:
        fou.writelines(data)
    fin.close()
    fou.close()


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It corrects the consecutives lines, which represent the same read, such that they look like read which has been split."""
    version = "%prog 0.11 beta"

    parser = optparse.OptionParser(usage = usage, description = description, version = version)

    parser.add_option("--input","-i",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file in SAM format.""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output file in SAM format.""")



    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)


    # running
    merge_sam(options.input_filename, options.output_filename)
    #
