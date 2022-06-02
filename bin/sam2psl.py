#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
It takes as input a file in SAM format and it converts into a PSL format file.



Author: Daniel Nicorici, Daniel.Nicorici@gmail.com

Copyright (c) 2009-2022 Daniel Nicorici

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


# info PSL
"""
========================================================
More about PSL format is here: http://genome.ucsc.edu/FAQ/FAQformat#format2


PSL format

PSL lines represent alignments, and are typically taken from files generated
by BLAT or psLayout. See the BLAT documentation for more details. All of the
following fields are required on each data line within a PSL file:

   1. matches - Number of bases that match that aren't repeats
   2. misMatches - Number of bases that don't match
   3. repMatches - Number of bases that match but are part of repeats
   4. nCount - Number of 'N' bases
   5. qNumInsert - Number of inserts in query
   6. qBaseInsert - Number of bases inserted in query
   7. tNumInsert - Number of inserts in target
   8. tBaseInsert - Number of bases inserted in target
   9. strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
  10. qName - Query sequence name
  11. qSize - Query sequence size
  12. qStart - Alignment start position in query
  13. qEnd - Alignment end position in query
  14. tName - Target sequence name
  15. tSize - Target sequence size
  16. tStart - Alignment start position in target
  17. tEnd - Alignment end position in target
  18. blockCount - Number of blocks in the alignment (a block contains no gaps)
  19. blockSizes - Comma-separated list of sizes of each block
  20. qStarts - Comma-separated list of starting positions of each block in query
  21. tStarts - Comma-separated list of starting positions of each block in target

Example:
Here is an example of an annotation track in PSL format. Note that line breaks have
been inserted into the PSL lines in this example for documentation display
purposes. Click here for a copy of this example that can be pasted into the
browser without editing.

track name=fishBlats description="Fish BLAT" useScore=1
59 9 0 0 1 823 1 96 +- FS_CONTIG_48080_1 1955 171 1062 chr22  47748585 13073589 13073753 2 48,20,  171,1042,  34674832,34674976,
59 7 0 0 1 55 1 55 +- FS_CONTIG_26780_1 2825 2456 2577 chr22  47748585 13073626 13073747 2 21,45,  2456,2532,  34674838,34674914,
59 7 0 0 1 55 1 55 -+ FS_CONTIG_26780_1 2825 2455 2676 chr22  47748585 13073727 13073848 2 45,21,  249,349,  13073727,13073827,

Be aware that the coordinates for a negative strand in a PSL line are handled in a special way.
In the qStart and qEnd fields, the coordinates indicate the position where the query matches
from the point of view of the forward strand, even when the match is on the reverse strand.
However, in the qStarts list, the coordinates are reversed.

Example:
Here is a 30-mer containing 2 blocks that align on the minus strand and 2 blocks that align on
the plus strand (this sometimes can happen in response to assembly errors):

0         1         2         3 tens position in query
0123456789012345678901234567890 ones position in query
            ++++          +++++ plus strand alignment on query
    --------    ----------      minus strand alignment on query
    12345678    1234567890

Plus strand:
     qStart=12
     qEnd=31
     blockSizes=4,5
     qStarts=12,26

Minus strand:
     qStart=4
     qEnd=26
     blockSizes=10,8
     qStarts=5,19

Essentially, the minus strand blockSizes and qStarts are what you would get if you
reverse-complemented the query. However, the qStart and qEnd are not reversed. To
convert one to the other:

     qStart = qSize - revQEnd
     qEnd = qSize - revQStart




Example:
Here is an example of an annotation track in PSL format. Note that line breaks have
been inserted into the PSL lines in this example for documentation display purposes.
This example can be pasted into the browser without editing.

browser position chr22:13073000-13074000
browser hide all
track name=fishBlats description="Fish BLAT" visibility=2
useScore=1
59 9 0 0 1 823 1 96 +- FS_CONTIG_48080_1 1955 171 1062 chr22  47748585 13073589 13073753 2 48,20,  171,1042,  34674832,34674976,
59 7 0 0 1 55 1 55 +- FS_CONTIG_26780_1 2825 2456 2577 chr22  47748585 13073626 13073747 2 21,45,  2456,2532,  34674838,34674914,
59 7 0 0 1 55 1 55 -+ FS_CONTIG_26780_1 2825 2455 2676 chr22  47748585 13073727 13073848 2 45,21,  249,349,  13073727,13073827,

Be aware that the coordinates for a negative strand in a PSL line are handled in a
special way. In the qStart and qEnd fields, the coordinates indicate the position
where the query matches from the point of view of the forward strand, even when the
match is on the reverse strand. However, in the qStarts list, the coordinates
are reversed.

Example:
Here is a 61-mer containing 2 blocks that align on the minus strand and 2 blocks
that align on the plus strand (this sometimes happens due to assembly errors):

0         1         2         3         4         5         6 tens position in query
0123456789012345678901234567890123456789012345678901234567890 ones position in query
                      ++++++++++++++                    +++++ plus strand alignment on query
    ------------------              --------------------      minus strand alignment on query
0987654321098765432109876543210987654321098765432109876543210 ones position in query negative strand coordinates
6         5         4         3         2         1         0 tens position in query negative strand coordinates
    123456789012345678              12345678901234567890

Plus strand:
     qStart=22
     qEnd=61
     blockSizes=14,5
     qStarts=22,56
     tStarts=122,156
     tSize=200

Minus strand:
     qStart=4
     qEnd=56
     blockSizes=20,18
     qStarts=5,39
     tStarts=45,79 (which is on positive strand 104,136)
     tSize=200

Essentially, the minus strand blockSizes and qStarts are what you would get if
you reverse-complemented the query. However, the qStart and qEnd are not
reversed. Use the following formulas to convert one to the other:
     Negative-strand-coordinate-qStart = qSize - qEnd   = 61 - 56 =  5
     Negative-strand-coordinate-qEnd   = qSize - qStart = 61 -  4 = 57


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

# PSL columns
psl_matches = 0
psl_misMatches = 1
psl_repMatches = 2
psl_nCount = 3
psl_qNumInsert = 4
psl_qBaseInsert = 5
psl_tNumInsert = 6
psl_tBaseInsert = 7
psl_strand = 8
psl_qName = 9
psl_qSize = 10
psl_qStart = 11
psl_qEnd = 12
psl_tName = 13
psl_tSize = 14
psl_tStart = 15
psl_tEnd = 16
psl_blockCount = 17
psl_blockSizes = 18
psl_qStarts = 19
psl_tStarts = 20
psl_seq = 21

#
psl_empty_line = [
'0',
'0',
'0',
'0',
'0',
'0',
'0',
'0',
'+',
's',
'0',
'0',
'0',
'r',
'0',
'0',
'0',
'0',
',',
',',
',']


#########################
def parse_cigar(c,toversion="1.3"):
    # parse CIGAR string
    # TOVERSION: if set to "1.3" the CIGAR is converted forcefully to CIGAR defined
    #          in SAM version 1.3 (i.e. neighbours "X" and "=" are joined and
    #          converted into one region of "M") if CIGAR is given as version 1.4
    #          If it is set to "1.4" not conversion to 1.3 is done is done!
    r = []
    d = ''
    mismatches_x = 0
    c = c.upper()
    for a in c:
        if a.isdigit():
            d = d + a
        elif a in cigar_set:
            dd = int(d)
            r.append((a,dd))
            if a == 'X':
                mismatches_x = mismatches_x + dd
            d = ''
        else:
            print >>sys.stderr,"ERROR: unknown CIGAR:",c
            sys.exit(1)
    if mismatches_x and toversion == '1.3':
        rr = []
        i = -1
        n = len(r)
        while True:
            i = i + 1
            if i == n:
                r = rr
                break
            elif r[i][0] in ('X', '=', 'M'):
                b = r[i][1]
                for j in xrange(i+1,n):
                    if r[j][0] in ('=','M','X'):
                        b = b + r[j][1]
                    else:
                        i = j - 1
                        break
                rr.append(('M',b))
            else:
                rr.append(r[i])
    return (r,mismatches_x)

#########################
def blocks(cigar, ig = 0, use_cigar_13 = True):
    # returns block of matches
    # input is from cigar()
    # NOTE: hard clipping is converted forecfully to soft clipping
    ir = 0 # index on read
    #ig = 0 # index on genome
    rr = [] # on read
    rg = [] # on genome
    match = 0
    mismatch = 0
    mismatch_x = 0
    mismatch_clip = 0
    insert_query = 0
    insert_query_count = 0
    insert_ref = 0
    insert_ref_count = 0
    seq_len = 0 # Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ
    #print >>sys.stderr,'cigar:',cigar
    (cig,mismatch_x) = parse_cigar(cigar, toversion = "1.3" if use_cigar_13 else "all")
    mismatch = mismatch_x
    #print >>sys.stderr,"parsed cigar:",cig
    for e in cig:
        if e[0] in ('S','H'):
            ir = ir + e[1] # read
            mismatch = mismatch + e[1]
            mismatch_clip = mismatch_clip + e[1]
            seq_len = seq_len + e[1]
        elif e[0] in ('I',):
            ir = ir + e[1] # read
            mismatch = mismatch + e[1]
            insert_query = insert_query + e[1]
            insert_query_count = insert_query_count + 1
            seq_len = seq_len + e[1]
        elif e[0] in ('X'):
            ir = ir + e[1] # read
            ig = ig + e[1] # reference/target seq
            mismatch = mismatch + e[1]
            mismatch_x = mismatch_x + e[1]
            seq_len = seq_len + e[1]
        elif e[0] in ('M','='):
            rr.append((ir,ir+e[1])) # read
            rg.append((ig,ig+e[1])) # reference/target seq
            ir = ir + e[1] # read
            ig = ig + e[1] # reference/target seq
            match = match + e[1]
            seq_len = seq_len + e[1]
        elif e[0] in ('D','N','P'):
            ig = ig + e[1] # reference/target seq
            insert_ref = insert_ref + e[1]
            insert_ref_count = insert_ref_count + 1
    #print >>sys.stderr,"cigar: query:",rr
    #print >>sys.stderr,"cigar: ref:",rg
    #print >>sys.stderr,"cigar: match:",match
    #print >>sys.stderr,"cigar: mismatch:",mismatch
    #print >>sys.stderr,"cigar: mismatch_clip:",mismatch_clip
    #print >>sys.stderr,"cigar: mismatch_x:",mismatch_x
    return (rr,rg,match,mismatch,mismatch_clip,mismatch_x,insert_ref,insert_ref_count,insert_query,insert_query_count,seq_len)


#########################
def get_psl(sam, lens, use_cigar_13=True , replace_string = '', read_sequence=False):
    # USE_CIGAR_13 - If True then the input CIGAR string is in format 1.4 then it will be converted into format 1.3
    #cig, qSize, tSize, tStart, strand):
    # returns PSL coordinates
    # input from blocks()
    #
    #  12. qStart - Alignment start position in query
    #  13. qEnd - Alignment end position in query
    #  18. blockCount - Number of blocks in the alignment (a block contains no gaps)
    #  19. blockSizes - Comma-separated list of sizes of each block
    #  20. qStarts - Comma-separated list of starting positions of each block in query
    #  15. tSize - Target sequence size
    #  16. tStart - Alignment start position in target
    #  17. tEnd - Alignment end position in target
    #  21. tStarts - Comma-separated list of starting positions of each block in target

    psl = None
    if sam and sam[sam_FLAG].isdigit():
        unmapped = True if int(sam[sam_FLAG]) & 0x4  else False
        if (not unmapped) and sam[sam_RNAME] != '*' and sam[sam_CIGAR] != '*' and sam[sam_QNAME] != '*':
            psl = psl_empty_line[:]

            # read sequence length
            psl[psl_tSize] = lens.get(sam[sam_RNAME],0)
            # reference name
            psl[psl_tName] = sam[sam_RNAME]
            # read name
            psl[psl_qName] = sam[sam_QNAME].replace(replace_string,'/',1) if replace_string else sam[sam_QNAME]

            # strand
            psl[psl_strand] = "-" if int(sam[sam_FLAG]) & 0x10  else '+'

            # start position
            psl[psl_tStart] = int(sam[sam_POS])-1

            (interval_query,interval_ref, match, mismatch, mismatch_clip, mismatch_x,insert_ref,insert_ref_count,insert_query,insert_query_count,seq_len) = blocks(sam[sam_CIGAR], ig = psl[psl_tStart], use_cigar_13 = use_cigar_13)

            # read sequence length
            if sam[sam_SEQ] != '*' and sam[sam_CIGAR].find('H') == -1:
                psl[psl_qSize] = len(sam[sam_SEQ])
            else:
                # • Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ
                psl[psl_qSize] = seq_len

            #
            psl[psl_qNumInsert] = insert_query_count
            psl[psl_qBaseInsert] = insert_query
            psl[psl_tNumInsert] = insert_ref_count
            psl[psl_tBaseInsert] = insert_ref

            # extract the mismatches from SAM (using tag NM:i)
            tag_nm_i = [e.partition("NM:i:")[2] for e in sam[sam_TAG:] if e.startswith('NM:i:')] # NM is mismatches per reads
            if not tag_nm_i:
                tag_nm_i = [e.partition("nM:i:")[2] for e in sam[sam_TAG:] if e.startswith('nM:i:')] # nM is not good because but is better than nothing because it is mismatches per fragment and not per read!
            tag_nm_i = int(tag_nm_i[0]) if tag_nm_i else 0
            if tag_nm_i > float(0.90)*seq_len:
                tag_nm_i = 0
            #print >>sys.stderr,"tag NM:i:",tag_nm_i

            # compute the matches and mismatches (include also the clipping as mismatches)
            mis = mismatch_clip + tag_nm_i
            #print >>sys.stderr,"mismatch_clip + tag_nm_i =",mis
            if mis >= mismatch:
                psl[psl_matches] = psl[psl_qSize] - mis
                psl[psl_misMatches] = mis
            else: # probably the tag NM:i are missing???!!!
                psl[psl_matches] = match
                psl[psl_misMatches] = mismatch

            if interval_query:
                psl[psl_qStart] = interval_query[0][0]
                psl[psl_qEnd] = interval_query[-1][1]
                #psl[tStart] = boxes[0][1][0]
                psl[psl_tEnd] = interval_ref[-1][1]
                psl[psl_blockCount] = len(interval_query)

                # this is how is the specification BUT BLAT does not follow the specification!!!
                # NOTE: BLAT _always_ gives the coordinates as everything is mapped on the forwward strand
                # even that it is mapped on the reverse strand
                #
                #if psl[psl_strand] == "+":
                #    psl[psl_blockSizes] = ','.join([str(e[1]-e[0]) for e in interval_query])+','
                #    psl[psl_qStarts] = ','.join([str(e[0]) for e in interval_query])+','
                #    psl[psl_tStarts] = ','.join([str(e[0]) for e in interval_ref])+','
                #elif psl[psl_strand] == "-":
                #    psl[psl_blockSizes] = ','.join([str(e[1]-e[0]) for e in interval_query[::-1]])+','
                #    psl[psl_qStarts] = ','.join([str(psl[psl_qSize]-e[1]) for e in interval_query[::-1]])+','
                #    psl[psl_tStarts] = ','.join([str(psl[psl_tSize]-e[0]-1) for e in interval_ref[::-1]])+','

                psl[psl_blockSizes] = ','.join([str(e[1]-e[0]) for e in interval_query])+','
                psl[psl_qStarts] = ','.join([str(e[0]) for e in interval_query])+','
                psl[psl_tStarts] = ','.join([str(e[0]) for e in interval_ref])+','

            if read_sequence:
                if len(psl) < psl_seq + 1:
                    psl.append(sam[sam_SEQ])

            psl = map(str,psl)

    return psl

#########################
def getlines(a_filename):
    # it gives chunks
    fin = None
    if a_filename == '-':
        fin = sys.stdin
    else:
        fin = open(a_filename,'r')
    header = dict()
    first = True
    while True:
        lines = fin.readlines(10**8)
        if not lines:
            break
        gc.disable()
        lines = [line.rstrip('\r\n').split('\t') for line in lines if line.rstrip('\r\n')]
        gc.enable()
        for line in lines:
            if line[0].startswith('@'):
                if line[0].startswith('@SQ') and line[1].startswith('SN:') and line[2].startswith('LN:'):
                    k = line[1][3:]
                    v = int(line[2][3:])
                    header[k] = v
                else:
                    pass
            else:
                if first:
                    first = False
                    yield header
                    header = None
                yield line
    if first and header:
        yield header
    fin.close()

############################
def sam2psl(file_in,file_ou, use_cigar_13 = True,replace_string = '',read_sequence=False):
    # It converts a SAM file to PSL file
    # USE_CIGAR_13 - If True then the input CIGAR string is in format 1.4 then it will be converted into format 1.3
    fou = None
    if file_ou == '-':
        fou = sys.stdout
    else:
        fou = open(file_ou,'w')

    # PSL data
    psl = []
    psl_empty_line = ['0']*21
    # processing
    i = 0
    size_lines = 10**6
    lengths = None
    for line in getlines(file_in):
        if i == 0:
            lengths = line
            i = i + 1
            continue
        i = i + 1
        temp = get_psl(line,lengths, use_cigar_13, replace_string, read_sequence)
        # saving
        if temp:
            psl.append('\t'.join(temp)+'\n')
            #print >>sys.stderr, '\t'.join(line)
            #print >>sys.stderr, '\t'.join(temp)
            #print >>sys.stderr, "-----------------------------------------------"
            if i > size_lines:
                fou.writelines(psl)
                psl = []
    if psl:
        fou.writelines(psl)
    fou.close()


if __name__ == '__main__':

    #command line parsing

    usage = "%prog [options]"
    description = """It takes as input a file in SAM format and it converts into a PSL format file."""
    version = "%prog 0.14 beta"

    parser = optparse.OptionParser(usage = usage, description = description, version = version)

    parser.add_option("--input","-i",
                      action="store",
                      type="string",
                      dest="input_filename",
                      help="""The input file in SAM format.""")

    parser.add_option("--skip-conversion-cigar-1.3","-4",
                      action = "store_true",
                      default = False,
                      dest = "skip_conversion_cigar_13",
                      help="By default if the CIGAR strings in the input SAM file are in the format defined in "+
                      "SAM version 1.4 (i.e. there are 'X' and '=') then the CIGAR string will be "+
                      "first converted into CIGAR string, which is described in SAM version 1.3, "+
                      "(i.e. there are no 'X' and '=' which are replaced with 'M') and afterwards into PSL format. "+
                      "Default is '%default'.")

    parser.add_option("--read-seq","-s",
                      action = "store_true",
                      default = False,
                      dest = "read_sequence",
                      help = """It adds to the PSL output as column 22, the sequence of the read. This is not anymore a valid PSL format.""")

    parser.add_option("--replace-read-ids","-r",
                      action = "store",
                      type = "string",
                      dest = "replace_reads_ids",
                      help = """In the reads ids (also known as query name in PSL) the string specified here will be replaced with '/' (which is used in Solexa for /1 and /2).""")

    parser.add_option("--output","-o",
                      action="store",
                      type="string",
                      dest="output_filename",
                      help="""The output file in PSL format.""")



    (options,args) = parser.parse_args()

    # validate options
    if not (options.input_filename and
            options.output_filename
            ):
        parser.print_help()
        sys.exit(1)

    t = options.replace_reads_ids if options.replace_reads_ids else ''

    # running
    sam2psl(
        options.input_filename,
        options.output_filename, 
        use_cigar_13 = (not options.skip_conversion_cigar_13),
        replace_string = t,
        read_sequence = options.read_sequence
        )
    #
