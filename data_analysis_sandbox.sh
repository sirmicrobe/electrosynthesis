-----------------------------------
#Using SPAdes assembly (Fangfang)
----------------------------------
wc -l ~/SPAdes/bw4g-feature.rna.cov.norm.txt
#88485 gene features in file

grep -c '>' bw4g.anno.feature.dna.fa
#151688

#need to add "cov_1" to header
awk '/>/{print $0("_cov_1")}!/>/' bw4g.anno.feature.dna.fa > bw4g.anno.feature.dna-mod_head.fa

sed -i 's/ /_/g' bw4g.anno.feature.dna-mod_head.fa && sed -i 's/\./_/g' bw4g.anno.feature.dna-mod_head.fa && sed -i 's/|/_/g' bw4g.anno.feature.dna-mod_head.fa

ruby /home/chrism/project/gscripts/fasta_contig_length-gc.rb /home/chrism/project/projects/cmarshall/metagenome/BW4gran/SPAdes/bw4g.anno.feature.dna-mod_head.fa > bw4g.anno.feature.dna-len.gc.fa

#Prodigal on contigs from spades
srun -p amd -N 1 --mem-per-cpu=8000 -c 32 -J prod -t 12:00:00 prodigal -i  /home/chrism/project/projects/cmarshall/metagenome/BW4gran/SPAdes/bw4g.anno.feature.dna-mod_head.fa -o spades-bw4g.genes -a spades-bw4g.genes.faa -d spades-bw4g.genes.fna -m -p meta 
#100.0% Searching, 94.2% matched

#ublast for spades
srun -p bigmem -N 1 --mem-per-cpu=16000 -c 16 -J ubl -t 36:00:00 /home/chrism/project/bin/usearch-64bit-sept2013/usearch64 -ublast spades-bw4g.genes.faa -db /home/chrism/project/databases/uniref90.udb -evalue 100 -accel 0.5 -maxhits 1 -blast6out spades-bw4g-uniref90.ublaste100.txt
#02:20:15 50.1Gb  100.0% Searching, 94.5% matched


#annotate blast output
#srun -p bigmem -N 1 --mem-per-cpu=16000 -c 16 -J ubl -t 36:00:00 ruby /home/chrism/project/gscripts/annot_blast_output.rb spades-bw4g-uniref90.ublaste100.txt /home/chrism/project/databases/uniref90.fasta spades-bw4g uniref -i bw4g.anno.feature.dna-len.gc.fa

#screen -r ubl2
srun -p amd -N 1 --mem-per-cpu=8000 -c 32 -J prod -t 12:00:00 ruby /home/chrism/project/gscripts/annot_blast_output.rb spades-bw4g-uniref90.ublaste100.txt /home/chrism/project/databases/uniref90.fasta spades-bw4g uniref -i bw4g.anno.feature.dna-len.gc.fa

wc -l annot.spades-bw4g-uniref90.ublaste100.txt
#142808 lines


#copy files to computer
scp -r chrism@midway.rcc.uchicago.edu:/home/chrism/project/projects/cmarshall/metagenome/BW4gran/SPAdes/annot.spades-bw4g-uniref90.ublaste100.txt /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/ARPA_Sequencing/metagenome/


#############################################
-------------------------------------
# ESOM steps - SPAdes
--------------------------------------
#############################################

#create file with two columns: 1st column contig name, 2nd column bin name

mkdir esom_spades

print_tetramer_freqs_deg_filter_esom.pl -s /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/BW4G.spades.all.contigs.fa -a /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/contig_bin_spades.txt -m 4000 -w 4000 -k 4

#make .cls file
#python ~/scripts/esom_make_cls_file.py -e <contig.names>  <2/3 column tab delimited bin file +/-colors> > col_filename.cls
python /home/chrism/project/gscripts/esom_make_cls_file.py -e /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/esom_spades/BW4G.spades.all.contigs.fa.names /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/esom_spades/contig_bin_spades.txt > spades-contig-bin-newcls.cls

# class file produced from BW4G.spades.all.contigs.fa.names and contig_bin_spades.txt on Tue Jul  8 10:49:59 2014
% 10809
% 0     NO_CLASS        255     255     255
% 1     bin     177     0       0
% 2     geo     232     0       0
% 3     sulfuro 255     40      0
% 4     sphaer  255     85      0
% 5     dsv-l   255     129     0
% 6     bacteroid-h     255     170     0
% 7     und-i   255     215     0
% 8     dsv-h   234     255     12
% 9     rhodo   199     255     47
% 10    aceto   160     255     86
% 11    bacteroid-l     121     255     125
% 12    methan  86      255     160
% 13    und-h   47      255     199
% 14    rhizo   12      244     234
% 15    phage   0       196     255
% 16    bacteroid       0       148     255
% 17    aceto-l 0       104     255
% 18    rhodo-l 0       56      255
% 19    phage-aceto     0       8       255
% 20    und     0       0       232
% 21    mobile  0       0       177
% 22    bacteriod-h     0       0       127



#Get coverage information
#ruby /home/chrism/project/gscripts/esom_cov_velvet.rb /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/BW4G.spades.all.contigs.fa /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/esom_spades/BW4G.spades.all.contigs.fa.names > spades-bw4g-names_cov.txt
#script didnt work properly because using spades instead of idba or velvet as input file
ruby /home/chrism/project/gscripts/esom_cov_idba.rb /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/BW4G.spades.all.contigs.fa /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/esom_spades/BW4G.spades.all.contigs.fa.names > spades-bw4g-names_cov.txt


#to start esom: Log in to nomachine and type the following:
/home/chrism/project/bin/esom/bin/esomana



#load lrn names and cls file
#did not normalize
#tools>training
#Training algorithm = K-batch
#for esom columns and rows: for 10809 data points (from cls file) =
# rows = (200*x)/12000
# columns = (328 * x) / 12000
# radius = (50 * x) / 12000
# rows:182 columns:298 radius:46 (added a couple of rows and columns to make the rows*columns >= five times total data points
#start value of radius = 50
#trainging epochs = 20

#after training:

#load .wts file
#load spades-contig-bin-newcls.cls


#projecting smaller fragments - 1kb
mkdir esom-1kb_projected

print_tetramer_freqs_deg_filter_esom.pl -s /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/BW4G.spades.all.contigs.fa -a /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/esom_spades/esom-projected/contig_bin_spades.txt -m 1000 -w 1000 -k 4

#make .cls file
#python ~/scripts/esom_make_cls_file.py -e <contig.names>  <2/3 column tab delimited bin file +/-colors> > col_filename.cls
python /home/chrism/project/gscripts/esom_make_cls_file.py -e /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/esom_spades/esom-1kb_projected/BW4G.spades.all.contigs.fa.proj.names /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/esom_spades/esom-1kb_projected/contig_bin_spades.txt > spades-bw4g-1kbproj-bin.cls
# class file produced from BW4G.spades.all.contigs.fa.proj.names and contig_bin_spades.txt on Thu Jul 17 11:55:51 2014
#% 61911


#Projecting

#moved .wts file into esom-projected folder with the .lrn file from the 1kb projection, renamed .wts file to match prefx of other projected files
#load .lrn
#did not do training
#load .wts
#tools>project

#load projected.bm file and 2kb.cls file
#draw class masks
#tools>classify
#fill in bin names of each class
#save new .cls file as sponge-new2kb-project.cls


#summarize the output of a new esom .cls file (after creating class masks). It will print a list of contigs, and the number of fragments classified by bin
ruby /home/chrism/project/gscripts/esom_merge_names-cls.rb	/home/chrism/scratch-midway/metagenome/BW4sup/SPAdes-bw4s/esom_spades/esom-1kb_projected/BW4G.spades.all.contigs.fa.proj.classmask.names /home/chrism/scratch-midway/metagenome/BW4sup/SPAdes-bw4s/esom_spades/esom-1kb_projected/BW4G.spades.all.contigs.fa.proj.classmask.cls > spades-bw4g-1kb-projected.txt

#not sure which names file is correct to use, testing both
ruby /home/chrism/project/gscripts/esom_merge_names-cls.rb	/home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/esom_spades/esom-1kb_projected/BW4G.spades.all.contigs.fa.proj.names /home/chrism/scratch-midway/metagenome/BW4gran/SPAdes/esom_spades/esom-1kb_projected/BW4G.spades.all.contigs.fa.proj.classmask.cls > spades-bw4g-1kb-projected2.txt
#both .names files give the same results



################################################
################################################

#checking mapping stats to SPAdes assembly, mine and Sebastien's

--------
#SPAdes
--------
#running SPAdes assembly
/space2/cmarshall/tools/SPAdes-3.1.0-Linux/bin/spades.py -o spades_bw4g --careful --pe1-1 /space2/cmarshall/work/metagenome/BW4G-DNA_S2_L001_R1_001.fastq --pe1-2 /space2/cmarshall/work/metagenome/BW4G-DNA_S2_L001_R2_001.fastq -t 8 -k 21,33,55,77,99,127

#checking spades assembly
python /space2/cmarshall/tools/quast-2.3/quast.py -o spades_quast_out -M 1000 -l spades21,spades33,spades77,spades55simple,spades55final,spades99,spades127,SPAdes /space2/cmarshall/work/metagenome/BW4granules/spades_bw4g/K21/simplified_contigs.fasta /space2/cmarshall/work/metagenome/BW4granules/spades_bw4g/K33/simplified_contigs.fasta /space2/cmarshall/work/metagenome/BW4granules/spades_bw4g/K77/simplified_contigs.fasta /space2/cmarshall/work/metagenome/BW4granules/spades_bw4g/K55/simplified_contigs.fasta /space2/cmarshall/work/metagenome/BW4granules/spades_bw4g/K55/final_contigs.fasta /space2/cmarshall/work/metagenome/BW4granules/spades_bw4g/K99/simplified_contigs.fasta /space2/cmarshall/work/metagenome/BW4granules/spades_bw4g/K127/final_contigs.fasta  /space2/cmarshall/work/metagenome/BW4granules/spades_bw4g/contigs.fasta

#All statistics are based on contigs of size >= 1000 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   spades21  spades33  spades77   spades55simple  spades55final  spades99  spades127  SPAdes  
# contigs (>= 0 bp)        127179    99775     81679      87428           87428          77050     68572      68379   
# contigs (>= 1000 bp)     22886     18199     15220      15710           15710          14598     13307      13307   
Total length (>= 0 bp)     92426844  99884209  100649740  100870837       100870837      99953391  98188301   98107633
Total length (>= 1000 bp)  50045623  62879433  66832472   65448803        65448803       67348247  67743412   67742829
# contigs                  22886     18199     15220      15710           15710          14598     13307      13307   
Largest contig             26318     141443    413444     257602          257602         413466    721097     721097  
Total length               50045623  62879433  66832472   65448803        65448803       67348247  67743412   67742829
GC (%)                     50.81     51.98     52.08      52.10           52.10          52.05     52.02      52.02   
N50                        2395      5452      11042      9023            9023           13313     21208      21208   
N75                        1517      2251      2830       2714            2714           2982      3375       3375    
L50                        5972      2569      1126       1461            1461           886       444        444     
L75                        12635     7208      4224       4867            4867           3716      2724       2724    
# N's per 100 kbp          0.00      0.00      0.00       0.00            0.00           0.00      0.00       0.01    
#best statistics: are final contigs.fa which is nearly identical to K127

--------
# Bowtie
--------

ruby /space2/cmarshall/gscripts/extract_contigs_by_length.rb /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/BW4G.spades.all.contigs.fa 1000 > spades_ff_bw4g_contigs_1kb.fa && grep '>' spades_ff_bw4g_contigs_1kb.fa | wc -l
#14227 contigs greater than 1kb

mkdir bowtie_spades_bw4g
cd bowtie_spades_bw4g
#fangfang assembly /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/BW4G.spades.all.contigs.fa 

#insert fragment length size and distribution
#creates an index file
bowtie2-build -f /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/spades_ff_bw4g_contigs_1kb.fa bowtie2-index-contigs

#Bowtie2 with >1kb reads
bowtie2 -q --phred33 --minins 0 --maxins 1500 -p 8 -x bowtie2-index-contigs -1 /space2/cmarshall/work/metagenome/BW4G-DNA_S2_L001_R1_001.fastq -2 /space2/cmarshall/work/metagenome/BW4G-DNA_S2_L001_R2_001.fastq -S bow2_spades-ff_bw4g.sam
#-------------------------------
4896577 reads; of these:
  4896577 (100.00%) were paired; of these:
    495850 (10.13%) aligned concordantly 0 times
    4203756 (85.85%) aligned concordantly exactly 1 time
    196971 (4.02%) aligned concordantly >1 times
    ----
    495850 pairs aligned concordantly 0 times; of these:
      59909 (12.08%) aligned discordantly 1 time
    ----
    435941 pairs aligned 0 times concordantly or discordantly; of these:
      871882 mates make up the pairs; of these:
        706134 (80.99%) aligned 0 times
        142544 (16.35%) aligned exactly 1 time
        23204 (2.66%) aligned >1 times
92.79% overall alignment rate
##93% of the raw metagenome reads mapped to the >1kb metagenome assembly
#---------------------------------

#now mapping the metatranscriptome to the SPAdes >1kb metagenome assembly
bowtie2 -q --phred33 --minins 0 --maxins 1500 -p 8 -x bowtie2-index-contigs -1 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.1 -2 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.2 -S bow2_spades-ff_bw4g-transcript.sam
#------------------------------------
14769238 reads; of these:
  14769238 (100.00%) were paired; of these:
    2398893 (16.24%) aligned concordantly 0 times
    11874536 (80.40%) aligned concordantly exactly 1 time
    495809 (3.36%) aligned concordantly >1 times
    ----
    2398893 pairs aligned concordantly 0 times; of these:
      181960 (7.59%) aligned discordantly 1 time
    ----
    2216933 pairs aligned 0 times concordantly or discordantly; of these:
      4433866 mates make up the pairs; of these:
        3905653 (88.09%) aligned 0 times
        450812 (10.17%) aligned exactly 1 time
        77401 (1.75%) aligned >1 times
86.78% overall alignment rate
#----------------------------------

#BWA 

#index -p output, -a is for database less than 2GB (bwtsw otherwise),
bwa index -p spades_bwa_index -a is /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/spades_ff_bw4g_contigs_1kb.fa

#bwa mem -t threads, -p interleaved pair, db_prefix, 
/space2/cmarshall/tools/bwa-0.7.9a/bwa mem -t 8 -p spades_bwa_index /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq > bwa-spades-bw4g-transcript-aln.sam

#convert sam to bam
samtools view -S bwa-spades-bw4g-transcript-aln.sam -b -o bwa-spades-bw4g-transcript-aln.bam

#get alignment stats from bwa
samtools flagstat bwa-spades-bw4g-transcript-aln.bam

#Alignment stats
#29892642 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
27473320 + 0 mapped (#91.91%:-nan%)
29892642 + 0 paired in sequencing
14945478 + 0 read1
14947164 + 0 read2
26830626 + 0 properly paired (89.76%:-nan%)
27409302 + 0 with itself and mate mapped
64018 + 0 singletons (0.21%:-nan%)
153020 + 0 with mate mapped to a different chr
117096 + 0 with mate mapped to a different chr (mapQ>=5)


#Dan assemblies
mkdir Dan_assemblies

#Dan's better assembly
/space2/cmarshall/work/metagenome/BW4granules/Dan_assemblies/BW4G-trimmed-assembly-single_reads.fa
#not quite as good
/space2/cmarshall/work/metagenome/BW4granules/Dan_assemblies/velvet_assembled_illumina_trimmed_contig-assembly.fa

ruby /space2/cmarshall/gscripts/extract_contigs_by_length.rb /space2/cmarshall/work/metagenome/BW4granules/Dan_assemblies/BW4G-trimmed-assembly-single_reads.fa 1000 > Dan1_bw4g_contigs_1kb.fa && grep '>' Dan1_bw4g_contigs_1kb.fa | wc -l
#1573 contigs greater than 1kb

mkdir bowtie_Dan1_bw4g
cd bowtie_Dan1_bw4g
#fangfang assembly /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/BW4G.spades.all.contigs.fa 

#insert fragment length size and distribution
#creates an index file
bowtie2-build -f /space2/cmarshall/work/metagenome/BW4granules/Dan_assemblies/Dan1_bw4g_contigs_1kb.fa bowtie2-index-contigs

#Bowtie2 with >1kb reads
bowtie2 -q --phred33 --minins 0 --maxins 1500 -p 8 -x bowtie2-index-contigs -1 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.1 -2 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.2 -S bow2_Dan1_bw4g.sam
-----------------------------
14769238 reads; of these:
  14769238 (100.00%) were paired; of these:
    9551313 (64.67%) aligned concordantly 0 times
    1713235 (11.60%) aligned concordantly exactly 1 time
    3504690 (23.73%) aligned concordantly >1 times
    ----
    9551313 pairs aligned concordantly 0 times; of these:
      23133 (0.24%) aligned discordantly 1 time
    ----
    9528180 pairs aligned 0 times concordantly or discordantly; of these:
      19056360 mates make up the pairs; of these:
        18632428 (97.78%) aligned 0 times
        196111 (1.03%) aligned exactly 1 time
        227821 (1.20%) aligned >1 times
36.92% overall alignment rate
------------------------------

#BWA 

#index -p output, -a <is> for database less than 2GB (bwtsw otherwise),
bwa index -p Dan1_bwa_index -a is /space2/cmarshall/work/metagenome/BW4granules/Dan_assemblies/Dan1_bw4g_contigs_1kb.fa

#bwa mem -t threads, -p interleaved pair, db_prefix, 
/space2/cmarshall/tools/bwa-0.7.9a/bwa mem -t 8 -p Dan1_bwa_index /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq > bwa-Dan1-bw4g-transcript-aln.sam

#convert sam to bam
samtools view -S bwa-Dan1-bw4g-transcript-aln.sam -b -o bwa-Dan1-bw4g-transcript-aln.bam

#get alignment stats from bwa
samtools flagstat bwa-Dan1-bw4g-transcript-aln.bam

#Alignment stats
#29695526 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
12100918 + 0 mapped (#40.75%:-nan%)
29695526 + 0 paired in sequencing
14847967 + 0 read1
14847559 + 0 read2
11719947 + 0 properly paired (39.47%:-nan%)
12038691 + 0 with itself and mate mapped
62227 + 0 singletons (0.21%:-nan%)
243664 + 0 with mate mapped to a different chr
30075 + 0 with mate mapped to a different chr (mapQ>=5)

#mapping Dan2 assembly
bowtie2-build -f /space2/cmarshall/work/metagenome/BW4granules/Dan_assemblies/velvet_assembled_illumina_trimmed_contig-assembly.fa bowtie2-index-contigs-dan2

bowtie2 -q --phred33 --minins 0 --maxins 1500 -p 8 -x bowtie2-index-contigs-dan2 -1 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.1 -2 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.2 -S bow2_Dan2_bw4g.sam
---------------------
14769238 reads; of these:
  14769238 (100.00%) were paired; of these:
    3841816 (26.01%) aligned concordantly 0 times
    9994983 (67.67%) aligned concordantly exactly 1 time
    932439 (6.31%) aligned concordantly >1 times
    ----
    3841816 pairs aligned concordantly 0 times; of these:
      146916 (3.82%) aligned discordantly 1 time
    ----
    3694900 pairs aligned 0 times concordantly or discordantly; of these:
      7389800 mates make up the pairs; of these:
        6768774 (91.60%) aligned 0 times
        542830 (7.35%) aligned exactly 1 time
        78196 (1.06%) aligned >1 times
77.08% overall alignment rate
----------------------
#forgot to cutoff above contigs at greater than 1kb


ruby /space2/cmarshall/gscripts/extract_contigs_by_length.rb /space2/cmarshall/work/metagenome/BW4granules/Dan_assemblies/velvet_assembled_illumina_trimmed_contig-assembly.fa 1000 > Dan2_bw4g_contigs_1kb.fa && grep '>' Dan2_bw4g_contigs_1kb.fa | wc -l
#2035

#BWA on Dan2 assembly

#index -p output, -a <is> for database less than 2GB (bwtsw otherwise),
bwa index -p Dan2_bwa_index -a is /space2/cmarshall/work/metagenome/BW4granules/Dan_assemblies/Dan2_bw4g_contigs_1kb.fa

#bwa mem -t threads, -p interleaved pair, db_prefix, 
/space2/cmarshall/tools/bwa-0.7.9a/bwa mem -t 8 -p Dan2_bwa_index /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq > bwa-Dan2-bw4g-transcript-aln.sam

#convert sam to bam
samtools view -S bwa-Dan2-bw4g-transcript-aln.sam -b -o bwa-Dan2-bw4g-transcript-aln.bam

#get alignment stats from bwa
samtools flagstat bwa-Dan2-bw4g-transcript-aln.bam

#Alignment stats
#29853540 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
23831502 + 0 mapped (#79.83%:-nan%)
29853540 + 0 paired in sequencing
14926103 + 0 read1
14927437 + 0 read2
23308592 + 0 properly paired (78.08%:-nan%)
23719153 + 0 with itself and mate mapped
112349 + 0 singletons (0.38%:-nan%)
71112 + 0 with mate mapped to a different chr
21335 + 0 with mate mapped to a different chr (mapQ>=5)

#Rename Dan's assemblies
mv BW4G-trimmed-assembly-single_reads.fa BW4G-trimmed-assembly-single_reads_Dan1.fa 

mv velvet_assembled_illumina_trimmed_contig-assembly.fa velvet_assembled_illumina_trimmed_contig-assembl_Dan2.fa 

#Final results:
Assembly    filename    bowtie2 bwa
Dan 1   (BW4G-trimmed-assembly-single_reads_Dan1.fa)    37% 41%
Dan 2   (velvet_assembled_illumina_trimmed_contig-assembl_Dan2.fa ) 77% 80%



#prodigal
srun -p bigmem -N 1 --mem-per-cpu=16000 -c 16 -J ubl -t 36:00:00 prodigal -i  /home/chrism/project/projects/cmarshall/metagenome/BW4gran/SPAdes/ BW4G.spades.all.contigs.fa -o spades-bw4g.genes -a spades-bw4g.genes.faa -d spades-bw4g.genes.fna -m -p meta -f gff

-----------------------------
#mannotator
/space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/mannotator
/space2/share/databases/Nr_mappings.txt
#mannotator -g rast.gff3,prodigal.gff3 -c contigs.fa -p <nr_path> -i <Nr_mappings_path>
mannotator -g bw4g.anno.gff -c bw4g.anno.fa -p /space2/share/databases/ncbi_nr -i /space2/share/databases/Nr_mappings.txt
#0 sequence(s) to BLAST in mannotator_86421415822075_unknowns.fna
#[blastall] FATAL ERROR: Input tag mismatch on Blast-def-line-set
#Error: Command 'blastall -p blastx -i mannotator_86421415822075_unknowns.fna -d /space2/share/databases/ncbi_nr -o mannotator_86421415822075_unknowns.blastx -m 8' failed with return status 16640

#try individual scripts to merge .gff files
scp scaffolds.sulfuro.fa scaffolds.dsv-h.fa scaffolds.aceto.fa scaffolds.bac-h.fa scaffolds.geo.fa scaffolds.rhodo.fa scaffolds.dsv-l.fa scaffolds.dsv-l-p.fa scaffolds.dsv-ul.fa scaffolds.sphaer.fa scaffolds.methan.fa scaffolds.rhizo.fa scaffolds.bac-l.fa scaffolds.und.fa cmarshall@sal:/space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin
scp sulfuro.93207.gff dsv-h.95583.gff aceto.95701.gff bacteroid-h.93188.gff geo.95582.gff rhodo.93205.gff dsv-l.93197.gff dsv-l-p.93195.gff dsv-ul.93201.gff sphaer.95581.gff methan.93203.gff rhizo.93204.gff bacteroid-l.93189.gff cmarshall@sal:/space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin

#smallbatch test
cat /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/scaffolds.sulfuro.fa /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/scaffolds.dsv-h.fa > sulfuro_dsv-h_combined.fa
combineGffOrfs -g /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/sulfuro.93207.gff,/space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/dsv-h.95583.gff -c sulfuro_dsv-h_combined.fa


combineGffOrfs -gffs|g gff_file1[,gff_file2[, ... ]] -contigs|c contigs_file
combineGffOrfs -gffs sulfuro.93207.gff,dsv-h.95583.gff,aceto.95701.gff,bacteroid-h.93188.gff,geo.95582.gff,rhodo.93205.gff,dsv-l.93197.gff,dsv-l-p.93195.gff,dsv-ul.93201.gff,sphaer.95581.gff,methan.93203.gff,rhizo.93204.gff,bacteroid-l.93189.gff -contigs 
#mannotator doesnt seem to work properly
-----------------------------------------------------------
#create a single fasta file
cat scaffolds.sulfuro.fa scaffolds.dsv-h.fa scaffolds.aceto.fa scaffolds.bac-h.fa scaffolds.geo.fa scaffolds.rhodo.fa scaffolds.dsv-l.fa scaffolds.dsv-l-p.fa scaffolds.dsv-ul.fa scaffolds.sphaer.fa scaffolds.methan.fa scaffolds.rhizo.fa scaffolds.bac-l.fa scaffolds.und.fa > all_RAST_scaffolds.fa
#remove first line form each .gff file
for i in *.gff
do sed -i -e "1d" $i
done
#create one gff file
cat sulfuro.93207.gff,dsv-h.95583.gff,aceto.95701.gff,bacteroid-h.93188.gff,geo.95582.gff,rhodo.93205.gff,dsv-l.93197.gff,dsv-l-p.93195.gff,dsv-ul.93201.gff,sphaer.95581.gff,methan.93203.gff,rhizo.93204.gff,bacteroid-l.93189.gff > all_RAST_annotat.gff
cat all_RAST_annotat.gff und.gff > all_annotated_rast_und.gff

--------------
#to get FPKM of dna (for RNA/DNA ratio analysis) need to first map reads to DNA
-------------

#index -p output, -a is for database less than 2GB (bwtsw otherwise),

cd /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin

#bwa index -p /space2/cmarshall/work/norman_demux/trim/bw4gran/spades-annot-map/map2rast/bw4all_RAST_genes -a is /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.fna

#bwa mem -t threads, -p interleaved pair, db_prefix, 
/space2/cmarshall/tools/bwa-0.7.9a/bwa mem -t 8 /space2/cmarshall/work/norman_demux/trim/bw4gran/spades-annot-map/map2rast/bw4all_RAST_genes /space2/cmarshall/work/metagenome/BW4granules/spades_cwm_bw4g/corrected/BW4G-DNA_S2_L001_R1_001.00.0_0.cor.fastq.gz /space2/cmarshall/work/metagenome/BW4granules/spades_cwm_bw4g/corrected/BW4G-DNA_S2_L001_R2_001.00.0_0.cor.fastq.gz > bw4gran_DNA_RAST_genes.sam


#convert sam to bam
samtools view -b bw4gran_DNA_RAST_genes.sam -o bw4gran_DNA_RAST_genes.bam -@ 8
samtools sort -T bw4gran_DNA_RAST_genes-pos-sort -o bw4gran_DNA_RAST_genes-pos-sort.bam -@ 8 bw4gran_DNA_RAST_genes.bam
samtools index bw4gran_DNA_RAST_genes-pos-sort.bam
samtools idxstats bw4gran_DNA_RAST_genes-pos-sort.bam > idxstat_bw4gran_DNA_RAST_genes-pos-sort.txt

samtools flagstat bw4gran_DNA_RAST_genes-pos-sort.bam > bw4gran_DNA_RAST_genes-pos-sort.log
10211092 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
449780 + 0 supplimentary
0 + 0 duplicates
9196673 + 0 mapped (90.07%:-nan%)
9761312 + 0 paired in sequencing
4880656 + 0 read1
4880656 + 0 read2
7625136 + 0 properly paired (78.12%:-nan%)
8400690 + 0 with itself and mate mapped
346203 + 0 singletons (3.55%:-nan%)
772720 + 0 with mate mapped to a different chr
685233 + 0 with mate mapped to a different chr (mapQ>=5)

scp cmarshall@sal:/space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/idxstat_bw4gran_DNA_RAST_genes-pos-sort.txt /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/ARPA_Sequencing/metatranscriptome/


#make a blast database of all annotated sequences
cd /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin//blastdb/

makeblastdb -in all_annotated_rast_und.fna -out bw -dbtype nucl

blastn -task blastn -db bw -query pceA.fa -out pceA.results.out  
blastn -task blastn -db bw -query hyf.fa -out hyf.results.out  
blastn -task blastn -db bw -query pfor_porA.fa -out pfor.results.out
blastn -task blastn -db bw -query pfor_porB.fa -out pforB.results.out
time blastn -task blastn -db bw -query atp_citrate_lyase_eproteobac.fa -out acl.results.out
time blastn -task blastn -db bw -query bcd.fa -out bcd.results.out
blastn -task blastn -db bw -query pha_synthase.fa -out pha.results.out 
#look at 16S qiime sequences matching to metagenome assembly
grep -A 1 '4supBacDNA' rep_set.fna > 4sup_rep_set.fasta
time blastn -task blastn -db bw -query 4sup_rep_set.fasta -out CCs.qiime16S.results.out -num_alignments 0 -num_descriptions 3
grep -A 1 '4bsupBacDNA' rep_set.fna > 4bsup_rep_set.fasta
time blastn -task blastn -db bw -query 4bsup_rep_set.fasta -out OCs.qiime16S.results.out -num_alignments 0 -num_descriptions 3
grep -A 1 '4bgranBacDNA' rep_set.fna > 4bgran_rep_set.fasta
time blastn -task blastn -db bw -query 4bgran_rep_set.fasta -out OCc.qiime16S.results.out -num_alignments 0 -num_descriptions 3
grep -A 1 '4granBacDNA' rep_set.fna > 4gran_rep_set.fasta
time blastn -task blastn -db bw -query 4gran_rep_set.fasta -out CCc.qiime16S.results.out -num_alignments 0 -num_descriptions 3

#top 20 otus in the list
#filter_fasta.py -f rep_set.fna -o rep_set_top20otu.fasta --sample_id_fp subsample_list_OC_CC.txt
for i in $(cat subsample_list_OC_CC.txt); do grep -A 1 "$i" rep_set.fna; done > rep_set_top20otu.fasta
time blastn -task blastn -db bw -query rep_set_top20otu.fasta -out top20qiime16S.results.out -num_alignments 0 -num_descriptions 3



#annotate with blast

# Testing out new binning software
srun -p bigmem -N 1 --mem-per-cpu=16000 -c 16 -J bin -t 36:00:00 time samtools view -b bow2_spades-ff_bw4g.sam -o bow2_spades-ff_bw4g.bam  -@ 8

#metabat binning tool installation
cd /home/chrism/project/build/berkeleylab-metabat-3a928aaf8e0a
#had to install the static binary in one folder + the most recent build containing the src/ and test/ folders and copy them into the berkeley folder
#also had to create a bam/ folder in the samtools-1.1 directory containing sam.h and bam.h
scons PREFIX=/home/chrism/project/bin/metabat/ SAMTOOLS_DIR=/home/chrism/project/build/samtools-1.1 HTSLIB_DIR=/home/chrism/project/build/samtools-1.1/htslib-1.2.1/htslib install

# --verysensitive = larger bins, more contamination, better for simpler communities
# --veryspecific = smaller bins, no contamination, few complete genomes
srun -p bigmem -N 1 --mem-per-cpu=16000 -c 16 -J bin -t 36:00:00 time /home/chrism/project/bin/metabat/bin/runMetaBat.sh /home/chrism/project/projects/cmarshall/spades_fangfang_bw4g/spades_ff_bw4g_contigs_1kb.fa /home/chrism/project/projects/cmarshall/spades_fangfang_bw4g/bowtie_spades_bw4g/bow2_spades-ff_bw4g.bam
#Number of clusters formed: 41 (bins)
# 30:27.81elapsed 101%CPU 

# taking the 41 bins and annotating
cat spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.1.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.2.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.3.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.4.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.5.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.6.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.7.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.8.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.9.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.10.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.11.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.12.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.13.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.14.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.15.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.16.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.17.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.18.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.19.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.20.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.21.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.22.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.23.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.24.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.25.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.26.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.27.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.28.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.29.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.30.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.31.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.32.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.33.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.34.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.35.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.36.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.37.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.38.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.39.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.40.fa spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.41.fa > metabat_41-bins_spades_ff_1kb.fa

#Prodigal on 41 metabat bins
srun -p bigmem -N 1 --mem-per-cpu=16000 -c 16 -J bin -t 36:00:00 time prodigal -i  metabat_41-bins_spades_ff_1kb.fa -o metabat_41-bins_spades_ff_1kb.genes -a metabat_41-bins_spades_ff_1kb.genes.faa -d metabat_41-bins_spades_ff_1kb.fna -m -p meta
#7:22.92elapsed

time ruby /home/chrism/project/gscripts/fasta_contig_length-gc.rb metabat_41-bins_spades_ff_1kb.fa > metabat_41-bins_spades_ff_1kb-len.gc.fa


#ublast for spades
srun -p bigmem -N 1 --mem-per-cpu=16000 -c 16 -J ubl -t 36:00:00 time /home/chrism/project/bin/usearch-64bit-sept2013/usearch64 -ublast metabat_41-bins_spades_ff_1kb.genes.faa -db /home/chrism/project/databases/uniref90.udb -evalue 100 -accel 0.5 -maxhits 1 -blast6out metabat_41-bins_spades_ff_1kb_bw4g-uniref90.ublaste100.txt
#01:07:41 45.7Gb  100.0% Searching, 97.4% matched


#annotate blast output
srun -p bigmem -N 1 --mem-per-cpu=16000 -c 16 -J ubl -t 36:00:00 time ruby /home/chrism/project/gscripts/annot_blast_output.rb metabat_41-bins_spades_ff_1kb_bw4g-uniref90.ublaste100.txt /home/chrism/project/databases/uniref90.fasta metabat_41-bins_spades_ff_1k uniref -i metabat_41-bins_spades_ff_1kb-len.gc.fa
#1154.37user 7.68system 19:23.78elapsed 99%CPU

#First line of each bin
for i in spades_ff_bw4g_contigs_1kb.fa.metabat-bins-.{1..41}.fa; do head -n 1 "$i" >> bin_line.txt; done

#running metaBAT again in --verysensitive mode
srun -p bigmem -N 1 --mem-per-cpu=16000 -c 16 -J bin -t 36:00:00 time /home/chrism/project/bin/metabat/bin/metabat --verysensitive -i /home/chrism/project/projects/cmarshall/spades_fangfang_bw4g/spades_ff_bw4g_contigs_1kb.fa -a spades_ff_bw4g_contigs_1kb.fa.depth.txt -o v_sens_bin


#BLASTING the metagenome to refseq  -outfmt 6 gives tabular file -num_threads increases processors 
cd /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/

#cd /space2/cmarshall/bin/ncbi-blast-2.2.29+/db
#makeblastdb –in nr –dbtype nucl –parse_seqids 

#test examples - with one gene
#blastx -db nr -query pfor_porA.fa -outfmt '6 qseqid sseqid sscinames pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -word_size 3 -num_threads 8 -out blastx-pfor_porA.txt -evalue 1e-3 -num_alignments 1 -culling_limit 1

#This took a REALLY long time - weeks!!
blastx -db nr -query all_annotated_rast_und.fna -outfmt '6 qseqid sseqid sscinames pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -word_size 3 -num_threads 8 -out blastx-bw4-RAST-metagenome.txt -evalue 1e-3 -num_alignments 1 -culling_limit 1

#curating data and putting it all into two folders for manuscript
#metagenome folder
mv blastdb/ /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/


#folder contains:
spades assembly over 1kb
EMIRGE files for CCc (BW4g) and CCs (BW4s)
#in diff_cov_bin subfolder:
binID.rastID.gff and fna = final binned annotation file and gene sequence file
bin_contigs = re-assembled scaffolds for each organism bin
blastx annotation



####################
#MUMmer - aceto
cd /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/ARPA_Sequencing/metagenome/diff_cov_bin/aceto/mummer_aceto

#aceto contigs /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/ARPA_Sequencing/metagenome/diff_cov_bin/diff_cov_bin/binned_contigs_for_RAST/scaffolds.aceto.fa
#aceto woodii fasta /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/Sequence_Scripts/ref_genomes/Acetobacterium_woodii_dsm_1030.GCA_000247605.1.28.dna.chromosome.Chromosome.fa

#align draft to reference
#promer - somewhat divergent sequences
#-p output prefix ref.fasta draft_contigs.fasta
promer -p promer_aceto /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/Sequence_Scripts/ref_genomes/Acetobacterium_woodii_dsm_1030.GCA_000247605.1.28.dna.chromosome.Chromosome.fa /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/ARPA_Sequencing/metagenome/diff_cov_bin/diff_cov_bin/binned_contigs_for_RAST/scaffolds.aceto.fa 

#view summary of alignment #-L minimum alignment length -I min % ID to display
show-coords -r -c -l -L 100 -I 50 promer_aceto.delta > promer_aceto.coords
show-aligns promer_aceto.delta "Acetobacterium_woodii_dsm_1030.GCA_000247605.1.28.dna.chromosome.Chromosome.fa" "scaffolds.aceto.fa" > promer_aceto.aligns

#run mapview - visualize alignment generated by promer
mapview -n 1 -p mapview promer_aceto.coords
xfig mapview_0.fig 
#to generate the above in pdf format: mapview -n 1 -f pdf -p mapview promer_aceto.coords

#run mummer 
mummer -mum -b -c /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/Sequence_Scripts/ref_genomes/Acetobacterium_woodii_dsm_1030.GCA_000247605.1.28.dna.chromosome.Chromosome.fa /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/ARPA_Sequencing/metagenome/diff_cov_bin/diff_cov_bin/binned_contigs_for_RAST/scaffolds.aceto.fa > mummer_aceto.mums
#view mummer output - cant get to work
mummerplot -t png -p mummer_aceto mummer_aceto.mums
mummerplot -t postscript -p mummer_aceto mummer_aceto.mums
mummerplot -t png -p promer_aceto promer_aceto.delta
#nucmer - closely related sequences
nucmer -maxmatch -p nucmer_aceto /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/Sequence_Scripts/ref_genomes/Acetobacterium_woodii_dsm_1030.GCA_000247605.1.28.dna.chromosome.Chromosome.fa /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/ARPA_Sequencing/metagenome/diff_cov_bin/diff_cov_bin/binned_contigs_for_RAST/scaffolds.aceto.fa

show-coords -r -c -l nucmer_aceto.delta > nucmer_aceto.coords
show-tiling nucmer_aceto.delta > nucmer_aceto.tiling
mapview -n 1 -p nuc_mapview nucmer_aceto.coords
xfig nuc_mapview_0.fig

delta-filter -q -r nucmer_aceto.delta > nuc_aceto.filter
mummerplot -t png -p nuc_aceto nuc_aceto.filter -R /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/Sequence_Scripts/ref_genomes/Acetobacterium_woodii_dsm_1030.GCA_000247605.1.28.dna.chromosome.Chromosome.fa -Q /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/ARPA_Sequencing/metagenome/diff_cov_bin/diff_cov_bin/binned_contigs_for_RAST/scaffolds.aceto.fa
##must delete comments in .gp referring to mouse
#nucmer did the best so far

show-tiling promer_aceto.delta > promer_aceto.tiling



# add characters at the beginning of each line using AWK
awk '{print "add_to_beginning"$0}' file
# add characters at the beginning of each line using SED
sed 's/^/add_to_beginning/' file
# add characters at the end of each line using AWK
awk '{print $0"append_to_end"}' file
# add characters at the end of each line using SED
sed 's/$/append_to_end/' file
# add /1 at end of each line beginning with @
awk '/^@/ {$0=$0"/1"}1' file
sed "/^@/ s/$/\/1/g" file



### Traitar = phenotype predictor
traitar phenotype /home/cwm47/electro sample_file_traitar.txt from_nucleotides /home/cwm47/electro/traitar_out -c 16

##################
#BW4 granules - metatranscriptome
#################

#159 BW4g metatranscriptome insert size

====================
#Quality filtering
====================


#adding /1 or /2 to reads so the below scripts can read pe
mv BW4gran_CAGATC_L005_R* trim/bw4gran/ 
ruby /space2/cmarshall/gscripts/fastq_add_pe_flag.rb BW4gran_CAGATC_L005_R1.fastq 1 && ruby /space2/cmarshall/gscripts/fastq_add_pe_flag.rb BW4gran_CAGATC_L005_R2.fastq 2

wc -l BW4gran_CAGATC_L005_R1.fastq
#109189060 #divided by 4 = 27297265

# Trimming sequences

java -jar /space2/cmarshall/bin/trimmomatic-0.27.jar PE r.BW4gran_CAGATC_L005_R1.fastq r.BW4gran_CAGATC_L005_R2.fastq bw4gran_r1 bw4gran_r1_se bw4gran_r2 bw4gran_r2_se ILLUMINACLIP:/space2/cmarshall/lib/illuminaClipping.fa:2:30:10
# Read Pairs: 27297265 Both Surviving: 26060368 (95.47%) Forward Only Surviving: 1236896 (4.53%) Reverse Only Surviving: 0 (0.00%) Dropped: 1 (0.00%)

wc -l bw4gran_r1_se
#4947584  
rm bw4gran_r2_se
wc -l bw4gran_r1
#104241472  #divide by 4 b/c each seq has 4 lines = reads
wc -l bw4gran_r2
#104,241,472   #divide by 4 b/c each seq has 4 lines = reads

/space2/cmarshall/gscripts/khmer/scripts/interleave_fastq.py -l bw4gran_r1 -r bw4gran_r2 -o bw4gran.pe.fq

wc -l bw4gran.pe.fq
#208482944 #  reads


# Quality filter paired end reads
fastq_quality_filter -Q33 -q 30 -p 75 -i bw4gran.pe.fq -o bw4gran.pe.qc.fq
wc -l bw4gran.pe.qc.fq
# 149369060  #divided by 4 = 37342265 (x100 =  billion bp)

fastq_quality_filter -Q33 -q 30 -p 75 -i bw4gran_r1_se -o bw4gran.se.qc
wc -l bw4gran.se.qc
# 4493208 # divided by 4 = 1123302 


#split out paired end reads and orphans
#/space2/cmarshall/tools/khmer/scripts/extract-paired-reads.py bw4bottle.pe.qc.fq
for i in *.pe.qc.fq
do
   /space2/cmarshall/tools/khmer/scripts/extract-paired-reads.py $i
done
#read 37342265 sequences, 14769238 pairs and 7803789 singletons

gzip bw4gran.pe.fq && gzip bw4gran_r1* && gzip bw4gran_r2 && gzip BW4gran_CAGATC_L005_R*


cat bw4gran.pe.qc.fq.se bw4gran.se.qc > bw4gran.se.qc.fq
wc -l bw4gran.se.qc.fq
# 35708364 #  reads

rm bw4gran.se.qc bw4gran.pe.qc.fq.se

#rename file
rm bw4gran.pe.qc.fq #remove because this file has singletons, is not just paired reads 
mv bw4gran.pe.qc.fq.pe bw4gran.pe.qc.fq

# Check quality - DID NOT RUN THIS for bw4gran
/space2/cmarshall/tools/FastQC/fastqc /space2/cmarshall/work/norman_demux//trim/bw4gran/bw4gran.pe.qc.fq

---------------
#Remove rRNA
----------------

# first pass to remove rRNA:
buildtrie --db /space2/share/databases/gg_13_5_otus/rep_set/85_otus.fasta

sortmerna --I /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq -n 7 --db /space2/share/databases/gg_13_5_otus/rep_set/97_otus.fasta /space2/share/databases/gg_13_5_otus/rep_set/85_otus.fasta /space2/cmarshall/tools/sortmerna-1.9/rRNA_databases/rfam-5s-database-id98.fasta /space2/cmarshall/tools/sortmerna-1.9/rRNA_databases/silva-bac-23s-database-id98.fasta /space2/cmarshall/tools/sortmerna-1.9/rRNA_databases/silva-arc-16s-database-id95.fasta /space2/cmarshall/tools/sortmerna-1.9/rRNA_databases/silva-arc-23s-database-id98.fasta /space2/cmarshall/tools/sortmerna-1.9/rRNA_databases/silva-bac-16s-database-id85.fasta --log bw4gran_log --accept rrna_bw4gran --other nonrrna_bw4gran --paired-in -m 67108864 -a 6
#----------------------------
##Results:
 total reads:   29538476
 non-rRNA:      24807340
 rRNA:          4731136
 % rRNA:        16.02%

/space2/share/databases/gg_13_5_otus/rep_set/97_otus.fasta              0.8255%
/space2/share/databases/gg_13_5_otus/rep_set/85_otus.fasta              0.01903%
/space2/cmarshall/tools/sortmerna-1.9/rRNA_databases/rfam-5s-database-id98.fasta                0.1089%
/space2/cmarshall/tools/sortmerna-1.9/rRNA_databases/silva-bac-23s-database-id98.fasta          2.416%
/space2/cmarshall/tools/sortmerna-1.9/rRNA_databases/silva-arc-16s-database-id95.fasta          0.3994%
/space2/cmarshall/tools/sortmerna-1.9/rRNA_databases/silva-arc-23s-database-id98.fasta          12.19%
/space2/cmarshall/tools/sortmerna-1.9/rRNA_databases/silva-bac-16s-database-id85.fasta          0.06055%
#--------------------------


fq2fa --paired /space2/cmarshall/work/norman_demux/trim/bw4gran/nonrrna_bw4gran.fq /space2/cmarshall/work/norman_demux/trim/bw4gran/nonrrna_bw4gran.fa && fq2fa bw4gran.se.qc.fq bw4gran.se.qc.fa

#add single end reads
cat nonrrna_bw4gran.fa bw4gran.se.qc.fa > bw4gran-singleandpaired.qc.fa
wc -l bw4gran-singleandpaired.qc.fa
#67,233,882  # divided by 2 = 33,616,941

==================
#BLAST
==================

#BLASTING the metagenome database  -outfmt 6 gives tabular file -num_threads increases processors 
blastn -db /space2/cmarshall/work/norman_demux/trim/bw4gran/blast/bw4g_db/bw4g_RAST_metagenome_header.fna  -task blastn -query /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran-singleandpaired.qc.fa -outfmt 6 -num_threads 8 -out blastn-bw4gran-RAST-metagenome.txt -evalue 1e-5 -num_alignments 1 -culling_limit 1
wc -l blastn-bw4gran-RAST-metagenome.txt
#11485881 lines



scp cmarshall@sal:/space2/cmarshall/work/norman_demux/trim/bw4gran/blastn-bw4gran-RAST-metagenome.txt  /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/ARPA_Sequencing/metatranscriptome/bw4gran


#columns in blast file: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore



/Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/ARPA_Sequencing/for_chrism_7apr_2014/transcript_blast_summary.rb

#Kim's script
#File preparation (file is assumed to by sorted by query and by highest to lowest bitscore):
#(a) Remove read duplicates (query column 1), then (b) sort by contig (db column 2): 
#awk 'a !~ $1; {a=$1}' file_in.txt | sort -V -k 2 > file_out.txt
awk 'a !~ $1; {a=$1}' blastn-bw4gran-RAST-metagenome.txt | sort -V -k 2 > sorted-blastn-bw4gran-RAST-metagenome.txt
#EXAMPLE:
#awk 'a !~ $1; {a=$1}' blastn-bw4g-RAST-metagenome.txt | sort -V -k 2 > sorted_blastn-bw4g-RAST-metagenome.txt


#ruby /home/chrism/project/gscripts/transcript_blast_summary.rb /home/chrism/project/projects/cmarshall/norman_demux/bw4bottle/sorted-blastn-bw4bottle-RAST-metagenome.txt


ruby /home/chrism/project/gscripts/transcript_blast_summary.rb sorted-blastn-bw4gran-RAST-metagenome.txt



#################
# Cufflinks
#################


mv bw4g_RAST_metagenome.fna bowtie2-bw4g-RAST-index.fa

#create bowtie index of metagenome
bowtie2-build -f bw4g_RAST_metagenome.fna bowtie2-bw4g-RAST-index

#extract paired end reads into r1 and r2 from quality filtered paired reads
python /space2/cmarshall/gscripts/khmer/sandbox/split-pe.py /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq
#DONE; split 29538476 sequences (14769238 left, 14769238 right)



# Tophat > Cufflinks > Cuffmerge > Cuffquant > Cuffnorm > Cuffdiff > CummeRbund

# Tophat mapping 
# -r inner distance between mate pair -o output directory <genome index base> <r1.fq> <r2.fq> 
/space2/cmarshall/bin/tophat-2.0.10.Linux_x86_64/tophat -r 50 -o tophat_bw4gran bowtie2-bw4g-RAST-index bw4gran.pe.qc.fq.1 bw4gran.pe.qc.fq.2
#tophat -r 50 -o tophat_brain /seqdata/indexes/hg19 brain_1.fq brain_2.fq tophat -r 50 -o tophat_liver /seqdata/indexes/hg19 liver_1.fq liver_2.fq tophat -r 50 -o tophat_heart /seqdata/indexes/hg19 heart_1.fq heart_2.fq

#tophat output to screen:
#left reads: min. length=101, max. length=101, 20957474 kept reads (967854 discarded)
#right reads: min. length=101, max. length=101, 20957474 kept reads (967854 discarded)
#A summary of the alignment counts can be found in tophat_bw4g/align_summary.txt
---------------------------------------
#Run complete: 05:30:50 elapsed
Left reads:
          Input     :  14769238
           Mapped   :   2966690 (20.1% of input)
            of these:      2259 ( 0.1%) have multiple alignments (0 have >20)
Right reads:
          Input     :  14769238
           Mapped   :   2951614 (20.0% of input)
            of these:      2204 ( 0.1%) have multiple alignments (0 have >20)
20.0% overall read mapping rate.

Aligned pairs:   2424795
     of these:      1712 ( 0.1%) have multiple alignments
                    2359 ( 0.1%) are discordant alignments
16.4% concordant pair alignment rate.
----------------------------------------

#concatenate all .gff files together to input in cufflinks -g flag
cat acetobacterium.62846.gff bacteroid1.62881.gff bacteroid2.62882.gff dsv1.62873.gff dsv2.62874.gff geo.62875.gff methan.62876.gff rhodo.62877.gff sphaer.62878.gff sulfuro.62879.gff > bw4g-RAST-annot.gff
#must first remove any line that contains trna gene
python remove.py
infile = bw4g-RAST-annot.gff
outfile = out.gff
python remove2.py
infile = out.gff
outfile = bw4g-RAST-annot-out.gff

#Cufflinks - assemble each sample independently

/space2/cmarshall/bin/cufflinks-2.2.0.Linux_x86_64/cufflinks -p 8 -g /space2/cmarshall/work/metagenome/BW4granules/extract_bins/bw4g-RAST-annot-out.gff -o cufflinks_bw4gran tophat_bw4gran/accepted_hits.bam 
#cufflinks -o cufflinks_brain tophat_brain/accepted_hits.bam -p 8
#cufflinks -o cufflinks_liver tophat_liver/accepted_hits.bam -p 8
#cufflinks -o cufflinks_heart tophat_liver/accepted_hits.bam -p 8
---------------------------------------
#Map Properties:
>       Normalized Map Mass: 3480105.50
>       Raw Map Mass: 3480105.50
>       Fragment Length Distribution: Empirical (learned)
>                     Estimated Mean: 170.57
>                  Estimated Std Dev: 53.29
[10:24:40] Assembling transcripts and estimating abundances.
> Processed 76639 loci
------------------------------------
#Merge the resulting assemblies - did not run
assemblies.txt:
    cufflinks_brain/transcripts.gtf
    cufflinks_liver/transcripts.gtf
    cufflinks_heart/transcripts.gtf
Now run the merge script:
cuffmerge -s /seqdata/fastafiles/hg19/hg19.fa assemblies.txt

#cuffquant
/space2/cmarshall/bin/cufflinks-2.2.0.Linux_x86_64/cuffquant -p 8 /space2/cmarshall/work/metagenome/BW4granules/extract_bins/bw4g-RAST-annot-out.gff tophat_bw4gran/accepted_hits.bam -o cuffquant-bw4gran  
--------------------------------
#Map Properties:
>       Normalized Map Mass: 0.00
>       Raw Map Mass: 0.00
>       Fragment Length Distribution: Truncated Gaussian (default)
>                     Default Mean: 200
>                  Default Std Dev: 80
[13:12:46] Calculating preliminary abundance estimates
[13:12:46] Quantifying expression levels in locus.
> Processed 28422 loci.
-------------------------------
#use output for cuffdiff

#cuffnorm
/space2/cmarshall/bin/cufflinks-2.2.0.Linux_x86_64/cuffnorm cufflinks_bw4gran/transcripts.gtf 
#cuffnorm transcripts.gff sample1.sam sample2.sam sample3.sam -o output_dir -p 8

#cuffdiff
#Take the annotated transcripts for your genome (as GFF or GTF) and provide them to cuffdiff along with the BAM files from TopHat for each replicate:
cuffdiff annotation.gtf mock_rep1.bam \ knockdown_rep1.bam

#cummeRbund


#compare known genome to merged assembly
cuffcompare -s /seqdata/fastafiles/hg19/hg19.fa -r known_annotation.gtf merged_asm/merged.gtf


#scripts path
/space2/cmarshall/bin/cufflinks-2.2.0.Linux_x86_64/cufflinks
/space2/cmarshall/bin/tophat-2.0.10.Linux_x86_64/tophat


##############
#Assemblies
##############

#velvet h
velveth bw4gran 41,73,2 -shortPaired -fastq /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq

#velvet g for haslengths = 41-73
for((n=41; n<=71; n=n+2)); do velvetg bw4gran_"$n" -read_trkg yes -amos_file yes; done

#bw4gran_41
Final graph has 787567 nodes and n50 of 110, max 3061, total 39669269, using 5657312/29538476 reads
#bwgran_57
Final graph has 245902 nodes and n50 of 208, max 5505, total 15931661, using 6496938/29538476 reads
#bwgran_61
Final graph has 176631 nodes and n50 of 258, max 6828, total 13094535, using 6617280/29538476 reads
#bwgran_65
Final graph has 129641 nodes and n50 of 327, max 7321, total 10832034, using 6760041/29538476 reads
#bw4gran_69
Final graph has 88537 nodes and n50 of 412, max 8442, total 8923729, using 6894585/29538476 reads
#bw4gran_71
Final graph has 75545 nodes and n50 of 447, max 8487, total 8157521, using 6959342/29538476 reads

##########
#Mapping
##########

#need to get insert size so mapping reads to contigs from bw4gran_59 directory
cd bw4gran_71
mkdir bowtie_contigs
cd bowtie_contigs
bowtie2-build -f /space2/cmarshall/work/norman_demux/trim/bw4gran/velvet/bw4gran_71/contigs.fa bowtie2-index-contigs



bowtie2 -q --phred33 --minins 0 --maxins 1000 -p 8 -x bowtie2-index-contigs -1 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.1 -2 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.2 -S bow2-bw4g71_r12.sam
================================
14769238 reads; of these:
  14769238 (100.00%) were paired; of these:
    8186616 (55.43%) aligned concordantly 0 times
    3413424 (23.11%) aligned concordantly exactly 1 time
    3169198 (21.46%) aligned concordantly >1 times
    ----
    8186616 pairs aligned concordantly 0 times; of these:
      207325 (2.53%) aligned discordantly 1 time
    ----
    7979291 pairs aligned 0 times concordantly or discordantly; of these:
      15958582 mates make up the pairs; of these:
        9128815 (57.20%) aligned 0 times
        1543296 (9.67%) aligned exactly 1 time
        5286471 (33.13%) aligned >1 times
69.10% overall alignment rate
================================

python /space2/cmarshall/gscripts/sam2insertsize.py bow2-bw4g71_r12.sam -v > bow2-bw4g_71-r12.txt
Total inserts read:..................... 20394075
First 10 insert sizes:.................. [  0   0   0   0 184 138 124   0   0   0]
Insert size mean:....................... 45.5
Insert size median:..................... 0
Insert stddev:.......................... 69.9
minimum insert size:.................... -1048
maximum insert size:.................... 1048
read length:............................ 101
pairs with insert size < read length:... 13855396 (67.9 %)
pairs with insert size < 2X read length: 19972041 (97.9 %)

#Do not know insert size so will be treating this data as single end reads from here on in
#Oases is then run on that directory, must include insert length -ins_length <int>:
oases bw4gran_79  

## And you can run Oases for repitition the same way:
for((n=59; n<=79; n=n+2)); do oases transcripts_"$n"; done


#now mapping the metatranscriptome to the SPAdes >1kb metagenome assembly
bowtie2 -q --phred33 --minins 0 --maxins 1500 -p 8 -x bowtie2-index-contigs -1 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.1 -2 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.2 -S bow2_spades-ff_bw4g-transcript.sam
#------------------------------------
14769238 reads; of these:
  14769238 (100.00%) were paired; of these:
    2398893 (16.24%) aligned concordantly 0 times
    11874536 (80.40%) aligned concordantly exactly 1 time
    495809 (3.36%) aligned concordantly >1 times
    ----
    2398893 pairs aligned concordantly 0 times; of these:
      181960 (7.59%) aligned discordantly 1 time
    ----
    2216933 pairs aligned 0 times concordantly or discordantly; of these:
      4433866 mates make up the pairs; of these:
        3905653 (88.09%) aligned 0 times
        450812 (10.17%) aligned exactly 1 time
        77401 (1.75%) aligned >1 times
86.78% overall alignment rate
#----------------------------------

#BWA 

#index -p output, -a is for database less than 2GB (bwtsw otherwise),
bwa index -p spades_bwa_index -a is /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/spades_ff_bw4g_contigs_1kb.fa

#bwa mem -t threads, -p interleaved pair, db_prefix, 
/space2/cmarshall/tools/bwa-0.7.9a/bwa mem -t 8 -p /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bwa_spades_bw4g/spades_bwa_index /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq > /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bwa_spades_bw4g/bwa-spades-bw4g-transcript-aln.sam

#convert sam to bam
samtools view -S bwa-spades-bw4g-transcript-aln.sam -b -o bwa-spades-bw4g-transcript-aln.bam

#get alignment stats from bwa
samtools flagstat bwa-spades-bw4g-transcript-aln.bam

#Alignment stats
#29892642 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
27473320 + 0 mapped (#91.91%:-nan%)
29892642 + 0 paired in sequencing
14945478 + 0 read1
14947164 + 0 read2
26830626 + 0 properly paired (89.76%:-nan%)
27409302 + 0 with itself and mate mapped
64018 + 0 singletons (0.21%:-nan%)
153020 + 0 with mate mapped to a different chr
117096 + 0 with mate mapped to a different chr (mapQ>=5)
----------------
----------------

#BWA - try 2 with annotated reads
#did not map well compared to full node mapping!!!!!

#index -p output, -a is for database less than 2GB (bwtsw otherwise),
bwa index -p spades_bwa_index2 -a is /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.feature.dna.fa

#bwa mem -t threads, -p interleaved pair, db_prefix, 
/space2/cmarshall/tools/bwa-0.7.9a/bwa mem -t 8 -p spades_bwa_index2 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq > bwa-spades-bw4g-transcript-aln2.sam

#convert sam to bam
samtools view -S bwa-spades-bw4g-transcript-aln2.sam -b -o bwa-spades-bw4g-transcript-aln2.bam

#get alignment stats from bwa
samtools flagstat bwa-spades-bw4g-transcript-aln2.bam

#29741646 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
9714736 + 0 mapped (32.66%:-nan%)
29741646 + 0 paired in sequencing
14868455 + 0 read1
14873191 + 0 read2
8988831 + 0 properly paired (30.22%:-nan%)
9378674 + 0 with itself and mate mapped
336062 + 0 singletons (1.13%:-nan%)
267932 + 0 with mate mapped to a different chr
256645 + 0 with mate mapped to a different chr (mapQ>=5)

#trying annotated mapping with bowtie to see if this works better than bwa
bowtie2-build -f /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.feature.dna.fa bowtie2_annotation_index

bowtie2 -q --phred33 --minins 0 --maxins 1500 -p 8 -x bowtie2_annotation_index -1 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.1 -2 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.2 -S bow2_spades-ff_bw4g-transcript-annot-aln.sam
-----------------------------
14769238 reads; of these:
  14769238 (100.00%) were paired; of these:
    11404848 (77.22%) aligned concordantly 0 times
    3243590 (21.96%) aligned concordantly exactly 1 time
    120800 (0.82%) aligned concordantly >1 times
    ----
    11404848 pairs aligned concordantly 0 times; of these:
      82882 (0.73%) aligned discordantly 1 time
    ----
    11321966 pairs aligned 0 times concordantly or discordantly; of these:
      22643932 mates make up the pairs; of these:
        21756539 (96.08%) aligned 0 times
        830394 (3.67%) aligned exactly 1 time
        56999 (0.25%) aligned >1 times
26.35% overall alignment rate
-------------------------------

#take home = do not use annotated file broken down by coding sequence

#looking at size of prodigal file versus assembly file
python /space2/cmarshall/gscripts/khmer/sandbox/assemstats.py 100 /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.feature.dna.fa /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/spades_ff_bw4g_contigs_1kb.fa
#------------------------------------------
filename sum n trim_n min med mean max n50 n50_len n90 n90_len
/space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.feature.dna.fa #sum:90489780 151688 146905 102 450 615 23760 34134 801 106552 327
/space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/spades_ff_bw4g_contigs_1kb.fa #sum:67472353 14227 14227 1000 1660 4742 574014 887 13129 8541 1446
#--------------------------------------------

#==================================================================
#creating small file to figure out mapping discrepency
fasta_formatter -i /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/spades_ff_bw4g_contigs_1kb.fa -o /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/spades_ff_bw4g_linewrap.fa
head -n 2 spades_ff_bw4g_linewrap.fa > spades_ff_bw4g_node1.fasta

prodigal -i  /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/spades_ff_bw4g_node1.fasta -o spades-bw4g-node1.genes -a spades-bw4g-node1.genes.faa -d spades-bw4g-node1.genes.fna -m -p meta 

#map transcript reads to prodigal and then to node 1
#index of node 1 contig
bowtie2-build -f /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test/spades_ff_bw4g_node1.fasta bowtie2_node1_index

#index of prodigal prediction file
bowtie2-build -f /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test/spades-bw4g-node1.genes.fna  bowtie2_node1_annotation_index

#mapping to node1 full contig
bowtie2 -q --phred33 --minins 0 --maxins 1500 -p 8 -x bowtie2_node1_index -1 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.1 -2 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.2 -S bow2_spades-ff_bw4g-node1-aln.sam --al-conc /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test/aligned-bowtie-node1.fa
-----------------------------------
14769238 reads; of these:
  14769238 (100.00%) were paired; of these:
    14767381 (99.99%) aligned concordantly 0 times
    1857 (0.01%) aligned concordantly exactly 1 time
    0 (0.00%) aligned concordantly >1 times
    ----
    14767381 pairs aligned concordantly 0 times; of these:
      27 (0.00%) aligned discordantly 1 time
    ----
    14767354 pairs aligned 0 times concordantly or discordantly; of these:
      29534708 mates make up the pairs; of these:
        29534648 (100.00%) aligned 0 times
        60 (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
0.01% overall alignment rate
----------------------------------
wc -l aligned-bowtie-node1.1.fa 
#7428 reads mapped in both forward and reverse aligned-bowtie-node1.1.fa


#mapping to prodigal prediction file
bowtie2 -q --phred33 --minins 0 --maxins 1500 -p 8 -x bowtie2_node1_annotation_index -1 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.1 -2 /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq.2 -S bow2_spades-ff_bw4g-node1-prodigal-aln.sam --al-conc /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test/aligned-prodigal-bowtie-node1.fa
--------------------------------------
14769238 reads; of these:
  14769238 (100.00%) were paired; of these:
    14767806 (99.99%) aligned concordantly 0 times
    1432 (0.01%) aligned concordantly exactly 1 time
    0 (0.00%) aligned concordantly >1 times
    ----
    14767806 pairs aligned concordantly 0 times; of these:
      23 (0.00%) aligned discordantly 1 time
    ----
    14767783 pairs aligned 0 times concordantly or discordantly; of these:
      29535566 mates make up the pairs; of these:
        29535320 (100.00%) aligned 0 times
        246 (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
0.01% overall alignment rate
------------------------------------

#repeat above steps with just the aligned files before moving on to the below steps

bowtie2 -q --phred33 --minins 0 --maxins 1500 -p 8 -x bowtie2_node1_index -1 /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test/aligned-bowtie-node1.1.fa  -2 /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test/aligned-bowtie-node1.2.fa -S nod1_aln-only.sam
---------------------------------------
1857 reads; of these:
  1857 (100.00%) were paired; of these:
    0 (0.00%) aligned concordantly 0 times
    1857 (100.00%) aligned concordantly exactly 1 time
    0 (0.00%) aligned concordantly >1 times
    ----
    0 pairs aligned concordantly 0 times; of these:
      0 (0.00%) aligned discordantly 1 time
    ----
    0 pairs aligned 0 times concordantly or discordantly; of these:
      0 mates make up the pairs; of these:
        0 (0.00%) aligned 0 times
        0 (0.00%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
100.00% overall alignment rate
----------------------------------
bowtie2 -q --phred33 --minins 0 --maxins 1500 -p 8 -x bowtie2_node1_annotation_index -1 /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test/aligned-bowtie-node1.1.fa  -2 /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test/aligned-bowtie-node1.2.fa -S nod1-prodigal_aln-only.sam
---------------------------------
1857 reads; of these:
  1857 (100.00%) were paired; of these:
    425 (22.89%) aligned concordantly 0 times
    1432 (77.11%) aligned concordantly exactly 1 time
    0 (0.00%) aligned concordantly >1 times
    ----
    425 pairs aligned concordantly 0 times; of these:
      4 (0.94%) aligned discordantly 1 time
    ----
    421 pairs aligned 0 times concordantly or discordantly; of these:
      842 mates make up the pairs; of these:
        642 (76.25%) aligned 0 times
        200 (23.75%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
82.71% overall alignment rate
----------------------------------

python /space2/cmarshall/gscripts/khmer/sandbox/assemstats.py 10 /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test/spades-bw4g-node1.genes.fna /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test/spades_ff_bw4g_node1.fasta
#-------------------------------
#filename sum n trim_n min med mean max n50 n50_len n90 n90_len
/space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test/spades-bw4g-node1.genes.fna #sum=528135 532 532 111 885 992 4533 147 1278 386 522
/space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test/spades_ff_bw4g_node1.fasta #sum=574014 1 1 574014 574014 574014 574014 1 574014 1 574014
#----------------------------------

#testing htseq on small subset of genes:
cd /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/small_file_mapping_test
htseq-count -r name -f bam nod1-prodigal_aln-only.sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.gff

htseq-count -r name -s reverse -i Name -t CDS -f bam nod1-prodigal_aln-only.sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.gff > count.txt

htseq-count -r name -i Name -t CDS nod1_aln-only.sam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.gff > count.txt
#this one might have worked:
1857 SAM alignment pairs processed.
__no_feature	963
__ambiguous	71
__too_low_aQual	12
__not_aligned	0
__alignment_not_unique	0

htseq-count -r name -i Name -t CDS nod1-prodigal_aln-only.sam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.gff > count.txt
__no_feature	1426
__ambiguous	0
__too_low_aQual	210
__not_aligned	221
__alignment_not_unique	0

#this below worked the best:
htseq-count -r name -i Name -s no -t CDS nod1_aln-only.sam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.gff > count.txt
__no_feature	50
__ambiguous	155
__too_low_aQual	12
__not_aligned	0
__alignment_not_unique	0

#change the gff file to number the genes 1-1xx,xxx
sed s/Name=/Name=_/ /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.gff > bw4g.anno_mod-temp.gff
awk -vRS=Name= '{$0=n$0;ORS=RT}++n' bw4g.anno_mod-temp.gff > bw4g.anno_mod.gff

htseq-count -r name -i Name -s no -t CDS nod1_aln-only.sam bw4g.anno_mod.gff > count.txt


#Text alignment viewer (based on the ncurses library). In the viewer, press ‘?’ for help and press ‘g’ to check the alignment start from a region in the format like ‘chr10:10,000,000’ or ‘=10,000,000’ when viewing the same reference sequence.
#view alignment for prodigal node1
samtools view -S nod1-prodigal_aln-only.sam -b -o nod1-prodigal_aln-only.bam
samtools sort -n nod1-prodigal_aln-only.bam nod1-prodigal_aln-only.name-sorted
samtools sort nod1-prodigal_aln-only.bam nod1-prodigal_aln-only.sorted


samtools tview nod1-prodigal_aln-only.name-sorted.bam spades-bw4g-node1.genes.fna
samtools tview nod1-prodigal_aln-only.sorted.bam bowtie2_node1_annotation_index.fa
#[bam_index_load] fail to load BAM index.




samtools tview bow2_spades-ff_bw4g-node1-prodigal-aln.sorted.bam aligned-prodigal-bowtie-node1.f

#samtools view -S bwa-spades-bw4g-aln2bw4s.sam -b -o bwa-spades-bw4g-aln2bw4s.bam && samtools flagstat bwa-spades-bw4g-aln2bw4s.bam
#output only aligned to separate sam file
=====================================================================================

#align to BW4s metagenome
#bwa mem -t threads, -p interleaved pair, db_prefix, 
/space2/cmarshall/tools/bwa-0.7.9a/bwa mem -t 8 -p /space2/cmarshall/work/metagenome/BW4supernatant/spades_bw4s/bwa_index/spades_bwa_index /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq > bwa-spades-bw4g-aln2bw4s.sam

samtools view -S bwa-spades-bw4g-aln2bw4s.sam -b -o bwa-spades-bw4g-aln2bw4s.bam && samtools flagstat bwa-spades-bw4g-aln2bw4s.bam

29855962 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
24392591 + 0 mapped (81.70%:-nan%)
29855962 + 0 paired in sequencing
14925312 + 0 read1
14930650 + 0 read2
23751322 + 0 properly paired (79.55%:-nan%)
24227955 + 0 with itself and mate mapped
164636 + 0 singletons (0.55%:-nan%)
90568 + 0 with mate mapped to a different chr
74179 + 0 with mate mapped to a different chr (mapQ>=5)










======================================  
#ht-seq for counting hits
-----------------------------
#for paired-end data, alignment must be sorted, can sort by name or alignment position using samtools sort function
samtools sort -n bwa-spades-bw4g-transcript-aln2.bam bwa-spades-bw4g-transcript-aln2_sorted.bam
#pay attention to stranded option with paired end: -s yes/no/reverse (try yes/default first but will most likely have to change it to -s no)
#htseq-count -r name -f bam/sam -s yes/no/reverse <alignment_file> <gff_file>
htseq-count -r name -f bam /space2/cmarshall/work/norman_demux/trim/bw4gran/spades-annot-map/bwa_map/bwa-spades-bw4g-transcript-aln2_sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.gff
# *aln2.bam file is the predicted gene file
----------------------------
#Warning: 190,592 reads with missing mate encountered.
14,972,390 SAM alignment pairs processed.
__no_feature    4,686,909
__ambiguous     0
__too_low_aQual 433,786
__not_aligned   9,851,695
__alignment_not_unique  0
----------------------------
samtools sort -n bwa-spades-bw4g-transcript-aln.bam bwa-spades-bw4g-transcript-aln-sorted.bam

htseq-count -r name -f bam bwa-spades-bw4g-transcript-aln-sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.gff
# *aln.bam is the total contig file
----------------------------
#Warning: 352,942 reads with missing mate encountered.
15,123,390 SAM alignment pairs processed.
__no_feature    13,575,097
__ambiguous     0
__too_low_aQual 370,043
__not_aligned   1,178,250
__alignment_not_unique  0
---------------------------
#trying with -s set to no
htseq-count -r name -f bam -s no /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bwa_spades_bw4g/bwa-spades-bw4g-transcript-aln-sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.gff > htseq_bw4g_transcript.txt
----------------------------
#Warning: 352942 reads with missing mate encountered.
15123390 SAM alignment pairs processed
__no_feature    13575097
__ambiguous     0
__too_low_aQual 370043
__not_aligned   1178250
__alignment_not_unique  0
----------------------------
htseq-count -r name -f bam -s reverse 
--------------------------



#GASSST - global alignment short sequence search tool
#For alignments with 90 % similarity minimum
#bank = DNA sequence of any length; query = very short DNA seqs
# -d bank.fasta -i query.fasta -o results -p 90
cd /space2/cmarshall/work/norman_demux/trim/bw4gran/spades-annot-map

Gassst -d /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.feature.dna.fa -i /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran-singleandpaired.qc.fa -o gassst_results -p 90 -n 8
#Forward and Reverse search :  9,438,061 Alignments found --> gassst_results

#Then call "gassst_to_sam" if you want to convert output to SAM format :
gassst_to_sam gassst_results gassst_results.sam && samtools view -bS gassst_results.sam > gassst_results.bam && samtools sort gassst_results.bam gassst_results_sorted

#using lower cutoff (80% matching)
Gassst -d /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/bw4g.anno.feature.dna.fa -i /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran-singleandpaired.qc.fa -o gassst_results2_80 -p 80 -n 8 #-s 5
#Forward and Reverse search :  10,468,758 Alignments found --> gassst_results2_80
#Execution time = 2 hour(s) 10 minute(s) 39 second(s)
#total = 33,616,941

#Then call "gassst_to_sam" if you want to convert output to SAM format :
gassst_to_sam gassst_results2_80 gassst_results2_80.sam && samtools view -bS gassst_results2_80.sam > gassst_results2_80.bam && samtools sort gassst_results2_80.bam gassst_results2_80_sorted


samtools index gassst_results2_80_sorted.bam
samtools idxstats gassst_results2_80_sorted.bam  >  gassst_result2_80.txt

#convert sam file into hits table
#sam to bam
samtools view -bS mated.sam > chris1.bam

./samtools sort chris1.bam chris2_sorted

./samtools index chris2_sorted

./samtools idxstats chris2_sorted.bam  >  result.txt

samtools sort bwa-spades-bw4g-transcript-aln2.bam bwa-spades-bw4g-transcript-aln2_sorted2
samtools index bwa-spades-bw4g-transcript-aln2_sorted2.bam && samtools idxstats bwa-spades-bw4g-transcript-aln2_sorted2.bam > bwa-spades-bw4g-transcript-aln2gene.txt


=========================
#mapping to RAST-annotated genomes and then creating transcript hit table
=========================

/space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_RAST_scaffolds.fa
/space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.gff

wc -l /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq
# 118,153,904 /4 = 29,538,476 /2 = 14,769,238

cd /space2/cmarshall/work/norman_demux/trim/bw4gran/spades-annot-map/map2rast
#BWA
#index -p output, -a is for database less than 2GB (bwtsw otherwise),
bwa index -p RAST_bwa_index -a is /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_RAST_scaffolds.fa

#bwa mem -t threads, -p interleaved pair, db_prefix, 
/space2/cmarshall/tools/bwa-0.7.9a/bwa mem -t 16 -p RAST_bwa_index /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq > bwa-RAST-bw4g-transcript-aln.sam

#convert sam to bam
samtools view -S bwa-RAST-bw4g-transcript-aln.sam -b -o bwa-RAST-bw4g-transcript-aln.bam

#get alignment stats from bwa
samtools flagstat bwa-RAST-bw4g-transcript-aln.bam
-------------------------------
29916265 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 duplicates
27751996 + 0 mapped (92.77%:-nan%)
29916265 + 0 paired in sequencing
14957013 + 0 read1
14959252 + 0 read2
27107274 + 0 properly paired (90.61%:-nan%)
27707239 + 0 with itself and mate mapped
44757 + 0 singletons (0.15%:-nan%)
200209 + 0 with mate mapped to a different chr
123594 + 0 with mate mapped to a different chr (mapQ>=5)
--------------------------------
#sort by name
samtools sort -n bwa-RAST-bw4g-transcript-aln.bam bwa-RAST-bw4g-transcript-aln-sorted

#htseq
htseq-count -f bam -r name -i ID -s no -t CDS bwa-RAST-bw4g-transcript-aln-sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.gff > count_bw4g_rast.txt
#Warning: 325046 reads with missing mate encountered.
#13,068,586 SAM alignment pairs processed.
__no_feature    7,430,834
__ambiguous     397,790
__too_low_aQual 597,831
__not_aligned   914,253
__alignment_not_unique  0
htseq-count -f bam -r name -i ID -s yes -t CDS bwa-RAST-bw4g-transcript-aln-sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.gff > count2_bw4g_rast.txt
#Warning: 325046 reads with missing mate encountered.
13068586 SAM alignment pairs processed.
__no_feature    9459210
__ambiguous     164164
__too_low_aQual 597831
__not_aligned   914253
__alignment_not_unique  0
htseq-count -f bam -r name -i ID -s reverse -t CDS bwa-RAST-bw4g-transcript-aln-sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.gff > count3_bw4g_rast.txt
#Warning: 325046 reads with missing mate encountered.
13068586 SAM alignment pairs processed.
__no_feature    9432486
__ambiguous     166413
__too_low_aQual 597831
__not_aligned   914253
__alignment_not_unique  0

htseq-count -f bam -r name -i ID -s no -t CDS bwa-RAST-bw4g-transcript-aln.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.gff > count4_bw4g_rast.txt
#Warning: 375707 reads with missing mate encountered.
15147014 SAM alignment pairs processed.
__no_feature    8613434
__ambiguous     461233
__too_low_aQual 693310
__not_aligned   1060784
__alignment_not_unique  0

htseq-count -f sam -r name -i ID -s no -t CDS bwa-RAST-bw4g-transcript-aln.sam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.gff > count5_bw4g_rast.txt
#Warning: 375707 reads with missing mate encountered.
15147014 SAM alignment pairs processed.
__no_feature    8613434
__ambiguous     461233
__too_low_aQual 693310
__not_aligned   1060784
__alignment_not_unique  0


#use new version of samtools:
samtools view -b bwa-RAST-bw4g-transcript-aln.sam -o bwa-RAST-bw4g-transcript-aln.bam -@ 8
samtools sort -n -T bwa-RAST-bw4g-transcript-aln-sorted -o bwa-RAST-bw4g-transcript-aln-sorted.bam -@ 8 bwa-RAST-bw4g-transcript-aln.bam

htseq-count -f bam -r name -i ID -s no -m intersection-nonempty -t CDS bwa-RAST-bw4g-transcript-aln-sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.gff > count6_bw4g_rast.txt
#Warning: Read HWI-ST222:238:D24U6ACXX:5:1101:1176:32467 claims to have an aligned mate which could not be found in an adjacent line.
#Warning: 375707 reads with missing mate encountered.
15147014 SAM alignment pairs processed.
__no_feature    9062269
__ambiguous     4
__too_low_aQual 693310
__not_aligned   1060784
__alignment_not_unique  0

htseq-count -f bam -r name -i ID -s no -m intersection-nonempty -t RNA bwa-RAST-bw4g-transcript-aln-sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.gff > count7_bw4g_rast.txt
#15147014 SAM alignment pairs processed.
__no_feature    13374383
__ambiguous     0
__too_low_aQual 693310
__not_aligned   1060784
__alignment_not_unique  0

htseq-count -f bam -r name -i ID -s no -m intersection-nonempty -t tRNA bwa-RAST-bw4g-transcript-aln-sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.gff > count8_bw4g_rast.txt
#Warning: 375707 reads with missing mate encountered.
15147014 SAM alignment pairs processed.
__no_feature    13374977
__ambiguous     0
__too_low_aQual 693310
__not_aligned   1060784
__alignment_not_unique  0

#comparing to see if samtools index is the same as htseq-count
samtools sort -T bwa-RAST-bw4g-transcript-aln-sorted -o bwa-RAST-bw4g-transcript-aln-sorted.bam -@ 8 bwa-RAST-bw4g-transcript-aln.bam
samtools index bwa-RAST-bw4g-transcript-aln-sorted.bam
samtools idxstats bwa-RAST-bw4g-transcript-aln-sorted.bam > idxstat_rast_bw4g_result.txt

#after indexing the sorted bam file
htseq-count -f bam -r pos -i ID -s no -m intersection-nonempty -t CDS bwa-RAST-bw4g-transcript-aln-sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.gff > count9_bw4g_rast.txt
#Warning: Mate records missing for 377788 records; first such record: <SAM_Alignment object: Paired-end read 'HWI-ST222:238:D24U6ACXX:5:1313:5844:61840' aligned to NODE_1_length_790048_cov_193.831_ID_1:[490245,490290)/->.
#Warning: Mate pairing was ambiguous for 11 records; mate key for first such record: ('HWI-ST222:238:D24U6ACXX:5:2201:3587:99922', 'second', 'NODE_1_length_790048_cov_193.831_ID_1', 490078, 'NODE_1_length_790048_cov_193.831_ID_1', 490154, 177).
15147027 SAM alignment pairs processed.
__no_feature    8819851
__ambiguous     4
__too_low_aQual 693312
__not_aligned   1060784
__alignment_not_unique  0

samtools sort -n -@ 8 -T bwa-RAST-bw4g-transcript-aln-name_sorted -o bwa-RAST-bw4g-transcript-aln-name_sorted.bam bwa-RAST-bw4g-transcript-aln.bam

htseq-count -f bam -r name -i ID -s no -m intersection-nonempty -t CDS bwa-RAST-bw4g-transcript-aln-name_sorted.bam /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.gff > count10_bw4g_rast.txt
#Warning: 375707 reads with missing mate encountered.
15147014 SAM alignment pairs processed.
__no_feature    8819848
__ambiguous     4
__too_low_aQual 693310
__not_aligned   1060784
__alignment_not_unique  0

######IMPORTANT#####
#the highest recovery of reads (htseq vs. idxstats) is with samtools idxstats
#with htseq - better to sort by position rather than name; and no diff between running on sam or bam file

scp cmarshall@sal:/space2/cmarshall/work/norman_demux/trim/bw4gran/spades-annot-map/map2rast/count6_bw4g_rast.txt cmarshall@sal:/space2/cmarshall/work/norman_demux/trim/bw4gran/spades-annot-map/map2rast/idxstat_rast_bw4g_result.txt /Users/chrismarshall/Documents/Chris_Argonne_Work_Mac/ARPA_Sequencing/metagenome/diff_cov_bin/diff_cov_bin/binned_contigs_for_RAST

rm -r map2rast #deleted whole folder and doing again:

#Using just idxstats on RAST fna files for mapping

cat sulfuro.93207.fna dsv-h.95583.fna aceto.95701.fna bacteroid-h.93188.fna geo.95582.fna rhodo.93205.fna dsv-l.93197.fna dsv-l-p.93195.fna dsv-ul.93201.fna sphaer.95581.fna methan.93203.fna rhizo.93204.fna bacteroid-l.93189.fna und.fna > all_annotated_rast_und.fna



#index -p output, -a is for database less than 2GB (bwtsw otherwise),
bwa index -p bw4all_RAST_genes -a is /space2/cmarshall/work/metagenome/BW4granules/spades_fangfang_bw4g/diff_cov_bin/all_annotated_rast_und.fna

#bwa mem -t threads, -p interleaved pair, db_prefix, 
/space2/cmarshall/tools/bwa-0.7.9a/bwa mem -t 16 -p bw4all_RAST_genes /space2/cmarshall/work/norman_demux/trim/bw4gran/bw4gran.pe.qc.fq > bw4all_RAST_genes.sam

#convert sam to bam
samtools view -b bw4all_RAST_genes.sam -o bw4all_RAST_genes.bam -@ 8
samtools sort -T bw4all_RAST_genes-pos-sort -o bw4all_RAST_genes-pos-sort.bam -@ 8 bw4all_RAST_genes.bam
samtools index bw4all_RAST_genes-pos-sort.bam
samtools idxstats bw4all_RAST_genes-pos-sort.bam > idxstat_bw4all_RAST_genes.txt

samtools flagstat bw4all_RAST_genes-pos-sort.bam 
29778162 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
239686 + 0 supplimentary
0 + 0 duplicates
14865500 + 0 mapped (49.92%:-nan%)
29538476 + 0 paired in sequencing
14769238 + 0 read1
14769238 + 0 read2
13824536 + 0 properly paired (46.80%:-nan%)
14236460 + 0 with itself and mate mapped
389354 + 0 singletons (1.32%:-nan%)
261450 + 0 with mate mapped to a different chr
219863 + 0 with mate mapped to a different chr (mapQ>=5)
