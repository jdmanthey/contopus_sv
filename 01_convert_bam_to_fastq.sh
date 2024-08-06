
bedtools bamtofastq -i m84046_240803_051617_s4.hifi_reads.bc2039.bam -fq C_sordidulus__KU39667.fastq

bedtools bamtofastq -i m84046_240803_051617_s4.hifi_reads.bc2040.bam -fq C_virens__KU39615.fastq

gzip C_sordidulus__KU39667.fastq

gzip C_virens__KU39615.fastq
