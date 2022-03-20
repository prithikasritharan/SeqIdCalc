# SeqIdCalc

Calculates the sequence identity scores from both the CIGAR and FAT-CIGAR string within the BAM file. The sequence identity score for the CIGAR string is written to the ZA custom tag within the BAM file and the sequence identity scores for the FAT-CIGAR string is written to the ZB custom tag.     
&nbsp;

&nbsp;


## Usage

It takes in an input BAM file and outputs the BAM file containing the sequence identity scores for all mapped reads. 

The script should be run with the following commands:
```python
python seq_id_calc.py bam_file outbam_file
```
