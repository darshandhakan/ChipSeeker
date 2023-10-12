# ChipSeeker Analysis ##
This analysis is used to identify the peaks and perform their annotations. It basically runs ChipSeeker tool on macs2 output and generates the plots.

The Steps to run this script are as follows:

-> singularity exec /storage/colddata/basesolve/tools/CHIPSEEKER/CHIPSEEKER_FINAL_JM.sif Rscript /storage/colddata/basesolve/tools/CHIPSEEKER/CHIPSEEKER.R -b macs2/ -s sample_table.txt -u TRUE

## The inputs are:
 -b macs2 output directory containing the individual folders and their outputs
 -s sample_metadata table containg the name of the sample and the macs2 output file name
 -u If the Human hg38 or hg19 needs to be used
