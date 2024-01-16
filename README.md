# KG Vejle TumorBoard DNA script

## Generel info:
Requires a samplesheet containing 5 columns in specific order (tab separated), without headerline:
1) caseID, 2) NPN normal WES, 3) NPN tumor WES, 4) NPN tumor RNA, 5) PCGR tumor value

Example samplesheet:

    johnDoe 112217976652	111184925465	111184925473    23

The script will automatically look for fastq or cram files in subfolders at /lnx01_data2/shared/dataArchive/. This location contains read-only access to the data archive, containing all FastQ files. Theres no need to copy or move any input data.

The user can point to a specific folder containing input data using the --fastq or --cram option 

This is only needed if input data only exists outside the data archive (e.g. if data are in personal folders or e.g. stored at other KG Vejle servers).

## Usage examples:

Analyze the cases in samplesheet.txt, starting with CRAM files, let the script find the relevant files:

        nextflow run KGVejle/tumorboardDNA -r main --samplesheet /path/to/samplesheet.txt

Analyze the cases in samplesheet.txt, starting with CRAM files, use only files in specific folder:

        nextflow run KGVejle/tumorboardDNA -r main --samplesheet /path/to/samplesheet.txt --cram /path/to/cramfolder/

Analyze the cases in samplesheet.txt, starting with Fastq files, let the script find the relevant fastqfiles at the archive:

        nextflow run KGVejle/tumorboardDNA -r main --samplesheet /path/to/samplesheet.txt --fastqInput

Analyze the cases in samplesheet.txt, starting with Fastq files, use only files in the specific folder:

        nextflow run KGVejle/tumorboardDNA -r main --samplesheet /path/to/samplesheet.txt --fastq /path/to/fastqfolder/



## Main options:

    --help                print this help message
    
    --genome              hg19 or hg38
                              Default: hg38
  
    --outdir              Select which folder to write output to.
                              Default: TN_WES_results
  
    --samplesheet         path to case samplesheet. Can contain multiple patients/cases (each case in a separate line). 
  
    --server              Select which server the analysis is performed on (kga01 or lnx01)
                              Default: lnx01
  
    --fastq               Path to dir with fastq files
                              Default: data storage dirs at lnx01 server or kga01 server
  
    --cram                Path to dir with CRAM files
  
    --fastqInput          Use Fastq as input, search the archive location for the relevant files for the samples in the samplesheet

    --skipQC              Do not run QC module
