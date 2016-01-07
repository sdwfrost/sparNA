# sparNA  

sparNA is a pipeline for assembling high depth paired-end viral sequencing reads. It is intended for use assembling highly diverse populations of viruses such as HIV and HCV.

Currently under active development. Do not expect things to work smoothly just yet. Feel free to get in touch via b|at|bede|dot|im

## Installing dependencies

### Mac OS X
*Using Homebrew and pip*  
- `brew tap homebrew/homebrew-science`
- `brew install python3 spades bowtie2 samtools vcftools bcftools seqtk`
- `pip3 install argh biopython khmer pandas plotly`  
*Finally, sign up for a plotly account and set an API key*
- `python -c "import plotly; plotly.tools.set_credentials_file(username='DemoAccount', api_key='lr1c37zw81')"`  

## Usage
N.B. Specify reads using absolute paths  
  
`./sparna.py --fwd-fq /path/to/reads_f.fastq --rev-fq /path/to/reads_r.fastq --norm-c-list 1,5,20 --norm-k-list 21,31 --asm-k-list 21,33,55,77 --min-len 1000 --threads 12`