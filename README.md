# sparNA  

sparNA is a pipeline for assembling high depth paired-end viral sequencing reads. It is intended for use assembling highly diverse populations of viruses such as HIV and HCV.

Currently under active development. Do not expect things to work smoothly just yet, but I'm happy to help. b|at|bede|dot|im or find me on [Twitter](https://twitter.com/beconstant) 

## Dependencies  
Java, Python 3.5, spades, bowtie2, samtools, vcftools, bcftools, seqtk, argh, biopython, khmer, pandas, plotly  
Trimmomatic (jar file) is bundled inside `res/`  

### Mac OS X
Using Homebrew and pip is by far the easiest approach
- `brew tap homebrew/homebrew-science`
- `brew install python3 spades bowtie2 samtools vcftools bcftools seqtk`
- `pip3 install argh biopython khmer pandas plotly`  
Finally, sign up for a Plotly account and set your [API key](https://plot.ly/settings/api)  
- `python3 -c "import plotly; plotly.tools.set_credentials_file(username='USERNAME', api_key='API_KEY')"`  

## Usage
N.B. Specify reads using absolute paths  
Encoutering issues? Set the logging level to `INFO` for verbose output  

`./sparna.py --fwd-fq /path/to/reads_f.fastq --rev-fq /path/to/reads_r.fastq --norm-c-list 1,5,20 --norm-k-list 21,31 --asm-k-list 21,33,55,77 --min-len 1000 --threads 12`