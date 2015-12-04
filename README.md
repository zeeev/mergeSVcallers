# mergeSVcallers
Creating an integrated SV is difficult.  The code associated with this project was designed to help merge SVs in consistent and straightforward way.  The inputs to mergeSVcallers are Tabix merged VCF files and the output is a merged VCF file.  MergeSVcallers can be re-run iteratively. 

## Downloading and installing:
```
git clone --recursive https://github.com/zeeev/mergeSVcallers.git
cd mergeSVcallers/
make
```

## Usage

```
Usage:
       mergeSVcallers -a ref.fasta -f a.vcf.gz,b.vcf.gz -t WHAM,LUMPY -s 500

 Required:
          -a - <STRING> - The samtools faidx indexed FASTA file
          -f - <STRING> - A comma separated list of tabix index VCF files
          -t - <STRING> - A comma seperated list of tags/identifiers for each file

 Optional:
          -s - <INT>   - Merge SVs with both breakpoints N BP away [100]
 Info:
          -This tool provides a simple set of operations to merge SVs.
          -Output is unsorted.
```
##TODO
- [ ] create a test suite
- [ ] merge by reciprocal overlap
- [ ] add a splitter function
- [ ] remap CN0->DEL
- [ ] add translocation functionality 
