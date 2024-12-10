#### Make VCF Pipeline
------
To run the pipeline for one script, edit the variables in the config.mk file. To run the pipeline for a set of samples, create a csv file following the format "Reference,Species(orsample),SRA". An example input file is included. To test the pipeline, run the following code:

```
cat input.csv | parallel --bar --eta --colsep , --header : make -f make_vcf/run.mk all ACC={Reference} SPECIES={Sample} SRA={SRA}
```

Note:
I've been able to run individual samples but have been encountering issues when trying to batch, the script is getting stuck on the snpEffect predictor. 