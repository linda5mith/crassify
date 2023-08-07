# crassify
Pipeline to taxonomically classify microbiome viruses in metagenomic data based on ICTV's VMR.

## 1. Generate matches.tsv using diamond v2 [https://github.com/bbuchfink/diamond]
Diamond input: your translated, unknown viral proteins
Diamond blastp your reference proteins against ICTV reference DB
For symmetrical best hit approach use --query-cover 50 --subject-cover 50 parameters else run with as recommended below:
```diamond blastp --query proteins.faa --db databases/ictv_reference.dmnd -o matches.tsv --ultra-sensitive --no-self-hits --evalue 0.00001```

## 2. Installing crassify
Install dependencies using environment.yml file:
```conda env create -f environment.yml```
```conda activate crassify```
```git clone https://github.com/linda5mith/crassify.git```

Set path to crassify installation
```PATH=$PATH:/home/user/programs/crassify```

## 3. Running crassify
Sample folders with example inputs are located:
python crassify.py -p /proteomes -m /diamond_hits/matches.tsv -db /metadata/ICTV_metadata.csv

## 4. Crassify output
protein_hits.csv - returns most homologous ICTV viral protein hit for each of your viral proteins that aligned to reference database.

genome_hits.csv - returns top 5 genome hits for each proteome. 

dists.csv - returns distances between all proteomes and related ICTV viral genomes.