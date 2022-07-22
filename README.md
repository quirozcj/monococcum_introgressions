# Analysis to detect Triticum monococcum introgressions in domesticated hexaploid wheat.

Script used in the publication " Title of the publication " 

In this analysis we used the 10 pangenome chromosome level assemblies from the publication "10 pangenome paper" plus Chinese Sprig genome assembly.
As a query samples, we used XX accesions of XX.

We employed a k-mer based approach to detect introgressions:

- Create k-mer databases:
we use kmc-3.0.1: ```script run_kmc.sh```

```sh
kmc -k31 -ci1 -m30 -t5 -fm assembly_id.fa assembly_id out_dir
```
- For raw reads:
```sh
accesion_id='accesion1'
in_dir=../${accesion_id}

kmc -k31 -ci1 -m30 -t5 -fq <(ls -d $in_dir/*.gz) accesion_id out_dir
```

- Detect variations:
We employed IBSpy (IBSpy-0.3.1) to quantify variaitons per 50kb windows: ``` script run_IBSpy.sh ```
For details about how IBSpy detects variaitons, please, read the documentation [here](https://github.com/Uauy-Lab/IBSpy)

```sh
IBSpy --kmer_size 31 --window_size 50000 --reference reference1.fa --database accesion1 \
--database_format kmc3 --output accesion1_vs_reference1_50000.tsv --compress
```

- Combine IBSpy output tables by window & reference:
Use the script: ```run_combine.sh```
This is to process a single table in downstream analysis. The script requres a metadata file with the name of the individual IBSpy output files with their corresponding names. See the ```metadata.tsv``` example. All the  indiviual files will be combined by reference and arbitrary windows size.
```sh
python3 combine_by_windows.py -i metadata.tsv -r reference_name -w 50000 -s variaitons -o reference_combined_queries_50000.tsv.gz
```

