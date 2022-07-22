# Analysis to detect Triticum monococcum introgressions in domesticated hexaploid wheat.

Scripts used in the publication " Title of the publication " 

In this analysis we used the 10 pangenome chromosome level assemblies from the publication "10 pangenome paper" plus Chinese Sprig genome assembly.
As a query samples, we used XX accesions of XX.

We employed a k-mer based approach to detect introgressions:

1. Build k-mer databases.\
We used kmc-3.0.1
	* ```script: run_kmc.sh```
- For genome assembly:

```sh
#example for genome assembly
kmc -k31 -ci1 -m30 -t5 -fm assembly_id.fa assembly_id out_dir
```
- For raw reads:
```sh
#example for raw reads
accesion_id='accesion1'
in_dir=../${accesion_id}

kmc -k31 -ci1 -m30 -t5 -fq <(ls -d $in_dir/*.gz) accesion_id out_dir
```

2. Detect variations.
We employed IBSpy (IBSpy-0.3.1) to quantify variaitons per 50kb windows.
	* ``` script: run_IBSpy.sh ```\
For details about how IBSpy detects variaitons, please, read the documentation [here](https://github.com/Uauy-Lab/IBSpy)

```sh
#example
IBSpy --kmer_size 31 --window_size 50000 --reference reference1.fa --database accesion1 \
--database_format kmc3 --output accesion1_vs_reference1_50000.tsv --compress
```

3. Combine IBSpy output tables by window & reference.\
	* ```script: run_combine.sh```\
	This is to proccess a single table in downstream analysis. The script requres a metadata file with the name of the individual IBSpy output files with their corresponding names. See the ```metadata.tsv``` example. All the  indiviual files will be combined by reference and arbitrary windows size.

```sh
#example
python3 combine_by_windows.py -i metadata.tsv -r reference_name -w 50000 -s variaitons -o reference_combined_queries_50000.tsv.gz
```

The list below of the variations tables combined by chrosmosome & reference are public available ``` here(link)```

	* arinalrfor_combined_queries_50000w.tsv.gz
	* chinese_combined_queries_50000w.tsv.gz
	* jagger_combined_queries_50000w.tsv.gz
	* julius_combined_queries_50000w.tsv.gz
	* lancer_combined_queries_50000w.tsv.gz
	* landmark_combined_queries_50000w.tsv.gz
	* mace_combined_queries_50000w.tsv.gz
	* norin61_combined_queries_50000w.tsv.gz
	* spelta_combined_queries_50000w.tsv.gz
	* stanley_combined_queries_50000w.tsv.gz
	* sy_mattis_combined_queries_50000w.tsv.gz


4. Define monococcum introgressions.\
To proces downstream data and define introgressions we use python & R notebook scripts.
Please, see ``` monococcum_min_variations.ipynb ``` jupyther notebooks and R scripts for detailed description.

