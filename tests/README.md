```bash
ilus WGS -c -C ilus_sys.yaml -L fastq.list -O ./tmp
```
```bash
ilus WGS -c -C ilus_sys.yaml -L fastq.list -O ./tmp -dr
```

```bash
ilus WGS -P gvcf,genotype -c -C ilus_sys.yaml -L fastq.list -O ./tmp

ilus WGS -P gvcf -c -C ilus_sys.yaml -L fastq.list -O ./tmp
ilus WGS -P genotype -c -C ilus_sys.yaml -L fastq.list -O ./tmp -f 
```

```bash
ilus genotype-joint-calling -C ilus_sys.yaml -L gvcf.list -O 03.genotype 
```

```bash
ilus genotype-joint-calling -C ilus_sys.yaml -L gvcf.list -O 03.genotype --as_pipe_shell_order -f
```

```bash
ilus VQSR -C ilus_sys.yaml -L vcf.list -O 03.genotype --as_pipe_shell_order -f 
```

```bash
python ../../scripts/yhbatch_slurm_jobs.py -I step1.bwa.sh -n 10 -t 3
```
python ../scripts/create_wgs_calling_intervals.py -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai -n GCA_000001405.15_GRCh38_no_alt_analysis_set.N.interval_list

python ../scripts/create_wgs_calling_intervals.py -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai -n GCA_000001405.15_GRCh38_no_alt_analysis_set.N.interval_list -w 5000000 | head -820 > wgs_calling_regions.GRCh38.5M.interval.bed
ilus WGS -c -C ilus_sys.yaml -L fastq.list -P align,BQSR,gvcf,VQSR -O ./
ilus WGS -c -n my_wgs -C ilus_sys.yaml -L fastq.list -O ./output
ilus split-jobs -I bigcs-II.step1.bwa.sh -n 53 -t 4
