```bash
ilus WGS -c -C ilus_sys.yaml -I sample_fastq.list -O ./tmp
```

```bash
ilus WGS -c -C ilus_sys.yaml -I sample_fastq.list -O ./tmp -dr
```

```bash
ilus WGS -P gvcf,genotype -c -C ilus_sys.yaml -I sample_fastq.list -O ./tmp
ilus WGS -P gvcf -c -C ilus_sys.yaml -I sample_fastq.list -O ./tmp
ilus WGS -P genotype -c -C ilus_sys.yaml -I sample_fastq.list -O ./tmp -f 
```

```bash
ilus genotype-joint-calling -C ilus_sys.yaml -I gvcf.list -O tmp/03.genotype 
```

```bash
ilus genotype-joint-calling -C ilus_sys.yaml -I gvcf.list -O tmp/03.genotype --as_pipe_shell_order -f
ilus genotype-joint-calling -C ilus_sys.yaml -I gvcf.list -O tmp/03.genotype --as_pipe_shell_order -f --use-sentieon
ilus genotype-joint-calling -C ilus_sys.yaml -I gvcf.list -O tmp/03.genotype --as_pipe_shell_order -f --use-sentieon --interval human_GRCh38.calling_interval.bed
```

```bash
ilus VQSR -C ilus_sys.yaml -I vcf.list -O 03.genotype --as_pipe_shell_order -f 
```

```bash
python ../../scripts/yhbatch_slurm_jobs.py -I step1.bwa.sh -n 10 -t 3
```

python ../scripts/create_wgs_calling_intervals.py -f human_GRCh38.fai -n human_GRCh38.N.interval_list -w 5000000 | head -820 > human_GRCh38.calling_interval.bed

ilus WGS -c -C ilus_sys.yaml -I sample_fastq.list -P align,BQSR,gvcf,VQSR -O ./tmp
ilus WGS -c -n my_wgs -C ilus_sys.yaml -I sample_fastq.list -O ./output

ilus WGS -c -n my_wgs -C ilus_sys.yaml -I sample_fastq.list -O ./tmp --interval chr2,chr3:4-10 -f
ilus WGS -c -n my_wgs -C ilus_sys.yaml -I sample_fastq.list -O ./tmp --use-sentieon --interval chr2,chr3:4-10 -f

ilus WGS -c -n my_wgs -C ilus_sys.yaml -I sample_fastq.list -O ./tmp --interval human_GRCh38.calling_interval.bed -f
ilus WGS -c -n my_wgs -C ilus_sys.yaml -I sample_fastq.list -O ./tmp --use-sentieon --interval human_GRCh38.calling_interval.bed -f

## Use sentieon to call variants
ilus WGS -c -n my_wgs -C ilus_sys.yaml -I sample_fastq.list -O ./tmp --use-sentieon -f
ilus WGS -c -n my_wgs -C ilus_sys.yaml -I sample_fastq.list -O ./tmp --use-sentieon -f --interval human_GRCh38.calling_interval.bed 
ilus WGS -c -n my_wgs -C ilus_sys.yaml -I sample_fastq.list -O ./tmp -f --interval human_GRCh38.calling_interval.bed

## WES
ilus WES -c -n my_wgs -C ilus_sys.yaml -I sample_fastq.list --capture-interval human_GRCh38.calling_interval.bed -O ./tmp -f
ilus WES -c -n my_wgs -C ilus_sys.yaml -I sample_fastq.list --capture-interval human_GRCh38.calling_interval.bed -O ./tmp -f --use-sentieon
