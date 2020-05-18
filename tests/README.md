```bash
ilus WGS -c -C ilus_sys.yaml -L fastq.list -O ./
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
