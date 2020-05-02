```bash
python ../ilus/runner.py WGS -C ilus_sys.yaml -L fastq.list -O ./
```

```bash
python ../ilus/runner.py genotype-joint-calling -C ilus_sys.yaml -L gvcf.list -O 03.genotype 
```

```bash
python ../ilus/runner.py genotype-joint-calling -C ilus_sys.yaml -L gvcf.list -O 03.genotype --as_pipe_shell_order -f
```
```bash
python ../../scripts/yhbatch_slurm_jobs.py -I step1.bwa.sh -n 10 -t 3
```

```bash
python ../ilus/runner.py VQSR -C ilus_sys.yaml -L vcf.list -O 03.genotype --as_pipe_shell_order -f 
```
