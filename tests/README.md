```bash
python ../ilus/runner.py pipeline -C ../config/ilus_sys.yaml -L fastq.list -O ./
```

```bash
python ../ilus/runner.py genotype -C ../config/ilus_sys.yaml -L gvcf.list -O 03.genotype 
```

```bash
python ../ilus/runner.py genotype -C ../config/ilus_sys.yaml -L gvcf.list -O 03.genotype --as_pipe_shell_order -f
```
