"""GATK modules

Author: Shujia Huang
Date: 2020-04-19 15:19:56
"""


class GATK(object):
    """A class for GATK
    """

    def __init__(self, config):
        """Constructor.
        """
        self.config = config
        if "gatk" not in self.config:
            raise ValueError("[Error] Missing `gatk` in config file.")

        self.gatk_options = self.config['gatk']
        self.gatk = self.gatk_options['gatk']  # gatk program

        self.reference_index = self.config["resources"]["reference"]
        self.reference_fasta = self.config["resources"]["reference"]
        self.resources_bundle = self.config["resources"]["bundle"]

    def markduplicates(self, input_bam, output_markdup_bam):
        java_options = "--java-options \"%s\"" % " ".join(self.gatk_options["markdup_java_options"]) \
            if "markdup_java_options" in self.gatk_options \
               and len(self.gatk_options["markdup_java_options"]) else ""

        out_metrics_fname = str(output_markdup_bam).replace(".bam", ".metrics.txt")
        return (f"time {self.gatk} {java_options} MarkDuplicates "
                f"-I {input_bam} "
                f"-M {out_metrics_fname} "
                f"-O {output_markdup_bam}")

    def baserecalibrator(self, input_bam, output_bqsr_bam, out_bqsr_recal_table):
        java_options = "--java-options \"%s\"" % " ".join(self.gatk_options["bqsr_java_options"]) \
            if "bqsr_java_options" in self.gatk_options \
               and len(self.gatk_options["bqsr_java_options"]) else ""

        known_site_1000G_indel = self.resources_bundle["1000G_known_indel"]
        known_site_mills_gold_indels = self.resources_bundle["mills"]
        known_site_dbsnp = self.resources_bundle["dbsnp"]

        # create recalibrate table file for BQSR
        recal_data_cmd = (f"time {self.gatk} {java_options} BaseRecalibrator "
                          f"-R {self.reference_fasta} "
                          f"--known-sites {known_site_dbsnp} "
                          f"--known-sites {known_site_1000G_indel} "
                          f"--known-sites {known_site_mills_gold_indels} "
                          f"-I {input_bam} "
                          f"-O {out_bqsr_recal_table}")

        # ApplyBQSR
        apply_bqsr_cmd = (f"time {self.gatk} {java_options} ApplyBQSR "
                          f"-R {self.reference_fasta} "
                          f"--bqsr-recal-file {out_bqsr_recal_table} "
                          f"-I {input_bam} "
                          f"-O {output_bqsr_bam}")

        if "capture_interval_file" in self.config:
            recal_data_cmd += f" -L {self.config['capture_interval_file']}"
            apply_bqsr_cmd += f" -L {self.config['capture_interval_file']}"

        return recal_data_cmd + " && " + apply_bqsr_cmd

    def haplotypecaller_gvcf(self, input_bam, output_gvcf_fname, interval=None):
        java_options = "--java-options \"%s\"" % " ".join(self.gatk_options["hc_gvcf_java_options"]) \
            if "hc_gvcf_java_options" in self.gatk_options \
               and len(self.gatk_options["hc_gvcf_java_options"]) else ""

        hc_options = " ".join(self.gatk_options["hc_gvcf_options"]) \
            if "hc_gvcf_options" in self.gatk_options else ""

        cmd = (f"time {self.gatk} {java_options} HaplotypeCaller "
               f"-R {self.reference_fasta} {hc_options} "
               f"--emit-ref-confidence GVCF "
               f"-I {input_bam} "
               f"-O {output_gvcf_fname}")

        if interval:
            cmd += f" -L {interval}"

        return cmd

    def combineGVCFs(self, input_sample_gvcfs, output_combineGVCF_fname, interval=None):
        """Combine GVCFs by GATK genomicsDBImport or CombineGVCFs module.
        """
        java_options = "--java-options \"%s\"" % " ".join(self.gatk_options["combineGVCFs_java_options"]) \
            if "combineGVCFs_java_options" in self.gatk_options \
               and len(self.gatk_options["combineGVCFs_java_options"]) else ""

        # set overwite existing genomicsdb workspace by default
        if "genomicsDBImport_options" not in self.gatk_options:
            self.gatk_options["genomicsDBImport_options"] = ["--overwrite-existing-genomicsdb-workspace true"]
        elif (("--overwrite-existing-genomicsdb-workspace false" not in self.gatk_options["genomicsDBImport_options"])
              and
              ("--overwrite-existing-genomicsdb-workspace true" not in self.gatk_options["genomicsDBImport_options"])):
            self.gatk_options["genomicsDBImport_options"].append("--overwrite-existing-genomicsdb-workspace true")

        genomicsDBImport_options = "%s" % " ".join(self.gatk_options["genomicsDBImport_options"])
        use_gDBI = self.gatk_options["use_genomicsDBImport"] \
            if "use_genomicsDBImport" in self.gatk_options else False

        # Create command line for GenomicsDBImport or CombineGVCFs
        if use_gDBI:
            # use GenomicsDBImport
            sample_name_map = input_sample_gvcfs[0]  # Only one file
            combine_gvcf_cmd = (f"time {self.gatk} {java_options} GenomicsDBImport {genomicsDBImport_options} "
                                f"-R {self.reference_fasta} "
                                f"--sample-name-map {sample_name_map} "
                                f"--genomicsdb-workspace-path {output_combineGVCF_fname}")
        else:
            gvcfs = " ".join(["-V %s" % s for s in input_sample_gvcfs])
            combine_gvcf_cmd = (f"time {self.gatk} {java_options} CombineGVCFs "
                                f"-R {self.reference_fasta} {gvcfs} "
                                f"-O {output_combineGVCF_fname}")
        if interval:
            combine_gvcf_cmd += f" -L {interval}"

        return combine_gvcf_cmd

    def genotypeGVCFs(self, input_combine_gvcf_fname, output_vcf_fname, interval=None):
        java_options = "--java-options \"%s\"" % " ".join(self.gatk_options["genotype_java_options"]) \
            if "genotype_java_options" in self.gatk_options \
               and len(self.gatk_options["genotype_java_options"]) else ""

        # location of the Single Nucleotide Polymorphism database (dbSNP) used to label known variants.
        use_gDBI = self.gatk_options["use_genomicsDBImport"] \
            if "use_genomicsDBImport" in self.gatk_options else False

        genotypeGVCFs_options = " ".join(self.gatk_options["genotypeGVCFs_options"]) \
            if "genotypeGVCFs_options" in self.gatk_options else ""

        if interval:
            genotypeGVCFs_options += f" -L {interval}"

        # Creat command line for genotypeGVCF
        if use_gDBI:
            genotype_cmd = (f"time {self.gatk} {java_options} GenotypeGVCFs "
                            f"-R {self.reference_fasta} {genotypeGVCFs_options} "
                            f"--dbsnp {self.resources_bundle['dbsnp']} "
                            f"-V gendb://{input_combine_gvcf_fname} "
                            f"-O {output_vcf_fname}")
        else:
            genotype_cmd = (f"time {self.gatk} {java_options} GenotypeGVCFs "
                            f"-R {self.reference_fasta} {genotypeGVCFs_options} "
                            f"--dbsnp {self.resources_bundle['dbsnp']} "
                            f"-V {input_combine_gvcf_fname} "
                            f"-O {output_vcf_fname}")

        # reindex with tabix for containing count metadata
        vcf_index_cmd = f"time {self.config['tabix']} -f -p vcf {output_vcf_fname}"

        return f"{genotype_cmd} && {vcf_index_cmd}"

    def variantrecalibrator(self, input_vcf, output_vcf_fname):
        java_options = "--java-options \"%s\"" % " ".join(self.gatk_options["vqsr_java_options"]) \
            if "vqsr_java_options" in self.gatk_options \
               and len(self.gatk_options["vqsr_java_options"]) else ""

        vqsr_snp_options = " ".join([str(x) for x in self.gatk_options.get("vqsr_snp_options", [])])
        vqsr_indel_options = " ".join([str(x) for x in self.gatk_options.get("vqsr_indel_options", [])])
        if "--max-gaussians" in vqsr_snp_options + vqsr_indel_options:
            raise ValueError("[ERROR] No need to set --max-gaussians for "
                             "`vqsr_snp_options` and `vqsr_indel_options`"
                             "in your configuration file (.yaml)")

        apply_snp_vqsr_options = " ".join([str(x) for x in self.gatk_options.get("apply_snp_vqsr_options", [])])
        apply_indel_vqsr_options = " ".join([str(x) for x in self.gatk_options.get("apply_indel_vqsr_options", [])])
        if '-mode' in vqsr_snp_options + vqsr_indel_options + \
                apply_snp_vqsr_options + apply_indel_vqsr_options:
            raise ValueError("[ERROR] Do not set '-mode' in `vqsr_snp_options` or "
                             "`vqsr_indel_options` or `apply_snp_vqsr_options` or "
                             "`apply_indel_vqsr_options` in configuration file")

        # Set name
        out_prefix = output_vcf_fname.replace(".gz", "").replace(".vcf", "")  # delete .vcf.gz
        out_snp_vqsr_fname = out_prefix + ".SNPs.vcf.gz"

        resource_hapmap = self.resources_bundle["hapmap"]
        resource_omni = self.resources_bundle["omni"]
        resource_1000G = self.resources_bundle["1000G"]
        resource_dbsnp = self.resources_bundle["dbsnp"]
        resource_mills_gold_indels = self.resources_bundle["mills"]
        resource_1000G_known_indel = self.resources_bundle["1000G_known_indel"]

        # SNP VQSR
        snp_vqsr_cmd = (
            f"time {self.gatk} {java_options} VariantRecalibrator {vqsr_snp_options} "
            f"-R {self.reference_fasta} "
            f"-V {input_vcf} "
            f"--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {resource_hapmap} "
            f"--resource:omini,known=false,training=true,truth=true,prior=12.0 {resource_omni} "
            f"--resource:1000G,known=false,training=true,truth=false,prior=10.0 {resource_1000G} "
            f"--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {resource_dbsnp} "
            f"-mode SNP "
            f"--max-gaussians {self.gatk_options['vqsr_snp_max_gaussians']} "
            f"--tranches-file {out_prefix}.SNPs.tranches.csv "
            f"-O {out_prefix}.SNPs.recal")
        apply_snp_vqsr_cmd = (f"time {self.gatk} {java_options} ApplyVQSR {apply_snp_vqsr_options} "
                              f"-R {self.reference_fasta} "
                              f"-V {input_vcf} "
                              f"--tranches-file {out_prefix}.SNPs.tranches.csv "
                              f"--recal-file {out_prefix}.SNPs.recal "
                              f"-mode SNP "
                              f"-O {out_snp_vqsr_fname}")

        # Indel VQSR after SNP
        indel_vqsr_cmd = (
            f"time {self.gatk} {java_options} VariantRecalibrator {vqsr_indel_options} "
            f"-R {self.reference_fasta} "
            f"-V {out_snp_vqsr_fname} "
            f"--resource:mills,known=false,training=true,truth=true,prior=12.0 {resource_mills_gold_indels} "
            f"--resource:1000G,known=false,training=true,truth=true,prior=10.0 {resource_1000G_known_indel} "
            f"--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {resource_dbsnp} "
            f"--tranches-file {out_prefix}.INDELs.tranches.csv "
            f"-mode INDEL "
            f"--max-gaussians {self.gatk_options['vqsr_indel_max_gaussians']} "
            f"-O {out_prefix}.INDELs.recal")
        apply_indel_vqsr_cmd = (
            f"time {self.gatk} {java_options} ApplyVQSR {apply_indel_vqsr_options} "
            f"-R {self.reference_fasta} "
            f"-V {out_snp_vqsr_fname} "
            f"--tranches-file {out_prefix}.INDELs.tranches.csv "
            f"--recal-file {out_prefix}.INDELs.recal "
            f"-mode INDEL "
            f"-O {output_vcf_fname} && rm -f {out_snp_vqsr_fname}* "
        )

        tabix = self.config["tabix"]
        vcf_index_cmd = f"time {tabix} -f -p vcf {output_vcf_fname}"
        return " && ".join([snp_vqsr_cmd,
                            apply_snp_vqsr_cmd,
                            indel_vqsr_cmd,
                            apply_indel_vqsr_cmd,
                            vcf_index_cmd])

    def collect_alignment_summary_metrics(self, input_bam, output_fname):
        java_options = f"--java-options {' '.join(self.gatk_options['CollectAlignmentSummaryMetrics_jave_options'])}" \
            if "CollectAlignmentSummaryMetrics_jave_options" in self.gatk_options \
               and len(self.gatk_options["CollectAlignmentSummaryMetrics_jave_options"]) else ""

        metric_options = " ".join(self.gatk_options["CollectAlignmentSummaryMetrics_options"]) \
            if "CollectAlignmentSummaryMetrics_options" in self.gatk_options \
               and len(self.gatk_options["CollectAlignmentSummaryMetrics_options"]) else ""

        # Add Adapter sequence
        for x in self.config["resources"].get("adapter_sequence", []):
            metric_options += f" --ADAPTER_SEQUENCE {x}"

        cmd = (f"time {self.gatk} {java_options} CollectAlignmentSummaryMetrics {metric_options} "
               f"-R {self.reference_fasta} "
               f"-I {input_bam} "
               f"-O {output_fname}")

        return cmd

    def mergevcfs(self, input_vcfs, output_fname):
        # Sometimes bcftools concat is better than MergeVcfs
        java_options = f"--java-options {' '.join(self.gatk_options['mergevcfs_java_options'])}" \
            if "mergevcfs_java_options" in self.gatk_options \
               and len(self.gatk_options["mergevcfs_java_options"]) else ""

        cmd = [f"time {self.gatk} {java_options} MergeVcfs -O {output_fname}"] + \
              [f"-I {f}" for f in input_vcfs]
        return " ".join(cmd)
