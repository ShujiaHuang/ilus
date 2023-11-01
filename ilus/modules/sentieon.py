""" NGS Pipeline and variant calling support for Sentieon tools

Sentieon provides optimized versions of standard tools like GATK HaplotypeCaller
and MuTect2 as well as their own developed versions. These require a license
from Sentieon for use:

http://sentieon.com/about/
https://support.sentieon.com/manual/usages/general/
https://www.goldenhelix.com/media/pdfs/whitepapers/Sentieon_Genomics_Tools.pdf

"""


def license_export(config):
    """Retrieve export statement for sentieon license server.
    """
    pls = "Please set environment variable SENTIEON_LICENSE to your Sentieon license file."
    return


class Sentieon(object):
    """A class for Sentieon.
    Sentieon doc: https://support.sentieon.com/manual
    """

    def __init__(self, config):
        """Constructor.
        """
        self.config = config
        if "sentieon" not in self.config:
            raise ValueError("[Error] Missing `sentieon` in config file.")

        self.sentieon = self.config['sentieon']['sentieon']
        self.reference_index = self.config["resources"]["reference"]
        self.reference_fasta = self.config["resources"]["reference"]
        self.driver_options = " ".join(
            [str(x) for x in self.config["sentieon"].get("sentieon_driver_options", [])])

    def bwamem(self, out_prefix, rgID, fastq1, fastq2=""):
        """Mapping reads with BWA-MEM and sorting BAM to coordinate.

        This function process the following processes:
            - sentieon bwa-mem alignment and output BAM
            - sentieon sort of BAM to coordinate
        """
        if fastq2 == ".":
            fastq2 = ""

        bwa_options = " ".join([str(x) for x in self.config["sentieon"].get("bwamem_options", [])])
        util_sort_options = " ".join(
            [str(x) for x in self.config["sentieon"].get("util_sort_options", [])])
        bwa_cmd = (f"time {self.sentieon} bwa mem {bwa_options} -R {rgID} {self.reference_index} "
                   f"{fastq1} {fastq2} | {self.sentieon} util sort {util_sort_options} --sam2bam "
                   f"-o {out_prefix}.sorted.bam -i -")

        return f"{out_prefix}.sorted.bam", bwa_cmd

    def alignment_metrics(self, input_bam, out_prefix):
        """ Alignment metrics. Call this function before `markduplicates`
        """
        adapter_seq = ",".join([str(x) for x in self.config["resources"].get("adapter_sequence", [])])
        if len(adapter_seq) == 0:
            adapter_seq = "' '"

        coverage_options = " ".join([
            str(x) for x in self.config["sentieon"].get("coverage_options", [])])
        align_metrics_cmd = (f"time {self.sentieon} driver {self.driver_options} "
                             f"-r {self.reference_fasta} "
                             f"-i {input_bam} "
                             f"--algo MeanQualityByCycle {out_prefix}.mq_metrics.txt "
                             f"--algo QualDistribution {out_prefix}.qd_metrics.txt "
                             f"--algo GCBias --summary {out_prefix}.gc_summary.txt {out_prefix}.gc_metrics.txt "
                             f"--algo AlignmentStat --adapter_seq {adapter_seq} {out_prefix}.aln_metrics.txt "
                             f"--algo InsertSizeMetricAlgo {out_prefix}.is_metrics.txt "
                             f"--algo CoverageMetrics {coverage_options} {out_prefix}.coverage_metrics")

        plot_cmd = (f"{self.sentieon} plot GCBias -o {out_prefix}.gc-report.pdf {out_prefix}.gc_metrics.txt && "
                    f"{self.sentieon} plot QualDistribution -o {out_prefix}.qd-report.pdf {out_prefix}.qd_metrics.txt && "
                    f"{self.sentieon} plot MeanQualityByCycle -o {out_prefix}.mq-report.pdf {out_prefix}.mq_metrics.txt && "
                    f"{self.sentieon} plot InsertSizeMetricAlgo -o {out_prefix}.is-report.pdf {out_prefix}.is_metrics.txt")

        return align_metrics_cmd + " && " + plot_cmd

    def markduplicates(self, input_bam, output_markdup_bam):
        """Markdup Duplicate Reads.
        """
        locus_collector_options = " ".join([
            str(x) for x in self.config["sentieon"].get("LocusCollector_options", [])
        ])
        dedup_options = " ".join([str(x) for x in self.config["sentieon"].get("dedup_options", [])])

        output_markdup_bam = str(output_markdup_bam)
        score_info_file = output_markdup_bam.replace(".bam", ".score.txt")
        dedup_metrics_f = output_markdup_bam.replace(".bam", ".metrics.txt")
        dedup_cmd = (f"time {self.sentieon} driver {self.driver_options} "
                     f"-i {input_bam} "
                     f"--algo LocusCollector {locus_collector_options} "
                     f"--fun score_info {score_info_file} && "
                     f"time {self.sentieon} driver {self.driver_options} "
                     f"-i {input_bam} "
                     f"--algo Dedup {dedup_options} --score_info {score_info_file} "
                     f"--metrics {dedup_metrics_f} {output_markdup_bam}")

        coverage_options = " ".join([
            str(x) for x in self.config["sentieon"].get("coverage_options", [])])

        cvg_metrics_fn = output_markdup_bam.replace(".bam", ".coverage_metrics")

        # '--omit_base_output' skip the output of the per locus coverage with no partition.
        # This option can be used when you do not use intervals to save space.
        cvg_cmd = (f"{self.sentieon} driver {self.driver_options} "
                   f"-r {self.reference_fasta} "
                   f"-i {output_markdup_bam} "
                   f"--algo CoverageMetrics {coverage_options} {cvg_metrics_fn}")

        return dedup_cmd + " && " + cvg_cmd

    def indelrealigner(self, input_bam, output_realig_bam):
        """Indel realignment for sample.
        """
        indel_realigner_options = " ".join([str(x) for x in self.config["sentieon"].get(
            "indel_realigner_options", [])])
        known_Mills_indels = self.config["resources"]["bundle"]["mills"]
        known_1000G_indels = self.config["resources"]["bundle"]["1000G_known_indel"]

        return (f"time {self.sentieon} driver {self.driver_options} "
                f"-r {self.reference_fasta} "
                f"-i {input_bam} "
                f"--algo Realigner {indel_realigner_options} "
                f"-k {known_Mills_indels} "
                f"-k {known_1000G_indels} {output_realig_bam}")

    def baserecalibrator(self, input_bam, out_bqsr_recal_table):
        """Base recalibration.

        Note:
            1. -k options to known SNP sites are optional, but recommended.
            2. Applying the ReadWriter algo is optional, the next step in
               HC will apply calibration table on the fly.
        """
        dbsnp = self.config["resources"]["bundle"]["dbsnp"]
        known_Mills_indels = self.config["resources"]["bundle"]["mills"]
        known_1000G_indels = self.config["resources"]["bundle"]["1000G_known_indel"]
        bqsr_recaltable_options = " ".join([str(x) for x in self.config["sentieon"].get(
            "bqsr_recaltable_options", [])])

        # Create recalibrate table file for BQSR
        recal_table_cmd = (f"time {self.sentieon} driver {self.driver_options} "
                           f"-r {self.reference_fasta} "
                           f"-i {input_bam} "
                           f"--algo QualCal {bqsr_recaltable_options} "
                           f"-k {dbsnp} "
                           f"-k {known_Mills_indels} "
                           f"-k {known_1000G_indels} {out_bqsr_recal_table}")

        # Apply BQSR
        recal_table_post_fn = f"{out_bqsr_recal_table}.post"
        apply_bqsr_cmd = (f"time {self.sentieon} driver {self.driver_options} "
                          f"-r {self.reference_fasta} "
                          f"-i {input_bam} "
                          f"-q {out_bqsr_recal_table} "
                          f"--algo QualCal {bqsr_recaltable_options} "
                          f"-k {dbsnp} "
                          f"-k {known_Mills_indels} "
                          f"-k {known_1000G_indels} {recal_table_post_fn} "

                          # 输出 BQSR 之后的 BAM 是不必要的，因为 HaplotypeCaller 可以通过读入
                          # Indel Realigner BAM 和 recal_table_fn 实现校正。如此既可以确保原
                          # 始的 base quality 不被修改，也能确保可以得到 BQSR 校正的 calling。
                          # f"--algo ReadWriter {output_bqsr_bam}"
                          )

        # Plot report
        recal_report_table = f"{out_bqsr_recal_table}.report.csv"
        bqsrreport_plot_fn = f"{out_bqsr_recal_table}.bqsrreport.pdf"
        plot_cmd = (f"time {self.sentieon} driver {self.driver_options} "
                    f"--algo QualCal --plot "
                    f"--before {out_bqsr_recal_table} "
                    f"--after {recal_table_post_fn} {recal_report_table} && "
                    f"{self.sentieon} plot QualCal "
                    f"-o {bqsrreport_plot_fn} {recal_report_table}")

        return recal_table_cmd + " && " + apply_bqsr_cmd + " && " + plot_cmd

    def haplotypecaller(self, input_bam, bqsr_recal_table, output_vcf_fname, interval=None):
        """HC Variant caller generating GVCF.

        :param input_bam: Input Indel Realigner BAM file
        :param bqsr_recal_table: BQSR recalibrator table.
        :param output_vcf_fname: gvcf
        :param interval: Interval string or file (BED/Picard)
        :return:
        """
        dbsnp = self.config["resources"]["bundle"]["dbsnp"]
        hc_options = " ".join([str(x) for x in self.config["sentieon"].get("hc_options", [])])

        if ("gvcf" in hc_options) and (".g.vcf" not in str(output_vcf_fname)):
            raise ValueError(f"[ERROR] {output_vcf_fname} missing .g.vcf in file name")

        driver_options = f"{self.driver_options} --interval {interval}" \
            if interval else self.driver_options

        return (f"time {self.sentieon} driver {driver_options} "
                f"-r {self.reference_fasta} "
                f"-i {input_bam} "
                f"-q {bqsr_recal_table} "
                f"--algo Haplotyper {hc_options} "
                f"-d {dbsnp} {output_vcf_fname}")

    def genotypeGVCFs(self, input_gvcfs_list, output_vcf_fname, interval=None):
        """Perform the joint calling by `GVCFtyper` module with input GVCFs.
        """
        # location of the Single Nucleotide Polymorphism database (dbSNP) used to
        # label known variants.
        dbsnp = self.config["resources"]["bundle"]["dbsnp"]
        gvcftyper_options = " ".join([
            str(x) for x in self.config["sentieon"].get("gvcftyper_options", [])])

        driver_options = f"{self.driver_options} --interval {interval}" \
            if interval else self.driver_options

        in_gvcf_arguments = " ".join([f"-v {gvcf_fn}" for gvcf_fn in input_gvcfs_list])
        return (f"time {self.sentieon} driver {driver_options} "
                f"-r {self.reference_fasta}  "
                f"--algo GVCFtyper {gvcftyper_options} "
                f"-d {dbsnp} "
                f"{in_gvcf_arguments} {output_vcf_fname}")

    def variantrecalibrator(self, input_vcf, output_vcf_fname):
        """Use VarCal and ApplyVarCal modules to do VQSR in sentieon.
        """
        resource_hapmap = self.config["resources"]["bundle"]["hapmap"]
        resource_omni = self.config["resources"]["bundle"]["omni"]
        resource_1000G = self.config["resources"]["bundle"]["1000G"]
        resource_dbsnp = self.config["resources"]["bundle"]["dbsnp"]
        resource_mills_gold_indels = self.config["resources"]["bundle"]["mills"]
        resource_1000G_known_indel = self.config["resources"]["bundle"]["1000G_known_indel"]

        vqsr_options = " ".join([str(x) for x in self.config["sentieon"].get("vqsr_options", [])])
        apply_snp_vqsr_options = " ".join([
            str(x) for x in self.config["sentieon"].get("apply_snp_vqsr_options", [])
        ])
        apply_indel_vqsr_options = " ".join([
            str(x) for x in self.config["sentieon"].get("apply_indel_vqsr_options", [])
        ])

        if '--var_type' in vqsr_options + apply_snp_vqsr_options + apply_indel_vqsr_options:
            raise ValueError("[ERROR] Do not set '--var_type' in `vqsr_options` or "
                             "`apply_snp_vqsr_options` or `apply_indel_vqsr_options` "
                             "in configuration file")

        # Set name
        out_prefix = str(output_vcf_fname).replace(".gz", "").replace(".vcf", "")  # delete .vcf.gz
        out_snp_vqsr_fname = out_prefix + ".SNPs.vcf.gz"
        snp_vqsr_cmd = (
            f"time {self.sentieon} driver {self.driver_options} "
            f"-r {self.reference_fasta} "
            f"--algo VarCal {vqsr_options} "
            f"-v {input_vcf} "
            f"--resource {resource_hapmap} --resource_param hapmap,known=false,training=true,truth=true,prior=15.0  "
            f"--resource {resource_omni} --resource_param omini,known=false,training=true,truth=truth,prior=12.0 "
            f"--resource {resource_1000G} --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 "
            f"--resource {resource_dbsnp} --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "
            f"--var_type SNP "
            f"--max-gaussians {self.config['sentieon']['vqsr_snp_max_gaussians']} "
            f"--tranches_file {out_prefix}.SNPs.tranches.csv "
            f"--plot_file {out_prefix}.plot.SNPs.csv "
            f"{out_prefix}.SNPs.recal"  # Output quality recalibrator data for SNPs
        )
        apply_snp_vqsr_cmd = (f"time {self.sentieon} driver {self.driver_options} "
                              f"-r {self.reference_fasta} "
                              f"--algo ApplyVarCal {apply_snp_vqsr_options} "
                              f"-v {input_vcf} "
                              f"--var_type SNP "
                              f"--recal {out_prefix}.SNPs.recal "
                              f"--tranches_file {out_prefix}.SNPs.tranches.csv "
                              f"{out_snp_vqsr_fname}")

        # Indel VQSR after SNP
        indel_vqsr_cmd = (
            f"time {self.sentieon} driver {self.driver_options} "
            f"-r {self.reference_fasta} "
            f"--algo VarCal {vqsr_options} "
            f"-v {out_snp_vqsr_fname} "
            f"--resource {resource_mills_gold_indels} --resource_param omini,known=false,training=true,truth=truth,prior=12.0 "
            f"--resource {resource_1000G_known_indel} --resource_param 1000G,known=false,training=true,truth=true,prior=10.0 "
            f"--resource {resource_dbsnp} --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 "
            f"--var_type INDEL "
            f"--max-gaussians {self.config['sentieon']['vqsr_indel_max_gaussians']} "
            f"--tranches_file {out_prefix}.INDELs.tranches.csv "
            f"--plot_file {out_prefix}.plot.INDELs.csv "
            f"{out_prefix}.INDELs.recal"  # Output quality recalibrator data for INDELs
        )
        apply_indel_vqsr_cmd = (f"time {self.sentieon} driver {self.driver_options} "
                                f"-r {self.reference_fasta} "
                                f"--algo ApplyVarCal {apply_indel_vqsr_options} "
                                f"-v {out_snp_vqsr_fname} "
                                f"--var_type INDEL "
                                f"--recal {out_prefix}.INDELs.recal "
                                f"--tranches_file {out_prefix}.INDELs.tranches.csv "
                                f"{output_vcf_fname} && rm -f {out_snp_vqsr_fname}")

        return " && ".join([snp_vqsr_cmd, apply_snp_vqsr_cmd, indel_vqsr_cmd, apply_indel_vqsr_cmd])
