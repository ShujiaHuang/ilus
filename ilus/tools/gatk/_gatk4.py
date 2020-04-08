"""Work with GATK commands
"""
import os
import sys
import re
import time

import multiprocessing

from ilus.launch import do
from ilus.utils import safe_makedir


def check_bundlefile_index(gatk, gatk_resource):
    """check GATK bundle is been indexed or not
    and create the index if not.

    Parameters
    ==========
        ``gatk``: String.
            A path to the GATK executor program

        ``gatk_resource``: A dict like.
            dataset's name => the path of bundle file
    """

    for k, datafile in gatk_resource.items():

        if not os.path.exists(datafile):
            sys.stderr.write(
                "[ERROR] %s is missing. Process Terminate!\n" %
                datafile)
            sys.exit(1)

        create_index = False
        cmd = [gatk, "IndexFeatureFile"]

        if datafile.endswith(".gz") and not os.path.exists(
            "%s.tbi" %
                datafile):
            create_index = True
            cmd.append("--feature-file %s" % datafile)

        elif datafile.endswith(".vcf") and not os.path.exists("%s.idx" % datafile):
            create_index = True
            cmd.append("--feature-file %s" % datafile)

        if create_index:

            cmd = " ".join(cmd)
            do.run(cmd, "GATK %s IndexFeatureFile done" % datafile)

    return


def _set_training_data(bundle_file):
    """Retrieve training data, returning an empty set of information if not avalable.

    Parameters:
        ``bundle_file``: A dict like.
            It's a dict that contain GATK bundle resource.
            Key is name of bundle, value is the file path.
    """

    if bundle_file is None:
        return {}

    bundle = {"SNP": [], "INDEL": []}
    # For SNP
    for name, resource in [("hapmap", "known=false,training=true,truth=true,prior=15.0"),
                           ("omni", "known=false,training=true,truth=false,prior=12.0"),
                           ("1000G", "known=false,training=true,truth=false,prior=10.0"),
                           ("dbsnp", "known=true,training=false,truth=false,prior=2.0")]:
        if name in bundle_file:
            bundle["SNP"].append([name, resource, bundle_file[name]])
        else:
            return {}

    # For Indel
    if "mills" in bundle_file:
        bundle["INDEL"].append(["mills",
                                "known=true,training=true,truth=true,prior=12.0",
                                bundle_file["mills"]])
    else:
        return {}

    return bundle


class GATKRunner(object):
    """Running GATK4 commandline tools.
    """

    def __init__(self, config, gatk_resource=None):
        """"""
        self.resources = _set_training_data(gatk_resource)
        self.reference = config["resources"]["reference"]
        self.gatk = config["variantcaller"]["gatk"]
        self.process_num = 1
        self.runable = ""

        if "process_num" in config["variantcaller"]:
            self.process_num = config["variantcaller"]["process_num"]

        if "runable_region" in config["variantcaller"]:
            self.runable = config["variantcaller"]["runable_region"]

        if self.process_num < 1:
            self.process_num = multiprocessing.cpu_count() / 4  # default detected

    def run(self, tools, options):
        """Run a GATK command with the provide options.
        """
        cmd = " ".join([self.gatk, tools] + ["%s %s" % (p, v)
                                             for p, v in options])
        if self.process_num == 1 or (
                tools != "HaplotypeCaller" and tools != "GenotypeGVCFs"):
            do.run(cmd, "GATK %s done" % tools)
            return

        runablebed = self.runable
        chrom = "All"  # prepared for wes/wgs of all chroms

        p = re.compile(r"\s+-L\s+")  # intervals
        if p.search(cmd):
            """
            eg1: chr1
            eg2: chr1:1-3 chr2:1-5
            eg3: chr1:1-2
            eg3: wes.bed
            """
            columns = p.split(cmd)
            interval = columns[1:-1] + [columns[-1].split()[0]]
            symbol = re.compile(r'[\:\-\,]+')
            if len(interval) == 1 and os.path.isfile(interval[0]):
                # wes bed file
                cmd = "%s %s" % (columns[0], " ".join(columns[-1].split()[1:]))
                runablebed = interval[0]

            elif len(interval) == 1 and (not symbol.search(interval[0])):
                # wgs -L 1
                cmd = "%s %s" % (columns[0], " ".join(columns[-1].split()[1:]))
                chrom = interval[0]

            else:
                sys.stderr.write(
                    "Because of the interval seted,to ensure the correct results,"
                    "the program will run in single process mode\n")
                do.run(cmd, "GATK %s done" % tools)
                return

        self.multiprocess(cmd, runablebed, self.process_num, chrom)
        return

    def multiprocess(self, cmd, bed, N, chr):

        sys.stdout.write("process_num is %s\nchr is %s\n" % (N, chr))

        p = re.compile(r"\s+-O\s+([\w\.\/]+)")  # output name
        output = p.search(cmd).group(1)
        name = os.path.basename(output)
        o_dir = os.path.abspath(os.path.dirname(output))
        tmp_dir = "%s/multiprocess" % o_dir

        safe_makedir(tmp_dir)

        regions_for_each_process = [[] for _ in range(N)]
        sub_vcf_files = []
        for n in range(N):
            new_name = "%s/%s_%s" % (tmp_dir, n, name)
            new_cmd = p.sub(" -O %s" % new_name, cmd)
            regions_for_each_process[n] = [new_cmd]
            sub_vcf_files.append(new_name)

        i = 0
        if bed != self.runable:
            '''
            1. wes 0.bed 1.bed ... 4.bed
            2. wgs multiple splited little bed 

            '''
            bed_region = [[] for _ in range(N)]
            with open(bed, "r") as b:
                for line in b:
                    chromid = line.split()[0]
                    start = int(line.split()[1])
                    end = int(line.split()[2])
                    bed_region[i % N].append([chromid, start, end])
                    i += 1

            if i < N:
                '''
                bed_region=[[['10',1000,5000]],[['10',118673526,134182548]],[['10',134182549,134185747]],[]]  i=3 N=4 test
                output:
                finaly bed_spit:[[['10', 1000, 1999], ['10', 118673526, 122550780]], [['10', 2000, 2999], ['10', 122550781, 126428035]], [['10', 3000, 3999], ['10', 126428036, 130305290], ['10', 134182549, 134185747]], [['10', 4000, 5000], ['10', 130305291, 134182548]]]

                '''
                bed_split = [[] for _ in range(N)]
                for en, N_list in enumerate(bed_region[:i]):
                    chrid, start, end = N_list[0]
                    size = end - start + 1
                    if size < 4000:
                        '''
                        Add the corresponding position that is suitble for
                        the bed of multiple short line
                        '''
                        bed_split[en].append([chrid, start, end])
                        continue
                    i = N+1  # if split the bed to N region; make i >N for all the regions_for_each_process

                    delta = int(size / N)
                    if delta == 0:
                        delta = 1

                    region_list = range(start - 1, end, delta)
                    for m, pos in enumerate(region_list):
                        s = pos + 1 if pos + 1 < end else end
                        e = pos + delta if pos + delta < end else end

                        if m + 1 == len(region_list) and e - s < 500:
                            '''
                            Just handle the last one circle
                            because:strong correlation with seted process num
                            4000/10=400 all the splited bed low than 1000
                            '''
                            bed_split[(m - 1) % N][-1][-1] = e
                        else:
                            bed_split[m % N].append([chrid, s, e])

                bed_region = bed_split

                '''
                    Handle it before run the multiproccess 
                    make program helthier and prevent from the few short region. 
                    if bed_region[-1] == []:  
                        bed_region = bed_region[:i]
                        regions_for_each_process = regions_for_each_process[:i]
                '''
            bed_name = os.path.basename(bed)
            for r, b in enumerate(bed_region):
                r_bed = "%s/%s_%s" % (tmp_dir, r, bed_name)
                with open(r_bed, "w") as w:
                    content_r = []
                    for listr in b:
                        content_r.append(
                            "%s\t%s\t%s" %
                            (listr[0], listr[1], listr[2]))

                    w.write("%s\n" % "\n".join(content_r))

                regions_for_each_process[r].append("-L %s" % r_bed)
        else:
            # wgs input the single runnable bed that must be multiple lines in yaml

            with open(bed, "r") as b:
                for line in b:
                    chromid, start, end = line.split()
                    if chr != "All":
                        # only run the specified chrom region
                        if chromid == chr:
                            regions_for_each_process[i % N].append(
                                "-L %s:%s-%s" % (chromid, start, end))
                        else:
                            continue
                    else:
                        regions_for_each_process[i % N].append(
                            "-L %s:%s-%s" % (chromid, start, end))
                    i += 1

        if i < N:  # lines is so few that can't satiesfy the multiprocess.Almost never
            regions_for_each_process = regions_for_each_process[:i]
            sys.stderr.write(
                "[WARNING] We suggest  checking your bed that is so short and few"
                "that can't satiesfy the multiprocess \n")
            sub_vcf_files = sub_vcf_files[:i]

        # general in wgs and wes
        jobs = list()
        for m in regions_for_each_process:
            commands_m = " ".join(m)
            p = multiprocessing.Process(
                target=do.run,
                kwargs={
                    "cmd": commands_m,
                    "descr": "%s done" %
                    commands_m})
            p.start()
            jobs.append(p)

        while True in [j.is_alive() for j in jobs]:
            try:
                time.sleep(1)

            except KeyboardInterrupt:
                sys.stderr.write('KeyboardInterrupt detected, terminating '
                                 'all processes...')
                for j in jobs:
                    j.terminate()


                sys.exit(1)

        for job in jobs:

            job.join()
            if job.exitcode != 0:
                sys.stderr.write("Terminal the multiprocess now !\n")
                job.terminate()
                return

        sys.stderr.write("GATK Multiprocess has finished\n")

        mergevcf = [self.gatk, "MergeVcfs"] + \
            ["-I %s" % f for f in sub_vcf_files]
        mergevcf.append("-O %s" % output)
        do.run(
            " ".join(mergevcf),
            "Merge the multiproccess's output has done\n")

        return
