"""Work with sambamba
"""
from ilus.launch import do

class SambambaRunner(object):
    def __init__(self, config):
        self.sambamba = config["sambamba"]["exe"]
        self.num_cores = config["sambamba"]["num_cores"]

    def run(self, tools, options):
        cmd = " ".join([self.sambamba, tools] + options)
        do.run(cmd, "sambamba %s done" % tools)

        return
