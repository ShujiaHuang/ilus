import sys
import yaml


def load_config(config_file):
    with open(config_file) as f:
        return yaml.safe_load(f)


if __name__ == "__main__":
    config = load_config(sys.argv[1])
    vqsr = " ".join([str(x) for x in config["sentieon"].get("vqsr_options", [])])
    adapter_seq = ",".join([str(x) for x in config["resources"].get("adapter_sequence", [])])
    print(config, "\n")
    print(vqsr)
    print(adapter_seq)
