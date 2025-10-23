import os
import argparse
import yaml
from oligo_designer_toolsuite.pipelines import GenomicRegionGenerator

def main():

    parser = argparse.ArgumentParser(
        prog="Real Dataset Reference DB Creation",
        usage="generate_real_dataset_reference_db [options]",
        description=main.__doc__,
    )
    parser.add_argument("-c", "--config", help="path to the configuration file", default="config/generate_real_dataset_reference_db.yaml")
    args = parser.parse_args()
    with open(args.config, "r") as handle:
        config = yaml.safe_load(handle)

    os.makedirs(config["dir_output"], exist_ok=True)

    genomic_region_genereator = GenomicRegionGenerator(dir_output=config["dir_output"])
    region_generator = genomic_region_genereator.load_annotations(source=config["source"], source_params=config["source_params"])
    files_fasta = genomic_region_genereator.generate_genomic_regions(
        region_generator=region_generator,
        genomic_regions=config["genomic_regions"],
        block_size=0,
    )
    
if __name__ == "__main__":
    main()