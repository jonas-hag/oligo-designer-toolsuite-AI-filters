import os
import argparse
import yaml
import time
import logging
from datetime import datetime

from oligo_designer_toolsuite.database import ReferenceDatabase
from oligo_designer_toolsuite.pipelines import GenomicRegionGenerator
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BowtieFilter,
    BlastNFilter,
)


def main():
    """
    Generate a bowtie index.
    """

    #########################
    # read in out arguments #
    #########################

    start = time.time()
    parser = argparse.ArgumentParser(
        prog="Real Dataset",
        usage="generate_real_dataset [options]",
        description=main.__doc__,
    )
    parser.add_argument("-c", "--config", help="path to the configuration file", default="config/generate_real_dataset_blastn.yaml")
    args = parser.parse_args()
    with open(args.config, "r") as handle:
        config = yaml.safe_load(handle)
    dataset_name = f"real_dataset_{config['alignment_method']}"
    

    ##############
    # set logger #
    ##############

    timestamp = datetime.now()
    file_logger = f"log_{dataset_name}_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt"
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.INFO,
        handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
    )
    logger = logging.getLogger("real_dataset_generation")

    ################################
    # generate the oligo sequences #
    ################################

    dir_output = "/home/haicu/jonas.hagenberg/projects/odt-ai/01_data/bowtie" + str(time.time())


    genomic_region_genereator = GenomicRegionGenerator(dir_output = dir_output)
    region_generator = genomic_region_genereator.load_annotations(source=config["source"], source_params=config["source_params"])
    files_fasta = genomic_region_genereator.generate_genomic_regions(
        region_generator = region_generator,
        genomic_regions  = config["genomic_regions"],
        block_size = 0,
    )

    ################################
    # generate the reference database #
    ################################

    logger.info("Generating reference database.")
    reference_database = ReferenceDatabase(dir_output=dir_output)
    reference_database.load_database_from_file(files=files_fasta, file_type="fasta", database_overwrite = True,)
    file_reference = reference_database.write_database_to_file(
            filename=f"db_reference",
        )
    logger.info("Generated reference database.")
    # log database information
    

    ################################################################
    # generate real off-targets and compute duplexing scores #
    ################################################################
    
    start_2 = time.time()
    if config["alignment_method"] == "blastn":
        alignment_method = BlastNFilter(
            search_parameters = config["search_parameters"],
            hit_parameters = config["hit_parameters"],
            names_search_output = config["names_search_output"],
            dir_output=dir_output
        )
    elif config["alignment_method"] == "bowtie":
        alignment_method = BowtieFilter(
            search_parameters = config["search_parameters"],
            dir_output=dir_output
        )
    else:
        raise ValueError("Unknown alignment method.")
    
    logger.info("Generating file index")

    alignment_method.set_reference_database(reference_database=reference_database)
    file_index = alignment_method.create_reference(n_jobs=config["n_jobs"])

    logger.info("Generated file index.")


if __name__ == "__main__":
    main()