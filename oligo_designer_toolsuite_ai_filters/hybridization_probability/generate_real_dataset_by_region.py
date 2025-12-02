import os
import shutil
import argparse
from pathlib import Path
import time
import yaml
import random
import numpy as np
import logging
from datetime import datetime, timedelta
import nupack
import joblib

from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BowtieFilter,
    BlastNFilter,
)

def generate_off_targets_region(
        oligo_fasta_file: str,
        config,
        file_reference,
        alignment_method: BlastNFilter,
    ):

    dir_output = "/localscratch/jonas.hagenberg/output_odt_real_blast_" + str(time.time())
    
    region_id = os.path.basename(oligo_fasta_file)
    region_id = region_id.replace("_filtered_sampled.fna", "")

    oligo_database = OligoDatabase(
        min_oligos_per_region=0,
        write_regions_with_insufficient_oligos=True,
        lru_db_max_in_memory=config["n_jobs"] + 1,
        database_name=f"oligo_database_{str(time.time())}",
        dir_output=dir_output,
    )

    oligo_database.load_database_from_fasta(
        files_fasta=oligo_fasta_file,
        sequence_type="oligo",
        region_ids=region_id,
        database_overwrite = True,
    )

    table_hits = alignment_method._run_filter(
        sequence_type='oligo',
        region_id=region_id,
        oligo_database=oligo_database,
        file_reference=file_reference,
        consider_hits_from_input_region=True,
        mode=2
    )

    # add the gaps
    references = alignment_method._get_references(table_hits, file_reference, region_id)
    queries = alignment_method._get_queries(oligo_database, table_hits, region_id, 'oligo')
    # unique_queries = list(set(queries))
    # align the references and queries by adding gaps
    gapped_queries, gapped_references = alignment_method._add_alignment_gaps(
        table_hits=table_hits, queries=queries, references=references
    )

    outfile = os.path.join(config["alignments_out_directory"], f"{region_id}_blast_results.csv")

    with open(outfile, "w") as f:
        f.write("gapped_query, gapped_referencen\n")
        for one_gapped_query, one_gapped_reference in zip(gapped_queries, gapped_references):
            f.write(f"{one_gapped_query}, {one_gapped_reference}\n")


    # create the output
    # targets = {}
    # on_targets = []
    # off_targets = []
    # sampled_temperatures = {}
    # with open(output_file, 'a') as file:
    #     file.write("start temp calculating for on-targets\n")
    # for query in unique_queries:
    #     temperatures = sample_temperatures(6)
    #     sampled_temperatures[query] = temperatures
    #     on_targets.extend(generate_datasamples(query, query, query, query, temperatures, 0))
    # with open(output_file, 'a') as file:
    #     file.write("start temp calculating for off-targets\n")
    # for query, reference, gapped_query, gapped_reference in zip(queries, references, gapped_queries, gapped_references):
    #     temperatures = sampled_temperatures[query]
    #     n_mismatches = sum(q != r for q, r in zip(gapped_query, gapped_reference))
    #     off_targets.extend(generate_datasamples(query, reference, gapped_query, gapped_reference, temperatures, n_mismatches))
    # return on_targets, off_targets

    shutil.rmtree(dir_output)
    return outfile


def main():
    # give every oligo a unique ID so blasted sequences can be traced back to the original oligo ?

    # 0. read in the oligos
    # 1. blast the oligos -> I tried to sample 2x oligos as needed per region, so there should be enough -> how to deal when there are not enough oligos?
    # 2. select 5 hits
    # 3. calculate NUPACK binding prop -> for 6 temperatures
    # 4. distribute to train/test/validation split

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
    # set random seed for reproducibility
    random.seed(config["seed"])
    rnd_gene_shuffling = np.random.RandomState(config["seed"] + 154872)
    # generate directories
    os.makedirs(config["alignments_out_directory"], exist_ok=True)
    # nupack run
    # nupack.config.cache = config["nupack_cache"]
    

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
    # generate the reference database #
    ################################

    logger.info("Generating reference database.")
    dir_output = "/localscratch/jonas.hagenberg/output_odt_reference_" + str(time.time())
    files_fasta = [
        os.path.join(config["annotation_path"], f'gene_{config["annotation_file"]}'),
        os.path.join(config["annotation_path"], f'intergenic_{config["annotation_file"]}'),
    ]
    reference_database = ReferenceDatabase(dir_output=dir_output)
    reference_database.load_database_from_file(files=files_fasta, file_type="fasta", database_overwrite = True,)
    file_reference = reference_database.write_database_to_file(
            filename=f"db_reference",
        )
    logger.info("Generated reference database.")

    ###########################
    # gather files with sampled oligos #
    ###########################

    oligo_files = [f for f in Path(config["oligo_files"]).rglob("*.fna") if f.is_file()]

    logger.info(f"Found {len(oligo_files)} files for processing.")

    # np.random.seed(config["seed"])
    # list_of_seeds = np.random.randint(1e8, size=len(oligo_files))

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


    logger.info("Start alignment.")
    blasted_oligos = joblib.Parallel(n_jobs=config["n_jobs"])(
        joblib.delayed(generate_off_targets_region)(
            oligo_fasta_file=one_oligo_file,
            config=config,
            file_reference=file_reference,
            alignment_method=alignment_method
        )
        for one_oligo_file in oligo_files
    )
    logger.info("Finish alignment.")

    shutil.rmtree(dir_output)


if __name__ == "__main__":
    main()