import os
import shutil
import argparse
from pathlib import Path
import time
import yaml
import logging
import logging.handlers
import multiprocessing
from datetime import datetime
import joblib
import re

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
        queue
    ):

    logger_ = logging.getLogger(__name__)
    if not logger_.hasHandlers():
        logger_.setLevel(logging.INFO)
        handler = logging.handlers.QueueHandler(queue)
        logger_.addHandler(handler)

    dir_output = "/localscratch/jonas.hagenberg/output_odt_real_blast_" + str(time.time())
    # for local testing
    # dir_output = "output_odt_real_blast_" + str(time.time())
    
    region_id = os.path.basename(oligo_fasta_file)
    region_id_match = re.search(r'^(.*?)(?=_filtered_sampled)', region_id)
    region_id = region_id_match.group(1)

    logger_.info(f"region: {region_id}")

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
    # align the references and queries by adding gaps
    gapped_queries, gapped_references = alignment_method._add_alignment_gaps(
        table_hits=table_hits, queries=queries, references=references
    )

    outfile = os.path.join(config["alignments_out_directory"], f"{region_id}_blast_results.csv")

    with open(outfile, "w") as f:
        f.write("query, gapped_query, gapped_reference\n")
        for one_query, one_gapped_query, one_gapped_reference in zip(queries, gapped_queries, gapped_references):
            f.write(f"{one_query}, {one_gapped_query}, {one_gapped_reference}\n")

    shutil.rmtree(dir_output)
    return outfile


def main():
    # give every oligo a unique ID so blasted sequences can be traced back to the original oligo ?

    # 0. read in the oligos
    # 1. blast the oligos -> I tried to sample 2x oligos as needed per region, so there should be enough -> how to deal when there are not enough oligos?
    #### the following steps will be performed in separate scripts
    # 2. select 5 hits
    # 3. calculate NUPACK binding prop -> for 6 temperatures
    # 4. distribute to train/test/validation split

    #########################
    # read in out arguments #
    #########################

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
    # generate directories
    os.makedirs(config["alignments_out_directory"], exist_ok=True)
    

    ##############
    # set logger #
    ##############

    timestamp = datetime.now()
    file_logger = f"log_{dataset_name}_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt"
    logging.basicConfig(
        format="%(asctime)s %(processName)s [%(levelname)s] %(message)s",
        level=logging.INFO,
        handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
    )

    ################################
    # generate the reference database #
    ################################

    dir_output = "/localscratch/jonas.hagenberg/output_odt_reference_" + str(time.time())
    # for local testing
    # dir_output = "output_odt_reference_" + str(time.time())
    files_fasta = [
        os.path.join(config["annotation_path"], f'gene_{config["annotation_file"]}'),
        os.path.join(config["annotation_path"], f'intergenic_{config["annotation_file"]}'),
    ]
    reference_database = ReferenceDatabase(dir_output=dir_output)
    reference_database.load_database_from_file(files=files_fasta, file_type="fasta", database_overwrite = True,)
    file_reference = reference_database.write_database_to_file(
            filename=f"db_reference",
        )

    ###########################
    # gather files with sampled oligos #
    ###########################

    oligo_files = [f for f in Path(config["oligo_files"]).glob(config["suffix"]) if f.is_file()]

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
    
    alignment_method.set_reference_database(reference_database=reference_database)
    file_index = alignment_method.create_reference(n_jobs=config["n_jobs"])

    queue = multiprocessing.Manager().Queue(-1) 
    listener = logging.handlers.QueueListener(queue, *logging.getLogger().handlers) 
    listener.start()

    blasted_oligos = joblib.Parallel(n_jobs=config["n_jobs"])(
        joblib.delayed(generate_off_targets_region)(
            oligo_fasta_file=one_oligo_file,
            config=config,
            file_reference=file_index,
            alignment_method=alignment_method,
            queue=queue
        )
        for one_oligo_file in oligo_files
    )
    listener.stop()

    shutil.rmtree(dir_output)


if __name__ == "__main__":
    main()