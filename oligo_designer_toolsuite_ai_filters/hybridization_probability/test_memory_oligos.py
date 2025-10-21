import os
import shutil
import argparse
import time
from typing import Union, List
import logging
from datetime import datetime, timedelta
import subprocess
import importlib.resources as resources
import pandas as pd
import numpy as np

from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_property_filter import PropertyFilter, HardMaskedSequenceFilter, SoftMaskedSequenceFilter

import psutil

def print_mem_usage():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 10**9  # in GB
    return f"Memory usage: {mem:.2f} GB\n"

def generate_oligos(n_jobs: int, dir_output: str, regions: list, oligo_fasta_file: Union[str, List[str]], logger):
    """Generate the oligo sequences.
    """

    ##### creating the oligo database #####
    # one database for train, test and validation is created
    oligo_database = OligoDatabase(
        min_oligos_per_region=0,
        write_regions_with_insufficient_oligos=True,
        lru_db_max_in_memory=n_jobs + 1,
        database_name=f"oligo_database_{str(time.time())}",
        dir_output=dir_output,
    )

    oligo_database.load_database_from_fasta(
        files_fasta=oligo_fasta_file,
        sequence_type="oligo",
        region_ids=regions,
        database_overwrite = True,
    )

    logger.info("Oligo database loaded into memory")
    logger.info(print_mem_usage())

    # Property filtering
    masked_seqeunces = HardMaskedSequenceFilter()
    soft_masked_seqeunces = SoftMaskedSequenceFilter()
    property_filter = PropertyFilter(filters=[masked_seqeunces, soft_masked_seqeunces])
    oligo_database = property_filter.apply(oligo_database=oligo_database, n_jobs=n_jobs, sequence_type="oligo")
    
    return oligo_database

def filter_oligos(oligo_fasta_file: str):
    """
    Implement the SoftMaskedSequenceFilter and HardMaskedSequenceFilter with awk/regex and directly work on the FASTA files,
    as ODT is currently not memory-efficient enough to work with such large files.
    """

    oligo_fasta_file_filtered = oligo_fasta_file.replace(".fna", "_filtered.fna")
    script_path = resources.files("oligo_designer_toolsuite_ai_filters") / "hybridization_probability" / "scripts" / "hard_soft_masked_sequence_filter.sh"
    try:
        completed = subprocess.run(
            ["bash", str(script_path), oligo_fasta_file, oligo_fasta_file_filtered],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=False,
        )
    except subprocess.CalledProcessError as e:
        out = e.stdout or ""
        err = e.stderr or ""
        
        print("Filter failed (rc={}):\nSTDOUT:\n{}\nSTDERR:\n{}".format(e.returncode, out, err))
        raise
    else:
        os.remove(oligo_fasta_file)
        return oligo_fasta_file_filtered
    
def determine_oligo_length(oligo_fasta_file: str):
    """
    Determine the oligo length by using awk.
    """

    oligo_fasta_file_length = oligo_fasta_file.replace(".fna", "_length.fna")
    script_path = resources.files("oligo_designer_toolsuite_ai_filters") / "hybridization_probability" / "scripts" / "determine_oligo_length.sh"
    try:
        completed = subprocess.run(
            ["bash", str(script_path), oligo_fasta_file, oligo_fasta_file_length],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=False,
        )
    except subprocess.CalledProcessError as e:
        out = e.stdout or ""
        err = e.stderr or ""
        
        print("Filter failed (rc={}):\nSTDOUT:\n{}\nSTDERR:\n{}".format(e.returncode, out, err))
        raise
    else:
        return oligo_fasta_file_length
    
def sample_oligos(oligo_fasta_file, oligo_fasta_file_length, sample_information, suffix, random_generator, logger):
    """
    Read in the oligo length file, sample according to the sample_information and create a new FASTA file.
    This is then read in as a oligo database. sample_information is a dict with the following structure:
        - interval: interval_name (arbitrary)
            - lower: lower bound (included in the interval)
            - upper: upper bound (included in the interval)
            - n: number of samples

    """
    length_information = pd.read_csv(oligo_fasta_file_length, names=["line_number", "length"])

    index_with_header = []
    for interval_name, interval_data in sample_information.items():
        # only check even line numbers as they contain the sequence
        temp_data = length_information[(length_information["line_number"] % 2 == 0) &
                                       (length_information["length"] >= interval_data["lower"]) &
                                       (length_information["length"] <= interval_data["upper"])]
        if temp_data.shape[0] < interval_data["n"]:
            logger.warning(f"The interval {interval_data['lower']} >= x <= {interval_data['upper']} "\
                           f"only has {temp_data.shape[0]} oligos instead of {interval_data['n']}")
            sampled_data = temp_data
        else:
            sampled_data = temp_data.sample(n=interval_data["n"], random_state=random_generator)

        # determine the line numbers of the sampled oligos and the line numbers
        # of the FASTA headers
        for row in sampled_data.itertuples():
            index_with_header.extend([row.line_number - 1, row.line_number])
    index_with_header.sort()

    # write the sampled lines into a new FASTA file
    oligo_fasta_file_sampled = oligo_fasta_file.replace(".fna", f"_{suffix}_sampled.fna")
    with open(oligo_fasta_file, "r") as infile:
        with open(oligo_fasta_file_sampled, "w") as outfile:
            for line_number, line in enumerate(infile):
                # line_number is 0-based, index_with_header 1-based
                if line_number + 1 in index_with_header:
                    outfile.write(line)

    os.remove(oligo_fasta_file)
    os.remove(oligo_fasta_file_length)
    return oligo_fasta_file_sampled


def main():
    """
    Test how much memory it needs to generate oligos in a sliding window fashion from a certain region with lengths 15-100.
    """

    #########################
    # read in out arguments #
    #########################

    start = time.time()
    parser = argparse.ArgumentParser(description="Test memory consumption of oligo databases")
    parser.add_argument("--region", action="store", dest="region", default="ABCA13")
    parser.add_argument("--n_jobs", action="store", dest="n_jobs", default=1, type=int)
    args = parser.parse_args()

    annotation_path = "/lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/output_odt_real_1757923915.7920792/annotation"
    annotation_file = "annotation_source-NCBI_species-Homo_sapiens_annotation_release-GCF_000001405.40-RS_2025_08_genome_assemly-unknown.fna"
    
    rng = np.random.default_rng(274390)

    interval_config_1 = {
        "interval_1": {"lower": 15, "upper": 20, "n": 40},
        "interval_2": {"lower": 21, "upper": 30, "n": 80},
        "interval_3": {"lower": 31, "upper": 40, "n": 80},
        "interval_4": {"lower": 41, "upper": 50, "n": 80},
        }
    
    interval_config_2 = {
        "interval_5": {"lower": 51, "upper": 60, "n": 24},
        "interval_6": {"lower": 61, "upper": 70, "n": 24},
        "interval_7": {"lower": 71, "upper": 80, "n": 24},
        "interval_8": {"lower": 81, "upper": 90, "n": 24},
        "interval_9": {"lower": 91, "upper": 100, "n": 24},
        }

    ##############
    # set logger #
    ##############

    timestamp = datetime.now()
    file_logger = f"log_memory_oligo_creation_{args.region}_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt"
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.INFO,
        handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
    )
    logger = logging.getLogger("memory test")
    logger.info("startup")
    logger.info(f"region: {args.region}")
    logger.info(f"n_jobs: {args.n_jobs}")
    logger.info(print_mem_usage())


    ################################
    # generate the oligo sequences #
    ################################

    dir_output_1 = "/localscratch/jonas.hagenberg/output_test_odt_memory_1_" + str(time.time())
    dir_output_2 = "/localscratch/jonas.hagenberg/output_test_odt_memory_2_" + str(time.time())
    # dir_output = "output_test_odt_memory_" + str(time.time())
    os.makedirs(dir_output_1, exist_ok=True)
    os.makedirs(dir_output_2, exist_ok=True)

    
    files_fasta = [
        os.path.join(annotation_path, f'gene_{annotation_file}'),
        os.path.join(annotation_path, f'intergenic_{annotation_file}'),
    ]

    ##### creating the oligo sequences #####
    logger.info("start generating oligo sequences")
    oligo_sequences_1 = OligoSequenceGenerator(dir_output=dir_output_1)
    oligo_sequences_2 = OligoSequenceGenerator(dir_output=dir_output_2)
    
    oligo_fasta_file_1 = oligo_sequences_1.create_sequences_sliding_window(
        files_fasta_in=files_fasta,
        length_interval_sequences=(15, 50),
        region_ids=args.region,
        stride=1,
        n_jobs=args.n_jobs,
    )
    logger.info(oligo_fasta_file_1)
    logger.info("sequences 15-50 generated")
    logger.info(print_mem_usage())
    # the second run appends the data to the same FASTA file
    oligo_fasta_file_2 = oligo_sequences_2.create_sequences_sliding_window(
        files_fasta_in=files_fasta,
        length_interval_sequences=(51, 100),
        region_ids=args.region,
        stride=1,
        n_jobs=args.n_jobs,
    )
    logger.info(oligo_fasta_file_2)
    logger.info("sequences 51-100 generated")
    logger.info(print_mem_usage())

    logger.info("filter oligos")
    oligo_fasta_file_filtered_1 = filter_oligos(oligo_fasta_file_1[0])
    logger.info(oligo_fasta_file_filtered_1)
    oligo_fasta_file_filtered_2 = filter_oligos(oligo_fasta_file_2[0])
    logger.info(oligo_fasta_file_filtered_2)
    logger.info("oligos filtered")
    logger.info(print_mem_usage())

    logger.info("determine oligo length")
    oligo_length_1 = determine_oligo_length(oligo_fasta_file_filtered_1)
    oligo_length_2 = determine_oligo_length(oligo_fasta_file_filtered_2)
    logger.info("oligo length determined")
    logger.info(print_mem_usage())

    logger.info("sample oligos")
    oligo_file_sampled_1 = sample_oligos(oligo_fasta_file_filtered_1, oligo_length_1, interval_config_1, "1", rng, logger)
    oligo_file_sampled_2 = sample_oligos(oligo_fasta_file_filtered_2, oligo_length_2, interval_config_2, "2", rng, logger)
    logger.info("oligos sampled")
    logger.info(print_mem_usage())

    logger.info("move sampled data")
    new_dir = "output_test_odt_memory_" + str(time.time())
    os.makedirs(new_dir)
    shutil.copy(oligo_file_sampled_1, os.path.join(new_dir, os.path.basename(oligo_file_sampled_1)))
    shutil.copy(oligo_file_sampled_2, os.path.join(new_dir, os.path.basename(oligo_file_sampled_2)))
    logger.info("sampled data moved")
    

    ################################
    # generate the reference database #
    ################################

    logger.info("Generating reference database.")
    reference_database = ReferenceDatabase(dir_output=dir_output_1)
    reference_database.load_database_from_file(files=files_fasta, file_type="fasta", database_overwrite = True,)
    file_reference = reference_database.write_database_to_file(
            filename=f"db_reference",
        )
    logger.info("Generated reference database.")
    logger.info(print_mem_usage())
    
    end = time.time()
    elapsed_total = end - start
    logger.info(f"Computational time: {str(timedelta(seconds=int(elapsed_total)))}")

    shutil.rmtree(dir_output_1) #remove oligo designer toolsuite output
    shutil.rmtree(dir_output_2) #remove oligo designer toolsuite output

# def main():
#     rng = np.random.default_rng(274390)
#     oligo_file = "ABCA13_small_mod.fna"
#     oligo_file_filtered = filter_oligos(oligo_file)
#     oligo_length = determine_oligo_length(oligo_file_filtered)
#     oligo_file_sampled = sample_oligos(oligo_file_filtered, oligo_length, {"interval_1": {"lower": 15, "upper": 15, "n": 10}}, rng)
#     print(oligo_file_sampled)

if __name__ == "__main__":
    main()