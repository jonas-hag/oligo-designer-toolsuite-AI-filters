import os
import shutil
import argparse
import time
from typing import Union, List
import logging
from datetime import datetime, timedelta

from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_property_filter import PropertyFilter, HardMaskedSequenceFilter, SoftMaskedSequenceFilter

import psutil

def print_mem_usage():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 10**9  # in GB
    return f"Memory usage: {mem:.2f} GB\n"

def generate_oligos(n_jobs: int, dir_output: str, regions: list, oligo_fasta_file: Union[str, List[str]]):
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

    # Property filtering
    masked_seqeunces = HardMaskedSequenceFilter()
    soft_masked_seqeunces = SoftMaskedSequenceFilter()
    property_filter = PropertyFilter(filters=[masked_seqeunces, soft_masked_seqeunces])
    oligo_database = property_filter.apply(oligo_database=oligo_database, n_jobs=n_jobs, sequence_type="oligo")
    
    return oligo_database

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
    parser.add_argument("--n_jobs", action="store", dest="n_jobs", default=1)
    args = parser.parse_args()

    annotation_path = "/lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/output_odt_real_1757923915.7920792/annotation"
    annotation_file = "annotation_source-NCBI_species-Homo_sapiens_annotation_release-GCF_000001405.40-RS_2025_08_genome_assemly-unknown.fna"
    

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

    # dir_output = "/localscratch/jonas.hagenberg/output_odt_real_" + str(time.time())
    dir_output = "output_test_odt_memory_" + str(time.time())
    os.makedirs(dir_output, exist_ok=True)

    
    files_fasta = [
        os.path.join(annotation_path, f'gene_{annotation_file}'),
        os.path.join(annotation_path, f'intergenic_{annotation_file}'),
    ]

    ##### creating the oligo sequences #####
    logger.info("start generating oligo sequences")
    oligo_sequences = OligoSequenceGenerator(dir_output=dir_output)
    
    oligo_fasta_file_1 = oligo_sequences.create_sequences_sliding_window(
        files_fasta_in=files_fasta,
        length_interval_sequences=(15, 50),
        region_ids=args.region,
        stride=1,
        n_jobs=args.n_jobs,
    )
    logger.info("sequences 15-50 generated")
    logger.info(print_mem_usage())
    oligo_fasta_file_2 = oligo_sequences.create_sequences_sliding_window(
        files_fasta_in=files_fasta,
        length_interval_sequences=(51, 100),
        region_ids=args.region,
        stride=2,
        overwrite=False,
        n_jobs=args.n_jobs,
    )
    oligo_fasta_files = oligo_fasta_file_1 + oligo_fasta_file_2
    logger.info("sequences 51-100 generated")
    logger.info(print_mem_usage())

    logger.info("Generating Oligo sequences.")
    oligo_database = generate_oligos(args.n_jobs, dir_output, args.region, oligo_fasta_files)
    logger.info("oligo data base generated")
    logger.info(print_mem_usage())
    

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
    logger.info(print_mem_usage())
    
    end = time.time()
    elapsed_total = end - start
    logger.info(f"Computational time: {str(timedelta(seconds=int(elapsed_total)))}")

    shutil.rmtree(dir_output) #remove oligo designer toolsuite output

if __name__ == "__main__":
    main()