import os
import shutil
import argparse
import yaml
import random
import copy
import time
from typing import Tuple, Union, List
import logging
from datetime import datetime
import iteration_utilities
from collections import Counter

from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator
from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_property_filter import PropertyFilter, HardMaskedSequenceFilter, SoftMaskedSequenceFilter
from oligo_designer_toolsuite.pipelines import GenomicRegionGenerator
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BowtieFilter,
    BlastNFilter,
)
from Bio.Seq import MutableSeq, Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import nupack
from math import log, ceil
import joblib
import numpy as np

from  oligo_designer_toolsuite_ai_filters.hybridization_probability.generate_artificial_dataset import split_genes_stratified, generate_datasamples, generate_dataset, sample_temperatures


base_pair = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'T', 't':'A', 'c':'G', 'g':'C'}


def generate_off_targets_region(
        oligo_database: OligoDatabase, 
        alignment_method: BlastNFilter, 
        file_index: str, region_id: str, 
        file_reference: str,
        sampled_oligos_per_region: int
    ):
    """Return a list of all the retrived off target sites in the following format:


    :param oligo_database: _description_
    :type oligo_database: OligoDatabase
    :param alignment_method: _description_
    :type alignment_method: BlastNFilter
    :param file_index: _description_
    :type file_index: str
    :param region_id: _description_
    :type region_id: str
    :param file_reference: _description_
    :type file_reference: str
    :param sampled_oligos_per_region: how many oligos to sample from each region before running the filter
    :type sampled_oligos_per_region: int
    """

    # downsample the oligo_database for better efficiency
    # assume 5 off-target hits per oligo
    number_regions = oligo_database.database.keys()

    # run the filter
    filtered_oligo_database = copy.deepcopy(oligo_database)
    oligo_ids = filtered_oligo_database.get_oligoid_list()
    oligo_id_sample = random.sample(population=oligo_ids, k=min(sampled_oligos_per_region, len(oligo_ids)))
    # filtered_oligo_database.filter_database_by_oligo(remove_region=False, oligo_ids=oligo_id_sample)

    table_hits = alignment_method._run_filter(
        sequence_type='oligo',
        region_id=region_id,
        oligo_database=filtered_oligo_database,
        file_reference=file_index,
        consider_hits_from_input_region=True,
        mode=2
    )

    # add the gaps
    references = alignment_method._get_references(table_hits, file_reference, region_id)
    queries = alignment_method._get_queries(filtered_oligo_database, table_hits, region_id, 'oligo')
    unique_queries = list(set(queries))
    # align the references and queries by adding gaps
    gapped_queries, gapped_references = alignment_method._add_alignment_gaps(
        table_hits=table_hits, queries=queries, references=references
    )

    # check how many off-target hits an oligo has
    # queries contains the oligo, one element for every off-target found
    print(region_id)
    print(Counter(queries))

    # create the output
    off_targets = []
    for query in unique_queries:
        off_targets.extend(generate_datasamples(query, query, query, query, sample_temperatures(6),0))
    for query, reference, gapped_query, gapped_reference in zip(queries, references, gapped_queries, gapped_references):
        n_mismatches = sum(q != r for q, r in zip(gapped_query, gapped_reference))
        off_targets.extend(generate_datasamples(query, reference, gapped_query, gapped_reference, sample_temperatures(2), n_mismatches))
    return off_targets


def generate_off_targets(
        oligo_database: OligoDatabase, 
        alignment_method: BlastNFilter, 
        file_index: str, 
        config: dict, 
        dataset_size: int, 
        file_reference: str
    ):

    # the oligo_database will be downsampled for better efficiency
    # for this, we need to calculate how many oligos we want to sample from each region
    # assume 5 off-target hits per oligo
    number_regions = oligo_database.database.keys()
    sampled_oligos_per_region = int(ceil(dataset_size / (5 * len(number_regions))))

    off_target_regions = joblib.Parallel(n_jobs=config["n_jobs"])(
        joblib.delayed(generate_off_targets_region)(
            oligo_database=oligo_database,
            alignment_method=alignment_method,
            file_index=file_index,
            region_id=region_id,
            file_reference=file_reference,
            sampled_oligos_per_region=sampled_oligos_per_region
        )
        for region_id in oligo_database.database.keys()
    )

    # off_target_regions = []
    # for region_id in oligo_database.database.keys():
    #     results = generate_off_targets_region(
    #         oligo_database=oligo_database,
    #         alignment_method=alignment_method,
    #         file_index=file_index,
    #         region_id=region_id,
    #         file_reference=file_reference,
    #         sampled_oligos_per_region=sampled_oligos_per_region
    #     )
    #     off_target_regions.append(results)

    # flatten the list
    off_target_regions = [off_target for region in off_target_regions for off_target in region]
    
    # sample
    if dataset_size > len(off_target_regions):
        logging.warning(f"Fewer oligos left to sample: {len(off_target_regions)} instead of {dataset_size}.")
    off_target_regions = random.sample(population=off_target_regions, k=min(dataset_size, len(off_target_regions)))

    # create dataset
    dataset = generate_dataset(alignments=off_target_regions)

    return dataset



def sample_oligos(oligo_database: OligoDatabase, oligos_per_region: int):
    for region in oligo_database.database.keys():
        oligo_ids = list(oligo_database.database[region].keys())
        if len(oligo_ids) > oligos_per_region:
            filtered_oligo_ids = random.sample(population=oligo_ids, k=len(oligo_ids) - oligos_per_region) # sample the ids to filter
            for oligo_id in filtered_oligo_ids:
                oligo_database.database[region].pop(oligo_id, None)
    return oligo_database


def generate_oligos(config: dict, dir_output: str, regions: list, oligo_fasta_file: Union[str, List[str]]):
    """Generate the oligo sequences.
    """

    ##### creating the oligo database #####
    # one database for train, test and validation is created
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
        region_ids=regions,
        database_overwrite = True,
    )

    # Property filtering
    masked_seqeunces = HardMaskedSequenceFilter()
    soft_masked_seqeunces = SoftMaskedSequenceFilter()
    property_filter = PropertyFilter(filters=[masked_seqeunces, soft_masked_seqeunces])
    oligo_database = property_filter.apply(oligo_database=oligo_database, n_jobs=config["n_jobs"], sequence_type="oligo")
    
    return oligo_database
    
    




def main():
    """
    Generate a real dataset containing oligos and real off-targets identified by a given alignment method (e.g. blastn, bowtie). 
    The oligos are extracted from a given list of genetic regions and sampled to match the desired dataset size.
    Oligos with lengths up to 50 nucleotides are created using a sliding window approach with a stride of 1,
    while longer oligos are created with a stride of 2.
    The genetic regions contain both genes and intergenic regions.
    These oligos are then aligned to their corresponding reference genome and real potential off-target regions identified.

    For every oligo, the free energy of binding to its on-target region and the off-target regions is computed using NUPACK.
    For the binding to its on-target region, 6 different temperatures are samples, with a focus on 35-65 degrees Celsius,
    the total range is 20-90 degrees Celsius. For the binding to the off-target regions, 2 temperatures are sampled,
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
    dataset_name = f"real_dataset_{config['alignment_method']}_{config['dataset_size']}_{config['oligo_length_min']}_{config['oligo_length_max']}_{config['dataset_size']}"
    # set random seed for reproducibility
    random.seed(config["seed"])
    rnd_gene_shuffling = np.random.RandomState(config["seed"] + 154872)
    # generate directories
    os.makedirs(config["dir_output"], exist_ok=True)
    plots_dir = os.path.join(config["dir_output"], f"{dataset_name}_plots")
    os.makedirs(plots_dir, exist_ok=True)
    # nupack run
    # nupack.config.threads = config["n_jobs"] # use all cores
    nupack.config.cache = config["nupack_cache"]
    

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

    dir_output = "output_odt_real_" + str(time.time())

    if config["precalculated_annotation_path"] is not None and config["precalculated_annotation_file"] is not None:
        files_fasta = [
            os.path.join(config["precalculated_annotation_path"], f'gene_{config["precalculated_annotation_file"]}'),
            os.path.join(config["precalculated_annotation_path"], f'intergenic_{config["precalculated_annotation_file"]}'),
        ]
    else:
        genomic_region_genereator = GenomicRegionGenerator(dir_output = dir_output)
        region_generator = genomic_region_genereator.load_annotations(source=config["source"], source_params=config["source_params"])
        files_fasta = genomic_region_genereator.generate_genomic_regions(
            region_generator = region_generator,
            genomic_regions  = config["genomic_regions"],
            block_size = 0,
        )

    

    # the gene list now also contains intergenic regions
    # genes_df must contain the columns region_id and region_type
    genes_df = pd.read_csv(config["file_genes"])
    genes = genes_df["region_id"].tolist()
    genes_train, genes_validation, genes_test = split_genes_stratified(genes_df, config["splits_size"], rnd_gene_shuffling)
    logger.info(f"Length of genes: {len(genes)}")
    logger.info(f"Length of genes train: {len(genes_train)}")
    logger.info(f"Length of genes validation: {len(genes_validation)}")
    logger.info(f"Length of genes test: {len(genes_test)}")

    ##### creating the oligo sequences #####
    oligo_sequences = OligoSequenceGenerator(dir_output=dir_output)
    if config["oligo_length_max"] > 50:
        oligo_fasta_file_1 = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta,
            length_interval_sequences=(config["oligo_length_min"], 50),
            region_ids=genes,
            stride=1,
            n_jobs=config["n_jobs"],
        )
        oligo_fasta_file_2 = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta,
            length_interval_sequences=(51, config["oligo_length_max"]),
            region_ids=genes,
            stride=2,
            overwrite=False,
            n_jobs=config["n_jobs"],
        )
        oligo_fasta_files = oligo_fasta_file_1 + oligo_fasta_file_2
    else:
        oligo_fasta_files = oligo_sequences.create_sequences_sliding_window(
            files_fasta_in=files_fasta,
            length_interval_sequences=(config["oligo_length_min"], config["oligo_length_max"]),
            region_ids=genes,
            stride=1,
            n_jobs=config["n_jobs"],
        )

    logger.info("Generating Oligo sequences.")
    oligo_database_train = generate_oligos(config, dir_output, genes_train, oligo_fasta_files)
    oligo_database_train = sample_oligos(oligo_database=oligo_database_train, oligos_per_region=config["oligos_per_region"])
    logger.info("Training set:")
    for gene in oligo_database_train.database.keys():
        logging.info(f"Gene {gene} has {len(oligo_database_train.database[gene].keys())} oligos.")
    oligo_database_validation = generate_oligos(config, dir_output, genes_validation, oligo_fasta_files)
    oligo_database_validation = sample_oligos(oligo_database=oligo_database_validation, oligos_per_region=config["oligos_per_region"])
    logger.info("Validation set:")
    for gene in oligo_database_validation.database.keys():
        logging.info(f"Gene {gene} has {len(oligo_database_validation.database[gene].keys())} oligos.")
    oligo_database_test = generate_oligos(config, dir_output, genes_test, oligo_fasta_files)
    oligo_database_test = sample_oligos(oligo_database=oligo_database_test, oligos_per_region=config["oligos_per_region"])
    logger.info("Test set:")
    for gene in oligo_database_test.database.keys():
        logger.info(f"Gene {gene} has {len(oligo_database_test.database[gene].keys())} oligos.")

    logger.info("Generated oligos.")

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
    if config["index_path"] is None:
        alignment_method.set_reference_database(reference_database=reference_database)
        file_index = alignment_method.create_reference(n_jobs=config["n_jobs"])
    else:
        file_index = config["index_path"]
    logger.info("Generated file index.")
    
    # sample the oligos
    sample_train = round(config["splits_size"][0]*config["dataset_size"])
    sample_validation = round(config["splits_size"][1]*config["dataset_size"])
    sample_test= config["dataset_size"] - sample_train - sample_validation

    # train
    logger.info("Generating real off-targets for the training set.")
    train_dataset = generate_off_targets(
        oligo_database= oligo_database_train, 
        alignment_method = alignment_method, 
        file_index = file_index, 
        config=config, 
        dataset_size=sample_train, 
        file_reference=file_reference
    )
    logger.info("Generated real off-targets for the training set.")

    # validation
    logger.info("Generating real off-targets for the validation set.")
    validation_dataset = generate_off_targets(
        oligo_database = oligo_database_validation, 
        alignment_method = alignment_method, 
        file_index = file_index, 
        config=config,
        dataset_size=sample_validation, 
        file_reference=file_reference
    )
    logger.info("Generated real off-targets for the validation set.")

    # test
    logger.info("Generating real off-targets for the test set.")
    test_dataset = generate_off_targets(
        oligo_database = oligo_database_test, 
        alignment_method = alignment_method, 
        file_index = file_index, 
        config=config,
        dataset_size=sample_test, 
        file_reference=file_reference
    )
    logger.info("Generated real off-targets for the test set.")

    ##################
    # write dataset #
    ##################

    file_train = os.path.join(config["dir_output"], f"{dataset_name}_train.csv")
    train_dataset.to_csv(file_train)
    file_validation = os.path.join(config["dir_output"], f"{dataset_name}_validation.csv")
    validation_dataset.to_csv(file_validation)
    file_test = os.path.join(config["dir_output"], f"{dataset_name}_test.csv")
    test_dataset.to_csv(file_test)
    logger.info(f"Dataset created and stored at: \n\t - {file_train},\n\t - {file_validation}, \n\t - {file_test}.")
    # plot distributions of the ground truths
    plt.figure(3)
    train_dataset["Source"] = "Train"
    validation_dataset["Source"] = "Validation"
    test_dataset["Source"] = "Test"
    dataset = pd.concat([train_dataset, validation_dataset, test_dataset])
    dataset["Free Energy"] = dataset["free_energy"] # rename for better understanding
    sns.boxplot(data=dataset, x="Source", y="Free Energy")
    plt.title("Free Energy distributions")
    plt.savefig(os.path.join(plots_dir,"Free_Energy_distribution.pdf"))
    
    # plot distributions of the n mismatches
    plt.figure(4)
    dataset["Number Mismatches"] = dataset["number_mismatches"]
    sns.boxplot(data=dataset, x="Source", y="Number Mismatches")
    plt.title("Number Mismatches distributions")
    plt.savefig(os.path.join(plots_dir,"Number_mismatches_distribution.pdf"))
    
    logger.info(f"Computational time: {time.time() - start} (off-targets generation: {time.time() - start_2})")
    del oligo_database_train
    del oligo_database_validation
    del oligo_database_test

    shutil.rmtree(dir_output) #remove oligo designer toolsuite output

if __name__ == "__main__":
    main()