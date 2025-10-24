import os
import shutil
import argparse
import yaml
import time
import logging
import logging.handlers
import multiprocessing
from datetime import datetime, timedelta
import subprocess
import importlib.resources as resources
import pandas as pd
import numpy as np
import joblib

from threading import Thread, Event

from psutil import Process

from oligo_designer_toolsuite.sequence_generator import OligoSequenceGenerator

class MemoryMonitor(Thread):
    """Threaded memory monitor that prints process+children RSS (MB) every interval seconds."""

    def __init__(self, interval: float = 10.0):
        super().__init__(daemon=True)
        self._stop_event = Event()
        self.interval = float(interval)
        self.start()

    def get_memory(self) -> int:
        "Return RSS of current process plus children in bytes."
        p = Process()
        memory = p.memory_info().rss
        for c in p.children(recursive=True):
            try:
                memory += c.memory_info().rss
            except Exception:
                pass
        return memory

    def run(self):
        memory_start = self.get_memory()
        while not self._stop_event.is_set():
            mem = self.get_memory() - memory_start
            print(f"[MemoryMonitor] {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {mem/1024**2:.1f} MB", flush=True)
            time.sleep(self.interval)

    def stop(self):
        self._stop_event.set()


def filter_oligos(oligo_fasta_file: str, logger):
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
        
        logger.warning("Filter failed (rc={}):\nSTDOUT:\n{}\nSTDERR:\n{}".format(e.returncode, out, err))
        raise
    else:
        os.remove(oligo_fasta_file)
        return oligo_fasta_file_filtered
    
def determine_oligo_length(oligo_fasta_file: str, logger):
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
        
        logger.warning("Filter failed (rc={}):\nSTDOUT:\n{}\nSTDERR:\n{}".format(e.returncode, out, err))
        raise
    else:
        return oligo_fasta_file_length
    
def sample_oligos(oligo_fasta_file, oligo_fasta_file_length, sample_information, logger, config, region, seed):
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
    information_fewer_oligos = []
    for interval_name, interval_data in sample_information.items():
        # only check even line numbers as they contain the sequence
        temp_data = length_information[(length_information["line_number"] % 2 == 0) &
                                       (length_information["length"] >= interval_data["lower"]) &
                                       (length_information["length"] <= interval_data["upper"])]
        if temp_data.shape[0] < interval_data["n"]:
            logger.warning(f"The interval {interval_data['lower']} >= x <= {interval_data['upper']} "\
                           f"only has {temp_data.shape[0]} oligos instead of {interval_data['n']}")
            information_fewer_oligos.append([interval_data['lower'], interval_data['upper'],
                                             interval_data['n'], temp_data.shape[0]])
            sampled_data = temp_data
        else:
            sampled_data = temp_data.sample(n=interval_data["n"], random_state=seed)

        # determine the line numbers of the sampled oligos and the line numbers
        # of the FASTA headers
        for row in sampled_data.itertuples():
            index_with_header.extend([row.line_number - 1, row.line_number])
    index_with_header.sort()

    # write the sampled lines into a new FASTA file
    oligo_fasta_file_sampled = oligo_fasta_file.replace(".fna", "_sampled.fna")
    with open(oligo_fasta_file, "r") as infile:
        with open(oligo_fasta_file_sampled, "w") as outfile:
            for line_number, line in enumerate(infile):
                # line_number is 0-based, index_with_header 1-based
                if line_number + 1 in index_with_header:
                    outfile.write(line)

    # write logging of too few oligos
    if len(information_fewer_oligos) > 0:
        information_fewer_oligos_df = pd.DataFrame(information_fewer_oligos,
                                                   columns=["lower", "upper", "n_expected", "n_actual"])
        information_fewer_oligos_df.to_csv(os.path.join(config["storage_dir"], "logs", f"{region}_problems_sampling.csv"), index=False)

    os.remove(oligo_fasta_file)
    os.remove(oligo_fasta_file_length)
    return oligo_fasta_file_sampled

def sample_oligos_one_region(region, config, queue, seed):

    start = time.time()
    logger_ = logging.getLogger(__name__)
    if not logger_.hasHandlers():
        logger_.setLevel(logging.INFO)
        handler = logging.handlers.QueueHandler(queue)
        logger_.addHandler(handler)

    ################################
    # generate the oligo sequences #
    ################################

    logger_.info(f"Start processing {region}")
    temp_dir_output = f"/localscratch/jonas.hagenberg/output_sample_oligos_{region}_" + str(time.time())
    os.makedirs(temp_dir_output, exist_ok=True)
    
    files_fasta = [
        os.path.join(config["annotation_path"], f'gene_{config["annotation_file"]}'),
        os.path.join(config["annotation_path"], f'intergenic_{config["annotation_file"]}'),
    ]

    ##### creating the oligo sequences #####
    logger_.info("start generating oligo sequences")
    oligo_sequences = OligoSequenceGenerator(dir_output=temp_dir_output)
    
    oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
        files_fasta_in=files_fasta,
        length_interval_sequences=(15, 100),
        region_ids=region,
        stride=1,
        n_jobs=1,
    )
    try:
        oligo_fasta_file_filtered = filter_oligos(oligo_fasta_file[0], logger_)
        oligo_length = determine_oligo_length(oligo_fasta_file_filtered, logger_)
        oligo_file_sampled = sample_oligos(oligo_fasta_file_filtered, oligo_length, config["interval_config"], logger_, config, region, seed)

        shutil.copy(oligo_file_sampled, os.path.join(config["storage_dir"], os.path.basename(oligo_file_sampled)))
    except Exception as e:
        logger_.exception(e)
        logger_.warning(f"Could not sample oligos for {region}")
    
    end = time.time()
    elapsed_total = end - start
    logger_.info(f"Computational time: {str(timedelta(seconds=int(elapsed_total)))}")

    shutil.rmtree(temp_dir_output) #remove oligo designer toolsuite output


def main():
    """
    Sample oligos in a sliding window fashion from a certain region with lengths 15-100.
    """

    ##############
    # set logger #
    ##############

    timestamp = datetime.now()
    file_logger = f"log_sample_oligos_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt"
    logging.basicConfig(
        format="%(asctime)s %(processName)s [%(levelname)s] %(message)s",
        level=logging.INFO,
        handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
    )

    parser = argparse.ArgumentParser(
        prog="Sample oligos",
        usage="generate_real_dataset_sample_oligos [options]",
        description=main.__doc__,
    )
    parser.add_argument("-c", "--config", help="path to the configuration file", default="config/generate_real_dataset_sample_oligos.yaml")
    args = parser.parse_args()
    with open(args.config, "r") as handle:
        config = yaml.safe_load(handle)

    regions_df = pd.read_csv(config["file_regions"])
    regions = regions_df["region_id"].tolist()

    np.random.seed(config["seed"])
    list_of_seeds = np.random.randint(1e8, size=len(regions))

    os.makedirs(config["storage_dir"], exist_ok=True)
    os.makedirs(os.path.join(config["storage_dir"], "logs"), exist_ok=True)

    monitor = MemoryMonitor()

    queue = multiprocessing.Manager().Queue(-1) 
    listener = logging.handlers.QueueListener(queue, *logging.getLogger().handlers) 
    listener.start()

    sampled_oligos = joblib.Parallel(n_jobs=config["n_jobs"])(
        joblib.delayed(sample_oligos_one_region)(
            region=region_id,
            config=config,
            queue=queue,
            seed=seed
        )
        for region_id, seed in zip(regions, list_of_seeds)
    )

    listener.stop()
    monitor.stop()
    monitor.join()

if __name__ == "__main__":
    main()