import pathlib
import os
import subprocess
import importlib.resources as resources

def main():
    # find all files
    base_path = "/lustre/groups/aiconsultants/projects/odt-ai/oligo-designer-toolsuite-AI-filters/data/hybridization_probability/real_dataset_sampled_oligos"
    new_out_dir = os.path.join(base_path, "split_fasta_files")
    os.makedirs(new_out_dir, exist_ok=True)
    all_files = list(pathlib.Path(base_path).glob("*.fna"))
    script_path = resources.files("oligo_designer_toolsuite_ai_filters") / "hybridization_probability" / "scripts" / "split_fasta_file_by_length.sh"

    for file in all_files:
        file_name = file.stem
        file_name_short = os.path.join(new_out_dir, f"{file_name}_short.fna")
        file_name_long = os.path.join(new_out_dir, f"{file_name}_long.fna")

        try:
            completed = subprocess.run(
                ["bash", str(script_path), file, file_name_short, file_name_long, "30"],
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


if __name__ == "__main__":
    main()