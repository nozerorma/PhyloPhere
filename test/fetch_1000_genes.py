import os
import subprocess
import random

def main():
    print("Fetching lists of files from marvin2...")
    
    # 1. Fetch remote alignments list
    cmd_align = [
        'ssh', 'marvin2',
        'ls -1 /homes/users/mramon/scratch/4.Generate_alignments_from_codons/results_ppga/20260615_200027/PROT/bmge/'
    ]
    res_align = subprocess.run(cmd_align, capture_output=True, text=True, check=True)
    all_alignments = [line.strip() for line in res_align.stdout.splitlines() if line.strip().endswith('.fa')]

    # 2. Fetch remote maps list
    cmd_map = [
        'ssh', 'marvin2',
        'ls -1 /homes/users/mramon/scratch/4.Generate_alignments_from_codons/results_ppga/20260615_200027/MAP/'
    ]
    res_map = subprocess.run(cmd_map, capture_output=True, text=True, check=True)
    all_maps = set(line.strip() for line in res_map.stdout.splitlines() if line.strip())

    # 3. Fetch remote ASR list
    cmd_asr = [
        'ssh', 'marvin2',
        'ls -1 /homes/users/mramon/scratch/2.Primates/1.Primates_data/asr/'
    ]
    res_asr = subprocess.run(cmd_asr, capture_output=True, text=True, check=True)
    all_asrs = set(line.strip() for line in res_asr.stdout.splitlines() if line.strip())

    print(f"Found {len(all_alignments)} alignments, {len(all_maps)} map files, and {len(all_asrs)} ASR directories.")

    # Filter to find the intersection where map and ASR both exist
    valid_filenames = []
    for name in all_alignments:
        prefix = name.rsplit('.', 1)[0]
        gene = prefix.split('.')[0]
        map_name = f"{prefix}.map.tsv"
        asr_name = f"asr_{gene}"
        if map_name in all_maps and asr_name in all_asrs:
            valid_filenames.append(name)

    print(f"Number of genes with valid alignments, maps, and ASR: {len(valid_filenames)}")

    if len(valid_filenames) < 1000:
        raise ValueError(f"Only {len(valid_filenames)} valid genes found, cannot sample 1000.")

    # Randomly sample 1000 genes
    random.seed(42)  # Seed for reproducibility
    filenames = random.sample(valid_filenames, 1000)
    print(f"Successfully sampled 1000 random genes (seed 42).")

    # Directories to create
    local_align_dir = '/home/miguel/IBE-UPF/PhD/PhyloPhere/test/inputs/alignments/Ali_1000'
    local_map_dir = '/home/miguel/IBE-UPF/PhD/PhyloPhere/test/inputs/map_1000'
    local_asr_dir = '/home/miguel/IBE-UPF/PhD/PhyloPhere/test/inputs/asr'

    os.makedirs(local_align_dir, exist_ok=True)
    os.makedirs(local_map_dir, exist_ok=True)
    os.makedirs(local_asr_dir, exist_ok=True)

    # Write files lists
    align_list = []
    map_list = []
    asr_list = []

    for name in filenames:
        # name = "A1BG.Homo_sapiens.fa" or "ABCF2.Papio_anubis.fa"
        prefix = name.rsplit('.', 1)[0]  # "A1BG.Homo_sapiens"
        gene = prefix.split('.')[0]       # "A1BG"
        align_list.append(name)
        map_list.append(f"{prefix}.map.tsv")
        asr_list.append(f"asr_{gene}")

    with open('/tmp/align_files.txt', 'w') as f:
        f.write('\n'.join(align_list) + '\n')

    with open('/tmp/map_files.txt', 'w') as f:
        f.write('\n'.join(map_list) + '\n')

    with open('/tmp/asr_files.txt', 'w') as f:
        f.write('\n'.join(asr_list) + '\n')

    print("Staged file lists in /tmp. Starting rsync transfers...")

    # Rsync command for alignments
    print("Transferring alignments...")
    subprocess.run([
        'rsync', '-avz', '--files-from=/tmp/align_files.txt',
        'marvin2:/homes/users/mramon/scratch/4.Generate_alignments_from_codons/results_ppga/20260615_200027/PROT/bmge/',
        local_align_dir + '/'
    ], check=True)

    # Rsync command for map files
    print("Transferring map files...")
    subprocess.run([
        'rsync', '-avz', '--files-from=/tmp/map_files.txt',
        'marvin2:/homes/users/mramon/scratch/4.Generate_alignments_from_codons/results_ppga/20260615_200027/MAP/',
        local_map_dir + '/'
    ], check=True)

    # Rsync command for ASR directories (requires recursive -r or -a)
    print("Transferring ASR directories...")
    subprocess.run([
        'rsync', '-avzr', '--files-from=/tmp/asr_files.txt',
        'marvin2:/homes/users/mramon/scratch/2.Primates/1.Primates_data/asr/',
        local_asr_dir + '/'
    ], check=True)

    print("Success! All transfers completed.")

if __name__ == '__main__':
    main()
