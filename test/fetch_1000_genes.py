import os
import subprocess

def main():
    genes_list_path = '/tmp/genes_list.txt'
    with open(genes_list_path, 'r') as f:
        filenames = [line.strip() for line in f if line.strip()]

    print(f"Read {len(filenames)} gene files from {genes_list_path}")

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
