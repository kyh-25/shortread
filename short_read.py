from Bio import SeqIO
import csv
import random


def mutate_sequence(seq, mutation_rate=0.03):
    """
    주어진 염기서열에서 mutation_rate 비율로 무작위 변이를 일으킴
    :return: (변이된 시퀀스, 변이 정보 리스트)
    """
    bases = ['A', 'C', 'G', 'T']
    seq_list = list(seq)
    num_mutations = int(len(seq) * mutation_rate)
    mutation_indices = random.sample(range(len(seq)), num_mutations)
    mutations = []

    for idx in mutation_indices:
        original = seq_list[idx]
        candidates = [b for b in bases if b != original]
        mutated = random.choice(candidates)
        seq_list[idx] = mutated
        mutations.append((idx, original, mutated))

    return ''.join(seq_list), mutations


def save_mutations_to_file(mutations, output_path):
    """ 변이 정보 CSV 저장 """
    with open(output_path, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Position", "Original", "Mutated"])
        for idx, orig, mut in mutations:
            writer.writerow([idx, orig, mut])


def save_fasta(sequence, output_path, header=">sequence"):
    """ FASTA 형식으로 염기서열 저장 """
    with open(output_path, "w") as f:
        f.write(f"{header}\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")


def read_limited_genome(fasta_path, genome_limit):
    """ FASTA 파일에서 genome_limit 만큼 유전체 읽기 """
    genome = []
    current_len = 0
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq = line.strip().upper().replace("N", "")
            needed = genome_limit - current_len
            if needed <= 0:
                break
            if len(seq) <= needed:
                genome.append(seq)
                current_len += len(seq)
            else:
                genome.append(seq[:needed])
                current_len += needed
                break
    return ''.join(genome)


def generate_short_reads(genome, read_length=100, num_reads=1000, overlap=False, step=1):
    """
    short read 생성 (입력은 시퀀스 문자열)
    """
    genome_length = len(genome)
    reads = []

    if overlap:
        for i in range(0, genome_length - read_length + step, step):
            if i + read_length < genome_length:
                read = genome[i:i + read_length]
            else:
                read = genome[i:genome_length]
            reads.append((read, i))
            if len(reads) >= num_reads:
                break
    else:
        for _ in range(num_reads):
            start = random.randint(0, genome_length - read_length)
            if start + read_length < genome_length:
                read = genome[start:start + read_length]
            else:
                read = genome[start:genome_length]
            reads.append((read, start))

    return reads


def save_reads_to_fasta(reads, output_path="short_reads.fasta"):
    """ short read 리스트를 FASTA로 저장 """
    with open(output_path, "w") as f:
        for i, (read, start_pos) in enumerate(reads):
            end_pos = start_pos + len(read)
            f.write(f">read_{i + 1}_{start_pos}_{end_pos}\n{read}\n")


def main():
    # 설정값
    fasta_file = "data/GCF_000001405.40_GRCh38.p14_genomic.fna"
    genome_limit = 300000
    num_reads = 100000
    read_length = 50
    mutation_rate = 0.1
    overlap = True
    step = 3  # overlap 사용 시 window 간격

    # 1. 유전체 읽기
    reference_seq = read_limited_genome(fasta_file, genome_limit)

    # 2. 원본 시퀀스 저장
    save_fasta(reference_seq, output_path=f"data/n{genome_limit}_m{num_reads}_l{read_length}_reference.fasta", header=">reference")

    # 3. 변이 적용
    mutated_seq, mutations = mutate_sequence(reference_seq, mutation_rate)

    # 4. 변이된 시퀀스 및 변이 정보 저장
    save_fasta(mutated_seq, output_path=f"data/n{genome_limit}_m{num_reads}_l{read_length}_mutated.fasta", header=">mutated")
    save_mutations_to_file(mutations,  output_path=f"data/n{genome_limit}_m{num_reads}_l{read_length}_mutations.csv")

    # 5. 변이된 유전체로부터 short read 생성 및 저장
    reads = generate_short_reads(mutated_seq, read_length, num_reads, overlap=overlap, step=step)
    save_reads_to_fasta(reads, output_path=f"data/n{genome_limit}_m{num_reads}_l{read_length}_short_reads.fasta")


if __name__ == '__main__':
    main()
