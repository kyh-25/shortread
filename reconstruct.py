import os
import time
from collections import defaultdict, Counter
from datetime import datetime

from stringMatch import *
from short_read import *

def reconstruct_sequence(reads, reference, read_length,max_mismatch=1):
    # 인덱싱
    sa = build_suffix_array(reference)
    bwt = build_bwt(reference, sa)
    occ = build_occ_table(bwt)
    c_table = build_c_table(bwt)

    # 각 위치에 복원할 base 후보 저장
    coverage = defaultdict(list)

    # 위치 탐색
    for read, _ in reads:
        matches = bwt_backward_search(bwt, c_table, occ, read, max_mismatch, max_return=1)
        if matches:
            sa_index = matches[0][0]
            ref_pos = sa[sa_index]
            # 삽입
            for i in range(len(read)):
                if ref_pos + i < len(reference):  # 범위 확인
                    coverage[ref_pos + i].append(read[i])

    # 복원 시퀀스 생성
    reconstructed = []
    for i in range(len(reference)):
        bases = coverage.get(i, [])
        if bases:
            most_common = Counter(bases).most_common(1)[0][0]
            reconstructed.append(most_common)
        else:
            reconstructed.append('N')  # 복원 불가능한 경우

    return ''.join(reconstructed)

def compare_sequences(seq1, seq2):

    total = 0
    match = 0
    for a, b in zip(seq1, seq2):
        total += 1
        if a == b:
            match += 1
    accuracy = match / total if total > 0 else 0
    return accuracy, match, total

def load_fasta(filepath):
    with open(filepath) as f:
        return ''.join(line.strip() for line in f if not line.startswith(">"))

def load_reads(filepath):
    reads = []
    with open(filepath) as f:
        read = ""
        for line in f:
            if line.startswith(">"):
                continue  # 헤더 무시
            else:
                read = line.strip()
                reads.append((read, 0))  # start 위치는 쓰지 않음
    return reads


def generate_data_if_not_exists(fasta_file, genome_limit, num_reads, read_length, mutation_rate, overlap, step):
    # 디렉토리 경로 생성
    case_dir = f"data/n{genome_limit}_m{num_reads}_l{read_length}_d{mutation_rate}"

    # 디렉토리가 이미 존재하면 모든 데이터가 생성된 것으로 판단
    if os.path.exists(case_dir):
        reference_seq = load_fasta(os.path.join(case_dir, "reference.fasta"))
        mutated_seq = load_fasta(os.path.join(case_dir, "mutated.fasta"))
        reads = load_reads(os.path.join(case_dir, "short_reads.fasta"))
        return reference_seq, mutated_seq, reads, case_dir

    # 디렉토리 생성
    os.makedirs(case_dir)
    # 데이터 생성
    # 1. 유전체 읽기
    reference_seq = read_limited_genome(fasta_file, genome_limit)
    # 2. 원본 시퀀스 저장
    save_fasta(reference_seq, os.path.join(case_dir, "reference.fasta"), header=">reference")
    # 3. 변이 적용
    mutated_seq, mutations = mutate_sequence(reference_seq, mutation_rate)
    # 4. 변이된 시퀀스 및 변이 정보 저장
    save_fasta(mutated_seq, os.path.join(case_dir, "mutated.fasta"), header=">mutated")
    save_mutations_to_file(mutations, os.path.join(case_dir, "mutations.csv"))
    # 5. 변이된 유전체로부터 short read 생성 및 저장
    reads = generate_short_reads(mutated_seq, read_length, num_reads, overlap=overlap, step=step)
    save_reads_to_fasta(reads, os.path.join(case_dir, "short_reads.fasta"))

    return reference_seq, mutated_seq, reads, case_dir


def run_case(genome_limit, num_reads, read_length, mutation_rate, overlap=True, step=3, max_mismatch=1):
    fasta_file = "data/GCF_000001405.40_GRCh38.p14_genomic.fna"

    #데이터 생성 혹은 로드
    reference, mutated, reads, case_dir = generate_data_if_not_exists(fasta_file, genome_limit, num_reads, read_length, mutation_rate, overlap, step)

    #복원
    start_time = time.time()
    reconstructed = reconstruct_sequence(reads, reference, read_length,max_mismatch)
    end_time = time.time()
    elapsed = end_time - start_time

    #결과 정확도
    accuracy, matched, total = compare_sequences(reconstructed, mutated)
    print(f"복원 정확도: {accuracy * 100:.2f}%")
    print(f"→ 소요 시간: {elapsed:.2f}초")

    #확인용
    save_fasta(reconstructed, os.path.join(case_dir, f"reconstructed_d{max_mismatch}.fasta"), header=">reconstructed")

    return {
        "genome_limit": genome_limit,
        "num_reads": num_reads,
        "read_length": read_length,
        "mutation_rate": mutation_rate,
        "max_mismatch":max_mismatch,
        "accuracy": accuracy,
        "time_sec": elapsed
    }

# 실행 메인
def main():
    # 설정값
    fasta_file = "data/GCF_000001405.40_GRCh38.p14_genomic.fna"

    cases = [
        # n      m     l    변이률 overlap 간격  d
        # L 변화
        (100000, 10000, 50, 0.03, True, 10, 2),
        (100000, 10000, 75, 0.03, True, 10, 2),
        (100000, 10000, 100, 0.03, True, 10, 2),
        (100000, 10000, 125, 0.03, True, 10, 2),
        (100000, 10000, 150, 0.03, True, 10, 2),
        # d 변화
        (100000, 10000, 50, 0.03, True, 10, 1),
        (100000, 10000, 50, 0.03, True, 10, 2),
        (100000, 10000, 50, 0.03, True, 10, 3),
        (100000, 10000, 50, 0.03, True, 10, 4),
        # m 변화
        (100000, 10000, 50, 0.03, True, 10, 2),
        (100000, 20000, 50, 0.03, True, 5, 2),
        (100000, 50000, 50, 0.03, True, 2, 2),
        #n 변화
        (100000, 10000, 50, 0.03, True, 10, 2),
        (300000,  30000, 50, 0.03,  True,   10,   2),
        (400000,  40000, 50, 0.03,   True,   10,   2),
        (700000,  70000, 50, 0.03,   True,   10,   2),
        (1000000, 100000, 50, 0.03, True, 10, 2),
        (3000000, 300000, 50, 0.03, True, 10, 2),
    ]

    results = []

    for case in cases:
        result = run_case(*case)
        results.append(result)

    # CSV로 저장
    now_str = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    with open(f"result/results_{now_str}", "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)



if __name__ == '__main__':
    main()
