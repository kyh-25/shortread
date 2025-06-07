# -*- coding: utf-8 -*-
import os
import time
from datetime import datetime

from BWTsearch import *
from short_read import *
from bruteForce import *

#원본 데이터 로드
def load_fasta(filepath):
    with open(filepath) as f:
        return ''.join(line.strip() for line in f if not line.startswith(">"))

#short read 데이터 로드
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

# ---------------------------------------------------------------------------------

#정확도 분석
def compare_sequences(seq1, seq2):
    total = 0
    match = 0
    for a, b in zip(seq1, seq2):
        total += 1
        if a == b:
            match += 1
    accuracy = match / total if total > 0 else 0
    return accuracy, match, total

# ---------------------------------------------------------------------------------

# 원본 데이터로 short_read 생성
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

# ---------------------------------------------------------------------------------

#복원
def run_case(genome_limit, num_reads, read_length, mutation_rate, overlap=True, step=3, max_mismatch=3):
    fasta_file = "data/GCF_000001405.40_GRCh38.p14_genomic.fna"

    #데이터 생성 혹은 로드
    reference, mutated, reads, case_dir = generate_data_if_not_exists(fasta_file, genome_limit, num_reads, read_length, mutation_rate, overlap, step)

    #BWT복원
    start_time = time.time()
    bwt_reconstructed = reconstruct_sequence_BWT(reads, reference,max_mismatch)
    end_time = time.time()
    #시간 계산
    bwt_elapsed = end_time - start_time

    #결과 정확도
    bwt_accuracy, matched, total = compare_sequences(bwt_reconstructed, mutated)
    print(f"복원 정확도: {bwt_accuracy * 100:.2f}%")
    print(f"→ 소요 시간: {bwt_elapsed:.2f}초")
    #확인용용
    save_fasta(bwt_reconstructed, os.path.join(case_dir, f"BWTreconstructed_d{max_mismatch}.fasta"), header=f">reconstructed_{bwt_accuracy}_{bwt_elapsed}")
   
    return {
        "genome_limit": genome_limit,
        "num_reads": num_reads,
        "read_length": read_length,
        "mutation_rate": mutation_rate,
        "max_mismatch":max_mismatch,
        "BWT_accuracy": bwt_accuracy,
        "BWT_time_sec": bwt_elapsed,
    }

# ---------------------------------------------------------------------------------

# 실행 메인
def main():
    # 설정값
    fasta_file = "data/GCF_000001405.40_GRCh38.p14_genomic.fna"

    mutate = 0.03
    max_mismatch = 3
    overlap = False

    #n, m, l, d 변수에 따른 케이스 설정
    cases = [
        # n      m     l  변이률 overlap(랜덤으로 할건지) 간격(랜덤 아닐시 간격을 기준으로 추출)  d
        #L 변화
        (100000, 20000, 10, mutate, overlap, 10, max_mismatch),
        (100000, 20000, 25, mutate, overlap, 10, max_mismatch),
        (100000, 20000, 50, mutate, overlap, 10, max_mismatch),
        (100000, 20000, 75, mutate, overlap, 10, max_mismatch),
        (100000, 20000, 100, mutate, overlap, 10, max_mismatch),
        (100000, 20000, 125, mutate, overlap, 10, max_mismatch),
        (100000, 20000, 150, mutate, overlap, 10, max_mismatch),
        (100000, 20000, 150, mutate, overlap, 10, 4),
        (100000, 20000, 150, mutate, overlap, 10, 5),
        #d 변화
        (100000, 20000, 50, mutate, overlap, 10, 1),
        (100000, 20000, 50, mutate, overlap, 10, 2),
        (100000, 20000, 50, mutate, overlap, 10, 3),
        (100000, 20000, 50, mutate, overlap, 10, 4),
        # m 변화
        (100000, 20000, 50, mutate, overlap, 10, max_mismatch),
        (100000, 40000, 50, mutate, overlap, 10, max_mismatch),
        (100000, 60000, 50, mutate, overlap, 10, max_mismatch),
        #n 변화
        (10000,  2000, 50, mutate, overlap, 10, max_mismatch),
        (30000,  6000, 50, mutate, overlap, 10, max_mismatch),
        (50000,  10000, 50, mutate, overlap, 10, max_mismatch),
        (70000,  14000, 50, mutate, overlap, 10, max_mismatch),
        (90000,  18000, 50, mutate, overlap, 10, max_mismatch),
        (120000,  24000, 50, mutate, overlap, 10, max_mismatch),
        (300000,  60000, 50, mutate,  overlap, 10, max_mismatch),
        (1000000, 200000, 50, mutate, overlap, 10, max_mismatch),
    ]

    results = []

    for case in cases:
        result = run_case(*case)
        results.append(result)

    # CSV로 저장
    now_str = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    with open(f"result/bwt/results_{now_str}", "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)



if __name__ == '__main__':
    main()
