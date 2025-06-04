from collections import defaultdict, Counter

def hamming_distance(s1, s2):
    """두 문자열의 불일치도 계산"""
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def reconstruct_sequence_brute_force(reads, reference, max_mismatch=3):
    coverage = defaultdict(list)
    ref_len = len(reference)

    for read, _ in reads:
        read_len = len(read)
        # reference 전체를 슬라이딩하며 read와 비교
        for pos in range(ref_len - read_len + 1):
            segment = reference[pos:pos+read_len]
            dist = hamming_distance(read, segment)
            if dist <= max_mismatch:
                # mismatch 이하인 경우에만 커버리지에 기록
                for i in range(read_len):
                    coverage[pos + i].append(read[i])
                # 여러 매치 중 첫 번째만 사용하려면 여기서 break 해도 됨
                # break

    # 복원 시퀀스 생성
    reconstructed = []
    for i in range(ref_len):
        bases = coverage.get(i, [])
        if bases:
            most_common = Counter(bases).most_common(1)[0][0]
            reconstructed.append(most_common)
        else:
            reconstructed.append('N')

    return ''.join(reconstructed)
