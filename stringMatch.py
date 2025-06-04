import functools
from collections import defaultdict,Counter

def build_suffix_array(text):
    """ 접미사 배열 생성 """
    text += '$'
    def compare(i, j):
        while i < len(text) and j < len(text):
            if text[i] != text[j]:
                return -1 if text[i] < text[j] else 1
            i += 1
            j += 1
        return -1 if i == len(text) and j != len(text) else (1 if j == len(text) and i != len(text) else 0)

    return sorted(range(len(text)), key=functools.cmp_to_key(compare))



def build_bwt(text,sa):
    return ''.join(text[i - 1] if i > 0 else '$' for i in sa)

def build_occ_table(bwt):
    """ Occurrence 테이블: 특정 위치까지 각 문자가 몇 번 나왔는지를 저장 """
    occ = defaultdict(list)
    count = defaultdict(int)
    for i, char in enumerate(bwt):
        count[char] += 1
        for c in set(bwt):
            occ[c].append(count[c])
    return occ

def build_c_table(bwt):
    """ C 테이블: 사전순으로 앞에 나오는 문자들의 개수 누적합 """
    counts = defaultdict(int)
    for c in bwt:
        counts[c] += 1
    c_table = {}
    total = 0
    for c in sorted(counts.keys()):
        c_table[c] = total
        total += counts[c]
    return c_table

def get_occ(occ_table, c, i):
    if i < 0:
        return 0
    return occ_table[c][i]


def bwt_backward_search(bwt, c_table, occ_table, pattern, max_mismatch=1, max_return=3):
    matches = []
    pattern_len = len(pattern)

    mismatch_priority = list(range(pattern_len - 1, -1, -1))

    visited = set()

    def search(i, l, r, mismatches, path, mismatch_positions):
        #너무 많은 매치
        if len(matches) >= max_return:
            return
        #오차 너무 많음
        if mismatches > max_mismatch:
            return
        #패턴 다 탐색
        if i < 0:
            for j in range(l, r):
                matches.append((j, mismatches, path[::-1], mismatch_positions))
            return
        # 중복 방지
        key = (i, l, r, mismatches)
        if key in visited:
            return
        visited.add(key)

        current_char = pattern[i]

        for a in c_table:
            l2 = c_table[a] + get_occ(occ_table, a, l - 1)
            r2 = c_table[a] + get_occ(occ_table, a, r - 1)
            if l2 < r2:
                if a == current_char:
                    search(i - 1, l2, r2, mismatches, path + a, mismatch_positions)
                elif (pattern_len - i - 1) in mismatch_priority:
                    search(i - 1, l2, r2, mismatches + 1, path + a, mismatch_positions + [pattern_len - i - 1])

    search(len(pattern) - 1, 0, len(bwt), 0, "", [])
    return matches

def reconstruct_sequence_BWT(reads, reference,max_mismatch=1):
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


def main():
    # 예시 참조 유전체
    reference = "ACGTCGACGTCG"
    suffix_array = build_suffix_array(reference)
    bwt = build_bwt(reference,suffix_array)
    occ_table = build_occ_table(bwt)
    c_table = build_c_table(bwt)

    # Short read
    read = "CGAC"

    matches = bwt_backward_search(
        bwt, c_table, occ_table,
        read,
        max_mismatch=1,
        max_return=1  # Bowtie는 하나만 찾고 끝냄
    )
    # 결과 출력
    print("Read:", read)
    for match in matches:
        sa_index, mismatches, aligned, mismatch_pos = match
        print(f"[hit] SA index: {sa_index}, mismatches: {mismatches}, aligned: {aligned}, mismatch pos: {mismatch_pos}")

if __name__ == '__main__':
    main()
