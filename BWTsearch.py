from collections import defaultdict,Counter

#멘버-마이어 알고리즘
def build_suffix_array(text):
    text += '$'  # 종료 문자 추가 (모든 접미사를 유일하게 만듦)
    n = len(text)
    k = 1  # 초기 접두사 길이 (1, 2, 4, 8, ...)
    
    # 초기 랭크: 문자 아스키 코드
    rank = [ord(c) for c in text]
    tmp = [0] * n  # 새 랭크를 임시 저장
    suffix_array = list(range(n))  # 초기 인덱스 배열

    while True:
        # 랭크 쌍을 기반으로 접미사 정렬
        #기본 sort = timsort
        suffix_array.sort(key=lambda i: (rank[i], rank[i + k] if i + k < n else -1))
        
        # 새 랭크 계산
        tmp[suffix_array[0]] = 0
        for i in range(1, n):
            prev = suffix_array[i - 1]
            curr = suffix_array[i]
            prev_rank = (rank[prev], rank[prev + k] if prev + k < n else -1)
            curr_rank = (rank[curr], rank[curr + k] if curr + k < n else -1)
            tmp[curr] = tmp[prev] + (1 if curr_rank != prev_rank else 0)
        
        rank = tmp[:]

        if max(rank) == n - 1:
            break  # 모든 접미사에 고유한 랭크 부여 완료

        k *= 2  # 접두사 길이 두 배 증가

    return suffix_array


def build_bwt(text,sa):
    return ''.join(text[i - 1] if i > 0 else '$' for i in sa)

def build_occ_table(bwt):
    """ Occurrence 테이블: 특정 위치까지 각 문자가 몇 번 나왔는지를 저장 """
    occ = defaultdict(list)
    count = defaultdict(int)
    for i, char in enumerate(bwt):
        count[char] += 1    #누적
        for c in set(bwt):
            occ[c].append(count[c]) #테이블 채워넣기기
    return occ

def build_c_table(bwt):
    """ C 테이블: 사전순으로 앞에 나오는 문자들의 개수 누적합 """
    #각 문자별 카운트
    counts = defaultdict(int)
    for c in bwt:
        counts[c] += 1
    #사전순으로 누적합
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
