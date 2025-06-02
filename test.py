from collections import defaultdict

def build_bwt(text,sa):
    # text += "$"
    # rotations = sorted([text[i:] + text[:i] for i in range(len(text))])
    # bwt = ''.join(row[-1] for row in rotations)
    # return bwt
    return ''.join(text[i - 1] if i > 0 else '$' for i in sa)

def build_suffix_array(text):
    """ 접미사 배열 생성 """
    text += "$"
    return sorted(range(len(text)), key=lambda i: text[i:])

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

# def bwt_backward_search(bwt, c_table, occ_table, pattern, max_mismatch=1):
#     matches = []
#     def search(i, l, r, mismatches, path):
#         if mismatches > max_mismatch:
#             return
#         if i < 0:
#             for j in range(l, r):
#                 matches.append((j, mismatches, path[::-1]))
#             return
#         for a in c_table:
#             l2 = c_table[a] + get_occ(occ_table, a, l - 1)
#             r2 = c_table[a] + get_occ(occ_table, a, r - 1)
#             if l2 < r2:
#                 mismatch_penalty = 0 if a == pattern[i] else 1
#                 search(i - 1, l2, r2, mismatches + mismatch_penalty, path + a)
#
#     search(len(pattern) - 1, 0, len(bwt), 0, "")
#     return matches

def bwt_backward_search(bwt, c_table, occ_table, pattern, max_mismatch=1, max_return=1):
    matches = []
    pattern_len = len(pattern)

    mismatch_priority = list(range(pattern_len - 1, -1, -1))

    visited = set()

    def search(i, l, r, mismatches, path, mismatch_positions):
        if len(matches) >= max_return:
            return
        if mismatches > max_mismatch:
            return
        if i < 0:
            for j in range(l, r):
                matches.append((j, mismatches, path[::-1], mismatch_positions))
            return

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

def main():
    # 예시 참조 유전체
    reference = "ACGTCGACGTCG"
    suffix_array = build_suffix_array(reference)
    bwt = build_bwt(reference,suffix_array)
    occ_table = build_occ_table(bwt)
    c_table = build_c_table(bwt)

    # Short read
    read = "CGAC"

    # # Soft alignment 수행
    # matches = bwt_backward_search(bwt, c_table, occ_table, read, max_mismatch=1)
    #
    # # 결과 출력
    # print("Read:", read)
    # for idx, mismatches, aligned in matches:
    #     pos = suffix_array[idx]
    #     print(f"Match at pos {pos} (mismatches: {mismatches}) → aligned: {aligned}")

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
