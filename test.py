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

def main():
    # 예시 참조 유전체
    reference = "banana"
    suffix_array = build_suffix_array(reference)
    bwt = build_bwt(reference,suffix_array)
    
    for sa in suffix_array:
        print(sa)
    for b in bwt:
        print(b)

 
if __name__ == '__main__':
    main()