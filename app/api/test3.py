import math

def micro(seq: str) -> float:

    res = 0
    l = len(seq)
    for i in range(l):
        if seq[i] in ['A', 'T']:
            res += 1
        else:
            res += 2
    return res

def calMicroScore(seq1: str, seq2: str) -> float:

    out_frame_score = 0.0
    frame_score = 0.0
    l = len(seq1)
    sad = set()

    for i in range(l):
        for j in range(2, l - i):
            tmp = seq1[i:i + j]
            for k in range(l - j):
                if tmp == seq2[k:k + j]:

                    ext = j
                    while (i + ext < l and k + ext < l and seq1[i + ext] == seq2[k + ext]):
                        ext += 1
                    tmp = seq2[k:k+ext]

                    key = (i, k, tmp)
                    if key in sad:
                        continue
                    sad.add(key)

                    delta = k - i + l
                    mh_score = math.exp(-delta / 20) * micro(tmp) * 100
                    frame_score += mh_score
                    if (delta) % 3 != 0:
                        out_frame_score += mh_score
                    print(i, k, delta,micro(tmp), tmp, mh_score)
    return out_frame_score / frame_score

print(calMicroScore("GCCACAGCAGCCTCTGA", "AGTTGGACAGCAAAACC"))