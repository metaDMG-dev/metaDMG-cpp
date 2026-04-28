#include "dust_score.h"

#include <htslib/sam.h>  // for bam_seqi
#include <string.h>      // for memset

static int decode_nt16_acgt(uint8_t nt16) {
    if (nt16 == 1) return 0; // A
    if (nt16 == 2) return 1; // C
    if (nt16 == 4) return 2; // G
    if (nt16 == 8) return 3; // T
    return -1;
}

double dust_score_nt16(const uint8_t *seq, int32_t l, int32_t window) {
    const int32_t WLEN = 3;
    const int32_t WTOT = 64;
    const int32_t WMASK = WTOT - 1;
    if (window < WLEN || window > 64)
        return -1.0;

    int32_t wCount[WTOT];
    int32_t wSeq[64];
    memset(wCount, 0, sizeof(wCount));
    memset(wSeq, 0, sizeof(wSeq));

    int64_t score = 0;
    int64_t maxScore = 0;
    int32_t t = 0;
    int32_t n = -WLEN;

    for (int32_t i = 0; i < l; ++i) {
        const int b = decode_nt16_acgt(bam_seqi(seq, i));
        if (b < 0)
            continue; // ignore N/ambiguous

        t = ((t << 2) | b) & WMASK;
        if (++n >= 0) {
            const int32_t k = n % window;
            if (n >= window) {
                const int32_t x = wSeq[k];
                if (wCount[x] > 0)
                    score -= --wCount[x];
                score += wCount[t]++;
                if (score > maxScore)
                    maxScore = score;
            } else {
                score += wCount[t]++;
            }
            wSeq[k] = t;
        }
    }

    if (n <= 0)
        return 0.0;
    if (n >= window)
        return (200.0 * (double)maxScore) / (window * (window - 1));
    return (200.0 * (double)score) / (n * (n + 1));
}
