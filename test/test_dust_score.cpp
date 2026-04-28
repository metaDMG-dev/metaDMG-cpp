#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "../dust_score.h"

static uint8_t base_to_nt16(char c) {
    switch (c) {
    case 'A': return 1;
    case 'C': return 2;
    case 'G': return 4;
    case 'T': return 8;
    default: return 15; // N/ambiguous
    }
}

static std::vector<uint8_t> pack_nt16(const std::string &s) {
    std::vector<uint8_t> out((s.size() + 1) / 2, 0);
    for (size_t i = 0; i < s.size(); ++i) {
        const uint8_t v = base_to_nt16(s[i]) & 0x0f;
        if ((i & 1) == 0)
            out[i / 2] |= (uint8_t)(v << 4); // even index -> high nibble
        else
            out[i / 2] |= v; // odd index -> low nibble
    }
    return out;
}

int main() {
    // Invalid window must return -1.
    const std::string tiny = "ACGTACGT";
    std::vector<uint8_t> tiny_packed = pack_nt16(tiny);
    const double invalid = dust_score_nt16(tiny_packed.data(), (int32_t)tiny.size(), 2);
    if (invalid != -1.0) {
        fprintf(stderr, "dust unit fail: expected -1.0 for invalid window, got %.6f\n", invalid);
        return 1;
    }

    // Low-complexity sequence should have higher dust score than high-complexity sequence.
    const std::string low(120, 'A');
    const std::string high = std::string("ACGTACGTACGTACGTACGTACGTACGTACGT"
                                         "TGCATGCATGCATGCATGCATGCATGCATGCA"
                                         "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT");
    std::vector<uint8_t> low_packed = pack_nt16(low);
    std::vector<uint8_t> high_packed = pack_nt16(high);

    const double low_score = dust_score_nt16(low_packed.data(), (int32_t)low.size(), 64);
    const double high_score = dust_score_nt16(high_packed.data(), (int32_t)high.size(), 64);

    if (!(low_score > high_score)) {
        fprintf(stderr, "dust unit fail: expected low complexity score > high complexity score (%.6f <= %.6f)\n", low_score, high_score);
        return 1;
    }
    if (!(low_score - high_score >= 5.0)) {
        fprintf(stderr, "dust unit fail: expected meaningful separation, got low-high=%.6f\n", low_score - high_score);
        return 1;
    }

    // Ns should be ignored; removing Ns should preserve score for same informative content.
    const std::string with_ns = "AAAAANNNNNAAAAANNNNNAAAAA";
    const std::string no_ns = "AAAAAAAAAAAAAAA";
    std::vector<uint8_t> with_ns_packed = pack_nt16(with_ns);
    std::vector<uint8_t> no_ns_packed = pack_nt16(no_ns);
    const double s_with_ns = dust_score_nt16(with_ns_packed.data(), (int32_t)with_ns.size(), 64);
    const double s_no_ns = dust_score_nt16(no_ns_packed.data(), (int32_t)no_ns.size(), 64);
    if (std::fabs(s_with_ns - s_no_ns) > 1e-9) {
        fprintf(stderr, "dust unit fail: expected Ns ignored (%.12f vs %.12f)\n", s_with_ns, s_no_ns);
        return 1;
    }

    return 0;
}
