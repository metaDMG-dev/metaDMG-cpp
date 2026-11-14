#include <cstdio>
#include <cstring>

#include <htslib/sam.h>

#include "../profile.h"

int main() {
    bam1_t *read = bam_init1();
    uint32_t cigar[] = {
        bam_cigar_gen(2, BAM_CSOFT_CLIP),
        bam_cigar_gen(4, BAM_CMATCH),
        bam_cigar_gen(2, BAM_CSOFT_CLIP),
    };
    const char qualities[] = "IIIIIIII";
    if (bam_set1(read, 3, "r1", 0, 0, 0, 60, 3, cigar, -1, -1, 0, 8,
                 "AACCCCTT", qualities, 0) < 0 ||
        bam_aux_append(read, "MD", 'Z', 2, reinterpret_cast<const uint8_t *>("4")) < 0) {
        std::fprintf(stderr, "Could not construct soft-clipped test read\n");
        bam_destroy1(read);
        return 1;
    }

    damage dmg(8, 1, 0);
    dmg.damage_analysis(read, 0, 1.0f);
    const triple &result = dmg.assoc.find(0)->second;

    float total = 0.0f;
    for (int position = 0; position < 8; ++position)
        for (int mismatch = 0; mismatch < 16; ++mismatch)
            total += result.mm5pF[position][mismatch] + result.mm3pF[position][mismatch];

    bam_destroy1(read);
    if (total != 8.0f) {
        std::fprintf(stderr, "Soft-clipped read contributed %g counts; expected 8\n", total);
        return 1;
    }
    std::puts("softclip profile test passed");
    return 0;
}
