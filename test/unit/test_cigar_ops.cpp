// Unit test: CIGAR operations that do not map 1:1 to 'M' must give the same
// damage profile as the equivalent plain 'M' alignment.
//
// The same 8 bp read (reference CCGTACGT, read TCGTACGT, i.e. a C->T at the
// first 5' position) is encoded with different CIGARs:
//   8M                      reference case
//   1X7=                    sequence match/mismatch ops (minimap2 --eqx)
//   4M<long>D4M             long deletion (used to overflow the reconstruction buffer)
//   4M<long>N4M             long reference skip (spliced alignment)
//   2H8M2H                  hard clips
// Every case must produce exactly the same 5' and 3' mismatch matrices.

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "profile.h"

static const int PRINT_LENGTH = 8;
static const int LONG_GAP = 1000000;  // far larger than the 512-byte initial buffer
static const char *SEQ = "TCGTACGT";
static const char *QUAL = "IIIIIIII";
static const int CT = 1 * 4 + 3;  // toIndex[C][T]

struct Profile {
    std::vector<float> mm5p;
    std::vector<float> mm3p;
};

static bool run_case(const char *name, const std::vector<uint32_t> &cigar, const std::string &md, Profile &out) {
    bam1_t *read = bam_init1();
    if (bam_set1(read, strlen(name), name, 0, 0, 100, 60, cigar.size(), cigar.data(), -1, -1, 0,
                 strlen(SEQ), SEQ, QUAL, 0) < 0 ||
        bam_aux_append(read, "MD", 'Z', md.size() + 1, reinterpret_cast<const uint8_t *>(md.c_str())) < 0) {
        std::fprintf(stderr, "[%s] could not construct test read\n", name);
        bam_destroy1(read);
        return false;
    }

    damage dmg(PRINT_LENGTH, 1, 0);
    dmg.damage_analysis(read, 0, 1.0f);
    const triple &result = dmg.assoc.find(0)->second;

    out.mm5p.clear();
    out.mm3p.clear();
    for (int pos = 0; pos < PRINT_LENGTH; ++pos)
        for (int mm = 0; mm < 16; ++mm) {
            out.mm5p.push_back(result.mm5pF[pos][mm]);
            out.mm3p.push_back(result.mm3pF[pos][mm]);
        }
    bam_destroy1(read);
    return true;
}

static float total(const Profile &p) {
    float sum = 0.0f;
    for (float v : p.mm5p) sum += v;
    for (float v : p.mm3p) sum += v;
    return sum;
}

int main() {
    int failures = 0;

    // Reference case: 8M with the mismatch at the first base.
    Profile ref;
    if (!run_case("plain_8M", {bam_cigar_gen(8, BAM_CMATCH)}, "0C7", ref))
        return 1;
    if (total(ref) != 2.0f * PRINT_LENGTH || ref.mm5p[0 * 16 + CT] != 1.0f) {
        std::fprintf(stderr, "[plain_8M] unexpected reference profile (total %g, 5' C->T at pos 0: %g)\n",
                     total(ref), ref.mm5p[0 * 16 + CT]);
        return 1;
    }

    std::string long_del_md = "0C3^" + std::string(LONG_GAP, 'A') + "4";

    struct Case {
        const char *name;
        std::vector<uint32_t> cigar;
        std::string md;
    } cases[] = {
        {"eqx_1X7=", {bam_cigar_gen(1, BAM_CDIFF), bam_cigar_gen(7, BAM_CEQUAL)}, "0C7"},
        {"long_deletion",
         {bam_cigar_gen(4, BAM_CMATCH), bam_cigar_gen(LONG_GAP, BAM_CDEL), bam_cigar_gen(4, BAM_CMATCH)},
         long_del_md},
        {"long_refskip",
         {bam_cigar_gen(4, BAM_CMATCH), bam_cigar_gen(LONG_GAP, BAM_CREF_SKIP), bam_cigar_gen(4, BAM_CMATCH)},
         "0C7"},
        {"hardclips",
         {bam_cigar_gen(2, BAM_CHARD_CLIP), bam_cigar_gen(8, BAM_CMATCH), bam_cigar_gen(2, BAM_CHARD_CLIP)},
         "0C7"},
    };

    for (const Case &c : cases) {
        Profile got;
        if (!run_case(c.name, c.cigar, c.md, got)) {
            ++failures;
            continue;
        }
        if (got.mm5p != ref.mm5p || got.mm3p != ref.mm3p) {
            std::fprintf(stderr, "[%s] profile differs from plain 8M (total %g, expected %g)\n", c.name,
                         total(got), total(ref));
            ++failures;
        } else {
            std::printf("[%s] ok\n", c.name);
        }
    }

    if (failures) {
        std::fprintf(stderr, "cigar ops test: %d case(s) failed\n", failures);
        return 1;
    }
    std::puts("cigar ops profile test passed");
    return 0;
}
