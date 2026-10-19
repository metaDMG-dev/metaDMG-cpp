#include <cmath>
#include <cstdio>
#include <cstring>

#include "../main_pmd.cpp"

static double reference_match(double damage, double total_error, double polymorphism) {
    const double base_error = total_error / 3.0;
    return (1.0 - damage) * (1.0 - base_error) * (1.0 - polymorphism)
         + damage * base_error * (1.0 - polymorphism)
         + base_error * polymorphism * (1.0 - damage);
}

static double reference_score(const char *seq, const char *ref, const unsigned char *qs, int len) {
    double ancient_log_likelihood = 0.0;
    double modern_log_likelihood = 0.0;
    for (int i = 0; i < len; ++i) {
        const double ancient_damage = ref[i] == 'C'
            ? 0.01 + 0.3 * std::pow(0.7, i)
            : 0.01 + 0.3 * std::pow(0.7, len - 1 - i);
        const double modern_damage = 0.001;
        const double ancient_match = reference_match(ancient_damage, std::pow(10.0, -qs[i] / 10.0), 0.001);
        const double modern_match = reference_match(modern_damage, std::pow(10.0, -qs[i] / 10.0), 0.001);
        if (ref[i] != 'C' && ref[i] != 'G')
            continue;
        const bool match = seq[i] == ref[i];
        ancient_log_likelihood += std::log(match ? ancient_match : 1.0 - ancient_match);
        modern_log_likelihood += std::log(match ? modern_match : 1.0 - modern_match);
    }
    return ancient_log_likelihood - modern_log_likelihood;
}

static char complement(char base) {
    return revTable[static_cast<unsigned char>(base)];
}

int main() {
    std::memset(revTable, 'n', sizeof(revTable));
    revTable['A'] = 'T';
    revTable['C'] = 'G';
    revTable['G'] = 'C';
    revTable['T'] = 'A';
    for (int i = 0; i < 256; ++i)
        qsToProp[i] = std::pow(10.0, -i / 10.0);

    const char seq[] = "TCCGAA";
    const char ref[] = "CCCGAA";
    const unsigned char qs[] = {30, 30, 20, 30, 30, 30};
    const int len = 6;

    char clean_ref[] = "CCCGAA";

    bam1_t forward = {};
    forward.core.flag = 0;
    const double expected_forward = reference_score(seq, clean_ref, qs, len);
    const double observed_forward = pmd_stat(&forward, const_cast<char *>(seq), clean_ref, len,
                                             const_cast<char *>(reinterpret_cast<const char *>(qs)));
    if (std::fabs(observed_forward - expected_forward) > 1e-12) {
        std::fprintf(stderr, "forward mismatch: %.17g != %.17g\n", observed_forward, expected_forward);
        return 1;
    }

    char reverse_seq[len + 1];
    char reverse_ref[len + 1];
    unsigned char reverse_qs[len];
    for (int i = 0; i < len; ++i) {
        reverse_seq[i] = complement(seq[len - 1 - i]);
        reverse_ref[i] = complement(clean_ref[len - 1 - i]);
        reverse_qs[i] = qs[len - 1 - i];
    }
    reverse_seq[len] = '\0';
    reverse_ref[len] = '\0';

    bam1_t reverse = {};
    reverse.core.flag = BAM_FREVERSE;
    const double observed_reverse = pmd_stat(&reverse, reverse_seq, reverse_ref, len,
                                             reinterpret_cast<char *>(reverse_qs));
    if (std::fabs(observed_reverse - expected_forward) > 1e-12) {
        std::fprintf(stderr, "reverse mismatch: %.17g != %.17g\n", observed_reverse, expected_forward);
        return 1;
    }

    std::printf("PMD model tests passed: %.12f\n", observed_forward);
    return 0;
}
