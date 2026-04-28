#ifndef METADMG_DUST_SCORE_H
#define METADMG_DUST_SCORE_H

#include <stdint.h>  // for int32_t, uint8_t

double dust_score_nt16(const uint8_t *seq, int32_t l, int32_t window);

#endif
