#ifndef SEQUENCE_ALIGNMENT_LOCAL_H
#define SEQUENCE_ALIGNMENT_LOCAL_H

#include <stdio.h>
#include <stdlib.h>

#include "alignment.h"

alignment_ops_t affine_gap_align_op_counts(const char *s1, const char *s2);
alignment_ops_t affine_gap_align_op_counts_options(const char *s1, size_t m, const char *s2, size_t n, alignment_options_t options);

#endif