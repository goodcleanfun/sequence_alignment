#ifndef SEQUENCE_ALIGNMENT_LOCAL_H
#define SEQUENCE_ALIGNMENT_LOCAL_H

#include <stdio.h>
#include <stdlib.h>

#include "alignment.h"
#include "num_arrays/uint32_array.h"

alignment_ops_t affine_gap_align_op_counts(const char *s1, const char *s2);
alignment_ops_t affine_gap_align_op_counts_options(const char *s1, const char *s2, alignment_options_t options);
alignment_ops_t affine_gap_align_op_counts_unicode(uint32_array *u1_array, uint32_array *u2_array);
alignment_ops_t affine_gap_align_op_counts_unicode_options(uint32_array *u1_array, uint32_array *u2_array, alignment_options_t options);

#endif