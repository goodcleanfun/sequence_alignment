#ifndef SEQUENCE_ALIGNMENT_H
#define SEQUENCE_ALIGNMENT_H

#include "vector/vector.h"

typedef struct {
    size_t gap_open_cost;
    size_t gap_extend_cost;
    size_t match_cost;
    size_t mismatch_cost;
    size_t transpose_cost;
    bool prefix_gap;
    bool suffix_gap;
    bool mismatches_alphanumeric_only;
} alignment_options_t;


static const alignment_options_t DEFAULT_ALIGNMENT_OPTIONS_AFFINE_GAP = {
    .gap_open_cost = 3,
    .gap_extend_cost = 2,
    .match_cost = 0,
    .mismatch_cost = 6,  // Should be greater than (gap_open_cost + gap_extend_cost)
    .transpose_cost = 4, // Should be less than 2x mismatch_cost
    .prefix_gap = false,
    .suffix_gap = false
};

typedef enum {
    ALIGN_MATCH,
    ALIGN_MISMATCH,
    ALIGN_TRANSPOSE,
    ALIGN_GAP_OPEN,
    ALIGN_GAP_EXTEND
} alignment_op;

typedef struct {
    size_t num_matches; // number of matched characters in the alignment
    size_t num_mismatches; // number of mismatched characters in the alignment
    size_t num_transpositions; // number of transposed characters
    size_t num_gap_opens; // total number of gaps
    size_t num_gap_extensions; // total number of characters across all gaps (insertions/deletions)
} alignment_ops_t;


#endif
