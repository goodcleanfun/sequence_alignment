#ifndef SEQUENCE_ALIGNMENT_H
#define SEQUENCE_ALIGNMENT_H

#include "vector/vector.h"

typedef struct {
    size_t insert_gap_open_cost;
    size_t insert_gap_extend_cost;
    size_t delete_gap_open_cost;
    size_t delete_gap_extend_cost;
    size_t match_cost;
    size_t mismatch_cost;
    size_t transpose_cost;
    bool prefix_gap;
    bool suffix_gap;
    bool ignore_non_alphanumeric;
} alignment_options_t;


static const alignment_options_t DEFAULT_ALIGNMENT_OPTIONS_AFFINE_GAP = {
    .insert_gap_open_cost = 3,
    .insert_gap_extend_cost = 2,
    .delete_gap_open_cost = 3,
    .delete_gap_extend_cost = 2,
    .match_cost = 0,
    .mismatch_cost = 6,  // Should be greater than or equal to (gap_open_cost + gap_extend_cost)
    .transpose_cost = 4, // Should be less than 2x mismatch_cost, ideally equal, should be less than (gap_open_cost + 2 x gap_extend_cost)
    .prefix_gap = false,
    .suffix_gap = false,
    .ignore_non_alphanumeric = true
};

typedef enum {
    ALIGN_NO_OP = 0,
    ALIGN_MATCH,
    ALIGN_MISMATCH,
    ALIGN_TRANSPOSE,
    ALIGN_DELETE_GAP_OPEN,
    ALIGN_DELETE_GAP_EXTEND,
    ALIGN_INSERT_GAP_OPEN,
    ALIGN_INSERT_GAP_EXTEND
} alignment_op;

typedef struct {
    size_t num_matches; // number of matched characters in the alignment
    size_t num_mismatches; // number of mismatched characters in the alignment
    size_t num_transpositions; // number of transposed characters
    size_t num_delete_gap_opens; // total number of gaps from deletes
    size_t num_delete_gap_extensions; // total number of characters across all delete gaps
    size_t num_insert_gap_opens; // total number of gaps from inserts
    size_t num_insert_gap_extensions; // total number of characters across all insert gaps
} alignment_ops_t;


#endif
