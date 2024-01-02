#include "local_alignment.h"
#include "utf8/utf8.h"
#include "utf8proc/utf8proc.h"

static alignment_ops_t NULL_ALIGNMENT_OPS = {
    .num_matches = 0,
    .num_mismatches = 0,
    .num_transpositions = 0,
    .num_delete_gap_opens = 0,
    .num_delete_gap_extensions = 0,
    .num_insert_gap_opens = 0,
    .num_insert_gap_extensions = 0
};

static inline bool utf8_is_non_alphanumeric(int32_t c) {
    int cat = utf8proc_category(c);
    return utf8_is_whitespace(c) || utf8_is_hyphen(c) || utf8_is_punctuation(cat);
}

typedef struct {
    size_t *costs;
    size_t *prev_costs;
    size_t *gap_costs;
    alignment_ops_t *ops;
    alignment_ops_t *prev_ops;
    alignment_ops_t *gap_ops;
} alignment_costs_t;

static bool affine_gap_costs_options(const char *s1, size_t m, const char *s2, size_t n, bool reverse, alignment_options_t options, alignment_ops_t *result) {
    bool swapped = false;
    // m should be the larger string, n should be the smaller
    if (m < n) {
        const char *tmp_str = s1;
        s1 = s2;
        s2 = tmp_str;
        size_t tmp_len = m;
        m = n;
        n = tmp_len;
        swapped = true;
    }

    if (m == n && utf8_compare_len(s1, s2, n) == 0) {
        alignment_ops_t exact_match_ops = NULL_ALIGNMENT_OPS;
        exact_match_ops.num_matches = n;
        *result = exact_match_ops;
        return true;
    }

    // Costs are a 3 x (m + 1) matrix, only uses 3 rows
    size_t *all_costs = malloc((m + 1) * 3 * sizeof(size_t));
    if (all_costs == NULL) {
        return false;
    }

    // Costs for the current row, CC in the Myers-Miller paper
    size_t *costs = all_costs;
    // Costs for previous row, used for transpositions as in Damerau-Levenshtein
    // This will track the cost at the i-1, j-2 position in order to revert back
    // to the cost at that position if the transpose is cheaper
    size_t *prev_costs = all_costs + m + 1;
    // Gap costs for inserts/deletes, DD in the Myers-Miller paper
    size_t *gap_costs = all_costs + 2 * (m + 1);

    // Edits are a 3 x (m + 1) matrix of alignment op counts, a breakdown of the cost by type
    alignment_ops_t *all_ops = malloc(3 * (m + 1) * sizeof(alignment_ops_t));
    if (all_ops == NULL) {
        free(all_costs);
        return false;
    }

    // Edits for the current row
    alignment_ops_t *ops = all_ops;
    // Edits for the previous row, used for transpositions in order to
    // revert back to the ops at the i - 1, j -2 at the end of the transpose
    // if the transpose is cheaper
    alignment_ops_t *prev_ops = all_ops + m + 1;
    // Edits for the current row, only for insertions/deletions
    alignment_ops_t *gap_ops = all_ops + 2 * (m + 1);

    size_t insert_gap_open_cost = options.insert_gap_open_cost;
    size_t insert_gap_extend_cost = options.insert_gap_extend_cost;
    size_t delete_gap_open_cost = options.delete_gap_open_cost;
    size_t delete_gap_extend_cost = options.delete_gap_extend_cost;
    size_t match_cost = options.match_cost;
    size_t mismatch_cost = options.mismatch_cost;
    size_t transpose_cost = options.transpose_cost;
    bool ignore_non_alphanumeric = options.ignore_non_alphanumeric;

    size_t e = 0, c = 0, s = 0;

    size_t max_cost = m + n;

    costs[0] = 0;
    gap_costs[0] = max_cost;
    ops[0] = NULL_ALIGNMENT_OPS;
    size_t t = delete_gap_open_cost;

    alignment_ops_t base_ops = NULL_ALIGNMENT_OPS;
    base_ops.num_delete_gap_opens++;

    for (size_t j = 1; j < m + 1; j++) {
        t += delete_gap_extend_cost;
        costs[j] = t;
        gap_costs[j] = t + delete_gap_open_cost;
        base_ops.num_delete_gap_extensions++;
        ops[j] = base_ops;
        gap_ops[j] = base_ops;
    }

    t = insert_gap_open_cost;
    base_ops = NULL_ALIGNMENT_OPS;
    base_ops.num_delete_gap_opens++;

    alignment_ops_t current_ops = NULL_ALIGNMENT_OPS;
    alignment_ops_t prev_char_ops = NULL_ALIGNMENT_OPS;
    alignment_ops_t prev_row_prev_char_ops = NULL_ALIGNMENT_OPS;
    alignment_ops_t prev_row_prev2_char_ops = NULL_ALIGNMENT_OPS;

    bool in_gap = false;
    size_t prev_c = 0;

    uint8_t *s2_ptr;
    if (!reverse) {
        s2_ptr = (uint8_t *)s2;
    } else {
        s2_ptr = (uint8_t *)s2 + n;
    }
    int32_t c1 = 0;
    int32_t c1_prev = 0;
    int32_t c2 = 0;
    int32_t c2_prev = 0;

    for (size_t i = 1; i < n + 1; i++) {
        // s = CC[0]
        s = costs[0];

        c2 = 0;
        ssize_t c2_len;
        if (!reverse) {
            c2_len = utf8proc_iterate(s2_ptr, -1, &c2);
        } else {
            c2_len = utf8proc_iterate_reversed(s2_ptr, -1, &c2);
        }
        if (c2 == 0) break;

        // CC[0] = c = t = t + h
        t += insert_gap_extend_cost;
        c = t;
        costs[0] = c;
        size_t prev_row_prev_cost = c;
        size_t prev_row_prev2_cost;  // cost for the previous row two characters ago for transposes
        
        prev_row_prev_char_ops = ops[0];
        base_ops.num_delete_gap_extensions++;
        prev_char_ops = base_ops;
        ops[0] = prev_char_ops;

        // e = t + g
        e = t + insert_gap_open_cost;

        alignment_op op = ALIGN_DELETE_GAP_OPEN;

        size_t min_cost = max_cost;
        uint8_t *s1_ptr; 
        if (!reverse) {
            s1_ptr = (uint8_t *)s1;
        } else {
            s1_ptr = (uint8_t *)s1 + m;
        }

        for (size_t j = 1; j < m + 1; j++) {
            // insertion
            // e = min(e, c + g) + h
            size_t min = e;

            c1 = 0;
            ssize_t c1_len = 0;
            if (!reverse) {
                c1_len = utf8proc_iterate(s1_ptr, -1, &c1);
            } else {
                c1_len = utf8proc_iterate_reversed(s1_ptr, -1, &c1);
            }
            if (c1 == 0) break;

            current_ops = ops[j];

            alignment_op insert_op = ALIGN_INSERT_GAP_OPEN;

            if ((c + insert_gap_open_cost) < min) {
                min = c + insert_gap_open_cost;
                insert_op = ALIGN_INSERT_GAP_OPEN;
            } else {
                insert_op = ALIGN_INSERT_GAP_EXTEND;
            }

            e = min + insert_gap_extend_cost;

            // deletion
            // DD[j] = min(DD[j], CC[j] + g) + h

            alignment_op delete_op = ALIGN_DELETE_GAP_OPEN;

            min = gap_costs[j];
            alignment_ops_t delete_ops = gap_ops[j];
            alignment_ops_t delete_ops_stored = delete_ops;
            delete_op = ALIGN_DELETE_GAP_OPEN;
            if (costs[j] + delete_gap_open_cost < min) {
                min = costs[j] + delete_gap_open_cost;

                delete_ops_stored = ops[j];
                delete_ops_stored.num_delete_gap_opens++;
            } else {
                delete_op = ALIGN_DELETE_GAP_EXTEND;
            }

            gap_costs[j] = min + delete_gap_extend_cost;
            delete_ops_stored.num_delete_gap_extensions++;
            gap_ops[j] = delete_ops_stored;

            // Cost
            // c = min(DD[j], e, s + w(a, b))

            alignment_op current_op = delete_op;


            min = gap_costs[j];

            if (e < min) {
                min = e;
                // Insert transition
                current_op = insert_op;
                current_ops = prev_char_ops;
            }

            bool c1_non_alphanumeric = utf8_is_non_alphanumeric((int32_t)c1);
            bool c2_non_alphanumeric = utf8_is_non_alphanumeric((int32_t)c2);
            bool both_non_alphanumeric = c1_non_alphanumeric && c2_non_alphanumeric;

            size_t w;
            if (c1 != c2) {
                w = mismatch_cost;
            } else {
                w = match_cost;
            }

            if (s + w < min) {
                min = s + w;

                // Match/mismatch transition
                current_ops = prev_row_prev_char_ops;

                if ((c1 == c2 || (ignore_non_alphanumeric && both_non_alphanumeric) )) {
                    current_op = ALIGN_MATCH;
                } else {
                    current_op = ALIGN_MISMATCH;
                }
            }

            if (j > 1 && i > 1 && c1 == c2_prev && c2 == c1_prev && c1 != c1_prev) {
                if (prev_row_prev2_cost + transpose_cost < min) {
                    s = prev_row_prev2_cost;
                    w = transpose_cost;
                    min = s + w;
                    current_op = ALIGN_TRANSPOSE;
                    current_ops = prev_row_prev2_char_ops;
                }
            }

            if (current_op == ALIGN_MATCH) {
                if (!both_non_alphanumeric || !ignore_non_alphanumeric) {
                    current_ops.num_matches++;
                }
            } else if (current_op == ALIGN_MISMATCH) {
                if (!both_non_alphanumeric && !ignore_non_alphanumeric) {
                    current_ops.num_mismatches++;
                }
            } else if (current_op == ALIGN_DELETE_GAP_EXTEND) {
                if (!c2_non_alphanumeric || !ignore_non_alphanumeric) {
                    current_ops.num_delete_gap_extensions++;
                }
            } else if (current_op == ALIGN_DELETE_GAP_OPEN) {
                if (!c2_non_alphanumeric || !ignore_non_alphanumeric) {
                    current_ops.num_delete_gap_opens++;
                    current_ops.num_delete_gap_extensions++;
                }
            } else if (current_op == ALIGN_INSERT_GAP_EXTEND) {
                if (!c1_non_alphanumeric || !ignore_non_alphanumeric) {
                    current_ops.num_insert_gap_extensions++;
                }
            } else if (current_op == ALIGN_INSERT_GAP_OPEN) {
                if (!c2_non_alphanumeric || !ignore_non_alphanumeric) {
                    current_ops.num_insert_gap_opens++;
                    current_ops.num_insert_gap_extensions++;
                }
            } else if (current_op == ALIGN_TRANSPOSE) {
                current_ops.num_transpositions++;
            }

            prev_row_prev2_cost = prev_costs[j];
            prev_costs[j] = prev_row_prev_cost;
            prev_row_prev_cost = costs[j];

            c = min;
            s = costs[j];
            costs[j] = c;

            prev_row_prev2_char_ops = prev_ops[j];
            prev_ops[j] = prev_row_prev_char_ops;
            prev_row_prev_char_ops = ops[j];

            prev_char_ops = current_ops;
            ops[j] = current_ops;

            if (!reverse) {
                s1_ptr += c1_len;
            } else {
                s1_ptr -= c1_len;
            }
            c1_prev = c1;
        }

        if (!reverse) {
            s2_ptr += c2_len;
        } else {
            s2_ptr -= c2_len;
        }
        c2_prev = c2;
    }

    alignment_ops_t best_ops = ops[m];
    if (swapped) {
        // swap the counts for insertions and deletions
        // so that they are in terms of the original strings
        size_t num_delete_gap_opens = best_ops.num_delete_gap_opens;
        size_t num_delete_gap_extensions = best_ops.num_delete_gap_extensions;
        size_t num_insert_gap_opens = best_ops.num_insert_gap_opens;
        size_t num_insert_gap_extensions = best_ops.num_insert_gap_extensions;

        best_ops.num_delete_gap_opens = num_insert_gap_opens;
        best_ops.num_delete_gap_extensions = num_insert_gap_extensions;
        best_ops.num_insert_gap_opens = num_delete_gap_opens;
        best_ops.num_insert_gap_extensions = num_delete_gap_extensions;
    }
    *result = best_ops;

    free(all_costs);
    free(all_ops);
    return true;
}

bool affine_gap_align_op_counts_options(const char *s1, size_t m, const char *s2, size_t n, alignment_options_t options, alignment_ops_t *result) {
    bool reverse = false;
    bool ret = affine_gap_costs_options(s1, m, s2, n, reverse, options, result);

    return ret;
}

bool affine_gap_align_op_counts(const char *s1, const char *s2, alignment_ops_t *result) {
    if (s1 == NULL || s2 == NULL) return false;
    size_t s1_size = strlen(s1) + 1;
    size_t m = utf8_len(s1, s1_size);
    size_t s2_size = strlen(s2) + 1;
    size_t n = utf8_len(s2, s2_size);
    return affine_gap_align_op_counts_options(s1, m, s2, n, DEFAULT_ALIGNMENT_OPTIONS_AFFINE_GAP, result);
}

