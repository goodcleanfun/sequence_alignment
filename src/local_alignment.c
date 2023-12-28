#include "local_alignment.h"
#include "utf8/utf8.h"


static alignment_ops_t NULL_ALIGNMENT_OPS = {
    .num_matches = 0,
    .num_mismatches = 0,
    .num_transpositions = 0,
    .num_gap_opens = 0,
    .num_gap_extensions = 0
};

static inline bool utf8_is_non_alphanumeric(int32_t c) {
    int cat = utf8proc_category(c);
    return utf8_is_whitespace(c) || utf8_is_hyphen(c) || utf8_is_punctuation(cat);
}

alignment_ops_t affine_gap_align_op_counts_unicode_options(uint32_array *u1_array, uint32_array *u2_array, alignment_options_t options) {
    if (u1_array->n < u2_array->n) {
        uint32_array *tmp_array = u1_array;
        u1_array = u2_array;
        u2_array = tmp_array;
    }

    size_t m = u1_array->n;
    size_t n = u2_array->n;

    uint32_t *u1 = u1_array->a;
    uint32_t *u2 = u2_array->a;

    alignment_ops_t result = NULL_ALIGNMENT_OPS;

    if (unicode_equals(u1_array, u2_array)) {
        result.num_matches = n;
        return result;
    }

    size_t num_bytes = (m + 1) * sizeof(size_t);

    size_t *costs = malloc(num_bytes);
    if (costs == NULL) {
        return NULL_ALIGNMENT_OPS;
    }

    size_t *prev_costs = malloc(num_bytes);
    if (prev_costs == NULL) {
        free(costs);
        return NULL_ALIGNMENT_OPS;
    }

    size_t *gap_costs = malloc(num_bytes);
    if (gap_costs == NULL) {
        free(prev_costs);
        free(costs);
        return NULL_ALIGNMENT_OPS;
    }

    alignment_ops_t *edits = malloc((m + 1) * sizeof(alignment_ops_t));
    if (edits == NULL) {
        free(gap_costs);
        free(prev_costs);
        free(costs);
        return NULL_ALIGNMENT_OPS;
    }


    alignment_ops_t *prev_edits = malloc((m + 1) * sizeof(alignment_ops_t));
    if (prev_edits == NULL) {
        free(edits);
        free(gap_costs);
        free(prev_costs);
        free(costs);
        return NULL_ALIGNMENT_OPS;
    }

    alignment_ops_t *gap_edits = malloc((m + 1) * sizeof(alignment_ops_t));
    if (gap_edits == NULL) {
        free(prev_edits);
        free(edits);
        free(gap_costs);
        free(prev_costs);
        free(costs);
        return NULL_ALIGNMENT_OPS;
    }

    size_t gap_open_cost = options.gap_open_cost;
    size_t gap_extend_cost = options.gap_extend_cost;
    size_t match_cost = options.match_cost;
    size_t mismatch_cost = options.mismatch_cost;
    size_t transpose_cost = options.transpose_cost;
    bool ignore_non_alphanumeric = options.ignore_non_alphanumeric;

    size_t e = 0, c = 0, s = 0;

    size_t max_cost = m + n;

    costs[0] = 0;
    gap_costs[0] = max_cost;
    edits[0] = NULL_ALIGNMENT_OPS;
    size_t t = gap_open_cost;

    alignment_ops_t base_edits = NULL_ALIGNMENT_OPS;
    base_edits.num_gap_opens++;

    for (size_t j = 1; j < m + 1; j++) {
        t += gap_extend_cost;
        costs[j] = t;
        gap_costs[j] = t + gap_open_cost;
        base_edits.num_gap_extensions++;
        edits[j] = base_edits;
        gap_edits[j] = base_edits;
    }

    t = gap_open_cost;
    base_edits = NULL_ALIGNMENT_OPS;
    base_edits.num_gap_opens++;

    alignment_ops_t current_edits = NULL_ALIGNMENT_OPS;
    alignment_ops_t prev_char_edits = NULL_ALIGNMENT_OPS;
    alignment_ops_t prev_row_prev_char_edits = NULL_ALIGNMENT_OPS;
    alignment_ops_t prev_row_prev2_char_edits = NULL_ALIGNMENT_OPS;

    bool in_gap = false;
    size_t prev_c = 0;

    for (size_t i = 1; i < n + 1; i++) {
        // s = CC[0]
        s = costs[0];
        uint32_t c2 = u2[i - 1];
        // CC[0] = c = t = t + h
        t += gap_extend_cost;
        c = t;
        costs[0] = c;
        size_t prev_row_prev_cost = c;
        size_t prev_row_prev2_cost;  // cost for the previous row two characters ago for transposes
        
        prev_row_prev_char_edits = edits[0];
        base_edits.num_gap_extensions++;
        prev_char_edits = base_edits;
        edits[0] = prev_char_edits;

        // e = t + g
        e = t + gap_open_cost;

        alignment_op op = ALIGN_GAP_OPEN;

        size_t min_cost = max_cost;

        for (size_t j = 1; j < m + 1; j++) {
            // insertion
            // e = min(e, c + g) + h
            size_t min = e;
            uint32_t c1 = u1[j - 1];

            current_edits = edits[j];

            alignment_op insert_op = ALIGN_GAP_OPEN;

            if ((c + gap_open_cost) < min) {
                min = c + gap_open_cost;
                insert_op = ALIGN_GAP_OPEN;
            } else {
                insert_op = ALIGN_GAP_EXTEND;
            }

            e = min + gap_extend_cost;

            // deletion
            // DD[j] = min(DD[j], CC[j] + g) + h

            alignment_op delete_op = ALIGN_GAP_OPEN;

            min = gap_costs[j];
            alignment_ops_t delete_edits = gap_edits[j];
            alignment_ops_t delete_edits_stored = delete_edits;
            delete_op = ALIGN_GAP_OPEN;
            if (costs[j] + gap_open_cost < min) {
                min = costs[j] + gap_open_cost;
                
                delete_edits_stored = edits[j];
                delete_edits_stored.num_gap_opens++;
            }

            gap_costs[j] = min + gap_extend_cost;
            delete_edits_stored.num_gap_extensions++;
            gap_edits[j] = delete_edits_stored;

            // Cost
            // c = min(DD[j], e, s + w(a, b))

            alignment_op current_op = delete_op;


            min = gap_costs[j];

            if (e < min) {
                min = e;
                // Insert transition
                current_op = insert_op;
                current_edits = prev_char_edits;
            }

            bool c1_non_alphanumeric = utf8_is_non_alphanumeric((int32_t)c1);
            bool c2_non_alphanumeric = utf8_is_non_alphanumeric((int32_t)c2);
            bool both_non_alphanumeric = c1_non_alphanumeric && c2_non_alphanumeric;

            size_t w = c1 != c2 && !both_non_alphanumeric ? mismatch_cost : match_cost;

            if (s + w < min) {
                min = s + w;

                // Match/mismatch/transpose transition
                current_edits = prev_row_prev_char_edits;

                if ((c1 == c2 || both_non_alphanumeric)) {
                    current_op = ALIGN_MATCH;
                } else {
                    current_op = ALIGN_MISMATCH;
                }
            }

            if (j > 1 && i > 1 && c1 == u2[i - 2] && c2 == u1[j - 2] && c1 != u1[j - 2]) {
                if (prev_row_prev2_cost + transpose_cost < min) {
                    s = prev_row_prev2_cost;
                    w = transpose_cost;
                    min = s + w;
                    current_op = ALIGN_TRANSPOSE;
                    current_edits = prev_row_prev2_char_edits;
                }
            }

            if (current_op == ALIGN_MATCH) {
                if (!both_non_alphanumeric) {
                    current_edits.num_matches++;
                }
            } else if (current_op == ALIGN_MISMATCH) {
                if (!c1_non_alphanumeric && !c2_non_alphanumeric) {
                    current_edits.num_mismatches++;
                }
            } else if (current_op == ALIGN_GAP_EXTEND) {
                current_edits.num_gap_extensions++;
            } else if (current_op == ALIGN_GAP_OPEN) {
                current_edits.num_gap_opens++;
                current_edits.num_gap_extensions++;
            } else if (current_op == ALIGN_TRANSPOSE) {
                current_edits.num_transpositions++;
            }

            prev_row_prev2_cost = prev_costs[j];
            prev_costs[j] = prev_row_prev_cost;
            prev_row_prev_cost = costs[j];

            c = min;
            s = costs[j];
            costs[j] = c;

            prev_row_prev2_char_edits = prev_edits[j];
            prev_edits[j] = prev_row_prev_char_edits;
            prev_row_prev_char_edits = edits[j];

            prev_char_edits = current_edits;
            edits[j] = current_edits;

        }
    }

    result = edits[m];
    free(costs);
    free(prev_costs);
    free(gap_costs);
    free(edits);
    free(prev_edits);
    free(gap_edits);

    return result;

}


alignment_ops_t affine_gap_align_op_counts_unicode(uint32_array *u1_array, uint32_array *u2_array) {
    return affine_gap_align_op_counts_unicode_options(u1_array, u2_array, DEFAULT_ALIGNMENT_OPTIONS_AFFINE_GAP);
}

alignment_ops_t affine_gap_align_op_counts_options(const char *s1, const char *s2, alignment_options_t options) {
    if (s1 == NULL || s2 == NULL) return NULL_ALIGNMENT_OPS;

    uint32_array *u1_array = unicode_codepoints(s1);
    if (u1_array == NULL) return NULL_ALIGNMENT_OPS;

    uint32_array *u2_array = unicode_codepoints(s2);

    if (u2_array == NULL) {
        uint32_array_destroy(u1_array);
        return NULL_ALIGNMENT_OPS;
    }

    alignment_ops_t edits = affine_gap_align_op_counts_unicode_options(u1_array, u2_array, options);

    uint32_array_destroy(u1_array);
    uint32_array_destroy(u2_array);

    return edits;
}


alignment_ops_t affine_gap_align_op_counts(const char *s1, const char *s2) {
    return affine_gap_align_op_counts_options(s1, s2, DEFAULT_ALIGNMENT_OPTIONS_AFFINE_GAP);
}
