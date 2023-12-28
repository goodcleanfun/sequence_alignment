#include "local_alignment.h"
#include "utf8/utf8.h"
#include "utf8proc/utf8proc.h"

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

alignment_ops_t affine_gap_align_op_counts_options(const char *s1, size_t m, const char *s2, size_t n, alignment_options_t options) {
    if (m < n) {
        const char *tmp_str = s1;
        s1 = s2;
        s2 = tmp_str;
        size_t tmp_len = m;
        m = n;
        n = tmp_len;
    }

    alignment_ops_t result = NULL_ALIGNMENT_OPS;

    if (m == n && utf8_compare_len(s1, s2, n) == 0) {
        result.num_matches = n;
        return result;
    }

    // Costs are a 3 x (m + 1) matrix, only uses 3 rows
    size_t *all_costs = malloc((m + 1) * 3 * sizeof(size_t));
    if (all_costs == NULL) {
        return NULL_ALIGNMENT_OPS;
    }

    // Costs for the current row, CC in the Myers-Miller paper
    size_t *costs = all_costs;
    // Costs for previous row, used for transpositions as in Damerau-Levenshtein
    // This will track the cost at the i-1, j-2 position in order to revert back
    // to the cost at that position if the transpose is cheaper
    size_t *prev_costs = all_costs + (m + 1);
    // Gap costs for inserts/deletes, DD in the Myers-Miller paper
    size_t *gap_costs = all_costs + 2 * (m + 1);

    // Edits are a 3 x (m + 1) matrix of alignment op counts, a breakdown of the cost by type
    alignment_ops_t *all_edits = malloc(3 * (m + 1) * sizeof(alignment_ops_t));
    if (all_edits == NULL) {
        free(all_costs);
        return NULL_ALIGNMENT_OPS;
    }

    // Edits for the current row
    alignment_ops_t *edits = all_edits;
    // Edits for the previous row, used for transpositions in order to
    // revert back to the edits at the i - 1, j -2 at the end of the transpose
    // if the transpose is cheaper
    alignment_ops_t *prev_edits = all_edits + (m + 1);
    // Edits for the current row, only for insertions/deletions
    alignment_ops_t *gap_edits = all_edits + 2 * (m + 1);

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

    uint8_t *s2_ptr = (uint8_t *)s2;
    int32_t c1 = 0;
    int32_t c1_prev = 0;
    int32_t c2 = 0;
    int32_t c2_prev = 0;

    for (size_t i = 1; i < n + 1; i++) {
        // s = CC[0]
        s = costs[0];

        c2 = 0;
        ssize_t c2_len = utf8proc_iterate(s2_ptr, -1, &c2);
        if (c2 == 0) break;

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
        uint8_t *s1_ptr = (uint8_t *)s1;

        for (size_t j = 1; j < m + 1; j++) {
            // insertion
            // e = min(e, c + g) + h
            size_t min = e;

            c1 = 0;
            ssize_t c1_len = utf8proc_iterate(s1_ptr, -1, &c1);
            if (c1 == 0) break;

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

            if (j > 1 && i > 1 && c1 == c2_prev && c2 == c1_prev && c1 != c1_prev) {
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

            s1_ptr += c1_len;
            c1_prev = c1;
        }

        s2_ptr += c2_len;
        c2_prev = c2;
    }

    result = edits[m];
    free(all_costs);
    free(all_edits);

    return result;

}

alignment_ops_t affine_gap_align_op_counts(const char *s1, const char *s2) {
    if (s1 == NULL || s2 == NULL) return NULL_ALIGNMENT_OPS;
    size_t s1_len = strlen(s1);
    size_t m = utf8_len(s1, s1_len);
    size_t s2_len = strlen(s2);
    size_t n = utf8_len(s2, s2_len);
    return affine_gap_align_op_counts_options(s1, m, s2, n, DEFAULT_ALIGNMENT_OPTIONS_AFFINE_GAP);
}
