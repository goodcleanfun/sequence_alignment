#include "greatest/greatest.h"
#include "local_alignment.h"

/*
Tests for local affine gap alignment with op counts
*/

TEST test_affine_gap_op_count_correctness(void) {
    const char *s1 = "bam 30 lafyette ave bk new yrok 11217";
    const char *s2 = "brooklyn academy of music 30 lafayette avenue brooklyn new york";
    alignment_ops_t ops;
    ASSERT(affine_gap_align_op_counts(s1, s2, &ops));
    // Test that it finds the optimal number of matches, transposes, and gap extensions
    ASSERT_EQ(ops.num_matches, 23);
    ASSERT_EQ(ops.num_mismatches, 0);
    ASSERT_EQ(ops.num_transpositions, 1);
    ASSERT_EQ(ops.num_insert_gap_extensions, 5);
    ASSERT_EQ(ops.num_delete_gap_extensions, 29);

    PASS();
}

// From Flouri, Kobert, Rognes and Stamatakis, 2015 "Are all global alignment algorithms and implementations correct?"
TEST test_sequence_alignment_symmetric(void) {
    const char *s1 = "AAATTTGC";
    const char *s2 = "CGCCTTAC";
    alignment_ops_t ops;
    ASSERT(affine_gap_align_op_counts(s1, s2, &ops));
    // Test that it finds the optimal number of matches, which is 3
    ASSERT_EQ(ops.num_matches, 3);
    ASSERT_EQ(ops.num_mismatches, 2);
    ASSERT_EQ(ops.num_transpositions, 0);
    ASSERT_EQ(ops.num_delete_gap_extensions, 3);
    ASSERT_EQ(ops.num_insert_gap_extensions, 3);

    // Test symmetry - note: these have to be the same length for this test
    ASSERT_EQ(strlen(s1), strlen(s2));
    alignment_ops_t ops2;
    ASSERT(affine_gap_align_op_counts(s2, s1, &ops2));
    ASSERT_EQ(ops.num_matches, ops2.num_matches);
    ASSERT_EQ(ops.num_mismatches, ops2.num_mismatches);
    ASSERT_EQ(ops.num_transpositions, ops2.num_transpositions);
    ASSERT_EQ(ops.num_insert_gap_opens, ops2.num_insert_gap_opens);
    ASSERT_EQ(ops.num_insert_gap_extensions, ops2.num_insert_gap_extensions);
    ASSERT_EQ(ops.num_delete_gap_opens, ops2.num_delete_gap_opens);
    ASSERT_EQ(ops.num_delete_gap_extensions, ops2.num_delete_gap_extensions);

    PASS();
}

SUITE(test_local_alignment_suite) {
    RUN_TEST(test_sequence_alignment_symmetric);
}


/* Add definitions that need to be in the test runner's main file. */
GREATEST_MAIN_DEFS();

int main(int argc, char **argv) {
    GREATEST_MAIN_BEGIN();      /* command-line options, initialization. */

    /* Tests can also be gathered into test suites. */
    RUN_SUITE(test_local_alignment_suite);

    GREATEST_MAIN_END();        /* display results */
}
