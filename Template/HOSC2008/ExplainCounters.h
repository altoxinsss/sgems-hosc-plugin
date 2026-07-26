#ifndef __EXPLAIN_COUNTERS_H__
#define __EXPLAIN_COUNTERS_H__

#include <stdio.h>

typedef struct
{
    unsigned long long template_attempts;
    unsigned long long template_success;
    unsigned long long template_fail;

    unsigned long long cumulant_eval_attempts;
    unsigned long long cumulant_eval_success;
    unsigned long long cumulant_eval_fail;

    unsigned long long moment_precompute_calls;
    unsigned long long moment_precompute_fail;
    unsigned long long moment_eval_calls;

    unsigned long long original_recursive_calls;
    unsigned long long original_recursive_node_visits;
    unsigned long long original_recursive_cache_queries;
    unsigned long long original_recursive_cache_hits;

    unsigned long long optimized_iterative_calls;

    unsigned long long avl_search_calls;
    unsigned long long avl_search_hits;
    unsigned long long avl_insert_calls;
} HOSCExplainCounters;

extern HOSCExplainCounters g_hosc_counters;

void HOSCResetCounters(void);
void HOSCDumpCounters(FILE* out);

#endif
