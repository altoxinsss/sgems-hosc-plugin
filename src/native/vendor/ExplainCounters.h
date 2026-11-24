#ifndef __EXPLAIN_COUNTERS_H__
#define __EXPLAIN_COUNTERS_H__

#include <stdio.h>

typedef struct
{
    unsigned __int64 template_attempts;
    unsigned __int64 template_success;
    unsigned __int64 template_fail;

    unsigned __int64 cumulant_eval_attempts;
    unsigned __int64 cumulant_eval_success;
    unsigned __int64 cumulant_eval_fail;

    unsigned __int64 moment_precompute_calls;
    unsigned __int64 moment_precompute_fail;
    unsigned __int64 moment_eval_calls;

    unsigned __int64 original_recursive_calls;
    unsigned __int64 original_recursive_node_visits;
    unsigned __int64 original_recursive_cache_queries;
    unsigned __int64 original_recursive_cache_hits;

    unsigned __int64 optimized_iterative_calls;

    unsigned __int64 avl_search_calls;
    unsigned __int64 avl_search_hits;
    unsigned __int64 avl_insert_calls;
} HOSCExplainCounters;

extern HOSCExplainCounters g_hosc_counters;

void HOSCResetCounters(void);
void HOSCDumpCounters(FILE* out);

#endif
