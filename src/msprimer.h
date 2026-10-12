#ifndef MSPRIMER_H
#define MSPRIMER_H

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "msprime.h"
#include "likelihood.h"

/* Error helpers: msp_strerror also understands tskit error codes. */
void msprimer_msp_error(int err, const char *context);
void msprimer_tsk_error(int err, const char *context);

/* Named list helpers */
SEXP msprimer_list_get(SEXP list, const char *name);
bool msprimer_has(SEXP list, const char *name);
double msprimer_get_double(SEXP list, const char *name, double dflt);
int msprimer_get_int(SEXP list, const char *name, int dflt);
bool msprimer_get_bool(SEXP list, const char *name, bool dflt);
SEXP msprimer_make_named_list(int n, const char **names);

/* Table collection <-> named list conversion (tskit asdict() layout) */
SEXP msprimer_tc_alloc_xp(tsk_table_collection_t **out);
void msprimer_tc_release(SEXP xp);
void msprimer_tables_from_list(tsk_table_collection_t *tables, SEXP d);
SEXP msprimer_tables_to_list(tsk_table_collection_t *tables);

/* .Call entry points */
SEXP msprimer_versions(void);
SEXP msprimer_tables_roundtrip(SEXP d);
SEXP msprimer_unpack_strings(SEXP data, SEXP offset);
SEXP msprimer_tables_dump(SEXP d, SEXP path);
SEXP msprimer_tables_load(SEXP path);
SEXP msprimer_ts_summary(SEXP d);
SEXP msprimer_tables_delete_intervals(SEXP d, SEXP left, SEXP right);
SEXP msprimer_tables_sort(SEXP d);

SEXP msprimer_sim_new(SEXP tables, SEXP seed, SEXP recomb_map, SEXP pop_config,
    SEXP migration_matrix, SEXP events, SEXP model, SEXP params);
SEXP msprimer_sim_run(SEXP xp, SEXP end_time, SEXP max_events);
SEXP msprimer_sim_set_model(SEXP xp, SEXP model);
SEXP msprimer_sim_get_model(SEXP xp);
SEXP msprimer_sim_finalise_tables(SEXP xp);
SEXP msprimer_sim_reset(SEXP xp);
SEXP msprimer_sim_get_tables(SEXP xp);
SEXP msprimer_sim_debug_demography(SEXP xp);
SEXP msprimer_sim_compute_population_size(SEXP xp, SEXP population, SEXP time);
SEXP msprimer_sim_state(SEXP xp);
SEXP msprimer_sim_ancestors(SEXP xp);
SEXP msprimer_sim_breakpoints(SEXP xp);
SEXP msprimer_sim_print_state(SEXP xp, SEXP path);
SEXP msprimer_sim_fenwick_drift(SEXP xp, SEXP label);
SEXP msprimer_sim_simple_mutate(SEXP xp, SEXP rate);
SEXP msprimer_sim_free(SEXP xp);

SEXP msprimer_sim_mutations(SEXP tables, SEXP seed, SEXP position, SEXP rate,
    SEXP model, SEXP discrete_genome, SEXP keep, SEXP start_time, SEXP end_time);
SEXP msprimer_check_matrix_model(SEXP alleles, SEXP root_distribution,
    SEXP transition_matrix);

SEXP msprimer_log_likelihood_arg(SEXP tables, SEXP Ne, SEXP recombination_rate);
SEXP msprimer_log_mutation_likelihood(SEXP tables, SEXP mutation_rate);

#endif
