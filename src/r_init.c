#include "msprimer.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_version.h>

void
msprimer_bug_assert_failed(const char *file, int line)
{
    Rf_error("msprimeR internal bug detected in %s at line %d. Please report this "
             "issue, ideally with a reproducible example.",
        file, line);
}

void
msprimer_msp_error(int err, const char *context)
{
    if (context != NULL) {
        Rf_error("Input error in %s: %s", context, msp_strerror(err));
    }
    Rf_error("%s", msp_strerror(err));
}

void
msprimer_tsk_error(int err, const char *context)
{
    if (context != NULL) {
        Rf_error("%s: %s", context, tsk_strerror(err));
    }
    Rf_error("%s", tsk_strerror(err));
}

SEXP
msprimer_list_get(SEXP list, const char *name)
{
    SEXP names;
    R_xlen_t j;

    if (list == R_NilValue || TYPEOF(list) != VECSXP) {
        return R_NilValue;
    }
    names = Rf_getAttrib(list, R_NamesSymbol);
    if (names == R_NilValue) {
        return R_NilValue;
    }
    for (j = 0; j < XLENGTH(list); j++) {
        if (strcmp(CHAR(STRING_ELT(names, j)), name) == 0) {
            return VECTOR_ELT(list, j);
        }
    }
    return R_NilValue;
}

bool
msprimer_has(SEXP list, const char *name)
{
    return msprimer_list_get(list, name) != R_NilValue;
}

double
msprimer_get_double(SEXP list, const char *name, double dflt)
{
    SEXP v = msprimer_list_get(list, name);
    if (v == R_NilValue) {
        return dflt;
    }
    if (!Rf_isNumeric(v) && !Rf_isLogical(v)) {
        Rf_error("'%s' is not a number", name);
    }
    if (XLENGTH(v) != 1) {
        Rf_error("'%s' must be a single number", name);
    }
    return Rf_asReal(v);
}

int
msprimer_get_int(SEXP list, const char *name, int dflt)
{
    double x = msprimer_get_double(list, name, (double) dflt);
    if (ISNAN(x) || x != (double) (int) x) {
        Rf_error("'%s' must be an integer", name);
    }
    return (int) x;
}

bool
msprimer_get_bool(SEXP list, const char *name, bool dflt)
{
    double x = msprimer_get_double(list, name, dflt ? 1.0 : 0.0);
    if (ISNAN(x)) {
        Rf_error("'%s' must be TRUE or FALSE", name);
    }
    return x != 0;
}

SEXP
msprimer_make_named_list(int n, const char **names)
{
    SEXP ret = PROTECT(Rf_allocVector(VECSXP, n));
    SEXP nms = PROTECT(Rf_allocVector(STRSXP, n));
    int j;
    for (j = 0; j < n; j++) {
        SET_STRING_ELT(nms, j, Rf_mkChar(names[j]));
    }
    Rf_setAttrib(ret, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ret;
}

SEXP
msprimer_versions(void)
{
    const char *names[] = { "gsl", "tskit", "kastore" };
    SEXP ret = PROTECT(msprimer_make_named_list(3, names));
    SEXP gsl = PROTECT(Rf_allocVector(INTSXP, 2));
    SEXP tsk = PROTECT(Rf_allocVector(INTSXP, 3));
    SEXP kas = PROTECT(Rf_allocVector(INTSXP, 3));

    INTEGER(gsl)[0] = GSL_MAJOR_VERSION;
    INTEGER(gsl)[1] = GSL_MINOR_VERSION;
    INTEGER(tsk)[0] = TSK_VERSION_MAJOR;
    INTEGER(tsk)[1] = TSK_VERSION_MINOR;
    INTEGER(tsk)[2] = TSK_VERSION_PATCH;
    INTEGER(kas)[0] = KAS_VERSION_MAJOR;
    INTEGER(kas)[1] = KAS_VERSION_MINOR;
    INTEGER(kas)[2] = KAS_VERSION_PATCH;
    SET_VECTOR_ELT(ret, 0, gsl);
    SET_VECTOR_ELT(ret, 1, tsk);
    SET_VECTOR_ELT(ret, 2, kas);
    UNPROTECT(4);
    return ret;
}

#define CALLDEF(name, n) { #name, (DL_FUNC) &name, n }

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(msprimer_versions, 0),
    CALLDEF(msprimer_tables_roundtrip, 1),
    CALLDEF(msprimer_unpack_strings, 2),
    CALLDEF(msprimer_tables_dump, 2),
    CALLDEF(msprimer_tables_load, 1),
    CALLDEF(msprimer_ts_summary, 1),
    CALLDEF(msprimer_tables_delete_intervals, 3),
    CALLDEF(msprimer_tables_sort, 1),
    CALLDEF(msprimer_sim_new, 8),
    CALLDEF(msprimer_sim_run, 3),
    CALLDEF(msprimer_sim_set_model, 2),
    CALLDEF(msprimer_sim_get_model, 1),
    CALLDEF(msprimer_sim_finalise_tables, 1),
    CALLDEF(msprimer_sim_reset, 1),
    CALLDEF(msprimer_sim_get_tables, 1),
    CALLDEF(msprimer_sim_debug_demography, 1),
    CALLDEF(msprimer_sim_compute_population_size, 3),
    CALLDEF(msprimer_sim_state, 1),
    CALLDEF(msprimer_sim_ancestors, 1),
    CALLDEF(msprimer_sim_breakpoints, 1),
    CALLDEF(msprimer_sim_print_state, 2),
    CALLDEF(msprimer_sim_fenwick_drift, 2),
    CALLDEF(msprimer_sim_simple_mutate, 2),
    CALLDEF(msprimer_sim_free, 1),
    CALLDEF(msprimer_sim_mutations, 9),
    CALLDEF(msprimer_check_matrix_model, 3),
    CALLDEF(msprimer_log_likelihood_arg, 3),
    CALLDEF(msprimer_log_mutation_likelihood, 2),
    { NULL, NULL, 0 },
};

void
R_init_msprimeR(DllInfo *dll)
{
    /* GSL's default error handler calls abort(); errors are reported through
     * return codes instead. */
    gsl_set_error_handler_off();
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
