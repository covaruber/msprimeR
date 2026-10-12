/* R interface to the msprime mutation generator. */
#include "msprimer.h"

#include <float.h>

typedef struct {
    mutation_model_t model;
    bool model_allocated;
    mutgen_t mutgen;
    bool mutgen_allocated;
    gsl_rng *rng;
} msprimer_mutgen_t;

static void
mutgen_finalizer(SEXP xp)
{
    msprimer_mutgen_t *self = (msprimer_mutgen_t *) R_ExternalPtrAddr(xp);
    if (self != NULL) {
        if (self->mutgen_allocated) {
            mutgen_free(&self->mutgen);
        }
        if (self->model_allocated) {
            mutation_model_free(&self->model);
        }
        if (self->rng != NULL) {
            gsl_rng_free(self->rng);
        }
        free(self);
        R_ClearExternalPtr(xp);
    }
}

static double *
doubles_of(SEXP v, const char *what, R_xlen_t *n)
{
    R_xlen_t j;
    double *out;
    if (!Rf_isNumeric(v)) {
        Rf_error("'%s' must be numeric", what);
    }
    *n = XLENGTH(v);
    out = (double *) R_alloc((size_t) (*n > 0 ? *n : 1), sizeof(double));
    for (j = 0; j < *n; j++) {
        out[j] = TYPEOF(v) == REALSXP ? REAL(v)[j] : (double) INTEGER(v)[j];
    }
    return out;
}

/* Allocates a matrix mutation model from R vectors. The transition matrix is
 * an R (column-major) num_alleles x num_alleles matrix. */
static int
alloc_matrix_model(mutation_model_t *model, SEXP alleles, SEXP root_distribution,
    SEXP transition_matrix)
{
    R_xlen_t num_alleles, n_root, n_tm, j, k;
    char **allele_ptrs;
    size_t *allele_length;
    double *root, *tm, *row_major;
    SEXP dim;

    if (TYPEOF(alleles) != STRSXP) {
        Rf_error("alleles must be character strings");
    }
    num_alleles = XLENGTH(alleles);
    root = doubles_of(root_distribution, "root_distribution", &n_root);
    if (n_root != num_alleles) {
        Rf_error("root distribution must have num_alleles elements");
    }
    dim = Rf_getAttrib(transition_matrix, R_DimSymbol);
    if (dim == R_NilValue || XLENGTH(dim) != 2) {
        Rf_error("transition_matrix must be a matrix");
    }
    if (INTEGER(dim)[0] != INTEGER(dim)[1]) {
        Rf_error("Square matrix required");
    }
    if (INTEGER(dim)[0] != num_alleles) {
        Rf_error("transition matrix must be a square matrix with num_alleles rows");
    }
    tm = doubles_of(transition_matrix, "transition_matrix", &n_tm);
    row_major = (double *) R_alloc((size_t) (n_tm > 0 ? n_tm : 1), sizeof(double));
    for (j = 0; j < num_alleles; j++) {
        for (k = 0; k < num_alleles; k++) {
            row_major[j * num_alleles + k] = tm[k * num_alleles + j];
        }
    }
    allele_ptrs = (char **) R_alloc((size_t) (num_alleles > 0 ? num_alleles : 1),
        sizeof(char *));
    allele_length = (size_t *) R_alloc(
        (size_t) (num_alleles > 0 ? num_alleles : 1), sizeof(size_t));
    for (j = 0; j < num_alleles; j++) {
        if (STRING_ELT(alleles, j) == NA_STRING) {
            Rf_error("alleles must not be NA");
        }
        allele_ptrs[j] = (char *) Rf_translateCharUTF8(STRING_ELT(alleles, j));
        allele_length[j] = strlen(allele_ptrs[j]);
    }
    return matrix_mutation_model_alloc(
        model, (size_t) num_alleles, allele_ptrs, allele_length, root, row_major);
}

SEXP
msprimer_check_matrix_model(SEXP alleles, SEXP root_distribution, SEXP transition_matrix)
{
    mutation_model_t model;
    int err;

    memset(&model, 0, sizeof(model));
    err = alloc_matrix_model(&model, alleles, root_distribution, transition_matrix);
    mutation_model_free(&model);
    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
    return R_NilValue;
}

static void
alloc_model(msprimer_mutgen_t *self, SEXP spec)
{
    SEXP type = msprimer_list_get(spec, "type");
    const char *t;
    double next;
    int err;

    if (type == R_NilValue || TYPEOF(type) != STRSXP) {
        Rf_error("mutation model type not specified");
    }
    t = CHAR(STRING_ELT(type, 0));
    if (strcmp(t, "matrix") == 0) {
        err = alloc_matrix_model(&self->model, msprimer_list_get(spec, "alleles"),
            msprimer_list_get(spec, "root_distribution"),
            msprimer_list_get(spec, "transition_matrix"));
    } else if (strcmp(t, "slim_v5") == 0) {
        err = slim_v5_mutation_model_alloc(&self->model,
            (int32_t) msprimer_get_double(spec, "type_id", 0),
            (int64_t) msprimer_get_double(spec, "next_id", 0),
            (int32_t) msprimer_get_double(spec, "slim_generation", 1),
            (size_t) msprimer_get_double(spec, "block_size", 0));
    } else if (strcmp(t, "slim_v6") == 0) {
        err = slim_v6_mutation_model_alloc(&self->model,
            (int64_t) msprimer_get_double(spec, "next_id", 0),
            (size_t) msprimer_get_double(spec, "block_size", 0));
    } else if (strcmp(t, "infinite_alleles") == 0) {
        next = msprimer_get_double(spec, "start_allele", 0);
        if (next < 0) {
            Rf_error("start_allele must be non-negative");
        }
        err = infinite_alleles_mutation_model_alloc(&self->model, (uint64_t) next, 0);
        if (err == 0 && msprimer_has(spec, "next_allele")) {
            self->model.params.infinite_alleles.next_allele
                = (uint64_t) msprimer_get_double(spec, "next_allele", next);
        }
    } else {
        Rf_error("Unknown mutation model type '%s'", t);
    }
    self->model_allocated = true;
    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
}

static SEXP
model_state(msprimer_mutgen_t *self, SEXP spec)
{
    const char *names[] = { "next_id" };
    const char *ia_names[] = { "next_allele" };
    const char *t = CHAR(STRING_ELT(msprimer_list_get(spec, "type"), 0));
    SEXP ret;

    if (strcmp(t, "slim_v5") == 0 || strcmp(t, "slim_v6") == 0) {
        ret = PROTECT(msprimer_make_named_list(1, names));
        SET_VECTOR_ELT(ret, 0,
            Rf_ScalarReal((double) self->model.params.slim_mutator.next_mutation_id));
    } else if (strcmp(t, "infinite_alleles") == 0) {
        ret = PROTECT(msprimer_make_named_list(1, ia_names));
        SET_VECTOR_ELT(ret, 0,
            Rf_ScalarReal((double) self->model.params.infinite_alleles.next_allele));
    } else {
        ret = PROTECT(Rf_allocVector(VECSXP, 0));
    }
    UNPROTECT(1);
    return ret;
}

SEXP
msprimer_sim_mutations(SEXP tables, SEXP seed, SEXP position, SEXP rate, SEXP model,
    SEXP discrete_genome, SEXP keep, SEXP start_time, SEXP end_time)
{
    const char *names[] = { "tables", "model_state" };
    msprimer_mutgen_t *self;
    tsk_table_collection_t *tc;
    SEXP xp, tc_xp, ret;
    double seed_value = Rf_asReal(seed);
    double t0 = Rf_asReal(start_time), t1 = Rf_asReal(end_time);
    double *pos, *r;
    R_xlen_t np, nr;
    int flags = 0;
    int err;

    if (ISNAN(seed_value) || seed_value <= 0 || seed_value >= 4294967296.0
        || seed_value != floor(seed_value)) {
        Rf_error("seeds must be greater than 0 and less than 2^32");
    }
    tc_xp = PROTECT(msprimer_tc_alloc_xp(&tc));
    msprimer_tables_from_list(tc, tables);

    self = calloc(1, sizeof(*self));
    if (self == NULL) {
        Rf_error("Out of memory");
    }
    xp = PROTECT(R_MakeExternalPtr(self, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(xp, mutgen_finalizer, TRUE);
    self->rng = gsl_rng_alloc(gsl_rng_default);
    if (self->rng == NULL) {
        Rf_error("Out of memory");
    }
    gsl_rng_set(self->rng, (unsigned long) seed_value);
    alloc_model(self, model);

    err = mutgen_alloc(&self->mutgen, self->rng, tc, &self->model, 0);
    self->mutgen_allocated = true;
    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
    pos = doubles_of(position, "position", &np);
    r = doubles_of(rate, "rate", &nr);
    if (np != nr + 1) {
        Rf_error("The position array must be one larger than rates");
    }
    err = mutgen_set_rate_map(&self->mutgen, (size_t) nr, pos, r);
    if (err != 0) {
        msprimer_msp_error(err, "mutation rate map");
    }
    if (ISNAN(t0)) {
        t0 = -DBL_MAX;
    }
    if (ISNAN(t1)) {
        t1 = DBL_MAX;
    }
    err = mutgen_set_time_interval(&self->mutgen, t0, t1);
    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
    if (Rf_asLogical(discrete_genome) == TRUE) {
        flags |= MSP_DISCRETE_SITES;
    }
    if (Rf_asLogical(keep) == TRUE) {
        flags |= MSP_KEEP_SITES;
    }
    err = mutgen_generate(&self->mutgen, flags);
    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
    ret = PROTECT(msprimer_make_named_list(2, names));
    SET_VECTOR_ELT(ret, 0, msprimer_tables_to_list(tc));
    SET_VECTOR_ELT(ret, 1, model_state(self, model));
    mutgen_finalizer(xp);
    msprimer_tc_release(tc_xp);
    UNPROTECT(3);
    return ret;
}
