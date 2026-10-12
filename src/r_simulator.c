/* R interface to the msprime simulator (msp_t), mirroring the Simulator class
 * of the Python low-level module. */
#include "msprimer.h"

#include <float.h>
#include <gsl/gsl_math.h>

typedef struct {
    msp_t sim;
    bool sim_allocated;
    gsl_rng *rng;
    tsk_table_collection_t tables;
    bool tables_initialised;
} msprimer_sim_t;

static void
sim_finalizer(SEXP xp)
{
    msprimer_sim_t *self = (msprimer_sim_t *) R_ExternalPtrAddr(xp);
    if (self != NULL) {
        if (self->sim_allocated) {
            msp_free(&self->sim);
        }
        if (self->tables_initialised) {
            tsk_table_collection_free(&self->tables);
        }
        if (self->rng != NULL) {
            gsl_rng_free(self->rng);
        }
        free(self);
        R_ClearExternalPtr(xp);
    }
}

/* Frees the simulator eagerly rather than waiting for garbage collection,
 * which R does not trigger based on memory allocated in C. */
SEXP
msprimer_sim_free(SEXP xp)
{
    if (TYPEOF(xp) == EXTPTRSXP) {
        sim_finalizer(xp);
    }
    return R_NilValue;
}

static msprimer_sim_t *
get_sim(SEXP xp)
{
    msprimer_sim_t *self;
    if (TYPEOF(xp) != EXTPTRSXP) {
        Rf_error("Not a simulator object");
    }
    self = (msprimer_sim_t *) R_ExternalPtrAddr(xp);
    if (self == NULL || !self->sim_allocated) {
        Rf_error("Simulator not initialised");
    }
    return self;
}

static double *
as_doubles(SEXP v, const char *what, R_xlen_t *n)
{
    R_xlen_t j;
    double *out;
    if (!Rf_isNumeric(v) && !Rf_isLogical(v)) {
        Rf_error("'%s' must be numeric", what);
    }
    *n = XLENGTH(v);
    out = (double *) R_alloc((size_t) (*n > 0 ? *n : 1), sizeof(double));
    for (j = 0; j < *n; j++) {
        out[j] = TYPEOF(v) == REALSXP ? REAL(v)[j]
                                      : (INTEGER(v)[j] == NA_INTEGER
                                                ? NAN
                                                : (double) INTEGER(v)[j]);
    }
    return out;
}

static int32_t *
as_int32s(SEXP v, const char *what, R_xlen_t *n)
{
    R_xlen_t j;
    double *x = as_doubles(v, what, n);
    int32_t *out = (int32_t *) R_alloc((size_t) (*n > 0 ? *n : 1), sizeof(int32_t));
    for (j = 0; j < *n; j++) {
        if (ISNAN(x[j]) || x[j] != floor(x[j]) || x[j] < INT32_MIN || x[j] > INT32_MAX) {
            Rf_error("'%s' must contain integers", what);
        }
        out[j] = (int32_t) x[j];
    }
    return out;
}

static const char *
model_name(SEXP model)
{
    SEXP name = msprimer_list_get(model, "name");
    if (name == R_NilValue || TYPEOF(name) != STRSXP || XLENGTH(name) != 1) {
        Rf_error("'name' not specified");
    }
    return CHAR(STRING_ELT(name, 0));
}

static double
required_number(SEXP list, const char *name)
{
    if (!msprimer_has(list, name)) {
        Rf_error("'%s' not specified", name);
    }
    return msprimer_get_double(list, name, 0);
}

static void
parse_simulation_model(msp_t *sim, SEXP model)
{
    const char *name = model_name(model);
    double psi, c, alpha, truncation_point, k;
    int err = 0;

    if (strcmp(name, "hudson") == 0) {
        err = msp_set_simulation_model_hudson(sim);
    } else if (strcmp(name, "dtwf") == 0) {
        err = msp_set_simulation_model_dtwf(sim);
    } else if (strcmp(name, "fixed_pedigree") == 0) {
        err = msp_set_simulation_model_fixed_pedigree(sim);
    } else if (strcmp(name, "smc") == 0) {
        err = msp_set_simulation_model_smc(sim);
    } else if (strcmp(name, "smc_prime") == 0) {
        err = msp_set_simulation_model_smc_prime(sim);
    } else if (strcmp(name, "smc_k") == 0) {
        k = required_number(model, "k");
        if (!(k >= 0.0)) {
            Rf_error("Must have k >= 0.0");
        }
        err = msp_set_simulation_model_smc_k(sim, k);
    } else if (strcmp(name, "dirac") == 0) {
        psi = required_number(model, "psi");
        c = required_number(model, "c");
        if (!(psi > 0 && psi <= 1.0)) {
            Rf_error("Must have 0 < psi <= 1");
        }
        if (!(c >= 0)) {
            Rf_error("c >= 0");
        }
        err = msp_set_simulation_model_dirac(sim, psi, c);
    } else if (strcmp(name, "beta") == 0) {
        alpha = required_number(model, "alpha");
        truncation_point = required_number(model, "truncation_point");
        err = msp_set_simulation_model_beta(sim, alpha, truncation_point);
    } else if (strcmp(name, "sweep_genic_selection") == 0) {
        err = msp_set_simulation_model_sweep_genic_selection(sim,
            required_number(model, "position"), required_number(model, "start_frequency"),
            required_number(model, "end_frequency"), required_number(model, "s"),
            required_number(model, "dt"));
        if (err != 0) {
            msprimer_msp_error(err, "sweep genic selection");
        }
    } else {
        Rf_error("Unknown simulation model");
    }
    if (err != 0) {
        msprimer_msp_error(err, "simulation model");
    }
}

static void
parse_population_configuration(msp_t *sim, SEXP pop_config)
{
    R_xlen_t j;
    SEXP item;
    int err;

    for (j = 0; j < XLENGTH(pop_config); j++) {
        item = VECTOR_ELT(pop_config, j);
        if (TYPEOF(item) != VECSXP) {
            Rf_error("population configuration must be a list of lists");
        }
        err = msp_set_population_configuration(sim, (int) j,
            required_number(item, "initial_size"), required_number(item, "growth_rate"),
            (bool) required_number(item, "initially_active"));
        if (err != 0) {
            msprimer_msp_error(err, "population configuration");
        }
    }
}

static void
parse_migration_matrix(msp_t *sim, SEXP migration_matrix)
{
    size_t N = msp_get_num_populations(sim);
    SEXP dim = Rf_getAttrib(migration_matrix, R_DimSymbol);
    R_xlen_t n;
    double *m, *row_major;
    size_t j, k;
    int err;

    if (dim == R_NilValue || XLENGTH(dim) != 2 || (size_t) INTEGER(dim)[0] != N
        || (size_t) INTEGER(dim)[1] != N) {
        Rf_error("migration matrix must be a N x N square matrix, where N is the "
                 "number of populations. The diagonal elements of this matrix must "
                 "be zero.");
    }
    m = as_doubles(migration_matrix, "migration_matrix", &n);
    row_major = (double *) R_alloc(N * N + 1, sizeof(double));
    for (j = 0; j < N; j++) {
        for (k = 0; k < N; k++) {
            row_major[j * N + k] = m[k * N + j];
        }
    }
    err = msp_set_migration_matrix(sim, N * N, row_major);
    if (err != 0) {
        msprimer_msp_error(err, "migration matrix");
    }
}

static double
optional_nan(SEXP item, const char *name)
{
    double x;
    if (!msprimer_has(item, name)) {
        return GSL_NAN;
    }
    x = msprimer_get_double(item, name, GSL_NAN);
    return ISNAN(x) ? GSL_NAN : x;
}

static void
parse_demographic_events(msp_t *sim, SEXP events)
{
    R_xlen_t j, N, a, b;
    SEXP item, value;
    const char *type;
    double time;
    int32_t *pops;
    double *props;
    R_xlen_t n2;
    int err;

    for (j = 0; j < XLENGTH(events); j++) {
        item = VECTOR_ELT(events, j);
        if (TYPEOF(item) != VECSXP) {
            Rf_error("demographic events must be a list of lists");
        }
        time = required_number(item, "time");
        if (time < 0) {
            Rf_error("negative times not valid");
        }
        value = msprimer_list_get(item, "type");
        if (value == R_NilValue || TYPEOF(value) != STRSXP) {
            Rf_error("'type' not specified");
        }
        type = CHAR(STRING_ELT(value, 0));
        err = 0;
        if (strcmp(type, "population_parameters_change") == 0) {
            err = msp_add_population_parameters_change(sim, time,
                (int) required_number(item, "population"),
                optional_nan(item, "initial_size"), optional_nan(item, "growth_rate"));
        } else if (strcmp(type, "migration_rate_change") == 0) {
            err = msp_add_migration_rate_change(sim, time,
                (int) required_number(item, "source"), (int) required_number(item, "dest"),
                required_number(item, "migration_rate"));
        } else if (strcmp(type, "mass_migration") == 0) {
            err = msp_add_mass_migration(sim, time, (int) required_number(item, "source"),
                (int) required_number(item, "dest"), required_number(item, "proportion"));
        } else if (strcmp(type, "activate_population_event") == 0) {
            err = msp_add_activate_population_event(
                sim, time, (int) required_number(item, "population"));
            if (err != 0) {
                msprimer_msp_error(err, "activate population event");
            }
        } else if (strcmp(type, "population_split") == 0) {
            value = msprimer_list_get(item, "derived");
            if (value == R_NilValue) {
                Rf_error("'derived' not specified");
            }
            pops = as_int32s(value, "derived", &N);
            if (N == 0) {
                Rf_error("Must provide at least one derived population");
            }
            err = msp_add_population_split(
                sim, time, (size_t) N, pops, (int) required_number(item, "ancestral"));
            if (err != 0) {
                msprimer_msp_error(err, "population split");
            }
        } else if (strcmp(type, "admixture") == 0) {
            value = msprimer_list_get(item, "ancestral");
            if (value == R_NilValue) {
                Rf_error("'ancestral' not specified");
            }
            pops = as_int32s(value, "ancestral", &N);
            if (N == 0) {
                Rf_error("Must provide at least one ancestral population");
            }
            value = msprimer_list_get(item, "proportions");
            if (value == R_NilValue) {
                Rf_error("'proportions' not specified");
            }
            props = as_doubles(value, "proportions", &n2);
            if (n2 != N) {
                Rf_error("proportions must be same size as ancestral");
            }
            err = msp_add_admixture(sim, time, (int) required_number(item, "derived"),
                (size_t) N, pops, props);
            if (err != 0) {
                msprimer_msp_error(err, "admixture");
            }
        } else if (strcmp(type, "symmetric_migration_rate_change") == 0) {
            value = msprimer_list_get(item, "populations");
            if (value == R_NilValue) {
                Rf_error("'populations' not specified");
            }
            pops = as_int32s(value, "populations", &N);
            if (N < 2) {
                Rf_error("Must provide at least two populations");
            }
            for (a = 0; a < N; a++) {
                for (b = 0; b < N; b++) {
                    if (a != b) {
                        err = msp_add_migration_rate_change(
                            sim, time, pops[a], pops[b], required_number(item, "rate"));
                        if (err != 0) {
                            msprimer_msp_error(err, "symmetric migration rate change");
                        }
                    }
                }
            }
        } else if (strcmp(type, "simple_bottleneck") == 0) {
            err = msp_add_simple_bottleneck(sim, time,
                (int) required_number(item, "population"),
                required_number(item, "proportion"));
        } else if (strcmp(type, "instantaneous_bottleneck") == 0) {
            err = msp_add_instantaneous_bottleneck(sim, time,
                (int) required_number(item, "population"),
                required_number(item, "strength"));
        } else if (strcmp(type, "census_event") == 0) {
            err = msp_add_census_event(sim, time);
        } else {
            Rf_error("Unknown demographic event type");
        }
        if (err != 0) {
            Rf_error("Input error in demographic_events[%d]: %s", (int) j,
                msp_strerror(err));
        }
    }
}

static void
parse_recombination_map(msp_t *sim, SEXP rate_map)
{
    R_xlen_t np, nr;
    double *position, *rate;
    int err;

    position = as_doubles(msprimer_list_get(rate_map, "position"), "position", &np);
    rate = as_doubles(msprimer_list_get(rate_map, "rate"), "rate", &nr);
    if (np != nr + 1) {
        Rf_error("The position array must be one larger than rates");
    }
    err = msp_set_recombination_map(sim, (size_t) nr, position, rate);
    if (err != 0) {
        msprimer_msp_error(err, "recombination map");
    }
}

SEXP
msprimer_sim_new(SEXP tables, SEXP seed, SEXP recomb_map, SEXP pop_config,
    SEXP migration_matrix, SEXP events, SEXP model, SEXP params)
{
    msprimer_sim_t *self;
    msp_t *sim;
    SEXP xp;
    double seed_value = Rf_asReal(seed);
    double start_time = msprimer_get_double(params, "start_time", -1);
    double gene_conversion_rate = msprimer_get_double(params, "gene_conversion_rate", 0);
    double gene_conversion_tract_length
        = msprimer_get_double(params, "gene_conversion_tract_length", 1.0);
    double additional_nodes = msprimer_get_double(params, "additional_nodes", 0);
    int err;

    if (ISNAN(seed_value) || seed_value <= 0 || seed_value >= 4294967296.0
        || seed_value != floor(seed_value)) {
        Rf_error("seeds must be greater than 0 and less than 2^32");
    }
    if (pop_config != R_NilValue && XLENGTH(pop_config) == 0) {
        Rf_error("Empty population configuration");
    }
    self = calloc(1, sizeof(*self));
    if (self == NULL) {
        Rf_error("Out of memory");
    }
    xp = PROTECT(R_MakeExternalPtr(self, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(xp, sim_finalizer, TRUE);

    self->rng = gsl_rng_alloc(gsl_rng_default);
    if (self->rng == NULL) {
        Rf_error("Out of memory");
    }
    gsl_rng_set(self->rng, (unsigned long) seed_value);
    err = tsk_table_collection_init(&self->tables, 0);
    if (err != 0) {
        msprimer_tsk_error(err, NULL);
    }
    self->tables_initialised = true;
    msprimer_tables_from_list(&self->tables, tables);

    sim = &self->sim;
    err = msp_alloc(sim, &self->tables, self->rng);
    self->sim_allocated = true;
    if (err != 0) {
        msprimer_msp_error(err, "simulator alloc");
    }
    if (model != R_NilValue) {
        parse_simulation_model(sim, model);
    }
    if (start_time >= 0) {
        err = msp_set_start_time(sim, start_time);
        if (err != 0) {
            msprimer_msp_error(err, "start time");
        }
    }
    err = msp_set_store_migrations(
        sim, msprimer_get_bool(params, "store_migrations", false));
    if (err != 0) {
        msprimer_msp_error(err, "store migrations");
    }
    err = msp_set_avl_node_block_size(
        sim, (size_t) msprimer_get_double(params, "avl_node_block_size", 10));
    if (err != 0) {
        msprimer_msp_error(err, "avl_node_block_size");
    }
    err = msp_set_segment_block_size(
        sim, (size_t) msprimer_get_double(params, "segment_block_size", 10));
    if (err != 0) {
        msprimer_msp_error(err, "segment_block_size");
    }
    err = msp_set_node_mapping_block_size(
        sim, (size_t) msprimer_get_double(params, "node_mapping_block_size", 10));
    if (err != 0) {
        msprimer_msp_error(err, "node_mapping_block_size");
    }
    msp_set_discrete_genome(sim, msprimer_get_bool(params, "discrete_genome", true));
    if (gene_conversion_rate != 0) {
        err = msp_set_gene_conversion_rate(sim, gene_conversion_rate);
        if (err != 0) {
            msprimer_msp_error(err, "set_gene_conversion_rate");
        }
        err = msp_set_gene_conversion_tract_length(sim, gene_conversion_tract_length);
        if (err != 0) {
            msprimer_msp_error(err, "set_gene_conversion_tract_length");
        }
    }
    if (recomb_map != R_NilValue) {
        parse_recombination_map(sim, recomb_map);
    }
    err = msp_set_ploidy(sim, msprimer_get_int(params, "ploidy", 2));
    if (err != 0) {
        msprimer_msp_error(err, "set_ploidy");
    }
    err = msp_set_num_labels(sim, (size_t) msprimer_get_int(params, "num_labels", 1));
    if (err != 0) {
        msprimer_msp_error(err, "set_num_labels");
    }
    if (pop_config != R_NilValue) {
        parse_population_configuration(sim, pop_config);
        if (migration_matrix == R_NilValue) {
            Rf_error("A migration matrix must be provided when a non-default "
                     "population configuration is used.");
        }
        parse_migration_matrix(sim, migration_matrix);
    } else if (migration_matrix != R_NilValue) {
        Rf_error("Cannot supply migration_matrix without population_configuration.");
    }
    if (events != R_NilValue) {
        parse_demographic_events(sim, events);
    }
    if (additional_nodes < 0 || additional_nodes >= 4294967296.0) {
        Rf_error("additional_nodes out of range");
    }
    msp_set_additional_nodes(sim, (uint32_t) additional_nodes);
    msp_set_coalescing_segments_only(
        sim, msprimer_get_bool(params, "coalescing_segments_only", true));
    msp_set_stop_at_local_mrca(sim, msprimer_get_bool(params, "stop_at_local_mrca", true));
    err = msp_initialise(sim);
    if (err != 0) {
        msprimer_msp_error(err, "initialise");
    }
    UNPROTECT(1);
    return xp;
}

SEXP
msprimer_sim_run(SEXP xp, SEXP r_end_time, SEXP r_max_events)
{
    msprimer_sim_t *self = get_sim(xp);
    double end_time = Rf_asReal(r_end_time);
    double max_events = Rf_asReal(r_max_events);
    int status;

    if (ISNAN(end_time) || end_time < 0) {
        Rf_error("end_time must be > 0");
    }
    if (!R_FINITE(end_time)) {
        end_time = DBL_MAX;
    }
    if (ISNAN(max_events) || max_events < 1) {
        Rf_error("max_events must be > 0");
    }
    if (max_events > (double) UINT32_MAX) {
        max_events = (double) UINT32_MAX;
    }
    status = msp_run(&self->sim, end_time, (unsigned long) max_events);
    if (status < 0) {
        msprimer_msp_error(status, NULL);
    }
    return Rf_ScalarInteger(status);
}

SEXP
msprimer_sim_set_model(SEXP xp, SEXP model)
{
    msprimer_sim_t *self = get_sim(xp);
    parse_simulation_model(&self->sim, model);
    return R_NilValue;
}

SEXP
msprimer_sim_get_model(SEXP xp)
{
    msprimer_sim_t *self = get_sim(xp);
    simulation_model_t *model = msp_get_model(&self->sim);
    const char *two[] = { "name", "", "" };
    SEXP ret;

    switch (model->type) {
        case MSP_MODEL_DIRAC:
            two[1] = "psi";
            two[2] = "c";
            ret = PROTECT(msprimer_make_named_list(3, two));
            SET_VECTOR_ELT(ret, 1, Rf_ScalarReal(model->params.dirac_coalescent.psi));
            SET_VECTOR_ELT(ret, 2, Rf_ScalarReal(model->params.dirac_coalescent.c));
            break;
        case MSP_MODEL_BETA:
            two[1] = "alpha";
            two[2] = "truncation_point";
            ret = PROTECT(msprimer_make_named_list(3, two));
            SET_VECTOR_ELT(ret, 1, Rf_ScalarReal(model->params.beta_coalescent.alpha));
            SET_VECTOR_ELT(
                ret, 2, Rf_ScalarReal(model->params.beta_coalescent.truncation_point));
            break;
        case MSP_MODEL_SWEEP:
            two[1] = "locus";
            ret = PROTECT(msprimer_make_named_list(2, two));
            SET_VECTOR_ELT(ret, 1, Rf_ScalarReal(model->params.sweep.position));
            break;
        case MSP_MODEL_SMC_K:
            two[1] = "k";
            ret = PROTECT(msprimer_make_named_list(2, two));
            SET_VECTOR_ELT(
                ret, 1, Rf_ScalarReal(model->params.smc_k_coalescent.hull_offset));
            break;
        default:
            ret = PROTECT(msprimer_make_named_list(1, two));
    }
    SET_VECTOR_ELT(ret, 0, Rf_mkString(msp_get_model_name(&self->sim)));
    UNPROTECT(1);
    return ret;
}

SEXP
msprimer_sim_finalise_tables(SEXP xp)
{
    msprimer_sim_t *self = get_sim(xp);
    int err = msp_finalise_tables(&self->sim);
    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
    return R_NilValue;
}

SEXP
msprimer_sim_reset(SEXP xp)
{
    msprimer_sim_t *self = get_sim(xp);
    int err = msp_reset(&self->sim);
    if (err < 0) {
        msprimer_msp_error(err, NULL);
    }
    return R_NilValue;
}

SEXP
msprimer_sim_get_tables(SEXP xp)
{
    msprimer_sim_t *self = get_sim(xp);
    return msprimer_tables_to_list(&self->tables);
}

SEXP
msprimer_sim_debug_demography(SEXP xp)
{
    msprimer_sim_t *self = get_sim(xp);
    double end_time;
    int err = msp_debug_demography(&self->sim, &end_time);
    if (err < 0) {
        msprimer_msp_error(err, NULL);
    }
    return Rf_ScalarReal(end_time);
}

SEXP
msprimer_sim_compute_population_size(SEXP xp, SEXP population, SEXP time)
{
    msprimer_sim_t *self = get_sim(xp);
    double size;
    int pop = Rf_asInteger(population);
    int err;

    if (pop == NA_INTEGER || pop < 0) {
        Rf_error("population out of bounds");
    }
    err = msp_compute_population_size(&self->sim, (size_t) pop, Rf_asReal(time), &size);
    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
    return Rf_ScalarReal(size);
}

static SEXP
size_t_matrix(const size_t *x, size_t N)
{
    SEXP ret = PROTECT(Rf_allocMatrix(REALSXP, (int) N, (int) N));
    size_t j, k;
    for (j = 0; j < N; j++) {
        for (k = 0; k < N; k++) {
            REAL(ret)[k * N + j] = (double) x[j * N + k];
        }
    }
    UNPROTECT(1);
    return ret;
}

SEXP
msprimer_sim_state(SEXP xp)
{
    const char *names[] = { "time", "start_time", "sequence_length", "ploidy",
        "discrete_genome", "record_migrations", "additional_nodes", "num_populations",
        "num_labels", "num_ancestors", "num_nodes", "num_edges", "num_migrations",
        "num_breakpoints", "num_common_ancestor_events",
        "num_rejected_common_ancestor_events", "num_recombination_events",
        "num_gene_conversion_events", "num_internal_gene_conversion_events",
        "num_noneffective_gene_conversion_events", "sum_internal_gc_tract_lengths",
        "num_multiple_recombination_events", "num_avl_node_blocks",
        "num_node_mapping_blocks", "num_segment_blocks", "num_fenwick_rebuilds",
        "segment_block_size", "avl_node_block_size", "node_mapping_block_size",
        "gene_conversion_tract_length", "migration_matrix", "num_migration_events",
        "population_configuration", "recombination_map" };
    const char *pc_names[] = { "initial_size", "growth_rate", "state" };
    const char *rm_names[] = { "position", "rate" };
    msprimer_sim_t *self = get_sim(xp);
    msp_t *sim = &self->sim;
    size_t N = msp_get_num_populations(sim);
    size_t *counts;
    double *mm;
    double initial_size, growth_rate;
    int state, err;
    size_t j, k;
    SEXP ret, m, pcs, pc, rm, pos, rate;

    ret = PROTECT(msprimer_make_named_list(34, names));
    SET_VECTOR_ELT(ret, 0, Rf_ScalarReal(msp_get_time(sim)));
    SET_VECTOR_ELT(ret, 1, Rf_ScalarReal(sim->start_time));
    SET_VECTOR_ELT(ret, 2, Rf_ScalarReal(sim->sequence_length));
    SET_VECTOR_ELT(ret, 3, Rf_ScalarInteger((int) sim->ploidy));
    SET_VECTOR_ELT(ret, 4, Rf_ScalarLogical(sim->discrete_genome));
    SET_VECTOR_ELT(ret, 5, Rf_ScalarLogical(msp_get_store_migrations(sim)));
    SET_VECTOR_ELT(ret, 6, Rf_ScalarReal((double) sim->additional_nodes));
    SET_VECTOR_ELT(ret, 7, Rf_ScalarReal((double) N));
    SET_VECTOR_ELT(ret, 8, Rf_ScalarReal((double) msp_get_num_labels(sim)));
    SET_VECTOR_ELT(ret, 9, Rf_ScalarReal((double) msp_get_num_ancestors(sim)));
    SET_VECTOR_ELT(ret, 10, Rf_ScalarReal((double) msp_get_num_nodes(sim)));
    SET_VECTOR_ELT(ret, 11, Rf_ScalarReal((double) msp_get_num_edges(sim)));
    SET_VECTOR_ELT(ret, 12, Rf_ScalarReal((double) msp_get_num_migrations(sim)));
    SET_VECTOR_ELT(ret, 13, Rf_ScalarReal((double) msp_get_num_breakpoints(sim)));
    SET_VECTOR_ELT(
        ret, 14, Rf_ScalarReal((double) msp_get_num_common_ancestor_events(sim)));
    SET_VECTOR_ELT(ret, 15,
        Rf_ScalarReal((double) msp_get_num_rejected_common_ancestor_events(sim)));
    SET_VECTOR_ELT(
        ret, 16, Rf_ScalarReal((double) msp_get_num_recombination_events(sim)));
    SET_VECTOR_ELT(
        ret, 17, Rf_ScalarReal((double) msp_get_num_gene_conversion_events(sim)));
    SET_VECTOR_ELT(ret, 18,
        Rf_ScalarReal((double) msp_get_num_internal_gene_conversion_events(sim)));
    SET_VECTOR_ELT(ret, 19,
        Rf_ScalarReal((double) msp_get_num_noneffective_gene_conversion_events(sim)));
    /* The Python module truncates this double to an integer; we keep it exact. */
    SET_VECTOR_ELT(ret, 20, Rf_ScalarReal(msp_get_sum_internal_gc_tract_lengths(sim)));
    SET_VECTOR_ELT(ret, 21, Rf_ScalarReal((double) sim->num_multiple_re_events));
    SET_VECTOR_ELT(ret, 22, Rf_ScalarReal((double) msp_get_num_avl_node_blocks(sim)));
    SET_VECTOR_ELT(
        ret, 23, Rf_ScalarReal((double) msp_get_num_node_mapping_blocks(sim)));
    SET_VECTOR_ELT(ret, 24, Rf_ScalarReal((double) msp_get_num_segment_blocks(sim)));
    SET_VECTOR_ELT(ret, 25, Rf_ScalarReal((double) sim->num_fenwick_rebuilds));
    SET_VECTOR_ELT(ret, 26, Rf_ScalarReal((double) sim->segment_block_size));
    SET_VECTOR_ELT(ret, 27, Rf_ScalarReal((double) sim->avl_node_block_size));
    SET_VECTOR_ELT(ret, 28, Rf_ScalarReal((double) sim->node_mapping_block_size));
    SET_VECTOR_ELT(ret, 29, Rf_ScalarReal(sim->gc_tract_length));

    mm = (double *) R_alloc(N * N + 1, sizeof(double));
    err = msp_get_migration_matrix(sim, mm);
    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
    m = PROTECT(Rf_allocMatrix(REALSXP, (int) N, (int) N));
    for (j = 0; j < N; j++) {
        for (k = 0; k < N; k++) {
            REAL(m)[k * N + j] = mm[j * N + k];
        }
    }
    SET_VECTOR_ELT(ret, 30, m);
    UNPROTECT(1);

    counts = (size_t *) R_alloc(N * N + 1, sizeof(size_t));
    err = msp_get_num_migration_events(sim, counts);
    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
    SET_VECTOR_ELT(ret, 31, size_t_matrix(counts, N));

    pcs = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) N));
    for (j = 0; j < N; j++) {
        msp_get_population_configuration(sim, j, &initial_size, &growth_rate, &state);
        pc = PROTECT(msprimer_make_named_list(3, pc_names));
        SET_VECTOR_ELT(pc, 0, Rf_ScalarReal(initial_size));
        SET_VECTOR_ELT(pc, 1, Rf_ScalarReal(growth_rate));
        SET_VECTOR_ELT(pc, 2, Rf_ScalarInteger(state));
        SET_VECTOR_ELT(pcs, (R_xlen_t) j, pc);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(ret, 32, pcs);
    UNPROTECT(1);

    rm = PROTECT(msprimer_make_named_list(2, rm_names));
    pos = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) sim->recomb_map.size + 1));
    rate = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) sim->recomb_map.size));
    memcpy(REAL(pos), sim->recomb_map.position,
        (sim->recomb_map.size + 1) * sizeof(double));
    if (sim->recomb_map.size > 0) {
        memcpy(REAL(rate), sim->recomb_map.rate, sim->recomb_map.size * sizeof(double));
    }
    SET_VECTOR_ELT(rm, 0, pos);
    SET_VECTOR_ELT(rm, 1, rate);
    SET_VECTOR_ELT(ret, 33, rm);
    UNPROTECT(3);

    UNPROTECT(1);
    return ret;
}

SEXP
msprimer_sim_ancestors(SEXP xp)
{
    msprimer_sim_t *self = get_sim(xp);
    size_t num_ancestors = msp_get_num_ancestors(&self->sim);
    segment_t **ancestors;
    segment_t *u;
    size_t j, num_segments, k;
    SEXP ret, m;
    int err;

    ancestors = (segment_t **) R_alloc(num_ancestors + 1, sizeof(segment_t *));
    err = msp_get_ancestors(&self->sim, ancestors);
    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
    ret = PROTECT(Rf_allocVector(VECSXP, (R_xlen_t) num_ancestors));
    for (j = 0; j < num_ancestors; j++) {
        num_segments = 0;
        for (u = ancestors[j]; u != NULL; u = u->next) {
            num_segments++;
        }
        m = PROTECT(Rf_allocMatrix(REALSXP, (int) num_segments, 4));
        k = 0;
        for (u = ancestors[j]; u != NULL; u = u->next) {
            REAL(m)[k] = u->left;
            REAL(m)[k + num_segments] = u->right;
            REAL(m)[k + 2 * num_segments] = (double) u->value;
            REAL(m)[k + 3 * num_segments] = (double) ancestors[j]->lineage->population;
            k++;
        }
        SET_VECTOR_ELT(ret, (R_xlen_t) j, m);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return ret;
}

SEXP
msprimer_sim_breakpoints(SEXP xp)
{
    msprimer_sim_t *self = get_sim(xp);
    size_t n = msp_get_num_breakpoints(&self->sim);
    size_t *bp = (size_t *) R_alloc(n + 1, sizeof(size_t));
    SEXP ret;
    size_t j;
    int err = msp_get_breakpoints(&self->sim, bp);

    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
    ret = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t) n));
    for (j = 0; j < n; j++) {
        REAL(ret)[j] = (double) bp[j];
    }
    UNPROTECT(1);
    return ret;
}

SEXP
msprimer_sim_print_state(SEXP xp, SEXP path)
{
    msprimer_sim_t *self = get_sim(xp);
    const char *filename = R_ExpandFileName(Rf_translateChar(STRING_ELT(path, 0)));
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        Rf_error("Cannot open file for writing");
    }
    msp_print_state(&self->sim, f);
    fclose(f);
    return R_NilValue;
}

SEXP
msprimer_sim_fenwick_drift(SEXP xp, SEXP r_label)
{
    msprimer_sim_t *self = get_sim(xp);
    int label = Rf_asInteger(r_label);
    double drift = 0;

    if (label == NA_INTEGER || label < 0 || label >= (int) msp_get_num_labels(&self->sim)) {
        Rf_error("bad label ID");
    }
    if (self->sim.recomb_mass_index != NULL) {
        drift = fenwick_get_numerical_drift(&self->sim.recomb_mass_index[label]);
    }
    return Rf_ScalarReal(drift);
}

/* Binary infinite-sites mutations on the simulator's own tables using its RNG,
 * as done by the legacy simulate(mutation_rate=...) interface. */
SEXP
msprimer_sim_simple_mutate(SEXP xp, SEXP r_rate)
{
    msprimer_sim_t *self = get_sim(xp);
    mutation_model_t model;
    mutgen_t mutgen;
    char *alleles[] = { "0", "1" };
    size_t allele_length[] = { 1, 1 };
    double root[] = { 1, 0 };
    double tm[] = { 0, 1, 1, 0 };
    double position[2];
    double rate = Rf_asReal(r_rate);
    int err;

    if (ISNAN(rate) || rate < 0) {
        Rf_error("mutation rate must be non-negative");
    }
    position[0] = 0;
    position[1] = self->tables.sequence_length;
    memset(&model, 0, sizeof(model));
    memset(&mutgen, 0, sizeof(mutgen));
    err = matrix_mutation_model_alloc(&model, 2, alleles, allele_length, root, tm);
    if (err == 0) {
        err = mutgen_alloc(&mutgen, self->rng, &self->tables, &model, 0);
    }
    if (err == 0) {
        err = mutgen_set_rate_map(&mutgen, 1, position, &rate);
    }
    if (err == 0) {
        err = mutgen_set_time_interval(&mutgen, -DBL_MAX, DBL_MAX);
    }
    if (err == 0) {
        err = mutgen_generate(&mutgen, 0);
    }
    mutgen_free(&mutgen);
    mutation_model_free(&model);
    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
    return R_NilValue;
}
