/* Conversion between tskit table collections and R named lists using the
 * same layout as tskit's (and tskitR's) TableCollection.asdict():
 *   - u32 columns (flags) are doubles, id columns are integers,
 *   - text/bytes columns are raw vectors with double offsets,
 *   - metadata schemas and time units are character scalars.
 */
#include "msprimer.h"

/*************************
 * Column readers
 *************************/

static SEXP
get_column(SEXP table, const char *table_name, const char *name, bool required)
{
    SEXP v = msprimer_list_get(table, name);
    if (v == R_NilValue && required) {
        Rf_error("%s table: required column '%s' not specified", table_name, name);
    }
    return v;
}

static void
check_length(const char *table_name, const char *name, R_xlen_t n, tsk_size_t expected)
{
    if ((tsk_size_t) n != expected) {
        Rf_error("%s table: column '%s' has length %ld but %ld rows expected",
            table_name, name, (long) n, (long) expected);
    }
}

static double *
read_double(SEXP v, const char *table_name, const char *name)
{
    R_xlen_t j, n = XLENGTH(v);
    double *out = (double *) R_alloc((size_t) (n > 0 ? n : 1), sizeof(double));

    if (TYPEOF(v) == REALSXP) {
        if (n > 0) {
            memcpy(out, REAL(v), (size_t) n * sizeof(double));
        }
    } else if (TYPEOF(v) == INTSXP || TYPEOF(v) == LGLSXP) {
        for (j = 0; j < n; j++) {
            out[j] = INTEGER(v)[j] == NA_INTEGER ? NAN : (double) INTEGER(v)[j];
        }
    } else {
        Rf_error("%s table: column '%s' must be numeric", table_name, name);
    }
    return out;
}

static tsk_id_t *
read_id(SEXP v, const char *table_name, const char *name)
{
    R_xlen_t j, n = XLENGTH(v);
    tsk_id_t *out = (tsk_id_t *) R_alloc((size_t) (n > 0 ? n : 1), sizeof(tsk_id_t));
    double x;

    if (TYPEOF(v) == INTSXP || TYPEOF(v) == LGLSXP) {
        for (j = 0; j < n; j++) {
            out[j] = INTEGER(v)[j] == NA_INTEGER ? TSK_NULL : INTEGER(v)[j];
        }
    } else if (TYPEOF(v) == REALSXP) {
        for (j = 0; j < n; j++) {
            x = REAL(v)[j];
            if (ISNAN(x)) {
                out[j] = TSK_NULL;
            } else {
                if (x != floor(x) || x < INT32_MIN || x > INT32_MAX) {
                    Rf_error("%s table: column '%s' must contain 32 bit integer IDs",
                        table_name, name);
                }
                out[j] = (tsk_id_t) x;
            }
        }
    } else {
        Rf_error("%s table: column '%s' must be integer", table_name, name);
    }
    return out;
}

static tsk_flags_t *
read_flags(SEXP v, const char *table_name, const char *name)
{
    R_xlen_t j, n = XLENGTH(v);
    tsk_flags_t *out
        = (tsk_flags_t *) R_alloc((size_t) (n > 0 ? n : 1), sizeof(tsk_flags_t));
    double x;

    for (j = 0; j < n; j++) {
        if (TYPEOF(v) == REALSXP) {
            x = REAL(v)[j];
        } else if (TYPEOF(v) == INTSXP || TYPEOF(v) == LGLSXP) {
            /* Negative integers are the bit patterns of large unsigned values. */
            x = INTEGER(v)[j] == NA_INTEGER ? NAN : (double) INTEGER(v)[j];
            if (x < 0) {
                x += 4294967296.0;
            }
        } else {
            Rf_error("%s table: column '%s' must be numeric", table_name, name);
        }
        if (ISNAN(x) || x < 0 || x >= 4294967296.0 || x != floor(x)) {
            Rf_error("%s table: column '%s' must contain unsigned 32 bit integers",
                table_name, name);
        }
        out[j] = (tsk_flags_t) x;
    }
    return out;
}

/* Reads a ragged column and its offset. Returns false if the column is absent. */
static bool
read_ragged(SEXP table, const char *table_name, const char *name, int type,
    tsk_size_t num_rows, void **data, tsk_size_t **offset)
{
    char offset_name[256];
    SEXP v, off;
    R_xlen_t j, n;
    double x;
    tsk_size_t *o;

    snprintf(offset_name, sizeof(offset_name), "%s_offset", name);
    v = msprimer_list_get(table, name);
    off = msprimer_list_get(table, offset_name);
    if (v == R_NilValue && off == R_NilValue) {
        *data = NULL;
        *offset = NULL;
        return false;
    }
    if (v == R_NilValue || off == R_NilValue) {
        Rf_error("%s table: %s and %s must be specified together", table_name, name,
            offset_name);
    }
    if (XLENGTH(off) != (R_xlen_t) num_rows + 1) {
        Rf_error("%s table: %s must have num_rows + 1 entries", table_name, offset_name);
    }
    n = XLENGTH(off);
    o = (tsk_size_t *) R_alloc((size_t) n, sizeof(tsk_size_t));
    for (j = 0; j < n; j++) {
        if (TYPEOF(off) == REALSXP) {
            x = REAL(off)[j];
        } else if (TYPEOF(off) == INTSXP) {
            x = (double) INTEGER(off)[j];
        } else {
            Rf_error("%s table: %s must be numeric", table_name, offset_name);
        }
        if (ISNAN(x) || x < 0 || x != floor(x)) {
            Rf_error("%s table: bad offset value in %s", table_name, offset_name);
        }
        o[j] = (tsk_size_t) x;
    }
    if (o[0] != 0) {
        Rf_error("%s table: %s must start with 0", table_name, offset_name);
    }
    for (j = 1; j < n; j++) {
        if (o[j] < o[j - 1]) {
            Rf_error("%s table: %s must be non-decreasing", table_name, offset_name);
        }
    }
    if (o[n - 1] != (tsk_size_t) XLENGTH(v)) {
        Rf_error("%s table: the last entry of %s must equal the length of %s",
            table_name, offset_name, name);
    }
    if (type == RAWSXP) {
        if (TYPEOF(v) != RAWSXP) {
            Rf_error("%s table: column '%s' must be a raw vector", table_name, name);
        }
        *data = (void *) RAW(v);
    } else if (type == REALSXP) {
        *data = (void *) read_double(v, table_name, name);
    } else {
        *data = (void *) read_id(v, table_name, name);
    }
    *offset = o;
    return true;
}

static const char *
read_string(SEXP v, const char *what, tsk_size_t *length)
{
    const char *s;
    if (TYPEOF(v) == RAWSXP) {
        *length = (tsk_size_t) XLENGTH(v);
        return (const char *) RAW(v);
    }
    if (TYPEOF(v) != STRSXP || XLENGTH(v) != 1 || STRING_ELT(v, 0) == NA_STRING) {
        Rf_error("'%s' must be a single string", what);
    }
    s = Rf_translateCharUTF8(STRING_ELT(v, 0));
    *length = (tsk_size_t) strlen(s);
    return s;
}

static void
set_schema(SEXP table, const char *table_name,
    int (*setter)(void *, const char *, tsk_size_t), void *t)
{
    SEXP v = msprimer_list_get(table, "metadata_schema");
    const char *s;
    tsk_size_t len;
    int err;

    if (v != R_NilValue) {
        s = read_string(v, "metadata_schema", &len);
        err = setter(t, s, len);
        if (err != 0) {
            msprimer_tsk_error(err, table_name);
        }
    }
}

#define SCHEMA_SETTER(type)                                                             \
    static int set_##type##_schema(void *t, const char *s, tsk_size_t len)              \
    {                                                                                   \
        return tsk_##type##_table_set_metadata_schema((tsk_##type##_table_t *) t, s, len); \
    }

SCHEMA_SETTER(individual)
SCHEMA_SETTER(node)
SCHEMA_SETTER(edge)
SCHEMA_SETTER(migration)
SCHEMA_SETTER(site)
SCHEMA_SETTER(mutation)
SCHEMA_SETTER(population)

static tsk_size_t
num_rows_of(SEXP table, const char *table_name, const char *col)
{
    SEXP v = get_column(table, table_name, col, true);
    return (tsk_size_t) XLENGTH(v);
}

static void
parse_individuals(tsk_individual_table_t *t, SEXP d)
{
    const char *tn = "individuals";
    tsk_size_t n = num_rows_of(d, tn, "flags");
    tsk_flags_t *flags = read_flags(get_column(d, tn, "flags", true), tn, "flags");
    void *location, *parents, *metadata;
    tsk_size_t *location_offset, *parents_offset, *metadata_offset;
    int err;

    read_ragged(d, tn, "location", REALSXP, n, &location, &location_offset);
    read_ragged(d, tn, "parents", INTSXP, n, &parents, &parents_offset);
    read_ragged(d, tn, "metadata", RAWSXP, n, &metadata, &metadata_offset);
    err = tsk_individual_table_clear(t);
    if (err == 0) {
        err = tsk_individual_table_append_columns(t, n, flags, (double *) location,
            location_offset, (tsk_id_t *) parents, parents_offset, (char *) metadata,
            metadata_offset);
    }
    if (err != 0) {
        msprimer_tsk_error(err, tn);
    }
    set_schema(d, tn, set_individual_schema, t);
}

static void
parse_nodes(tsk_node_table_t *t, SEXP d)
{
    const char *tn = "nodes";
    tsk_size_t n = num_rows_of(d, tn, "flags");
    SEXP v;
    tsk_flags_t *flags = read_flags(get_column(d, tn, "flags", true), tn, "flags");
    double *time;
    tsk_id_t *population = NULL, *individual = NULL;
    void *metadata;
    tsk_size_t *metadata_offset;
    int err;

    v = get_column(d, tn, "time", true);
    check_length(tn, "time", XLENGTH(v), n);
    time = read_double(v, tn, "time");
    v = get_column(d, tn, "population", false);
    if (v != R_NilValue) {
        check_length(tn, "population", XLENGTH(v), n);
        population = read_id(v, tn, "population");
    }
    v = get_column(d, tn, "individual", false);
    if (v != R_NilValue) {
        check_length(tn, "individual", XLENGTH(v), n);
        individual = read_id(v, tn, "individual");
    }
    read_ragged(d, tn, "metadata", RAWSXP, n, &metadata, &metadata_offset);
    err = tsk_node_table_clear(t);
    if (err == 0) {
        err = tsk_node_table_append_columns(t, n, flags, time, population, individual,
            (char *) metadata, metadata_offset);
    }
    if (err != 0) {
        msprimer_tsk_error(err, tn);
    }
    set_schema(d, tn, set_node_schema, t);
}

static void
parse_edges(tsk_edge_table_t *t, SEXP d)
{
    const char *tn = "edges";
    tsk_size_t n = num_rows_of(d, tn, "left");
    const char *cols[] = { "right", "parent", "child" };
    double *left, *right;
    tsk_id_t *parent, *child;
    void *metadata;
    tsk_size_t *metadata_offset;
    int j, err;

    for (j = 0; j < 3; j++) {
        check_length(tn, cols[j], XLENGTH(get_column(d, tn, cols[j], true)), n);
    }
    left = read_double(get_column(d, tn, "left", true), tn, "left");
    right = read_double(get_column(d, tn, "right", true), tn, "right");
    parent = read_id(get_column(d, tn, "parent", true), tn, "parent");
    child = read_id(get_column(d, tn, "child", true), tn, "child");
    read_ragged(d, tn, "metadata", RAWSXP, n, &metadata, &metadata_offset);
    err = tsk_edge_table_clear(t);
    if (err == 0) {
        err = tsk_edge_table_append_columns(
            t, n, left, right, parent, child, (char *) metadata, metadata_offset);
    }
    if (err != 0) {
        msprimer_tsk_error(err, tn);
    }
    set_schema(d, tn, set_edge_schema, t);
}

static void
parse_migrations(tsk_migration_table_t *t, SEXP d)
{
    const char *tn = "migrations";
    tsk_size_t n = num_rows_of(d, tn, "left");
    const char *cols[] = { "right", "node", "source", "dest", "time" };
    double *left, *right, *time;
    tsk_id_t *node, *source, *dest;
    void *metadata;
    tsk_size_t *metadata_offset;
    int j, err;

    for (j = 0; j < 5; j++) {
        check_length(tn, cols[j], XLENGTH(get_column(d, tn, cols[j], true)), n);
    }
    left = read_double(get_column(d, tn, "left", true), tn, "left");
    right = read_double(get_column(d, tn, "right", true), tn, "right");
    node = read_id(get_column(d, tn, "node", true), tn, "node");
    source = read_id(get_column(d, tn, "source", true), tn, "source");
    dest = read_id(get_column(d, tn, "dest", true), tn, "dest");
    time = read_double(get_column(d, tn, "time", true), tn, "time");
    read_ragged(d, tn, "metadata", RAWSXP, n, &metadata, &metadata_offset);
    err = tsk_migration_table_clear(t);
    if (err == 0) {
        err = tsk_migration_table_append_columns(t, n, left, right, node, source, dest,
            time, (char *) metadata, metadata_offset);
    }
    if (err != 0) {
        msprimer_tsk_error(err, tn);
    }
    set_schema(d, tn, set_migration_schema, t);
}

static void
parse_sites(tsk_site_table_t *t, SEXP d)
{
    const char *tn = "sites";
    tsk_size_t n = num_rows_of(d, tn, "position");
    double *position = read_double(get_column(d, tn, "position", true), tn, "position");
    void *ancestral_state, *metadata;
    tsk_size_t *ancestral_state_offset, *metadata_offset;
    int err;

    if (!read_ragged(d, tn, "ancestral_state", RAWSXP, n, &ancestral_state,
            &ancestral_state_offset)) {
        Rf_error("sites table: required column 'ancestral_state' not specified");
    }
    read_ragged(d, tn, "metadata", RAWSXP, n, &metadata, &metadata_offset);
    err = tsk_site_table_clear(t);
    if (err == 0) {
        err = tsk_site_table_append_columns(t, n, position, (char *) ancestral_state,
            ancestral_state_offset, (char *) metadata, metadata_offset);
    }
    if (err != 0) {
        msprimer_tsk_error(err, tn);
    }
    set_schema(d, tn, set_site_schema, t);
}

static void
parse_mutations(tsk_mutation_table_t *t, SEXP d)
{
    const char *tn = "mutations";
    tsk_size_t n = num_rows_of(d, tn, "site");
    SEXP v;
    tsk_id_t *site, *node, *parent = NULL;
    double *time = NULL;
    void *derived_state, *metadata;
    tsk_size_t *derived_state_offset, *metadata_offset;
    int err;

    check_length(tn, "node", XLENGTH(get_column(d, tn, "node", true)), n);
    site = read_id(get_column(d, tn, "site", true), tn, "site");
    node = read_id(get_column(d, tn, "node", true), tn, "node");
    v = get_column(d, tn, "parent", false);
    if (v != R_NilValue) {
        check_length(tn, "parent", XLENGTH(v), n);
        parent = read_id(v, tn, "parent");
    }
    v = get_column(d, tn, "time", false);
    if (v != R_NilValue) {
        check_length(tn, "time", XLENGTH(v), n);
        time = read_double(v, tn, "time");
    }
    if (!read_ragged(
            d, tn, "derived_state", RAWSXP, n, &derived_state, &derived_state_offset)) {
        Rf_error("mutations table: required column 'derived_state' not specified");
    }
    read_ragged(d, tn, "metadata", RAWSXP, n, &metadata, &metadata_offset);
    err = tsk_mutation_table_clear(t);
    if (err == 0) {
        err = tsk_mutation_table_append_columns(t, n, site, node, parent, time,
            (char *) derived_state, derived_state_offset, (char *) metadata,
            metadata_offset);
    }
    if (err != 0) {
        msprimer_tsk_error(err, tn);
    }
    set_schema(d, tn, set_mutation_schema, t);
}

static void
parse_populations(tsk_population_table_t *t, SEXP d)
{
    const char *tn = "populations";
    SEXP off = msprimer_list_get(d, "metadata_offset");
    tsk_size_t n = off == R_NilValue ? 0 : (tsk_size_t) XLENGTH(off) - 1;
    void *metadata;
    tsk_size_t *metadata_offset;
    int err;

    if (off != R_NilValue && XLENGTH(off) == 0) {
        Rf_error("populations table: metadata_offset must have num_rows + 1 entries");
    }
    read_ragged(d, tn, "metadata", RAWSXP, n, &metadata, &metadata_offset);
    err = tsk_population_table_clear(t);
    if (err == 0 && n > 0) {
        err = tsk_population_table_append_columns(
            t, n, (char *) metadata, metadata_offset);
    }
    if (err != 0) {
        msprimer_tsk_error(err, tn);
    }
    set_schema(d, tn, set_population_schema, t);
}

static void
parse_provenances(tsk_provenance_table_t *t, SEXP d)
{
    const char *tn = "provenances";
    SEXP off = msprimer_list_get(d, "timestamp_offset");
    tsk_size_t n = off == R_NilValue ? 0 : (tsk_size_t) XLENGTH(off) - 1;
    void *timestamp, *record;
    tsk_size_t *timestamp_offset, *record_offset;
    int err;

    if (off != R_NilValue && XLENGTH(off) == 0) {
        Rf_error("provenances table: timestamp_offset must have num_rows + 1 entries");
    }
    read_ragged(d, tn, "timestamp", RAWSXP, n, &timestamp, &timestamp_offset);
    read_ragged(d, tn, "record", RAWSXP, n, &record, &record_offset);
    if ((timestamp == NULL) != (record == NULL)) {
        Rf_error("provenances table: timestamp and record must be specified together");
    }
    err = tsk_provenance_table_clear(t);
    if (err == 0 && n > 0) {
        err = tsk_provenance_table_append_columns(t, n, (char *) timestamp,
            timestamp_offset, (char *) record, record_offset);
    }
    if (err != 0) {
        msprimer_tsk_error(err, tn);
    }
}

static SEXP
required_table(SEXP d, const char *name)
{
    SEXP v = msprimer_list_get(d, name);
    if (v == R_NilValue || TYPEOF(v) != VECSXP) {
        Rf_error("table collection: '%s' table must be a named list", name);
    }
    return v;
}

static void
parse_reference_sequence(tsk_reference_sequence_t *ref, SEXP d)
{
    SEXP v;
    const char *s;
    tsk_size_t len;
    int err = 0;

    v = msprimer_list_get(d, "metadata_schema");
    if (v != R_NilValue) {
        s = read_string(v, "reference_sequence metadata_schema", &len);
        err = tsk_reference_sequence_set_metadata_schema(ref, s, len);
    }
    v = msprimer_list_get(d, "metadata");
    if (err == 0 && v != R_NilValue) {
        s = read_string(v, "reference_sequence metadata", &len);
        err = tsk_reference_sequence_set_metadata(ref, s, len);
    }
    v = msprimer_list_get(d, "data");
    if (err == 0 && v != R_NilValue) {
        s = read_string(v, "reference_sequence data", &len);
        err = tsk_reference_sequence_set_data(ref, s, len);
    }
    v = msprimer_list_get(d, "url");
    if (err == 0 && v != R_NilValue) {
        s = read_string(v, "reference_sequence url", &len);
        err = tsk_reference_sequence_set_url(ref, s, len);
    }
    if (err != 0) {
        msprimer_tsk_error(err, "reference_sequence");
    }
}

void
msprimer_tables_from_list(tsk_table_collection_t *tables, SEXP d)
{
    SEXP v, ins, rem;
    const char *s;
    tsk_size_t len;
    tsk_id_t *ins_ids, *rem_ids;
    int err;

    if (TYPEOF(d) != VECSXP) {
        Rf_error("table collection must be a named list");
    }
    v = msprimer_list_get(d, "sequence_length");
    if (v == R_NilValue || !Rf_isNumeric(v) || XLENGTH(v) != 1) {
        Rf_error("table collection: 'sequence_length' must be a number");
    }
    tables->sequence_length = Rf_asReal(v);

    v = msprimer_list_get(d, "metadata_schema");
    if (v != R_NilValue) {
        s = read_string(v, "metadata_schema", &len);
        err = tsk_table_collection_set_metadata_schema(tables, s, len);
        if (err != 0) {
            msprimer_tsk_error(err, "metadata_schema");
        }
    }
    v = msprimer_list_get(d, "metadata");
    if (v != R_NilValue) {
        s = read_string(v, "metadata", &len);
        err = tsk_table_collection_set_metadata(tables, s, len);
        if (err != 0) {
            msprimer_tsk_error(err, "metadata");
        }
    }
    v = msprimer_list_get(d, "time_units");
    if (v != R_NilValue) {
        s = read_string(v, "time_units", &len);
        err = tsk_table_collection_set_time_units(tables, s, len);
        if (err != 0) {
            msprimer_tsk_error(err, "time_units");
        }
    }
    parse_individuals(&tables->individuals, required_table(d, "individuals"));
    parse_nodes(&tables->nodes, required_table(d, "nodes"));
    parse_edges(&tables->edges, required_table(d, "edges"));
    parse_migrations(&tables->migrations, required_table(d, "migrations"));
    parse_sites(&tables->sites, required_table(d, "sites"));
    parse_mutations(&tables->mutations, required_table(d, "mutations"));
    parse_populations(&tables->populations, required_table(d, "populations"));
    parse_provenances(&tables->provenances, required_table(d, "provenances"));

    v = msprimer_list_get(d, "indexes");
    if (v != R_NilValue) {
        ins = msprimer_list_get(v, "edge_insertion_order");
        rem = msprimer_list_get(v, "edge_removal_order");
        if ((ins == R_NilValue) != (rem == R_NilValue)) {
            Rf_error("edge_insertion_order and edge_removal_order must be specified "
                     "together");
        }
        if (ins != R_NilValue) {
            if (XLENGTH(ins) != XLENGTH(rem)
                || (tsk_size_t) XLENGTH(ins) != tables->edges.num_rows) {
                Rf_error("edge_insertion_order and edge_removal_order must be the "
                         "same length as the number of edges");
            }
            ins_ids = read_id(ins, "indexes", "edge_insertion_order");
            rem_ids = read_id(rem, "indexes", "edge_removal_order");
            err = tsk_table_collection_set_indexes(tables, ins_ids, rem_ids);
            if (err != 0) {
                msprimer_tsk_error(err, "indexes");
            }
        }
    }
    v = msprimer_list_get(d, "reference_sequence");
    if (v != R_NilValue) {
        parse_reference_sequence(&tables->reference_sequence, v);
    }
}

/*************************
 * Writers
 *************************/

static SEXP
make_doubles(const double *x, tsk_size_t n)
{
    SEXP ret = Rf_allocVector(REALSXP, (R_xlen_t) n);
    if (n > 0) {
        memcpy(REAL(ret), x, n * sizeof(double));
    }
    return ret;
}

static SEXP
make_ids(const tsk_id_t *x, tsk_size_t n)
{
    SEXP ret = Rf_allocVector(INTSXP, (R_xlen_t) n);
    if (n > 0) {
        memcpy(INTEGER(ret), x, n * sizeof(tsk_id_t));
    }
    return ret;
}

static SEXP
make_flags(const tsk_flags_t *x, tsk_size_t n)
{
    SEXP ret = Rf_allocVector(REALSXP, (R_xlen_t) n);
    tsk_size_t j;
    for (j = 0; j < n; j++) {
        REAL(ret)[j] = (double) x[j];
    }
    return ret;
}

static SEXP
make_offsets(const tsk_size_t *x, tsk_size_t n)
{
    SEXP ret = Rf_allocVector(REALSXP, (R_xlen_t) (n + 1));
    tsk_size_t j;
    for (j = 0; j <= n; j++) {
        REAL(ret)[j] = (double) x[j];
    }
    return ret;
}

static SEXP
make_raw(const char *x, tsk_size_t n)
{
    SEXP ret = Rf_allocVector(RAWSXP, (R_xlen_t) n);
    if (n > 0) {
        memcpy(RAW(ret), x, n);
    }
    return ret;
}

static SEXP
make_string(const char *x, tsk_size_t n)
{
    SEXP ret = PROTECT(Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(ret, 0, Rf_mkCharLenCE(x == NULL ? "" : x, (int) n, CE_UTF8));
    UNPROTECT(1);
    return ret;
}

/* Fills the list element-by-element; elements are protected by the list. */
#define SET(list, i, value) SET_VECTOR_ELT(list, i, value)

static SEXP
write_individuals(const tsk_individual_table_t *t)
{
    const char *names[] = { "flags", "location", "location_offset", "parents",
        "parents_offset", "metadata", "metadata_offset", "metadata_schema" };
    SEXP ret = PROTECT(msprimer_make_named_list(8, names));
    tsk_size_t n = t->num_rows;
    SET(ret, 0, make_flags(t->flags, n));
    SET(ret, 1, make_doubles(t->location, t->location_length));
    SET(ret, 2, make_offsets(t->location_offset, n));
    SET(ret, 3, make_ids(t->parents, t->parents_length));
    SET(ret, 4, make_offsets(t->parents_offset, n));
    SET(ret, 5, make_raw(t->metadata, t->metadata_length));
    SET(ret, 6, make_offsets(t->metadata_offset, n));
    SET(ret, 7, make_string(t->metadata_schema, t->metadata_schema_length));
    UNPROTECT(1);
    return ret;
}

static SEXP
write_nodes(const tsk_node_table_t *t)
{
    const char *names[] = { "flags", "time", "population", "individual", "metadata",
        "metadata_offset", "metadata_schema" };
    SEXP ret = PROTECT(msprimer_make_named_list(7, names));
    tsk_size_t n = t->num_rows;
    SET(ret, 0, make_flags(t->flags, n));
    SET(ret, 1, make_doubles(t->time, n));
    SET(ret, 2, make_ids(t->population, n));
    SET(ret, 3, make_ids(t->individual, n));
    SET(ret, 4, make_raw(t->metadata, t->metadata_length));
    SET(ret, 5, make_offsets(t->metadata_offset, n));
    SET(ret, 6, make_string(t->metadata_schema, t->metadata_schema_length));
    UNPROTECT(1);
    return ret;
}

static SEXP
write_edges(const tsk_edge_table_t *t)
{
    const char *names[] = { "left", "right", "parent", "child", "metadata",
        "metadata_offset", "metadata_schema" };
    SEXP ret = PROTECT(msprimer_make_named_list(7, names));
    tsk_size_t n = t->num_rows;
    SET(ret, 0, make_doubles(t->left, n));
    SET(ret, 1, make_doubles(t->right, n));
    SET(ret, 2, make_ids(t->parent, n));
    SET(ret, 3, make_ids(t->child, n));
    SET(ret, 4, make_raw(t->metadata, t->metadata_length));
    SET(ret, 5, make_offsets(t->metadata_offset, n));
    SET(ret, 6, make_string(t->metadata_schema, t->metadata_schema_length));
    UNPROTECT(1);
    return ret;
}

static SEXP
write_migrations(const tsk_migration_table_t *t)
{
    const char *names[] = { "left", "right", "node", "source", "dest", "time",
        "metadata", "metadata_offset", "metadata_schema" };
    SEXP ret = PROTECT(msprimer_make_named_list(9, names));
    tsk_size_t n = t->num_rows;
    SET(ret, 0, make_doubles(t->left, n));
    SET(ret, 1, make_doubles(t->right, n));
    SET(ret, 2, make_ids(t->node, n));
    SET(ret, 3, make_ids(t->source, n));
    SET(ret, 4, make_ids(t->dest, n));
    SET(ret, 5, make_doubles(t->time, n));
    SET(ret, 6, make_raw(t->metadata, t->metadata_length));
    SET(ret, 7, make_offsets(t->metadata_offset, n));
    SET(ret, 8, make_string(t->metadata_schema, t->metadata_schema_length));
    UNPROTECT(1);
    return ret;
}

static SEXP
write_sites(const tsk_site_table_t *t)
{
    const char *names[] = { "position", "ancestral_state", "ancestral_state_offset",
        "metadata", "metadata_offset", "metadata_schema" };
    SEXP ret = PROTECT(msprimer_make_named_list(6, names));
    tsk_size_t n = t->num_rows;
    SET(ret, 0, make_doubles(t->position, n));
    SET(ret, 1, make_raw(t->ancestral_state, t->ancestral_state_length));
    SET(ret, 2, make_offsets(t->ancestral_state_offset, n));
    SET(ret, 3, make_raw(t->metadata, t->metadata_length));
    SET(ret, 4, make_offsets(t->metadata_offset, n));
    SET(ret, 5, make_string(t->metadata_schema, t->metadata_schema_length));
    UNPROTECT(1);
    return ret;
}

static SEXP
write_mutations(const tsk_mutation_table_t *t)
{
    const char *names[] = { "site", "node", "time", "derived_state",
        "derived_state_offset", "parent", "metadata", "metadata_offset",
        "metadata_schema" };
    SEXP ret = PROTECT(msprimer_make_named_list(9, names));
    tsk_size_t n = t->num_rows;
    SET(ret, 0, make_ids(t->site, n));
    SET(ret, 1, make_ids(t->node, n));
    SET(ret, 2, make_doubles(t->time, n));
    SET(ret, 3, make_raw(t->derived_state, t->derived_state_length));
    SET(ret, 4, make_offsets(t->derived_state_offset, n));
    SET(ret, 5, make_ids(t->parent, n));
    SET(ret, 6, make_raw(t->metadata, t->metadata_length));
    SET(ret, 7, make_offsets(t->metadata_offset, n));
    SET(ret, 8, make_string(t->metadata_schema, t->metadata_schema_length));
    UNPROTECT(1);
    return ret;
}

static SEXP
write_populations(const tsk_population_table_t *t)
{
    const char *names[] = { "metadata", "metadata_offset", "metadata_schema" };
    SEXP ret = PROTECT(msprimer_make_named_list(3, names));
    SET(ret, 0, make_raw(t->metadata, t->metadata_length));
    SET(ret, 1, make_offsets(t->metadata_offset, t->num_rows));
    SET(ret, 2, make_string(t->metadata_schema, t->metadata_schema_length));
    UNPROTECT(1);
    return ret;
}

static SEXP
write_provenances(const tsk_provenance_table_t *t)
{
    const char *names[] = { "timestamp", "timestamp_offset", "record", "record_offset" };
    SEXP ret = PROTECT(msprimer_make_named_list(4, names));
    SET(ret, 0, make_raw(t->timestamp, t->timestamp_length));
    SET(ret, 1, make_offsets(t->timestamp_offset, t->num_rows));
    SET(ret, 2, make_raw(t->record, t->record_length));
    SET(ret, 3, make_offsets(t->record_offset, t->num_rows));
    UNPROTECT(1);
    return ret;
}

SEXP
msprimer_tables_to_list(tsk_table_collection_t *tables)
{
    const char *names[] = { "encoding_version", "sequence_length", "metadata_schema",
        "metadata", "time_units", "individuals", "nodes", "edges", "migrations", "sites",
        "mutations", "populations", "provenances", "indexes", "reference_sequence" };
    const char *index_names[] = { "edge_insertion_order", "edge_removal_order" };
    const char *ref_names[] = { "metadata_schema", "metadata", "data", "url" };
    int num_items = 13;
    bool has_index = tsk_table_collection_has_index(tables, 0);
    bool has_ref = tsk_table_collection_has_reference_sequence(tables);
    SEXP ret, ev, idx, ref, nms;
    tsk_reference_sequence_t *rs = &tables->reference_sequence;

    ret = PROTECT(msprimer_make_named_list(15, names));
    ev = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(ev)[0] = 1;
    INTEGER(ev)[1] = 6;
    SET(ret, 0, ev);
    SET(ret, 1, Rf_ScalarReal(tables->sequence_length));
    SET(ret, 2, make_string(tables->metadata_schema, tables->metadata_schema_length));
    SET(ret, 3, make_raw(tables->metadata, tables->metadata_length));
    SET(ret, 4, make_string(tables->time_units, tables->time_units_length));
    SET(ret, 5, write_individuals(&tables->individuals));
    SET(ret, 6, write_nodes(&tables->nodes));
    SET(ret, 7, write_edges(&tables->edges));
    SET(ret, 8, write_migrations(&tables->migrations));
    SET(ret, 9, write_sites(&tables->sites));
    SET(ret, 10, write_mutations(&tables->mutations));
    SET(ret, 11, write_populations(&tables->populations));
    SET(ret, 12, write_provenances(&tables->provenances));
    if (has_index) {
        idx = PROTECT(msprimer_make_named_list(2, index_names));
        SET(idx, 0, make_ids(tables->indexes.edge_insertion_order,
                        tables->indexes.num_edges));
        SET(idx, 1, make_ids(tables->indexes.edge_removal_order,
                        tables->indexes.num_edges));
        SET(ret, num_items, idx);
        num_items++;
        UNPROTECT(1);
    }
    if (has_ref) {
        ref = PROTECT(msprimer_make_named_list(4, ref_names));
        SET(ref, 0, make_string(rs->metadata_schema, rs->metadata_schema_length));
        SET(ref, 1, make_raw(rs->metadata, rs->metadata_length));
        SET(ref, 2, make_string(rs->data, rs->data_length));
        SET(ref, 3, make_string(rs->url, rs->url_length));
        SET(ret, num_items, ref);
        if (!has_index) {
            /* Keep names aligned when the index slot is skipped. */
            nms = Rf_getAttrib(ret, R_NamesSymbol);
            SET_STRING_ELT(nms, num_items, Rf_mkChar("reference_sequence"));
        }
        num_items++;
        UNPROTECT(1);
    }
    ret = Rf_lengthgets(ret, num_items);
    UNPROTECT(2);
    return ret;
}

/*************************
 * Table collection external pointers
 *************************/

static void
tc_finalizer(SEXP xp)
{
    tsk_table_collection_t *tc = (tsk_table_collection_t *) R_ExternalPtrAddr(xp);
    if (tc != NULL) {
        tsk_table_collection_free(tc);
        free(tc);
        R_ClearExternalPtr(xp);
    }
}

/* Frees the table collection now rather than at garbage collection. */
void
msprimer_tc_release(SEXP xp)
{
    tc_finalizer(xp);
}

/* Allocates an empty table collection owned by the returned (unprotected)
 * external pointer, so that any subsequent R error does not leak memory. */
SEXP
msprimer_tc_alloc_xp(tsk_table_collection_t **out)
{
    tsk_table_collection_t *tc = calloc(1, sizeof(*tc));
    SEXP xp;
    int err;

    if (tc == NULL) {
        Rf_error("Out of memory");
    }
    xp = PROTECT(R_MakeExternalPtr(tc, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(xp, tc_finalizer, TRUE);
    err = tsk_table_collection_init(tc, 0);
    if (err != 0) {
        msprimer_tsk_error(err, NULL);
    }
    *out = tc;
    UNPROTECT(1);
    return xp;
}

static void
build_index_if_needed(tsk_table_collection_t *tc)
{
    int err;
    if (!tsk_table_collection_has_index(tc, 0)) {
        err = tsk_table_collection_build_index(tc, 0);
        if (err != 0) {
            msprimer_tsk_error(err, "building table indexes");
        }
    }
}

SEXP
msprimer_tables_roundtrip(SEXP d)
{
    tsk_table_collection_t *tc;
    SEXP ret, xp;

    xp = PROTECT(msprimer_tc_alloc_xp(&tc));

    msprimer_tables_from_list(tc, d);
    ret = PROTECT(msprimer_tables_to_list(tc));
    msprimer_tc_release(xp);
    UNPROTECT(2);
    return ret;
}

/* Splits a ragged raw column into a character vector using its offsets. */
SEXP
msprimer_unpack_strings(SEXP data, SEXP offset)
{
    R_xlen_t j, n = XLENGTH(offset) - 1;
    double start, end;
    SEXP ret;

    if (TYPEOF(data) != RAWSXP || TYPEOF(offset) != REALSXP || n < 0) {
        Rf_error("unpack_strings requires a raw vector and double offsets");
    }
    ret = PROTECT(Rf_allocVector(STRSXP, n));
    for (j = 0; j < n; j++) {
        start = REAL(offset)[j];
        end = REAL(offset)[j + 1];
        if (start < 0 || end < start || end > (double) XLENGTH(data)) {
            Rf_error("bad offsets");
        }
        SET_STRING_ELT(ret, j,
            Rf_mkCharLenCE((const char *) RAW(data) + (R_xlen_t) start,
                (int) (end - start), CE_UTF8));
    }
    UNPROTECT(1);
    return ret;
}

static const char *
path_arg(SEXP path)
{
    if (TYPEOF(path) != STRSXP || XLENGTH(path) != 1 || STRING_ELT(path, 0) == NA_STRING) {
        Rf_error("file must be a single path");
    }
    return R_ExpandFileName(Rf_translateChar(STRING_ELT(path, 0)));
}

SEXP
msprimer_tables_dump(SEXP d, SEXP path)
{
    tsk_table_collection_t *tc;
    const char *filename;
    SEXP xp;
    int err;

    xp = PROTECT(msprimer_tc_alloc_xp(&tc));
    filename = path_arg(path);
    msprimer_tables_from_list(tc, d);
    build_index_if_needed(tc);
    err = tsk_table_collection_dump(tc, filename, 0);
    if (err != 0) {
        msprimer_tsk_error(err, "writing tree sequence file");
    }
    msprimer_tc_release(xp);
    UNPROTECT(1);
    return R_NilValue;
}

SEXP
msprimer_tables_load(SEXP path)
{
    tsk_table_collection_t *tc;
    const char *filename;
    SEXP ret, xp;
    int err;

    xp = PROTECT(msprimer_tc_alloc_xp(&tc));
    filename = path_arg(path);
    /* load initialises the collection itself. */
    tsk_table_collection_free(tc);
    memset(tc, 0, sizeof(*tc));
    err = tsk_table_collection_load(tc, filename, 0);
    if (err != 0) {
        msprimer_tsk_error(err, "reading tree sequence file");
    }
    ret = PROTECT(msprimer_tables_to_list(tc));
    msprimer_tc_release(xp);
    UNPROTECT(2);
    return ret;
}

/* Equivalent of tskit's TableCollection.tree_sequence(): builds the indexes
 * if needed, validates that the tables define a tree sequence and returns
 * list(tables = <indexed tables>, summary = <basic properties>). */
SEXP
msprimer_ts_summary(SEXP d)
{
    const char *names[] = { "num_trees", "num_samples", "num_nodes", "num_edges",
        "num_sites", "num_mutations", "num_individuals", "num_populations",
        "num_migrations", "num_provenances", "sequence_length", "discrete_genome",
        "discrete_time", "min_time", "max_time", "time_units" };
    const char *out_names[] = { "tables", "summary" };
    tsk_table_collection_t *tc;
    tsk_treeseq_t ts;
    SEXP ret, out, xp;
    int err;

    xp = PROTECT(msprimer_tc_alloc_xp(&tc));
    msprimer_tables_from_list(tc, d);
    build_index_if_needed(tc);
    memset(&ts, 0, sizeof(ts));
    err = tsk_treeseq_init(&ts, tc, 0);
    if (err != 0) {
        tsk_treeseq_free(&ts);
        msprimer_tsk_error(err, "invalid tree sequence");
    }
    ret = PROTECT(msprimer_make_named_list(16, names));
    SET(ret, 0, Rf_ScalarReal((double) tsk_treeseq_get_num_trees(&ts)));
    SET(ret, 1, Rf_ScalarReal((double) tsk_treeseq_get_num_samples(&ts)));
    SET(ret, 2, Rf_ScalarReal((double) tsk_treeseq_get_num_nodes(&ts)));
    SET(ret, 3, Rf_ScalarReal((double) tsk_treeseq_get_num_edges(&ts)));
    SET(ret, 4, Rf_ScalarReal((double) tsk_treeseq_get_num_sites(&ts)));
    SET(ret, 5, Rf_ScalarReal((double) tsk_treeseq_get_num_mutations(&ts)));
    SET(ret, 6, Rf_ScalarReal((double) tsk_treeseq_get_num_individuals(&ts)));
    SET(ret, 7, Rf_ScalarReal((double) tsk_treeseq_get_num_populations(&ts)));
    SET(ret, 8, Rf_ScalarReal((double) tsk_treeseq_get_num_migrations(&ts)));
    SET(ret, 9, Rf_ScalarReal((double) tsk_treeseq_get_num_provenances(&ts)));
    SET(ret, 10, Rf_ScalarReal(tsk_treeseq_get_sequence_length(&ts)));
    SET(ret, 11, Rf_ScalarLogical(tsk_treeseq_get_discrete_genome(&ts)));
    SET(ret, 12, Rf_ScalarLogical(tsk_treeseq_get_discrete_time(&ts)));
    SET(ret, 13, Rf_ScalarReal(tsk_treeseq_get_min_time(&ts)));
    SET(ret, 14, Rf_ScalarReal(tsk_treeseq_get_max_time(&ts)));
    SET(ret, 15, make_string(tc->time_units, tc->time_units_length));
    tsk_treeseq_free(&ts);
    out = PROTECT(msprimer_make_named_list(2, out_names));
    SET(out, 0, msprimer_tables_to_list(tc));
    SET(out, 1, ret);
    msprimer_tc_release(xp);
    UNPROTECT(3);
    return out;
}

SEXP
msprimer_tables_sort(SEXP d)
{
    tsk_table_collection_t *tc;
    SEXP ret, xp;
    int err;

    xp = PROTECT(msprimer_tc_alloc_xp(&tc));
    msprimer_tables_from_list(tc, d);
    err = tsk_table_collection_sort(tc, NULL, 0);
    if (err != 0) {
        msprimer_tsk_error(err, "sorting tables");
    }
    ret = PROTECT(msprimer_tables_to_list(tc));
    msprimer_tc_release(xp);
    UNPROTECT(2);
    return ret;
}

/* Port of tskit's TableCollection.keep_intervals(simplify=False) applied to the
 * complement of the specified intervals (i.e., delete_intervals). */
SEXP
msprimer_tables_delete_intervals(SEXP d, SEXP r_left, SEXP r_right)
{
    tsk_table_collection_t *tc, *copy;
    R_xlen_t num_del = XLENGTH(r_left);
    double *del_left = REAL(r_left), *del_right = REAL(r_right);
    double *keep_left, *keep_right, s, e, L;
    tsk_size_t num_keep = 0, j, k;
    tsk_bool_t *keep_sites, *keep_mutations;
    tsk_id_t *site_map;
    const tsk_edge_table_t *edges;
    const tsk_migration_table_t *migs;
    tsk_edge_t edge;
    tsk_migration_t mig;
    double last;
    int err = 0;
    SEXP ret, xp, xp_copy;

    xp = PROTECT(msprimer_tc_alloc_xp(&tc));
    xp_copy = PROTECT(msprimer_tc_alloc_xp(&copy));
    msprimer_tables_from_list(tc, d);
    L = tc->sequence_length;
    keep_left = (double *) R_alloc((size_t) num_del + 1, sizeof(double));
    keep_right = (double *) R_alloc((size_t) num_del + 1, sizeof(double));
    last = 0;
    for (k = 0; k < (tsk_size_t) num_del; k++) {
        if (del_left[k] < last || del_right[k] <= del_left[k] || del_right[k] > L) {
            Rf_error("Intervals must be non-overlapping, sorted and within the sequence");
        }
        if (del_left[k] > last) {
            keep_left[num_keep] = last;
            keep_right[num_keep] = del_left[k];
            num_keep++;
        }
        last = del_right[k];
    }
    if (last < L) {
        keep_left[num_keep] = last;
        keep_right[num_keep] = L;
        num_keep++;
    }
    err = tsk_table_collection_copy(tc, copy, TSK_NO_INIT);
    if (err != 0) {
        msprimer_tsk_error(err, NULL);
    }
    edges = &copy->edges;
    migs = &copy->migrations;
    tsk_edge_table_clear(&tc->edges);
    tsk_migration_table_clear(&tc->migrations);
    for (k = 0; k < num_keep; k++) {
        s = keep_left[k];
        e = keep_right[k];
        for (j = 0; j < edges->num_rows; j++) {
            tsk_edge_table_get_row(edges, (tsk_id_t) j, &edge);
            if (!(edge.right <= s || edge.left >= e)) {
                err = tsk_edge_table_add_row(&tc->edges, TSK_MAX(s, edge.left),
                    TSK_MIN(e, edge.right), edge.parent, edge.child, edge.metadata,
                    edge.metadata_length);
                if (err < 0) {
                    msprimer_tsk_error(err, NULL);
                }
            }
        }
        for (j = 0; j < migs->num_rows; j++) {
            tsk_migration_table_get_row(migs, (tsk_id_t) j, &mig);
            if (!(mig.right <= s || mig.left >= e)) {
                err = tsk_migration_table_add_row(&tc->migrations, TSK_MAX(s, mig.left),
                    TSK_MIN(e, mig.right), mig.node, mig.source, mig.dest, mig.time,
                    mig.metadata, mig.metadata_length);
                if (err < 0) {
                    msprimer_tsk_error(err, NULL);
                }
            }
        }
    }
    keep_sites = (tsk_bool_t *) R_alloc(tc->sites.num_rows + 1, sizeof(tsk_bool_t));
    site_map = (tsk_id_t *) R_alloc(tc->sites.num_rows + 1, sizeof(tsk_id_t));
    for (j = 0; j < tc->sites.num_rows; j++) {
        keep_sites[j] = false;
        for (k = 0; k < num_keep; k++) {
            if (tc->sites.position[j] >= keep_left[k]
                && tc->sites.position[j] < keep_right[k]) {
                keep_sites[j] = true;
            }
        }
    }
    keep_mutations
        = (tsk_bool_t *) R_alloc(tc->mutations.num_rows + 1, sizeof(tsk_bool_t));
    for (j = 0; j < tc->mutations.num_rows; j++) {
        keep_mutations[j] = keep_sites[tc->mutations.site[j]];
    }
    err = tsk_site_table_keep_rows(&tc->sites, keep_sites, 0, site_map);
    if (err == 0) {
        err = tsk_mutation_table_keep_rows(&tc->mutations, keep_mutations, 0, NULL);
    }
    if (err != 0) {
        msprimer_tsk_error(err, NULL);
    }
    for (j = 0; j < tc->mutations.num_rows; j++) {
        tc->mutations.site[j] = site_map[tc->mutations.site[j]];
    }
    err = tsk_table_collection_sort(tc, NULL, 0);
    if (err != 0) {
        msprimer_tsk_error(err, "sorting tables");
    }
    ret = PROTECT(msprimer_tables_to_list(tc));
    msprimer_tc_release(xp);
    msprimer_tc_release(xp_copy);
    UNPROTECT(3);
    return ret;
}
