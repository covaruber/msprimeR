/* Likelihood computations: the Hudson ARG likelihood (C library) and the
 * infinite sites mutation likelihood (ported from the Python implementation). */
#include "msprimer.h"

#include <math.h>
#include <float.h>

typedef struct {
    tsk_treeseq_t ts;
    tsk_tree_t tree;
    bool ts_init;
    bool tree_init;
} msprimer_ts_holder_t;

static void
holder_finalizer(SEXP xp)
{
    msprimer_ts_holder_t *self = (msprimer_ts_holder_t *) R_ExternalPtrAddr(xp);
    if (self != NULL) {
        if (self->tree_init) {
            tsk_tree_free(&self->tree);
        }
        if (self->ts_init) {
            tsk_treeseq_free(&self->ts);
        }
        free(self);
        R_ClearExternalPtr(xp);
    }
}

/* Frees the tree sequence and its source tables now. */
static void
release_treeseq(SEXP xp)
{
    SEXP tc_xp = R_ExternalPtrProtected(xp);
    holder_finalizer(xp);
    msprimer_tc_release(tc_xp);
}

/* Builds a tree sequence (with indexes) owned by an external pointer that is
 * returned protected; the caller must UNPROTECT 2 (table xp and holder xp). */
static msprimer_ts_holder_t *
make_treeseq(SEXP tables, SEXP *xp_out)
{
    tsk_table_collection_t *tc;
    msprimer_ts_holder_t *self;
    SEXP tc_xp, xp;
    int err;

    tc_xp = PROTECT(msprimer_tc_alloc_xp(&tc));
    msprimer_tables_from_list(tc, tables);
    self = calloc(1, sizeof(*self));
    if (self == NULL) {
        Rf_error("Out of memory");
    }
    xp = PROTECT(R_MakeExternalPtr(self, R_NilValue, tc_xp));
    R_RegisterCFinalizerEx(xp, holder_finalizer, TRUE);
    err = tsk_treeseq_init(&self->ts, tc, TSK_TS_INIT_BUILD_INDEXES);
    self->ts_init = true;
    if (err != 0) {
        msprimer_tsk_error(err, "invalid tree sequence");
    }
    err = tsk_tree_init(&self->tree, &self->ts, 0);
    self->tree_init = true;
    if (err != 0) {
        msprimer_tsk_error(err, NULL);
    }
    *xp_out = xp;
    return self;
}

SEXP
msprimer_log_likelihood_arg(SEXP tables, SEXP r_Ne, SEXP r_recombination_rate)
{
    SEXP xp;
    msprimer_ts_holder_t *h;
    tsk_tree_t *tree;
    tsk_size_t num_nodes, j;
    double Ne = Rf_asReal(r_Ne);
    double rho = Rf_asReal(r_recombination_rate);
    double ret = 0;
    bool has_re = false;
    int err;

    if (ISNAN(rho) || rho < 0) {
        Rf_error("recombination_rate must be >= 0");
    }
    h = make_treeseq(tables, &xp);
    tree = &h->tree;
    num_nodes = tsk_treeseq_get_num_nodes(&h->ts);
    for (err = tsk_tree_first(tree); err == TSK_TREE_OK; err = tsk_tree_next(tree)) {
        for (j = 0; j < num_nodes; j++) {
            if (tree->num_children[j] > 2) {
                Rf_error("ARG likelihood encountered a polytomy. Tree sequences must "
                         "contain binary mergers only for valid likelihood evaluation.");
            }
        }
        if (tree->num_children[tree->virtual_root] > 1
            && tsk_treeseq_get_num_edges(&h->ts) > 1) {
            Rf_error("ARG likelihood encountered a tree with multiple roots. All "
                     "local trees must have a single mrca for valid likelihood "
                     "evaluation.");
        }
    }
    if (err < 0) {
        msprimer_tsk_error(err, NULL);
    }
    for (j = 0; j < num_nodes; j++) {
        if (h->ts.tables->nodes.flags[j] & MSP_NODE_IS_RE_EVENT) {
            has_re = true;
        }
    }
    if (tsk_treeseq_get_num_trees(&h->ts) > 1 && !has_re) {
        Rf_error("ARG likelihood only valid for tree sequences where recombinations "
                 "have been recorded via nodes flagged with `NODE_IS_RE_EVENT`.");
    }
    err = msp_log_likelihood_arg(&h->ts, rho, Ne, &ret);
    if (err != 0) {
        msprimer_msp_error(err, NULL);
    }
    release_treeseq(xp);
    UNPROTECT(2);
    return Rf_ScalarReal(ret);
}

static double
branch_length(const tsk_tree_t *tree, const double *time, tsk_id_t u)
{
    tsk_id_t p = tree->parent[u];
    return p == TSK_NULL ? 0.0 : time[p] - time[u];
}

SEXP
msprimer_log_mutation_likelihood(SEXP tables, SEXP r_mutation_rate)
{
    SEXP xp;
    msprimer_ts_holder_t *h;
    tsk_tree_t *tree;
    const tsk_table_collection_t *tc;
    const double *time;
    double mu = Rf_asReal(r_mutation_rate);
    double total_material = 0, ret, potential_branch_length;
    const tsk_site_t *sites;
    tsk_size_t num_sites, j, num_mutations;
    tsk_id_t child, parent;
    int err;

    if (ISNAN(mu) || mu < 0) {
        Rf_error("mutation_rate must be >= 0");
    }
    h = make_treeseq(tables, &xp);
    tree = &h->tree;
    tc = h->ts.tables;
    time = tc->nodes.time;
    for (j = 0; j < tc->edges.num_rows; j++) {
        total_material += (tc->edges.right[j] - tc->edges.left[j])
                          * (time[tc->edges.parent[j]] - time[tc->edges.child[j]]);
    }
    num_mutations = tc->mutations.num_rows;
    if (mu == 0) {
        ret = num_mutations == 0 ? 0 : R_NegInf;
    } else {
        ret = (double) num_mutations * log(total_material * mu) - total_material * mu;
        for (err = tsk_tree_first(tree); err == TSK_TREE_OK; err = tsk_tree_next(tree)) {
            err = tsk_tree_get_sites(tree, &sites, &num_sites);
            if (err != 0) {
                msprimer_tsk_error(err, NULL);
            }
            for (j = 0; j < num_sites; j++) {
                if (sites[j].mutations_length == 0) {
                    continue;
                }
                child = sites[j].mutations[0].node;
                parent = tree->parent[child];
                potential_branch_length = branch_length(tree, time, child);
                /* Extend upwards through unary nodes, stopping below the root. */
                while (parent != TSK_NULL && tree->parent[parent] != TSK_NULL
                       && tree->num_children[parent] == 1) {
                    child = parent;
                    parent = tree->parent[child];
                    potential_branch_length += branch_length(tree, time, child);
                }
                child = sites[j].mutations[0].node;
                while (tree->num_children[child] == 1) {
                    child = tree->left_child[child];
                    potential_branch_length += branch_length(tree, time, child);
                }
                ret += log(potential_branch_length / total_material);
            }
        }
        if (err < 0) {
            msprimer_tsk_error(err, NULL);
        }
    }
    release_treeseq(xp);
    UNPROTECT(2);
    return Rf_ScalarReal(ret);
}
