"""Generate reference outputs from Python msprime for msprimeR parity tests."""
import json
import os
import sys

import msprime
import numpy as np

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
values = {}


def dump(name, ts):
    ts.dump(os.path.join(out, name + ".trees"))


def complex_demography():
    d = msprime.Demography()
    d.add_population(name="A", initial_size=1000, growth_rate=0.01)
    d.add_population(name="B", initial_size=2000)
    d.add_population(name="C", initial_size=500)
    d.add_population(name="ADMIX", initial_size=800)
    d.add_population(name="ANC", initial_size=5000)
    d.set_symmetric_migration_rate(["A", "B"], 1e-3)
    d.set_migration_rate("C", "B", 5e-4)
    d.add_admixture(50, derived="ADMIX", ancestral=["A", "B"], proportions=[0.3, 0.7])
    d.add_population_parameters_change(100, initial_size=300, population="A", growth_rate=0)
    d.add_mass_migration(150, source="C", dest="B", proportion=0.2)
    d.add_migration_rate_change(200, rate=2e-3, source="A", dest="C")
    d.add_symmetric_migration_rate_change(250, ["A", "B"], 0)
    d.add_population_split(300, derived=["A", "B", "C"], ancestral="ANC")
    d.add_census(400)
    return d


# 1 basic
dump("basic", msprime.sim_ancestry(5, sequence_length=1e5, recombination_rate=1e-8,
                                   population_size=1e4, random_seed=42))
# 2 island model with migrations recorded, continuous genome
d = msprime.Demography.island_model([1000, 1000, 500], 1e-3)
dump("island", msprime.sim_ancestry({0: 3, 1: 2, 2: 2}, demography=d, sequence_length=5e4,
                                    recombination_rate=1e-8, record_migrations=True,
                                    discrete_genome=False, random_seed=7))
# 3 complex demography, samples by name & historical sample sets
d = complex_demography()
dump("complex", msprime.sim_ancestry(
    [msprime.SampleSet(3, "A"), msprime.SampleSet(2, "B", time=10),
     msprime.SampleSet(2, "C", ploidy=1), msprime.SampleSet(2, "ADMIX")],
    demography=d, sequence_length=1e5, recombination_rate=1e-8, random_seed=11))
# 4 models
dump("dtwf", msprime.sim_ancestry(4, population_size=100, model="dtwf",
                                  sequence_length=1e4, recombination_rate=1e-6,
                                  random_seed=3))
dump("multi_model", msprime.sim_ancestry(
    4, population_size=100, sequence_length=1e4, recombination_rate=1e-6,
    model=[msprime.DiscreteTimeWrightFisher(duration=20), msprime.StandardCoalescent()],
    random_seed=4))
dump("smc", msprime.sim_ancestry(5, population_size=1e4, model="smc", sequence_length=1e5,
                                 recombination_rate=1e-8, random_seed=5))
dump("smc_prime", msprime.sim_ancestry(5, population_size=1e4, model="smc_prime",
                                       sequence_length=1e5, recombination_rate=1e-8,
                                       random_seed=5))
dump("smck", msprime.sim_ancestry(5, population_size=1e4, model=msprime.SMCK(k=100),
                                  sequence_length=1e5, recombination_rate=1e-8,
                                  random_seed=6))
dump("beta", msprime.sim_ancestry(6, population_size=1e4,
                                  model=msprime.BetaCoalescent(alpha=1.5),
                                  sequence_length=1e4, recombination_rate=1e-8,
                                  random_seed=8))
dump("dirac", msprime.sim_ancestry(6, population_size=1e4,
                                   model=msprime.DiracCoalescent(psi=0.5, c=1),
                                   sequence_length=1e4, recombination_rate=1e-8,
                                   random_seed=9))
dump("sweep", msprime.sim_ancestry(
    5, population_size=1e3, sequence_length=1e4, recombination_rate=1e-7,
    model=[msprime.SweepGenicSelection(position=5000, start_frequency=1 / 2000,
                                       end_frequency=0.99, s=0.25, dt=1e-6),
           msprime.StandardCoalescent()], random_seed=10))
dump("gc", msprime.sim_ancestry(5, population_size=1e4, sequence_length=1e5,
                                recombination_rate=1e-8, gene_conversion_rate=1e-8,
                                gene_conversion_tract_length=100, random_seed=12))
dump("full_arg", msprime.sim_ancestry(4, population_size=1e4, sequence_length=1e4,
                                      recombination_rate=1e-7, record_full_arg=True,
                                      discrete_genome=False, random_seed=13))
dump("additional_nodes", msprime.sim_ancestry(
    4, population_size=1e4, sequence_length=1e4, recombination_rate=1e-7,
    coalescing_segments_only=False,
    additional_nodes=msprime.NodeType.RECOMBINANT | msprime.NodeType.COMMON_ANCESTOR,
    random_seed=14))
dump("end_time", msprime.sim_ancestry(6, population_size=1e4, end_time=5000,
                                      sequence_length=1e4, random_seed=15))
dump("ploidy1", msprime.sim_ancestry(6, ploidy=1, population_size=100, random_seed=16))
dump("ploidy4", msprime.sim_ancestry(3, ploidy=4, population_size=100,
                                     sequence_length=1000, recombination_rate=1e-4,
                                     random_seed=17))
dump("no_local_mrca", msprime.sim_ancestry(4, population_size=100, sequence_length=1000,
                                           recombination_rate=1e-4, random_seed=18,
                                           stop_at_local_mrca=False))
# rate map with missing flanks
rm = msprime.RateMap(position=[0, 1000, 5000, 9000, 10000],
                     rate=[np.nan, 1e-7, 5e-8, np.nan])
dump("ratemap", msprime.sim_ancestry(4, population_size=1e4, recombination_rate=rm,
                                     random_seed=19))
# replicates
for j, ts in enumerate(msprime.sim_ancestry(3, population_size=100, num_replicates=3,
                                            sequence_length=100, recombination_rate=1e-3,
                                            random_seed=20)):
    dump(f"rep_{j}", ts)
dump("rep_index", msprime.sim_ancestry(3, population_size=100, replicate_index=2,
                                       sequence_length=100, recombination_rate=1e-3,
                                       random_seed=20))
# initial state (recapitation)
ts0 = msprime.sim_ancestry(4, population_size=100, end_time=10, sequence_length=1000,
                           recombination_rate=1e-4, random_seed=21)
dump("init_base", ts0)
dump("recap", msprime.sim_ancestry(initial_state=ts0, population_size=200,
                                   recombination_rate=1e-4, random_seed=22))
# fixed pedigree
pb = msprime.PedigreeBuilder()
founders = [pb.add_individual(time=3) for _ in range(4)]
g2 = [pb.add_individual(time=2, parents=[founders[i], founders[(i + 1) % 4]])
      for i in range(4)]
g1 = [pb.add_individual(time=1, parents=[g2[i], g2[(i + 2) % 4]]) for i in range(4)]
for i in range(3):
    pb.add_individual(time=0, parents=[g1[i], g1[i + 1]])
ped = pb.finalise(sequence_length=1000)
dump("pedigree_tables", ped.tree_sequence())
dump("pedigree_sim", msprime.sim_ancestry(initial_state=ped, model="fixed_pedigree",
                                          recombination_rate=1e-3, random_seed=23))
# legacy simulate
dump("simulate_basic", msprime.simulate(10, Ne=1000, length=1e4, recombination_rate=1e-7,
                                        mutation_rate=1e-7, random_seed=24))
dump("simulate_popconf", msprime.simulate(
    population_configurations=[msprime.PopulationConfiguration(sample_size=4, initial_size=500),
                               msprime.PopulationConfiguration(sample_size=3, initial_size=500,
                                                               growth_rate=0.01)],
    migration_matrix=[[0, 1e-3], [1e-3, 0]],
    demographic_events=[msprime.MassMigration(time=200, source=1, dest=0, proportion=1),
                        msprime.CensusEvent(time=300)],
    length=1e4, recombination_rate=1e-7, random_seed=25))
# mutations
base = msprime.sim_ancestry(6, population_size=1e4, sequence_length=1e4,
                            recombination_rate=1e-8, random_seed=30)
dump("mut_base", base)
mut_models = {
    "jc69": msprime.JC69(), "hky": msprime.HKY(kappa=2.0, equilibrium_frequencies=[0.1, 0.2, 0.3, 0.4]),
    "f84": msprime.F84(kappa=3.0, equilibrium_frequencies=[0.3, 0.2, 0.2, 0.3]),
    "gtr": msprime.GTR(relative_rates=[1, 2, 3, 4, 5, 6]),
    "blosum62": msprime.BLOSUM62(), "pam": msprime.PAM(), "binary": msprime.BinaryMutationModel(),
    "jc69_si": msprime.JC69(state_independent=True),
    "hky_si": msprime.HKY(kappa=1.5, state_independent=True),
    "smm": msprime.SMM(lo=2, hi=10), "tpm": msprime.TPM(p=0.5, m=0.3, lo=2, hi=10),
    "el2": msprime.EL2(m=0.4, u=0.6, v=0.01, lo=2, hi=10),
    "microsat": msprime.MicrosatMutationModel(s=0.1, u=0.5, v=0.02, p=0.8, m=0.5, lo=1, hi=20),
    "infinite_alleles": msprime.InfiniteAlleles(start_allele=5),
    "slim_v5": msprime.SLiMv5MutationModel(type=2, next_id=10),
    "slim_v6": msprime.SLiMv6MutationModel(next_id=3),
}
for name, model in mut_models.items():
    dump("mut_" + name, msprime.sim_mutations(base, rate=1e-7, model=model, random_seed=31))
values["microsat_root"] = list(msprime.MicrosatMutationModel(s=0.1, u=0.5, v=0.02, p=0.8, m=0.5, lo=1, hi=20).root_distribution)
values["hky_matrix"] = msprime.HKY(kappa=2.0, equilibrium_frequencies=[0.1, 0.2, 0.3, 0.4]).transition_matrix.tolist()
ia = msprime.InfiniteAlleles()
m1 = msprime.sim_mutations(base, rate=1e-7, model=ia, random_seed=32)
m2 = msprime.sim_mutations(base, rate=1e-7, model=ia, random_seed=33)
dump("mut_ia_chain", m2)
values["ia_next"] = ia.next_allele
dump("mut_cont", msprime.sim_mutations(base, rate=1e-7, discrete_genome=False, random_seed=34))
mrm = msprime.RateMap(position=[0, 3000, 10000], rate=[1e-7, 5e-7])
dump("mut_ratemap", msprime.sim_mutations(base, rate=mrm, start_time=100, end_time=20000,
                                          random_seed=35))
dump("mut_keep", msprime.sim_mutations(msprime.sim_mutations(base, rate=1e-7, random_seed=36),
                                       rate=1e-7, random_seed=37))
dump("mutate_legacy", msprime.mutate(base, rate=1e-7, random_seed=38))
# likelihoods
arg = msprime.sim_ancestry(4, population_size=1, sequence_length=10, recombination_rate=0.1,
                           record_full_arg=True, discrete_genome=False, random_seed=40)
dump("lik_arg", arg)
values["log_arg_likelihood"] = msprime.log_arg_likelihood(arg, recombination_rate=0.1, Ne=1)
mts = msprime.sim_mutations(arg, rate=0.5, discrete_genome=False, random_seed=41)
dump("lik_mut", mts)
values["log_mutation_likelihood"] = msprime.log_mutation_likelihood(mts, 0.5)
# demography debugger
dbg = complex_demography().debug()
values["dbg_epoch_start"] = list(dbg.epoch_start_time)
values["dbg_size_history"] = dbg.population_size_history.tolist()
values["dbg_traj"] = dbg.population_size_trajectory([0, 10, 120, 500]).tolist()
steps = np.linspace(1, 1000, 30)
r, p = dbg.coalescence_rate_trajectory(steps, {"A": 2, "B": 2}, double_step_validation=False)
values["dbg_rates"] = r.tolist()
values["dbg_p"] = p.tolist()
P = dbg.lineage_probabilities([10, 160, 350])
values["dbg_lineage_probs"] = P.tolist()
dbg2 = msprime.Demography.island_model([100, 200], 1e-2).debug()
values["mean_coal"] = dbg2.mean_coalescence_time({0: 2})
# species tree & demes
st = msprime.Demography.from_species_tree("((A:10.0,B:10.0)AB:20.0,C:30.0)", initial_size=100)
dump("species_tree", msprime.sim_ancestry({"A": 2, "B": 2, "C": 2}, demography=st,
                                          random_seed=50))
import demes
yaml = """
time_units: years
generation_time: 25
demes:
  - name: ANC
    epochs: [{end_time: 5000, start_size: 2000}]
  - name: X
    ancestors: [ANC]
    epochs: [{start_size: 1000, end_size: 4000}]
  - name: Y
    ancestors: [ANC]
    epochs: [{end_time: 1000, start_size: 500}, {start_size: 800}]
migrations:
  - demes: [X, Y]
    rate: 1e-4
pulses:
  - {sources: [X], dest: Y, time: 2000, proportions: [0.1]}
"""
with open(os.path.join(out, "model.yaml"), "w") as f:
    f.write(yaml)
dm = msprime.Demography.from_demes(demes.loads(yaml))
dump("demes", msprime.sim_ancestry({"X": 2, "Y": 2}, demography=dm, sequence_length=1e4,
                                   recombination_rate=1e-8, random_seed=51))
values["demes_events"] = [str(e) for e in dm.events]
with open(os.path.join(out, "values.json"), "w") as f:
    json.dump(values, f)
print("done")
