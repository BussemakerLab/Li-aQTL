import time
start_time = time.time()

import pandas as pd
import numpy as np
import tensorflow as tf
import tensorflow_probability as tfp
import sys
import os


### Load data

TF = sys.argv[1]
tis = sys.argv[2]

counts_DF = pd.read_csv("example_data/" + tis + ".counts.downsampled.csv", index_col=0)
samples = counts_DF.columns[1:]

pairs_DF = pd.read_csv("example_data/gene_pairs_neighboring.csv", index_col=0)


def makeDataSet(LFC_DF, TF):

    lfc = LFC_DF[TF].loc[np.invert(np.isnan(LFC_DF[TF]))]
    pairs_fil = pairs_DF[pairs_DF['g1'].isin(lfc.index) & pairs_DF['g2'].isin(lfc.index)]

    g1_counts = counts_DF[samples].loc[pairs_fil['g1']]
    g2_counts = counts_DF[samples].loc[pairs_fil['g2']]
    g1_sum = g1_counts.sum(axis=1)
    g2_sum = g2_counts.sum(axis=1)
    if_sum = (g1_counts.add(g2_counts.values) != 0)
    pairs_fil = pairs_fil[(g1_sum.values != 0) & (g2_sum.values != 0) & if_sum.all(axis=1).values]

    diff_lfc = lfc.loc[pairs_fil['g1']].values - lfc.loc[pairs_fil['g2']].values
    logr_sum = np.log(g1_sum.loc[pairs_fil['g1']].values / g2_sum.loc[pairs_fil['g2']].values)

    def makeSampleDF(s):
        return pd.DataFrame(
            {'g1_name': pairs_fil['g1'].values,
             'g2_name': pairs_fil['g2'].values,
             'g1': counts_DF[s].loc[pairs_fil['g1']].values,
             'g2': counts_DF[s].loc[pairs_fil['g2']].values,
             'diff_sig': diff_lfc,
             'logr_sum': logr_sum},
            index=pairs_fil['g1'].values)

    return [makeSampleDF(s=s) for s in samples]


LFC_DF = pd.read_csv("example_data/ENCODE_CRISPRi_shrinkLFC.csv", index_col=0)


df_ls = makeDataSet(LFC_DF=LFC_DF, TF=TF)

df_all = pd.concat(df_ls, keys=[("s" + str(s)) for s in range(len(df_ls))])


### Run joint fit

def loss(ag, lrhog, bs):
    av = tf.concat([ag for i in range(num_samples)], axis=0)
    lrhov = tf.concat([lrhog for i in range(num_samples)], axis=0)
    bv = tf.repeat(bs, repeats=num_pairs)

    rho = 1 / (1 + tf.math.exp(-lrhov))
    k = tf.cast(df_all["g1"], dtype=tf.float64)
    n = tf.cast(df_all["g1"] +df_all["g2"], dtype=tf.float64)
    x = tf.cast(df_all["diff_sig"], dtype=tf.float64)

    alpha = 1 / (1 + tf.math.exp(-(av + bv * x))) * (1 - rho) / rho
    beta = 1 / (1 + tf.math.exp(av + bv  * x)) * (1 - rho) / rho
    dist = tfp.distributions.BetaBinomial(total_count=n, concentration1=alpha, concentration0=beta)
    ll = dist.log_prob(k)
    if np.isinf(ll).sum() != 0:
      print(f"Inf: {np.where(np.isinf(ll))[0]}")
    return -tf.reduce_mean(ll)


def make_init(num_samples, num_pairs):
    init_a = tf.constant(np.repeat(0.0, num_pairs), dtype=tf.float64)
    init_lrho = tf.constant(np.repeat(np.log(0.1 / 0.9), num_pairs), dtype=tf.float64)
    init_b = tf.constant(np.repeat(0.01, num_samples), dtype=tf.float64)
    return tf.concat([init_a, init_lrho, init_b], axis=0)


def loss_and_gradient(params):
    return tfp.math.value_and_gradient(
        lambda params: loss(ag=params[0:num_pairs],
                            lrhog=params[num_pairs:(2*num_pairs)],
                            bs=params[(2*num_pairs):(2*num_pairs + num_samples)]),
        params)



num_samples = len(df_ls)
num_pairs = len(df_ls[0])
init_v = make_init(num_samples=num_samples, num_pairs=num_pairs)

pair_fit = tfp.optimizer.bfgs_minimize(loss_and_gradient, initial_position=init_v, max_iterations=10000)
print(f"converged: {pair_fit.converged}\niterations: {pair_fit.num_iterations}")

optim_a = pair_fit.position[0:num_pairs]
optim_rho = 1 / (1 + np.exp(-pair_fit.position[num_pairs:(2*num_pairs)]))
optim_b = pair_fit.position[(2*num_pairs):(2*num_pairs + num_samples)]


### Save coefficients (inferred activities)

with tf.io.gfile.GFile('example_data/coefs_{}_{}.csv'.format(TF, tis), 'w') as f:
     np.savetxt(f, optim_b, delimiter=',')


print("Run: %s hours ---" % round((time.time() - start_time)/3600, 2))


