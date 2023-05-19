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

LFC_DF = pd.read_csv("example_data/ENCODE_CRISPRi_shrinkLFC.csv", index_col=0)


lfc = LFC_DF[TF].loc[np.invert(np.isnan(LFC_DF[TF]))]
counts_fil = counts_DF[counts_DF["gene"].isin(lfc.index)]
g_sum = counts_fil[samples].sum(axis=1)  #filter if (all genes=0)
counts_fil = counts_fil[(g_sum.values != 0)]
del g_sum


def makeSampleDF(s):
    return pd.DataFrame(
        {'count': counts_fil[s].values,
          'sig': lfc.loc[counts_fil['gene']].values},
        index=counts_fil['gene'].values)
            

df_ls = [makeSampleDF(s=s) for s in samples]
df_all = pd.concat(df_ls, keys=[("s" + str(s)) for s in range(len(df_ls))])

num_samples = len(df_ls)
num_genes = len(df_ls[0])


### Run fit (with fixed a) 

g_mean = counts_fil[samples].mean(axis=1)
g_var = counts_fil[samples].var(axis=1)
ag = np.log(g_mean)
ag = tf.cast(ag, dtype=tf.float64)

est_rg = np.square(g_mean)/(g_var - g_mean)
est_rg[est_rg <= 0] = 0.01
est_rg[np.isinf(est_rg)] = 0.01

init_rlog = tf.cast(np.log(est_rg), dtype=tf.float64)


def loss(rlog, bs):
    av = tf.concat([ag for i in range(num_samples)], axis=0)
    rg = tf.math.exp(rlog)
    rv = tf.concat([rg for i in range(num_samples)], axis=0)
    bv = tf.repeat(bs, repeats=num_genes)
    k = tf.cast(df_all["count"], dtype=tf.float64)
    x = tf.cast(df_all["sig"], dtype=tf.float64)
    mu = tf.math.exp(av + bv * x)
    dist = tfp.distributions.NegativeBinomial(total_count=rv, probs=1/(1+rv/mu))
    ll = dist.log_prob(k)
    if np.isinf(ll).sum() != 0:
      ind = np.where(np.isinf(ll))[0]
      print(f"Inf: {ind}")
    return -tf.reduce_mean(ll) 
#-tf.reduce_mean(tf.boolean_mask(ll, tf.math.is_finite(ll)))


def loss_and_gradient(params):
    return tfp.math.value_and_gradient(
        lambda params: loss(rlog=params[0:num_genes],
                            bs=params[num_genes:(num_genes + num_samples)]),
        params)


init_b = tf.constant(np.repeat(0.0, num_samples), dtype=tf.float64)
init_v = tf.concat([init_rlog, init_b], axis=0)

fit_rb = tfp.optimizer.bfgs_minimize(loss_and_gradient, initial_position=init_v, max_iterations=10000)
print(f"converged: {fit_rb.converged}\niterations: {fit_rb.num_iterations}")


optim_r = fit_rb.position[0:num_genes]
optim_b = fit_rb.position[(num_genes):(num_genes + num_samples)]


with tf.io.gfile.GFile('example_data/GeneModel_coefs.{}.{}.b.csv'.format(TF, tis), 'w') as f:
     np.savetxt(f, optim_b, delimiter=',')

print("Run: %s hours ---" % round((time.time() - start_time)/3600, 2))




