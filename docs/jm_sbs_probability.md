# Joint-Multiplicity SBS Spectrum

Let $g=(i,j)$ represent a joint-multiplicity (JM) group and let $m$
represent an SBS signature. The target probability is

$$
\theta_{g,m}
=
P(\mathrm{SBS}=m \mid \mathrm{JM}=g).
$$

## Hard JM Assignment

For mutation $n$, SigProfilerAssignment provides the SBS probabilities

$$
p_{n,m}
=
P(\mathrm{SBS}_n=m \mid X_n),
$$

where $X_n$ contains the mutation context and information used by the fitted
signature model.

## Preferred Input: Shared SSNVs as One Sample

For this analysis, the biological unit is the shared SSNV. A shared SSNV
occurred once, even though it is observed in both tumor samples. Therefore, the
cleanest approach is to build one mutation table where each shared SSNV appears
once, then run SigProfilerAssignment on this shared-SSNV set as one sample.

This produces one SBS probability vector per shared SSNV:

$$
p_{n,m}
=
P(\mathrm{SBS}_n=m \mid X_n).
$$

Here, $X_n$ is the mutation context and fitted signature model information used
to assign SBS probabilities for mutation $n$. In the simplest interpretation,
$X_n$ is the trinucleotide context of mutation $n$. If the assignment uses a
fitted exposure profile, then $X_n$ also includes that fitted shared-SSNV
signature background.

This avoids assigning two different SBS probability vectors to the same shared
mutation. The resulting target is

$$
P(\mathrm{SBS}=m \mid \mathrm{shared\ SSNV}, \mathrm{JM}=g).
$$

The interpretation is specific to the selected shared-SSNV subset, not the whole
tumor sample.

## Shared SNVs With Sample-Specific SBS Probabilities

If SigProfilerAssignment is instead run separately on the two tumor samples, the
same shared SSNV can have two SBS probability vectors:

$$
p^{(1)}_{n,m}
=
P(\mathrm{SBS}_n=m \mid X^{(1)}_n)
$$

and

$$
p^{(2)}_{n,m}
=
P(\mathrm{SBS}_n=m \mid X^{(2)}_n).
$$

Before estimating the JM-specific SBS spectrum, collapse these into one
mutation-level SBS probability vector:

$$
\bar{p}_{n,m}
=
P(\mathrm{SBS}_n=m \mid X^{(1)}_n, X^{(2)}_n).
$$

A simple default is the average of the two sample-specific probabilities:

$$
\bar{p}_{n,m}
=
\frac{p^{(1)}_{n,m}+p^{(2)}_{n,m}}{2}.
$$

If the two samples have different reliability for mutation $n$, a weighted
average can be used:

$$
\bar{p}_{n,m}
=
\frac{
w_{n,1}p^{(1)}_{n,m}
+
w_{n,2}p^{(2)}_{n,m}
}{
w_{n,1}+w_{n,2}
}.
$$

The weights can represent read depth, alternate-read support, or another
sample-level confidence measure. A more aggressive alternative is a
product-of-experts combination:

$$
\bar{p}_{n,m}
=
\frac{
p^{(1)}_{n,m}p^{(2)}_{n,m}
}{
\sum_k p^{(1)}_{n,k}p^{(2)}_{n,k}
}.
$$

This emphasizes SBS signatures supported by both samples, but can become too
sharp when one sample has noisy probabilities. This sample-specific approach is
less clean than fitting the shared SSNVs as one sample, because the same mutation
is forced to reconcile two sample-level assignment backgrounds.

For mutations assigned to JM group $g$, the expected number of mutations
attributed to SBS signature $m$ is

$$
\widehat{C}_{g,m}
=
\sum_{n:\mathrm{JM}_n=g} p_{n,m}.
$$

The estimated SBS spectrum of JM group $g$ is

$$
\widehat{\theta}_{g,m}
=
\frac{\widehat{C}_{g,m}}
     {\sum_k \widehat{C}_{g,k}}.
$$

If the signature probabilities for every mutation sum to one, then

$$
\sum_k \widehat{C}_{g,k}=N_g,
$$

where $N_g$ is the number of mutations assigned to JM group $g$. Therefore,

$$
\widehat{\theta}_{g,m}
=
\frac{1}{N_g}
\sum_{n:\mathrm{JM}_n=g}p_{n,m}.
$$

If sample-specific SBS probabilities are used instead, replace $p_{n,m}$ with
the combined mutation-level probability $\bar{p}_{n,m}$ in the formulas above.

For the 3:1 versus 2:1 pure-gain analysis, the two spectra of interest are

$$
P(\mathrm{SBS}=m \mid \mathrm{JM}=(1,1))
$$

and

$$
P(\mathrm{SBS}=m \mid \mathrm{JM}=(2,1)).
$$

Because these mutations belong to the truncal subroute and accumulated in the
diploid state, the JM groups do not necessarily represent different mutation
times. Their joint multiplicities result from subsequent copy-number gains and
copy fates.

## Partial Maximum-Likelihood JM Assignment

For shared SSNV $n$, let $g=(i,j)$ be a candidate JM group. Let $r_{n,s}$ be
the alternate-read count and $d_{n,s}$ be the total read depth in sample $s$.
Given purity, copy number, and other model parameters collected as $\phi$, the
expected VAF in sample $s$ under JM group $g$ is

$$
v_s(g,\phi).
$$

The read-count likelihood for mutation $n$ under JM group $g$ can be written as

$$
L_{n,g}
=
P(r_n \mid d_n,\phi,\mathrm{JM}_n=g).
$$

For two samples, assuming the read counts are conditionally independent across
samples given the JM group and model parameters,

$$
L_{n,g}
=
\prod_{s=1}^{2}
\mathrm{Binomial}
\left(
r_{n,s};
d_{n,s},
v_s(g,\phi)
\right).
$$

Equivalently, the log-likelihood is

$$
\ell_{n,g}
=
\sum_{s=1}^{2}
\log
\mathrm{Binomial}
\left(
r_{n,s};
d_{n,s},
v_s(g,\phi)
\right).
$$

The hard partial maximum-likelihood JM assignment is

$$
\widehat{\mathrm{JM}}_n
=
\arg\max_g L_{n,g}
=
\arg\max_g \ell_{n,g}.
$$

This is a maximum-likelihood assignment, not automatically a probability. It
uses only the read-count likelihood under the candidate JM groups. If priors
over JM groups are not included, this is equivalent to treating candidate JM
groups as equally likely before observing the read counts.

The hard-assignment SBS spectrum can then be estimated as

$$
\widehat{C}_{g,m}
=
\sum_{n:\widehat{\mathrm{JM}}_n=g}p_{n,m}
$$

and

$$
\widehat{\theta}_{g,m}
=
\frac{
\sum_{n:\widehat{\mathrm{JM}}_n=g}p_{n,m}
}{
\sum_{n:\widehat{\mathrm{JM}}_n=g}1
}.
$$

## Probabilistic JM Assignment

If SPUR provides a probability that mutation $n$ belongs to JM group $g$,
define

$$
q_{n,g}
=
P(\mathrm{JM}_n=g \mid D_n).
$$

If only partial likelihoods are available, they can be normalized into soft JM
weights:

$$
q_{n,g}
=
\frac{L_{n,g}\pi_g}
     {\sum_h L_{n,h}\pi_h},
$$

where $\pi_g$ is the prior probability of JM group $g$. With equal priors,

$$
q_{n,g}
=
\frac{L_{n,g}}
     {\sum_h L_{n,h}}.
$$

In log space, this is

$$
q_{n,g}
=
\exp
\left(
\ell_{n,g}
-
\log\sum_h \exp(\ell_{n,h})
\right)
$$

for equal priors.

An approximate expected joint count is then

$$
\widehat{C}_{g,m}
=
\sum_n q_{n,g}p_{n,m},
$$

and the normalized spectrum is

$$
\widehat{\theta}_{g,m}
=
\frac{\sum_n q_{n,g}p_{n,m}}
     {\sum_n q_{n,g}}.
$$

This calculation approximates

$$
P(\mathrm{SBS}_n=m,\mathrm{JM}_n=g)
\approx
p_{n,m}q_{n,g}.
$$

It assumes that the evidence used for SBS assignment and the evidence used for
JM assignment are conditionally independent under the fitted models.

Again, if sample-specific SBS probabilities are used, replace $p_{n,m}$ with
$\bar{p}_{n,m}$ after combining the sample-specific vectors.

## Viewing in VS Code

1. Open this Markdown file in VS Code.
2. Press `Ctrl+Shift+V` to open the rendered preview.
3. Alternatively, press `Ctrl+K` and then `V` to open the preview beside the
   Markdown source.
4. If equations do not render, enable the `Markdown > Math: Enabled` setting
   and verify that VS Code's built-in `Markdown Language Features` extension is
   enabled.
