# ipeval: an R package for validating predictions under interventions

## Abstract

We present the R package ipeval, which facilitates evaluation of
predictive performance under hypothetical intervention settings using
observational data. The package currently supports binary outcomes and
time-to-event outcomes under binary (point) interventions. It implements
methods to assess counterfactual predictive performance using inverse
probability of treatment weighting (IPTW).

## Introduction

Certain prediction models aim to estimate an individual’s risk under
hypothetical intervention scenarios. These are referred to as
interventional prediction models or models for prediction under
interventions (where intervention can e.g. be a certain medical
treatment but could also be a certain behavioral change or health
policy).¹ For example, a model may estimate a patient’s risk if
untreated, if undergoing surgery, or if initiating medication.^(2–4)
Models may also facilitate risk predictions under several treatment
options.⁵

Because such predictions under interventions are often intended to
inform medical decision-making, rigorous validation is essential. In
standard prediction settings, validation involves comparing estimated
risks with observed outcomes in a validation dataset. This approach is
not directly applicable to predictions under interventions in
observational data. Each individual typically receives only one
intervention, and outcomes that they would have had under alternative
interventions remain unobserved. These are referred to as counterfactual
outcomes in causal inference terminology.⁶

When models are used to guide treatment decisions, their performance
must be evaluated across all relevant intervention options for each
individual, not only the observed intervention option.^(1,7,8) Because
clinicians generally assign treatments based on patient characteristics,
treatment groups may not be comparable. Because of this, a model that
accurately estimates the risk of a patient under treatment may not
perform well estimating their counterfactual risk under no treatment,
even if the model also performs well in the group of patients that were
actually untreated.

Evaluating performance of the predictions under only the observed
treatments reflects predictive performance under the historical
treatment assignment mechanism present in the observed data. If such
models are used to inform treatment decisions in new patients, the
treatment assignment mechanism may change, rendering conventional
performance measures less relevant. In fact, a model that performs well
under historical treatment assignment mechanism can lead to harmful
decisions when applied to new patients.⁹

The appropriate target is counterfactual predictive performance: the
agreement between estimated risks and the outcomes that would be
observed if all individuals were assigned the specific intervention
option of interest. For example, how well do estimated risks align with
outcomes in a hypothetical scenario where all patients receive
treatment? Despite its importance, a recent review indicates that this
type of validation is rarely conducted.¹⁰ This package addresses this
gap by providing tools to perform such evaluations for binary and
time-to-event outcomes, based on the work of Keogh and Van Geloven
(2024).¹

## Methods

The implementation uses the observed validation data to approximate a
counterfactual dataset representing a population in which all
individuals receive a specified treatment option. This is achieved via
inverse probability of treatment weighting (IPTW). Each individual is
weighted by the inverse of the probability of receiving their observed
treatment conditional on confounders. This reweighting allows
individuals who received a given treatment to represent similar
individuals who did not receive this treatment. Validity of this
approach relies on standard causal assumptions: conditional
exchangeability, consistency, positivity, and correct specification of
the treatment model. More details are given elsewhere.^(1,6,7)

Predictive performance under a given intervention is then evaluated by
comparing estimated risks with observed outcomes in the weighted
(pseudo-)population. The package supports the following performance
metrics: area under the receiver operating characteristic curve (AUC),
Brier score, observed-expected (O/E) ratio, and calibration plots.

## Illustration

To demonstrate package functionality, we require a validation dataset
and interventional predictions to validate. We simulate data with binary
treatment $A$, binary outcome $Y$, continuous confounder $L$ and an
additional predictor $P$ (see Figure 1). Treatment assignment depends on
$L$, and treatment has a protective effect on the outcome (if $A = 1$
the probability of $Y = 1$ is lower than if $A = 0$). Table 1 shows the
first few rows of the data.

![Figure 1](dag.png)

Figure 1

## References

1\.

Keogh, R. H. & Van Geloven, N. [Prediction Under Interventions:
Evaluation of Counterfactual Performance Using Longitudinal
Observational Data](https://doi.org/10.1097/EDE.0000000000001713).
*Epidemiology (Cambridge, Mass.)* **35**, 329–339 (2024).

2\.

SCORE2 working group and ESC Cardiovascular risk collaboration. [SCORE2
risk prediction algorithms: New models to estimate 10-year risk of
cardiovascular disease in
Europe](https://doi.org/10.1093/eurheartj/ehab309). *European Heart
Journal* **42**, 2439–2454 (2021).

3\.

Nashef, S. A. M. *et al.* [EuroSCORE
II](https://doi.org/10.1093/ejcts/ezs043). *European Journal of
Cardio-Thoracic Surgery: Official Journal of the European Association
for Cardio-Thoracic Surgery* **41**, 734-744; discussion 744-745 (2012).

4\.

O’Brien, E. C. *et al.* [The ORBIT bleeding score: A simple bedside
score to assess bleeding risk in atrial
fibrillation](https://doi.org/10.1093/eurheartj/ehv476). *European Heart
Journal* **36**, 3258–3264 (2015).

5\.

Hageman, S. H. J. *et al.* [Estimation of recurrent atherosclerotic
cardiovascular event risk in patients with established cardiovascular
disease: The updated SMART2
algorithm](https://doi.org/10.1093/eurheartj/ehac056). *European Heart
Journal* **43**, 1715–1727 (2022).

6\.

Hernan, M. A. & Robins, J. M. *Causal Inference: What If*. (Chapman &
Hall/CRC, Boca Raton, 2020).

7\.

Boyer, C. B., Dahabreh, I. J. & Steingrimsson, J. A. [Estimating and
Evaluating Counterfactual Prediction
Models](https://doi.org/10.1002/sim.70287). *Statistics in Medicine*
**44**, e70287 (2025).

8\.

Pajouheshnia, R., Peelen, L. M., Moons, K. G. M., Reitsma, J. B. &
Groenwold, R. H. H. [Accounting for treatment use when validating a
prognostic model: A simulation
study](https://doi.org/10.1186/s12874-017-0375-8). *BMC Med Res
Methodol* **17**, 103 (2017).

9\.

Amsterdam, W. A. C. van, Geloven, N. van, Krijthe, J. H., Ranganath, R.
& Cinà, G. [When accurate prediction models yield harmful
self-fulfilling
prophecies](https://doi.org/10.1016/j.patter.2025.101229). *Patterns*
**6**, 101229 (2025).

10\.

Lin, L., Sperrin, M., Jenkins, D. A., Martin, G. P. & Peek, N. [A
scoping review of causal methods enabling predictions under hypothetical
interventions](https://doi.org/10.1186/s41512-021-00092-9). *Diagnostic
and Prognostic Research* **5**, 3 (2021).
