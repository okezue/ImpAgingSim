# Response to Reviewers

*Independent control of sequence-correlation amplitude tunes post-quench composition fluctuations in A/B heteropolymer melts*

Alexander Okezue Bell

## Reviewer 2

### Comment 11: RPA calculation (follow-up)

> *Summary of the comment.* The RPA treatment is unclear. A matrix treatment along the lines of *Scaling Law for Sequence-Induced Demixing of Compositionally Identical Copolymers* was expected, rather than an unexplained scalar equation or an apparent dilute limit in which chain repulsions do not matter. The paragraph preceding Eq. (11) is also unclear. Since the RPA was optional, it may be corrected or removed.

**Response.** Thank you for identifying the ambiguity in the RPA presentation and for pointing me to Rumyantsev and Gavrilov, *Scaling Law for Sequence-Induced Demixing of Compositionally Identical Copolymers*. I have corrected this part of the manuscript and clarified the purpose and limits of the calculation.

I now begin with the full A/B structure-factor matrix for a compressible melt. The inverse response contains both the intramolecular correlation matrix and a finite reference contribution describing the common repulsive interactions. I explicitly derive the total-density and composition modes. Under A/B exchange symmetry, the common repulsive term remains in the density stiffness and cancels from the explicit composition stiffness. The resulting scalar composition expression is therefore a symmetry reduction of the matrix theory; no dilute limit is taken. The added Supplementary Information derivation also gives the density–composition coupling term that must be retained away from this symmetry.

I have clarified the normalization: the scalar composition spectrum is twice the symmetric contrast used to interpret the balanced-composition peak estimator. The paragraph preceding Eq. (11) now explains how the one-chain reference statistic is evaluated and averaged, and states explicitly that the plotted regressions are empirical comparisons, not fits or validations of the RPA. The fitted slope and intercept are not interpreted as physical RPA parameters.

I have also corrected an overinterpretation of the interaction-kernel diagnostic. The unweighted transform of the attractive force-equivalent kernel provides a bare interaction parameter, not an independently calibrated effective interaction for the dense melt. A negative inverse stiffness in that closure diagnoses its failure as a stable homogeneous Gaussian reference; it does not establish the stability or a phase boundary of the simulated post-quench system. I have revised the abstract, Results, Discussion, Methods and Fig. 4d consistently, and identify the zero-wavevector value as a formal long-wavelength limit rather than an observed canonical mode.

The cited work provides the appropriate matrix framework and proceeds to a Gaussian fluctuation free energy for demixing between distinct, compositionally identical copolymer populations. The observable in this study is the monomer-composition spectrum of a quenched stochastic sequence ensemble that permits chain-to-chain composition variation. I now explain this distinction and the relationship between the two formulations without applying the two-population demixing law to the present data.

These corrections retain the RPA framework with its assumptions made explicit. They do not require new molecular-dynamics trajectories or change the simulation results, empirical regression coefficients, or non-RPA conclusions.

**Location in revision:** Abstract; paragraph preceding Eq. (11) and the matrix-RPA formulation in Results; Fig. 4d and its caption; corresponding Discussion and Methods passages; Supplementary Information section *Two-component melt response and the symmetric projection*.
