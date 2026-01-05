# Amyloid-Fibril-to-Ionic-microenvironments

This repository contains the computational workflows, analysis scripts, and supporting data for the study:

Thermodynamics of Amyloid Peptide Transfer from Crystalline Fibrils into Ionic Microenvironments

The project investigates the thermodynamic stability of amyloid-forming peptides as they transition from an ordered crystalline fibril state into solvated ionic environments, using a combination of molecular dynamics simulations, nonequilibrium free-energy methods, and quantum chemical calculations.

----------------------------------
In this work, we introduce a free-energy framework that explicitly tracks peptide detachment from amyloid crystals into solution. The approach combines:

All-atom steered molecular dynamics (SMD)

Jarzynski-based free-energy reconstruction

Symmetry-adapted perturbation theory (SAPT) for ion–peptide interaction analysis

The framework is applied to the NNQQ peptide, a minimal amyloidogenic motif derived from the yeast prion Sup35, in different ionic environments (water, NaCl, and KCl).

The overall computational workflow is illustrated below.

<p align="center">
  <img src="Jar-Sch.png" alt="Schematic of the computational workflow" width="700"/>
</p>


Figure: Schematic representation of the computational framework used to construct amyloid crystal models, perform nonequilibrium pulling simulations, calculate free-energy profiles, and analyze ion–peptide interactions using quantum chemical methods.
