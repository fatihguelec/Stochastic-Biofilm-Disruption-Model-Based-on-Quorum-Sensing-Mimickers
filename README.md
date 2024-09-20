# **Stochastic Biofilm Disruption Model based on Quorum Sensing Mimickers**

This repository contains the MATLAB code developed for the paper "**A Stochastic Biofilm Disruption Model based on Quorum Sensing Mimickers**," published in IEEE Transactions on Molecular, Biological and Multi-Scale Communications. The code implements a stochastic simulation method to model the disruption of biofilms using quorum sensing (QS) mimickers. The results generated from this code correspond to Figures 2-6 in the paper.

## **Abstract**

Quorum sensing (QS) mimickers can be used as an effective tool to disrupt biofilms, which consist of communicating bacteria and extracellular polymeric substances (EPS). This paper presents a stochastic biofilm disruption model using QS mimickers. The model is based on a chemical reaction network (CRN) with four different states, capturing the processes of biofilm formation and its disruption via QS mimickers. A state-based stochastic simulation algorithm is proposed to simulate this CRN. The model is validated using experimental results of *Pseudomonas aeruginosa* biofilm disrupted by rosmarinic acid, a QS mimicker. The results highlight the inherent randomness of the CRN, revealing uncertainties in state transitions. Besides the QS activation threshold, two additional thresholds are identified for EPS and bacterial disruption, providing a realistic model for biofilm disruption through QS mimickers.

## **Code Overview**

This repository includes MATLAB code files named after Figures 2-6 from the paper. The code is based on a modified version of the direct Gillespie algorithm:

- **D. T. Gillespie, “Exact stochastic simulation of coupled chemical reactions,” J. Phys. Chem., vol. 81, no. 25, pp. 2340–2361, 1977.**

The modifications incorporate a four-state mechanism representing different stages of biofilm formation and disruption. The first two states model biofilm formation before (downregulation) and after (upregulation) QS activation, depending on the concentrations of autoinducer and QS mimickers. The last two states represent EPS disruption and bacterial killing, both triggered by QS mimickers through two distinct thresholds. The stochastic simulation captures randomness in state transitions and changes in bacterial and EPS concentrations within the biofilm.

## **How to Run the Code**

- The code files ending with "Fig. 2" and "Fig. 3" validate the results extracted from **A. Corral-Lugo, A. Daddaoua, A. Ortega, M. Espinosa-Urgel, and T. Krell, “Rosmarinic acid is a homoserine lactone mimic produced by plants that activates a bacterial quorum-sensing regulator,” Sci. Signal., vol. 9, no. 409, p. RA1, 2016**. These files require `.xlsx` data files, which should be placed in the same directory as the MATLAB code files.
- Additionally, make sure the `gillespie_direct.m` file is also in the same directory as the code files.
- Simply run the provided MATLAB code files to generate the results corresponding to Figures 2-6 of the paper.

## **Citation Requirement**

If you use or build upon this code in your research, or if this code is used for any academic work, publication, or research, proper attribution and citation of the paper is **required**. Please cite the paper in any related publications, presentations, or derived works.

**Gulec, F., & Eckford, A. W., "A Stochastic Biofilm Disruption Model based on Quorum Sensing Mimickers," IEEE Transactions on Molecular, Biological and Multi-Scale Communications, vol. 9, no. 3, pp. 346-350, 2023.**

