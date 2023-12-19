CDDC README File
=================================================
1. Running the Code
To reproduce the experiments in [1], you need to perform three steps: 
(a) Run mexC.m file. Note that this step is optional and only used for the One-Hidden-Layer Networks Problem.
(b) To obtain the accuracy for the 6 applications, run demo_accuracy1_*.m / demo_accuracy2_*.m / demo_accuracy3_*.m / demo_accuracy4_*.m / demo_accuracy5_*.m / demo_accuracy6_*.m 
(c) To obtain the efficiency for the 6 applications, run demo_cpu1_*.m / demo_cpu2_*.m / demo_cpu3_*.m / demo_cpu4_*.m / demo_cpu5_*.m / demo_cpu6_*.m 

2. Directory tree
We demonstrate the directory tree of our code, including the names of some main functions and their descriptions.

------ solvers
  |      |
  |      |------ L1PCA_MSCR.m                   : Implementation of Multi-Stage Convex Relaxation for the L1 Norm Generalized Eigenvalue Problem
  |      |------ L1PCA_TolandDual.m             : Implementation of Toland-Dual Method for the L1 Norm Generalized Eigenvalue Problem
  |      |------ L1PCA_SubGrad.m                : Implementation of SubGradient Method for the L1 Norm Generalized Eigenvalue Problem
  |      |------ L1PCA_CD_SCA.m                 : Implementation of CD-SCA for the L1 Norm Generalized Eigenvalue Problem
  |      |------ L1PCA_CD_SNCA.m                : Implementation of CD-SNCA for the L1 Norm Generalized Eigenvalue Problem
  |      |------ L1PCA_ComputeTrueObj.m         : The objective function used in the L1 Norm Generalized Eigenvalue Problem
  |      |
  |      |
  |      |------ OneNet_MSCR.m                  : Implementation of Multi-Stage Convex Relaxation for the One-Hidden-Layer Networks Problem
  |      |------ OneNet_PDCA.m                  : Implementation of Proximal DC Algorithm for the One-Hidden-Layer Networks Problem
  |      |------ OneNet_SubGrad.m               : Implementation of SubGradient Method for the One-Hidden-Layer Networks Problem
  |      |------ OneNet_CD_SCA.m                : Implementation of CD-SCA for the One-Hidden-Layer Networks Problem
  |      |------ OneNet_CD_SNCA.m               : Implementation of CD-SNCA for the One-Hidden-Layer Networks Problem
  |      |------ OneNet_ComputeTrueObj.m        : The objective function used in the One-Hidden-Layer Networks Problem
  |      |
  |      |
  |      |------ SparseOpt_MSCR.m               : Implementation of Multi-Stage Convex Relaxation for the DC Penalized Sparsity Constrained Optimization Problem
  |      |------ SparseOpt_PDCA.m               : Implementation of Proximal DC Algorithm for the DC Penalized Sparsity Constrained Optimization Problem
  |      |------ SparseOpt_SubGrad.m            : Implementation of SubGradient Method for the DC Penalized Sparsity Constrained Optimization Problem
  |      |------ SparseOpt_CD_SCA.m             : Implementation of CD-SCA for the DC Penalized Sparsity Constrained Optimization Problem
  |      |------ SparseOpt_CD_SNCA.m            : Implementation of CD-SNCA for the DC Penalized Sparsity Constrained Optimization Problem
  |      |------ SparseOpt_ComputeTrueObj.m     : The objective function used in the DC Penalized Sparsity Constrained Optimization Problem
  |      |
  |      |
  |      |------ L4EIG_MSCR.m                   : Implementation of Multi-Stage Convex Relaxation for the L4 Norm Regularized Eigenvalue Problem
  |      |------ L4EIG_PDCA.m                   : Implementation of Proximal DC Algorithm for the L4 Norm Regularized Eigenvalue Problem
  |      |------ L4EIG_GProj.m                  : Implementation of Gradient Projection Method for the L4 Norm Regularized Eigenvalue Problem
  |      |------ L4EIG_SCA.m                    : Implementation of CD-SCA for the L4 Norm Regularized Eigenvalue Problem
  |      |------ L4EIG_SNCA.m                   : Implementation of CD-SNCA for the L4 Norm Regularized Eigenvalue Problem
  |      |------ L4EIG_ComputeTrueObj.m         : The objective function used in the L4 Norm Regularized Eigenvalue Problem
  |      |
  |      |
  |      |------ BinaryOpt_MSCR.m               : Implementation of Multi-Stage Convex Relaxation for the DC Penalized Binary Optimization Problem
  |      |------ BinaryOpt_PDCA.m               : Implementation of Proximal DC Algorithm for the DC Penalized Binary Optimization Problem
  |      |------ BinaryOpt_SubGrad.m            : Implementation of SubGradient Method for the DC Penalized Binary Optimization Problem
  |      |------ BinaryOpt_CD_SCA.m             : Implementation of CD-SCA for the DC Penalized Binary Optimization Problem
  |      |------ BinaryOpt_CD_SNCA.m            : Implementation of CD-SNCA for the DC Penalized Binary Optimization Problem
  |      |------ BinaryOpt_ComputeTrueObj.m     : The objective function used in the DC Penalized Binary Optimization Problem
  |      |
  |      |
  |      |------ PhaseRetrivalAltMin.m          : Implementation of Alternating Minimization Method for the Phase Retrieval Problem
  |      |------ PhaseRetrivalGradProj.m        : Implementation of Gradient Projection Method for the Phase Retrieval Problem
  |      |------ PhaseRetrivalGPM.m             : Implementation of Generalized Power Method for the Phase Retrieval Problem
  |      |------ PhaseRetrivalCD_SCA.m          : Implementation of CD-SCA for the Phase Retrieval Problem
  |      |------ PhaseRetrivalCD_SNCA.m         : Implementation of CD-SNCA for the Phase Retrieval Problem
  |      |------ PhaseRetrival_ComputeTrueObj.m : The objective function used in the Phase Retrieval Problem
  |
  |
  |
------- util
  |      |
  |      |------ APG.m                          : Accelerated Proximal Gradient Method for solving convex composite optimization problems
  |      |------ PhaseRetrivalComputeQ.m        : Compute the auxiliary matrix for the Phase Retrieval Problem
  |      |------ ProjBinary.m                   : Euclidean projection for the binary constraints: min_v ||v - a||_2^2, s.t. v \in {-1,+1}^n
  |      |------ prox_l1.m                      : Proximal operator for the L1 norm: min_{x} 0.5||x-a||^2 + lambda*||x||_1
  |      |------ SpectralNorm.m                 : Estimating the spectral norm of the matrix G'G using the Lanczos algorithm
  |      |------ top_k_subgrad.m                : Compute the subgradient of the top k norm function (in absolute value)
  |      |------ GetAvgCell.m                   : Compute the average vector for multiple variable-length vectors. Note: we truncate the vectors to make them the same length.
  |
  |
  |
----- proximal
  |      |
  |      |------ ncvx_prox_l1.m                 : Nonconvex proximal operator for the l1 norm function
  |      |------ ncvx_prox_linf.m               : Nonconvex proximal operator for the linf norm function
  |      |------ ncvx_prox_relu.m               : Nonconvex proximal operator for the relu function
  |      |------ ncvx_prox_l1topk.m             : Nonconvex proximal operator for the l1-top DC function
  |      |------ ncvx_prox_poly4.m              : Nonconvex proximal operator for the 4th-order polynomial function
  |      |------ ncvx_prox_linfl2.m             : Nonconvex proximal operator for the linf-l2 function
  |      |------ ncvx_prox_quad.m               : Nonconvex proximal operator for the quadratic function
  |
  |
  |
------ data
  |      |
  |      |------ GetDataSet.m                   : Retrieve the data and perform data preprocessing
  |      |------ NormData.m                     : Normalize the data
  |      |------ GetDataMeasDCBinary.m          : Generate the measurement vector for the DC Penalized Binary Optimization Problem
  |      |------ GetDataMeasDCSparse.m          : Generate the measurement vector for the DC Penalized Sparsity Constrained Optimization Problem
  |      |------ GetDataMeasOneNet.m            : Generate the measurement vector for the One-Hidden-Layer Networks Problem
  |      |------ GetDataMeasPhaseRetrival.m     : Generate the measurement vector for the Phase Retrieval Problem
  |      |------ *.mat                          : five matlab datasets
  |
  |
  |
----TableFigure
         |
         |------ DrawFigurei_*.m                : Generate the figures for the i-th application
         |------ ShowLatexTablei.m              : Generate the latex source codes for the i-th application
	 
	 
3. REFERENCES:
[1] Ganzhao Yuan. Coordinate Descent Methods for DC Minimization: Optimality Conditions and Global Convergence. Proceedings of the AAAI Conference on Artificial Intelligence (AAAI), 2023.
[2] Ganzhao Yuan. CD Methods for DC Minimization. Submitted for publication.


	 
