# Pomona College Summer Undergraduate Research Program (SURP) 2019

Adam Guo, Pomona College '22  
Topic: Inversion of Laplace transform of a linear combination of point masses, with applications in medical imaging  
Research supervisor: Hrushikesh Mhaskar, Claremont Graduate University  
Faculty mentor: Adolfo Rumbos, Pomona College

Presented at Claremont Center for the Mathematical Sciences (CCMS) [Poster Session Fall 2019](https://colleges.claremont.edu/ccms/event/poster-session-fall-2019/) and Pomona College [2019 Intensive Summer Experience Poster Conference](https://www.pomona.edu/events/2019-intensive-summer-experience-poster-conference).

[poster](https://github.com/adamguos/surp-2019/blob/master/poster.pdf)

## Reference papers

- [Chui and Mhaskar - A Fourier-invariant method for locating point-masses and computing their attributes](https://arxiv.org/abs/1707.09319)
- [Sabett et al. - L1, Lp, L2, and elastic net penalties for regularization of Gaussian component distributions in magnetic resonance relaxometry](https://www.researchgate.net/publication/327692674_L1_Lp_L2_and_elastic_net_penalties_for_regularization_of_Gaussian_component_distributions_in_magnetic_resonance_relaxometry)
- [Townsend et al. - Fast computation of Gauss quadrature nodes and weights on the whole real line](https://arxiv.org/abs/1410.5286)
- [Tygert - Recurrence relations and fast algorithms](https://arxiv.org/abs/cs/0609081)

## Directories

Main project resides in inverse-laplace. Other directories contain supporting code/projects.

- practice: implementations of algorithms including Kaczmarz method, conjugate residual method, and Gauss-Legendre quadrature rules
- nfft-hermite: implementations of NFFT transforms based on the orthonormal Hermite polynomials
	- adc: [arrowhead divide and conquer algorithm](https://zenodo.org/record/1236142/files/article.pdf)
- gauss-quadrature: Gauss quadrature rules based on orthonormal Hermite polynomials
	- gausshermite: computes quadrature nodes and weights at the roots of the Hermite polynomials
	- gausshermitesc: computes quadrature weights at scattered (random) nodes
- inverse-laplace: Main project, centered around inverse_laplace.jl
	- wolframscript: Legacy WolframScript implementation
