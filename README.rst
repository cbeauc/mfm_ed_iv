.. REFERENCES
.. |qtab24| replace:: Quirouette et al. 2024
.. _qtab24: https://arxiv.org/abs/2412.12960
.. |psimon| replace:: Simon et al. 2016
.. _psimon: https://doi.org/10.1038/srep24154
.. |midSIN| replace:: midSIN
.. _midSIN: https://midSIN.roadcake.org
.. |phymcmc| replace:: phymcmc
.. _phymcmc: https://github.com/cbeauc/phymcmc

mfm_ed_iv
=========

Simulates the course of an in vitro virus infection using an ODE model, a stochastic model (SM), or a mean-field version of the SM (MFSM) where random variables are replaced by their mean. The code includes a script showing how to estimate model parameters from viral titer measurements taken over the course of a virus infection using the raw outcomes (number of infected wells at each dilution) of an endpoint dilution (TCID50) assay. This project includes the code and data associated with |qtab24|_.

.. image:: https://img.shields.io/badge/GitHub-cbeauc%2Fmfm_ed_iv-blue.svg?style=flat
    :target: https://github.com/cbeauc/mfm_ed_iv
.. image:: https://img.shields.io/badge/license-GPL-blue.svg?style=flat
    :target: https://github.com/cbeauc/mfm_ed_iv/blob/master/LICENSE


Download
--------

To get the code, scripts, and data (with git)

	$ git clone https://github.com/cbeauc/mfm_ed_iv.git


Key files and directories
-------------------------

``models.py``
	function ``infection_simulator`` reproduces the single-cycle and multiple-cycle infections as described in |qtab24|_. New work for a different data set or experiment should write a replacement for this function, and make use of any one of the ODE model (``solve_ode_ivp``), stochastic model (SM, ``solve_stoch_ivp``) or its mean-field counterpart (MFSM, ``solve_meanfield_ivp``).


``estimator.py``
	``base_params_dict``
		base values of all fixed and estimated parameters, as described in Table 3 of |qtab24|_

	``get_PVestED``
		estimates P\ :sub:`V->Establishment`\  in the ED assay

	``params``
		whose method ``validate`` imposes the constraints on the model parameters

	``data_packer``
		parses the data files and makes some pre-computations from the data to ensure the posterior probability can be computed faster for each parameter set attempt, and is specific to the data and experiments in |qtab24|_

	``estimator``
		computes the posterior probability of a certain parameter set based on either the SSR (``ln_LikeSSR``) or ED (``ln_LikeED``) likelihood function.


``data/``
	directory containing the complete dataset from |psimon|_ where the first column is always time in hours. The ``SIN`` files contain the log\ :sub:`10`\  highest likelihood SIN concentration (log\ :sub:`10`\  SIN/ml) in the sample as estimated by |midSIN|_; ``RNA`` is the viral RNA concentration (vRNA/ml); ``parsed`` contains the log\ :sub:`10`\  fold-dilution (-1 = 0.1 or 10-fold, -2 = 0.01 or 100-fold, etc.) for all 8 columns of the assay, and the remaining 8 columns are a the number of infected wells out of 4 replicates at each dilution in the ED assay.

``chains/``
	directory containing the final (thinned, post burn-in) chain files used in |qtab24|_, whose name indicates the variant, e.g. model choice of MFM (ODE) vs SM; likelihood choice of SSR (Gaussian) vs ED (based on endpoint dilution); and units choice of SIN vs IV. The chain files are in the |phymcmc|_ hdf5 format.


Usage
-----

To create a chain file to collect MCMC samples that will make-up your posterior distribution, you'll need to install |phymcmc|_ and then you can run the estimator

	$ ./run_mcmc [choice]

where ``[choice]`` is of the form ``[model]_[likelihood]_[units]`` where ``[model]`` is either ``mfm`` or ``sm``; ``[likelihood]`` is either ``ssr`` or ``ed``; and ``[units]`` is either ``sin`` or ``iv``. In |qtab24|_, the analysis suggests ``mfm_ed_iv`` is your best option, whereas ``sm_ed_iv`` is the most accurate, but most computationally costly.


Attribution
-----------

If you make use of this code, make sure to cite our paper.

The BibTeX entry is::

	@ARTICLE{quirouette24,
		AUTHOR = "Christian Quirouette and Risavarshni Thevakumaran and Kyosuke Adachi and Catherine A. A. Beauchemin",
		TITLE = "Does the random nature of cell-virus interactions during in vitro infections affect TCID50 measurements and parameter estimation by mathematical models?",
		JOURNAL = "arXiv",
		MONTH = "December",
		YEAR = "2024",
		DOI = "10.48550/arXiv.2412.12960",
		EPRINT = "arXiv:2412.12960",
	}


License
-------

mfm_ed_iv is free software made available under the GNU General Public License Version 3. For details see the LICENSE file.
