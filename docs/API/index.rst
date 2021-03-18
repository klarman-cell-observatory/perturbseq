
API
===


You can import perturbseq as::

   import perturbseq as perturb                                                                                            

.. note::
   Testnote

Preprocessing: `pp`
-------------------

Perturbation assignments and queries

.. module:: perturbseq.pp
.. currentmodule:: perturbseq

.. autosummary::
   :toctree: .

   pp.compute_TPT
   pp.get_perturbations
   pp.perturb_overlap_obs
   pp.annotate_controls
   pp.subset_singly_perturbed
   pp.perturbs_per_cell
   pp.cells_per_perturb
   pp.delete_guides_from_varnames

Sub/down-sampling

.. module:: perturbseq.pp
.. currentmodule:: perturbseq

.. autosummary::
   :toctree: .

   pp.subsample_cells
   pp.downsample_counts

Summarization across perturbations

.. module:: perturbseq.pp
.. currentmodule:: perturbseq

.. autosummary::
   :toctree: .   

   pp.obs_mean

Preparing for modeling
   
.. module:: perturbseq.pp
.. currentmodule:: perturbseq

.. autosummary::
   :toctree: .

   pp.obs_to_design_matrix
   pp.split_train_valid_test


Tools: `tl`
-------------------

.. module:: perturbseq.tl
.. currentmodule:: perturbseq

.. autosummary::
   :toctree: .

   tl.moi

