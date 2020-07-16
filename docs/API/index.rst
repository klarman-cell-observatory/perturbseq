
API
===

Import perturbseq as::

   import perturbseq as perturb                                                     

                                       
Read in data: `io`
------------------

Read in guides.

.. module:: perturbseq.io
.. currentmodule:: perturbseq

.. autosummary::
   :toctree: .

   io.read_perturbations_csv
 

Preprocessing: `pp`
-------------------

.. module:: perturbseq.pp
.. currentmodule:: perturbseq

.. autosummary::
   :toctree: .

   pp.annotate_controls
   pp.remove_guides_from_gene_names
   pp.perturbations_per_cell
   pp.subsample_cells
   pp.cells_per_perturbation
   pp.perturb2obs
   pp.moi



