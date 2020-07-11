
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

   perturbseq.pp.annotate_controls
   perturbseq.pp.remove_guides_from_gene_names
   perturbseq.pp.perturbations_per_cell
   perturbseq.pp.perturbations_per_cell
   perturbseq.pp.subsample_cells
   perturbseq.pp.cells_per_perturbation
   perturbseq.pp.perturb2obs
   perturbseq.pp.moi

Tools: `tl`
-------------------

.. module:: perturbseq.tl
.. currentmodule:: perturbseq

.. autosummary::
   :toctree: .

