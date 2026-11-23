Python reference
================

Generated from the type stubs ``bitvf/__init__.pyi`` and ``bitvf/io.pyi``,
whose docstrings are the ones an IDE shows. Types of otswap in the
signatures are documented in `otswap's reference
<https://simosartori.github.io/otswap/python.html>`_.

bitvf
-----

Errors
~~~~~~

.. autoapiexception:: bitvf.Error

Configuration
~~~~~~~~~~~~~

.. autoapiclass:: bitvf.LengthUnit
   :members:

.. autoapiclass:: bitvf.Smoothing
   :members:

.. autoapiclass:: bitvf.Config
   :members:

.. autoapiclass:: bitvf.BoxConfig
   :members:

.. autoapiclass:: bitvf.BoxRsd
   :members:

.. autoapiclass:: bitvf.LightconeConfig
   :members:

.. autoapiclass:: bitvf.LightconeRsd
   :members:

.. autoapiclass:: bitvf.SecondPassRandoms
   :members:

Void finder, box
~~~~~~~~~~~~~~~~

.. autoapifunction:: bitvf.find_voids_box

.. autoapifunction:: bitvf.find_voids_box_from_divergence

.. autoapifunction:: bitvf.find_voids_box_from_field

.. autoapifunction:: bitvf.box_grid

.. autoapifunction:: bitvf.watershed_box

Void finder, lightcone
~~~~~~~~~~~~~~~~~~~~~~

.. autoapifunction:: bitvf.find_voids_lightcone

.. autoapifunction:: bitvf.find_voids_lightcone_from_cartesian

.. autoapifunction:: bitvf.lightcone_grid

Pipeline, box
~~~~~~~~~~~~~

.. autoapifunction:: bitvf.run_box

Pipeline, lightcone
~~~~~~~~~~~~~~~~~~~

.. autoapifunction:: bitvf.run_lightcone

Results
~~~~~~~

.. autoapiclass:: bitvf.Result
   :members:

.. autoapiclass:: bitvf.BasinCatalog
   :members:

.. autoapiclass:: bitvf.DivergenceGrid
   :members:

.. autoapiclass:: bitvf.TracerAssignment
   :members:

.. autoapiclass:: bitvf.Counts
   :members:

.. autoapiclass:: bitvf.Partition
   :members:

.. autoapiclass:: bitvf.BitVFResult
   :members:

Aliases
~~~~~~~

.. autoapidata:: bitvf.Verbosity

.. autoapidata:: bitvf.AngleUnit

.. autoapidata:: bitvf.RedshiftCutArg

.. autoapidata:: bitvf.Box

bitvf.io
--------

The columns and keywords of each file are on the :doc:`lightcone/formats` and
:doc:`box/formats` pages.

Writers
~~~~~~~

.. autoapiclass:: bitvf.io.Product
   :members:

.. autoapifunction:: bitvf.io.file_name

.. autoapifunction:: bitvf.io.write_voids

.. autoapifunction:: bitvf.io.write_overdense

.. autoapifunction:: bitvf.io.write_divergence

.. autoapifunction:: bitvf.io.write_galaxies

Readers
~~~~~~~

.. autoapifunction:: bitvf.io.read_divergence

.. autoapifunction:: bitvf.io.read_field

.. autoapiclass:: bitvf.io.Field
   :members:
