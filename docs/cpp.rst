C++ reference
=============

Generated from the public headers ``bitvf/VoidFinder.h``, ``bitvf/BitVF.h``
and ``bitvf/io.h`` by Doxygen, with only what carries a comment included.
``bitvf/BitVF.h`` includes ``bitvf/VoidFinder.h``; ``bitvf/io.h`` is included
separately. Types of otswap in the signatures are documented in `otswap's
reference <https://simosartori.github.io/otswap/cpp.html>`_.

bitvf/VoidFinder.h
------------------

Errors
~~~~~~

.. doxygenclass:: bitvf::Error
   :members:

Array layout
~~~~~~~~~~~~

.. doxygengroup:: array_layout
   :content-only:

Configuration
~~~~~~~~~~~~~

.. doxygenenum:: bitvf::LengthUnit

.. doxygenstruct:: bitvf::Smoothing
   :members:

.. doxygenstruct:: bitvf::Config
   :members:

Results
~~~~~~~

.. doxygenstruct:: bitvf::Result
   :members:

.. doxygenstruct:: bitvf::BasinCatalog
   :members:

.. doxygenstruct:: bitvf::DivergenceGrid
   :members:

.. doxygenstruct:: bitvf::TracerAssignment
   :members:

.. doxygenstruct:: bitvf::Counts
   :members:

Void finder, box
~~~~~~~~~~~~~~~~

.. doxygenfunction:: bitvf::findVoidsBox(const std::vector<double>&, const std::vector<double>&, const std::vector<std::uint8_t>&, double, const Config&, const std::optional<otswap::BoxLimits>&, const std::vector<std::int64_t>&, const DivergenceGrid*)

.. doxygenfunction:: bitvf::findVoidsBox(const otswap::Result&, const Config&, const DivergenceGrid*)

.. doxygenfunction:: bitvf::findVoidsBoxFromDivergence

.. doxygenfunction:: bitvf::findVoidsBoxFromField

Void finder, lightcone
~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: bitvf::findVoidsLightcone(const std::vector<double>&, const std::vector<double>&, const std::vector<std::uint8_t>&, const std::vector<double>&, const otswap::Mask&, const otswap::DistanceTable&, const otswap::MpsProfile&, const Config&, const otswap::RedshiftCut&, const std::vector<std::int64_t>&, const DivergenceGrid*)

.. doxygenfunction:: bitvf::findVoidsLightcone(const otswap::Result&, const std::vector<double>&, const otswap::Mask&, const Config&, const DivergenceGrid*)

.. doxygenfunction:: bitvf::findVoidsLightconeFromCartesian

Grid, lightcone
~~~~~~~~~~~~~~~

.. doxygenfunction:: bitvf::lightconeGrid(const otswap::Result&, const std::vector<double>&, const otswap::Mask&, std::optional<double>)

.. doxygenfunction:: bitvf::lightconeGrid(const std::vector<double>&, const std::vector<double>&, const otswap::Mask&, const otswap::DistanceTable&, const otswap::MpsProfile&, double, const otswap::RedshiftCut&)

Grid and watershed, box
~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: bitvf::boxGrid(const otswap::BoxLimits&, double, double)

.. doxygenfunction:: bitvf::boxGrid(const otswap::Result&, std::optional<double>)

.. doxygenfunction:: bitvf::watershedBox

.. doxygenstruct:: bitvf::Partition
   :members:

bitvf/BitVF.h
-------------

Configuration
~~~~~~~~~~~~~

.. doxygenstruct:: bitvf::BoxConfig
   :members:

.. doxygenstruct:: bitvf::BoxRsd
   :members:

.. doxygenstruct:: bitvf::LightconeConfig
   :members:

.. doxygenstruct:: bitvf::LightconeRsd
   :members:

.. doxygenenum:: bitvf::SecondPassRandoms

Pipeline, box
~~~~~~~~~~~~~

.. doxygenfunction:: bitvf::runBox(const std::vector<double>&, const std::vector<double>&, const BoxConfig&, const std::vector<std::int64_t>&)

.. doxygenfunction:: bitvf::runBox(const std::vector<double>&, const std::vector<double>&, const otswap::DistanceTable&, const BoxConfig&, const std::vector<std::int64_t>&)

Pipeline, lightcone
~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: bitvf::runLightcone

.. doxygenstruct:: bitvf::BitVFResult
   :members:

bitvf/io.h
----------

The columns and keywords of each file are on the :doc:`lightcone/formats` and
:doc:`box/formats` pages.

Writers
~~~~~~~

.. doxygenenum:: bitvf::io::Product

.. doxygenfunction:: bitvf::io::fileName

.. doxygenfunction:: bitvf::io::writeVoids

.. doxygenfunction:: bitvf::io::writeOverdense

.. doxygenfunction:: bitvf::io::writeDivergence

.. doxygenfunction:: bitvf::io::writeGalaxies

Readers
~~~~~~~

.. doxygenfunction:: bitvf::io::readDivergence

.. doxygenfunction:: bitvf::io::readField

.. doxygenstruct:: bitvf::io::Field
   :members:
