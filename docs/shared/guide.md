# Shared parts of the guides

Not a page of its own: the box and lightcone guides include the parts
between the markers.

<!-- results-start -->
Each search returns a result holding the catalogues of the voids and of the
overdense basins, with what it counted. Two larger products are kept only
when asked for, since they take memory: the divergence grid, with
`keepDivergence` (Python `keep_divergence`), and where each tracer ended up,
its cell, its basin and its distance from the basin's centre, with
`keepAssignment` (Python `keep_assignment`).
<!-- results-end -->

<!-- writing-start -->
A search writes no file. The writers of `bitvf::io` (Python `bitvf.io`)
write the four products of a result, each to its own file: `writeVoids` and
`writeOverdense`, one row per basin; `writeDivergence`, one row per cell,
which needs `keepDivergence`; and `writeGalaxies`, one row per tracer, which
needs `keepAssignment`. `fileName` (Python `file_name`) gives the
conventional names, such as `voids_<name>` and, for a redshift-space pass,
`voids_zspace_<name>`. The extension chooses the format: FITS for `.fits`,
`.fit` and `.fits.gz`,
ASCII otherwise. Each file records the main parameters of the search that
made it as keywords of its header, such as the smoothing and the cell size.
<!-- writing-end -->

<!-- reproducibility-start -->
Given the same inputs and a seed other than 0, every result is the same bit
for bit whatever the number of threads, on every supported platform and
compiler; only the elapsed times differ. A seed of 0, otswap's default, is
replaced by one drawn at the start of the run and recorded in otswap's
result (`otswap::Result::config`), so that any run can be repeated; setting
a seed of your own saves reading it back. The number of threads follows
`OMP_NUM_THREADS`.
<!-- reproducibility-end -->

<!-- prints-start -->
`Config::verbosity` (Python `verbosity`) sets the report, at otswap's
levels: `Silent` prints nothing; `Normal`, the default, prints what the
selection left out and one line per call with its time; `Detailed` also
prints one line per step of the search. Each line starts with `bitvf: ` and
is written as soon as its step finishes, to `Config::log` in C++
(`std::clog` by default), to `sys.stdout` in Python. Every count printed can
also be read from the result's `counts`. In Python, Ctrl+C stops a search
at its next report line.
<!-- prints-end -->

<!-- errors-start -->
Every failure raises `bitvf::Error` (Python `bitvf.Error`), with a message
that names the quantity and the rule broken; an error raised by otswap
inside a bitvf call has a message starting `otswap: `. Arguments are checked
before any work, so that a wrong value fails at once rather than after the
reconstruction. `bitvf::Error` is unrelated to `otswap::Error`: catching one
does not catch the other.
<!-- errors-end -->
