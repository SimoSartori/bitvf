# The method

BitVF, the Back-in-time Void Finder (Sartori et al. 2026), is a dynamical void
finder that identifies cosmic voids as non-overlapping regions of coherent mass
outflow. It reconstructs the Lagrangian displacement field, which moves each
tracer back in time from its observed position to its initial position in a
nearly uniform distribution, and measures the outflow through the divergence
of this field, which is negative in regions that matter has flowed out of,
since going back in time the displacements converge into them. BitVF computes
the divergence on a grid and divides it into basins around its local minima:
the basins whose minimum is negative are the voids, and the remaining ones are
overdense regions.

The displacement field is reconstructed by solving a discrete optimal
transport problem between the tracers and a uniform distribution of points,
with the [otswap](https://simosartori.github.io/otswap/) package. Apart from
the cosmological parameters needed to convert redshifts into distances, the
reconstruction makes no cosmological assumption, and it accounts for the bias
of the tracers by construction. Because the displacements are shaped by the
motion of all the tracers rather than by local tracer counts, the voids found
by BitVF have been shown to be less sensitive to shot noise, tracer
subsampling and redshift-space distortions than those of topological void
finders. The residual redshift-space distortions can be further corrected
with a fiducial cosmology and a value of the tracer bias, as described below.

BitVF works with two geometries:

- **Box geometry**, for simulation snapshots: the tracers fill a rectangular
  volume in comoving Cartesian coordinates, and the uniform distribution used
  by the reconstruction is drawn inside it.
- **Lightcone geometry**, for galaxy surveys and mock catalogues: it handles
  the observational setup of real data, with tracers in sky coordinates and
  redshift, a random catalogue sampling the survey volume and angular
  systematics, and an angular mask describing the footprint.

The method runs in six steps:

1. [Reconstruction](#1-reconstruction): the displacement of every tracer is
   computed by optimal transport.
2. [Divergence on a grid](#2-divergence-on-a-grid): the divergence of the
   displacements is estimated in each cell of a grid.
3. [Smoothing](#3-smoothing): the divergence is smoothed with a Gaussian.
4. [Basins](#4-basins): the grid is divided into basins around the local
   minima of the divergence.
5. [Voids and overdense regions](#5-voids-and-overdense-regions): the centre,
   size and density of each basin are measured.
6. [Correction of redshift-space distortions](#6-correction-of-redshift-space-distortions):
   optionally, the tracers are moved to real space and reconstructed again.

Each step ends with the parameters that control it, which are fields of the
configuration structures of the [Python](python.rst) and [C++](cpp.rst)
libraries ([names side by side](naming.md) maps one to the other), and with
a paragraph on what changes in the lightcone geometry.

## 1. Reconstruction

The reconstruction, performed by otswap, pairs each tracer one to one with a
point of a uniform distribution, a *random*, so that the sum of the squared
comoving distances between paired points is as small as possible; the
vector from a tracer to its random is the tracer's displacement. The
pairing is repeated on independent sets of randoms, the *realizations*, and
BitVF uses all of them, so that more realizations give a less noisy
divergence. The reconstruction also yields the *mean particle separation*
(mps), the cube root of the volume per tracer, which sets the scale of every
later step. A search can also start from displacements computed earlier or
by another code.

In the box geometry the randoms are optional: when none are given, as many
randoms as tracers are drawn uniformly in the box for each realization. The
box is the one given by the user or, otherwise, the smallest box containing
the tracers and the randoms. A tracer outside a given box keeps its row in
the outputs, with a flag, but takes no part in the search.

**Parameters**

- {py:attr}`BoxConfig.ot <bitvf.BoxConfig.ot>`,
  {py:attr}`LightconeConfig.ot <bitvf.LightconeConfig.ot>`
  ({cpp:member}`BoxConfig::ot <bitvf::BoxConfig::ot>`,
  {cpp:member}`LightconeConfig::ot <bitvf::LightconeConfig::ot>`): the
  settings of the reconstruction, among them the number of realizations,
  the convergence threshold and the seed. otswap's defaults are a reasonable
  start; its
  [guide](https://simosartori.github.io/otswap/guide.html#the-parameters-that-matter)
  explains them.
- {py:attr}`BoxConfig.limits <bitvf.BoxConfig.limits>`
  ({cpp:member}`BoxConfig::limits <bitvf::BoxConfig::limits>`): the box;
  for a simulation snapshot, give the box of the simulation.
- {py:attr}`BoxConfig.mps <bitvf.BoxConfig.mps>`
  ({cpp:member}`BoxConfig::mps <bitvf::BoxConfig::mps>`): the mean particle
  separation. When it is left empty, the default, otswap computes it.
- {py:class}`otswap.DistanceTable` ({cpp:class}`otswap::DistanceTable`),
  lightcone only: the cosmological parameters that convert redshifts into
  comoving distances, given as a distance table; use the fiducial cosmology
  of the analysis.
- {py:attr}`LightconeConfig.n_bins <bitvf.LightconeConfig.n_bins>`
  ({cpp:member}`LightconeConfig::nBins <bitvf::LightconeConfig::nBins>`),
  lightcone only: the number of redshift bins of the mps(z) profile. It has
  no default; 30 is a reasonable start, and a small sample needs fewer bins,
  since each must hold enough tracers.
- {py:attr}`LightconeConfig.cut <bitvf.LightconeConfig.cut>`
  ({cpp:member}`LightconeConfig::cut <bitvf::LightconeConfig::cut>`),
  lightcone only: the range of observed redshift; by default nothing is
  cut.
- {py:attr}`otswap.Config.reject_crossings`
  ({cpp:member}`otswap::Config::rejectCrossings`) and
  {py:attr}`otswap.Config.max_unobserved_pixels_crossed`
  ({cpp:member}`otswap::Config::maxUnobservedPixelsCrossed`), lightcone
  only: whether to reject the displacements whose path on the sky crosses
  unobserved pixels of the mask, and how many consecutive unobserved pixels
  a path may cross. The defaults, rejection with no pixel allowed, are a
  reasonable start.

**In the lightcone geometry**, the randoms are required, because they are
the only description of the survey volume and its systematics. They must
number at least the tracers times the realizations; 30 randoms per tracer,
the minimum with otswap's default of 30 realizations, is a reasonable
start. The mps depends on redshift and is measured
as a profile mps(z). The tracers and randoms on unobserved pixels of the
mask, or with an observed redshift outside the cut, are left out; the
tracers keep their rows, with a flag.

## 2. Divergence on a grid

Let $\boldsymbol{\Psi}$ be the Lagrangian displacement, which carries
matter forward in time, so that the reconstruction measures the
back-in-time displacement $\mathbf{d} = -\boldsymbol{\Psi}$. At linear
order, the density contrast is

$$
\delta \simeq -\nabla \cdot \boldsymbol{\Psi} = \nabla \cdot \mathbf{d} ,
$$

so voids, where $\delta < 0$, are regions where the divergence of
$\mathbf{d}$ is negative. BitVF estimates this divergence on a grid whose
cell side, the *cell size*, is a multiple of the mps. Smaller cells resolve
smaller voids but hold fewer displacements, so the divergence is noisier.
In the box geometry the cells are cubes covering the box, and the grid is
not periodic, even for a periodic simulation.

A displacement that leaves the cell of its tracer contributes to two cells.
To the cell it leaves, it adds a positive amount: its component along the
outward normal $\hat{\mathbf{n}}$ of the face through which it exits,
divided by the cell side $\ell$. To the cell where it ends, it adds a
negative amount, computed in the same way at the face through which it
enters. The divergence of a cell $c$ is the mean of these contributions over
the $W_c$ displacement ends counted in the cell, for all the realizations,

$$
\Theta_c = \frac{1}{W_c} \sum_{e \in c} \frac{\hat{\mathbf{n}}_e \cdot \mathbf{d}_e}{\ell_e} .
$$

A displacement that ends in its starting cell is counted in $W_c$ but adds
nothing to the sum. The *weight* of the cell, $W_c/R$ with $R$ the number
of realizations, measures how many displacements support its value. A cell
that no displacement reaches, a *filled cell*, takes its value from its
neighbours and keeps a weight of 0.

**Parameters**

- {py:attr}`Config.divergence_cell_size <bitvf.Config.divergence_cell_size>`
  ({cpp:member}`Config::divergenceCellSize <bitvf::Config::divergenceCellSize>`):
  the cell size, in units of the mps. The default,
  $\min\left(1, (20/R)^{1/3}\right)$, shrinks the cells when more
  realizations fill them, and is a reasonable start.

**In the lightcone geometry**, the cells are shells in comoving distance
from the observer, divided into rings of declination and sectors of right
ascension. Each shell is one cell size times mps(z) wide, so the cells grow
with redshift. A cell whose centre falls on an unobserved pixel of the
mask is *masked*: it has no divergence and belongs to no basin, so a smaller
cell size follows the footprint more closely. A displacement that starts or
ends in a masked cell, or ends outside the grid, is not counted. A cell that
the mask isolates can remain empty, and then has a divergence of 0.

## 3. Smoothing

Because each cell holds a finite number of displacements, the divergence
is noisy, so it is smoothed with a Gaussian of width $\sigma$: each cell takes
the mean of the divergence of the cells around it, weighted by their weight
times $\exp(-r^2/2\sigma^2)$, where $r$ is the distance between the cell
centres. Filled cells, of weight 0, do not contribute. A larger $\sigma$ merges small basins into larger ones.
The smoothed weight is kept beside the smoothed divergence: where it is
close to 0, the smoothed value comes mostly from filled cells and should be
read with care.

**Parameters**

- {py:attr}`Smoothing.sigma <bitvf.Smoothing.sigma>`
  ({cpp:member}`Smoothing::sigma <bitvf::Smoothing::sigma>`): the width of
  the Gaussian. 1 mps is a reasonable start, and a width of 0 turns the
  smoothing off.
- {py:attr}`Smoothing.unit <bitvf.Smoothing.unit>`
  ({cpp:member}`Smoothing::unit <bitvf::Smoothing::unit>`): the unit of the
  width, either the mps, the default, or the units of the coordinates.

**In the lightcone geometry**, masked cells are left out of every mean, and
a width given in mps is converted with the value of mps(z) at the shell of
the cell being smoothed.

## 4. Basins

The smoothed divergence, or the divergence itself when the smoothing is
off, is divided into basins by a *watershed*. Each cell is linked to the
neighbouring cell with the lowest divergence, provided that value is lower
than its own. A cell with no lower neighbour is a local minimum and becomes
the *seed* of a basin, which contains every cell whose chain of links leads
to that seed. Every unmasked cell belongs to exactly one basin. A basin whose seed
has a negative divergence is a void; the others are overdense basins. In
the box geometry, the neighbours of a cell are the 26 cells that share a
face, an edge or a corner with it.

**Parameters**: none. The number and size of the basins follow from the
cell size and the smoothing.

**In the lightcone geometry**, the cells are not aligned, so neighbours are
found from the distance between cell centres, a rule that selects the same
26 cells on a lattice of equal cubes. Masked cells belong to no basin.

## 5. Voids and overdense regions

The *core* of a void is the part of its basin where the divergence is
negative, and the *extended void* is the whole basin. For an overdense
basin, the core and the extended basin coincide. The centre of a basin
starts at the centre $\mathbf{x}_0$ of its seed cell and is moved according
to the divergence $\Theta_n$ of the neighbouring cells $n$,

$$
\mathbf{x}_c = \mathbf{x}_0 + s \, \frac{\sum_n \Theta_n (\mathbf{x}_n - \mathbf{x}_0)}{\sum_n |\Theta_n|} ,
$$

where $s = -1$ for a void, whose centre moves towards lower divergence, and
$s = +1$ for an overdense basin, whose centre moves towards higher
divergence. The centre is kept inside the seed cell. Each basin is then
described by:

- `radius` and `radiusExt`, the radii of the spheres with the same volume
  as the core and as the extended basin;
- `deltaAv`, the mean density contrast of the core: the number of tracers
  that start in it (`nTracers`) divided by its volume and by the mean
  density $\mathrm{mps}^{-3}$, minus 1;
- `delta2Mps`, the density contrast within a sphere of radius 2 mps around
  the centre;
- `edgeFlag`, which is 0 when the basin does not touch the boundary of the
  grid, 1 when only the extended basin touches it, and 2 when the core
  does;
- `contamination`, the fraction of the basin's surface on the boundary.

Each tracer is assigned to the basin of the cell in which it starts.

**Parameters**: none.

**In the lightcone geometry**, the centre is also given in sky coordinates
and redshift, and the densities use mps(z) at the centre. Because a sphere
near the edge of the survey lies partly outside it, `delta2Mps` also counts
the randoms in the sphere: it is $n_t/(n_r N_t/N_r) - 1$, where $n_t$
and $n_r$ are the numbers of tracers and randoms in the sphere and
$N_t/N_r$ is the ratio of tracers to randoms over the whole search. When no
random lies in the sphere, `delta2Mps` is undefined and set to NaN. The
masked cells are part of the boundary of the grid, so a basin that touches
the edge of the footprint, the redshift limits or a hole in the mask is
flagged. A tracer that starts in a masked cell belongs to no basin.

## 6. Correction of redshift-space distortions

This step is optional. The redshift of a tracer includes its peculiar
velocity along the line of sight, so its *redshift-space* position differs
from its real-space one. otswap estimates the velocity of each tracer from
its displacements, averaged over its neighbours, using linear theory with
the growth rate of a fiducial cosmology and the bias of the tracers, and
moves the tracer to real space; otswap's
[method page](https://simosartori.github.io/otswap/method.html#redshift-space-correction)
derives the correction. The corrected tracers are then reconstructed a
second time with the same seed, and the voids are searched on this second
reconstruction, on the grid of the first and with the same smoothing and
mps. The voids of the first, redshift-space, reconstruction can be searched
as well. A tracer that the correction does not move, or moves outside the
grid, keeps its row with a flag.

In the box geometry, the tracers are moved along one axis, the line of
sight, with the growth rate at the redshift of the box and a constant bias.
The second reconstruction always reuses the randoms of the first, whether
they were given or drawn, since in a box the uniform distribution does not
change.

**Parameters**

- {py:class}`bitvf.BoxRsd` ({cpp:struct}`bitvf::BoxRsd`): the axis of the
  line of sight, the redshift of the box and the bias, all set by the
  simulation, and the width of the average over neighbouring tracers, for
  which 10 Mpc/h is a reasonable start, to be checked in each analysis. The
  growth rate comes from the distance table given to the run.
- {py:class}`bitvf.LightconeRsd` ({cpp:struct}`bitvf::LightconeRsd`): the
  bias of the tracers as a table in redshift, and the width of the average
  over neighbouring tracers, chosen as in the box.
- {py:attr}`LightconeRsd.second_pass_randoms <bitvf.LightconeRsd.second_pass_randoms>`
  ({cpp:member}`LightconeRsd::secondPassRandoms <bitvf::LightconeRsd::secondPassRandoms>`):
  the randoms of the second reconstruction, described below; the default is
  a reasonable start.
- {py:attr}`BoxRsd.zspace_voids <bitvf.BoxRsd.zspace_voids>`,
  {py:attr}`LightconeRsd.zspace_voids <bitvf.LightconeRsd.zspace_voids>`
  ({cpp:member}`BoxRsd::zspaceVoids <bitvf::BoxRsd::zspaceVoids>`,
  {cpp:member}`LightconeRsd::zspaceVoids <bitvf::LightconeRsd::zspaceVoids>`):
  whether to search the redshift-space reconstruction as well; it is off by
  default.

**In the lightcone geometry**, each tracer is moved along its own line of
sight, in redshift only, with a bias that depends on redshift. A tracer
that the correction moves outside the redshift cut is flagged as well. The
second reconstruction needs randoms that follow the real-space sample, and
three choices are available. By default, the randoms of the first
reconstruction keep their sky positions and take new redshifts, drawn from
those of the corrected tracers. Alternatively, the randoms of the first
reconstruction can be reused unchanged, keeping the observed redshift
distribution, or the user can provide randoms built for the real-space
sample.

## Products

A search returns the catalogues of the voids and of the overdense basins,
with the properties of step 5, and counts of the tracers used and left out.
On request it also returns the divergence grid and the assignment of every
tracer to its cell and basin; a run with the correction can also return the
redshift-space voids, both reconstructions and the corrected tracers. Their
files are described on the format pages of the
[lightcone](lightcone/formats.md) and [box](box/formats.md) geometries.

## References

- Sartori, S., et al. 2026, *The Back-in-time Void Finder: dynamical
  identification of cosmic voids through optimal transport reconstruction*,
  [arXiv:2601.15378](https://arxiv.org/abs/2601.15378)
- otswap, the optimal transport reconstruction used by BitVF: its
  [documentation](https://simosartori.github.io/otswap/) and the
  [references](https://simosartori.github.io/otswap/method.html#references)
  of its method page.
