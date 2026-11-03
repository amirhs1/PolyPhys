# Density: local number density and volume fraction

This page explains what {mod}`polyphys.analyze.density` computes and why it
computes it the way it does: the local number density of a species from a
histogram of its bead centres, and the local volume fraction, which has to
account for a bead of finite size reaching into the bins next to the one
that holds its centre. Every routine is a numpy or scipy kernel on arrays;
reading frames, binning positions, averaging over frames, and placing the
result in the artifact lineage are the caller's job.

## Notation

$x$
: The binned coordinate: the axial position $z$ of a cylinder, or the
  distance $r$ from the axis of a cylinder or from a point in a cube.
  Lengths are in the monomer diameter.

$e_0 < e_1 < \dots < e_B$
: The bin edges, `edges`; bin $i$ is $[e_i, e_{i+1})$ and has width
  $h_i$. There are $B$ bins, `n_bins`.

$w(x)\,dx$, $V_i$
: The measure of the coordinate, $A\,dx$ for slabs of cross-section $A$,
  $2\pi L x\,dx$ for cylindrical shells of length $L$, $4\pi x^2\,dx$ for
  spherical shells; and the bin volume
  $V_i = \int_{e_i}^{e_{i+1}} w(x)\,dx$, `volumes`.

$H_i$, $F$
: The bead centres counted in bin $i$ over $F$ frames, `counts` and
  `n_frames`.

$\rho_i$
: The local number density of bin $i$, `rho`, in beads per cubic length
  unit.

$a$, $R_b = a/2$, $v = \pi a^3 / 6$
: The bead diameter, `bead_diameter`, its radius, and its volume.

$s(x \to j)$
: The *share kernel*: the volume of one bead whose centre is at coordinate
  $x$ that lies inside bin $j$. Summed over $j$ it is $v$ whenever the bead
  is inside the domain.

$S_{ij}$
: The *share matrix*, `share`, from {func}`~polyphys.analyze.density.share_matrix`:
  the volume of a test bead centred in bin $i$ that lies in bin $j$.

$\phi_i$
: The local volume fraction of bin $i$, `phi`, dimensionless.

## The number density

Counting centres per bin volume,

$$
\rho_i = \frac{H_i}{F\, V_i},
$$

is {func}`~polyphys.analyze.density.number_density`. Summed against the
volumes it returns the number of beads per frame,
$\sum_i \rho_i V_i = N$, which is the first thing to check on a new profile;
pass `n_particles` and the function checks that no bead was lost when the
positions were binned. The volumes come from
{func}`~polyphys.analyze.density.bin_volumes` in closed form, where the
research code integrated a polynomial numerically:

| Scheme | `(geometry, direction)` | $V_i$ |
| --- | --- | --- |
| spherical shell | `('cubic', 'r')` | $\tfrac{4\pi}{3}(e_{i+1}^3 - e_i^3)$ |
| cylindrical shell | `('cylindrical', 'r')` | $\pi L\,(e_{i+1}^2 - e_i^2)$ |
| slab | `('cylindrical', 'z')` | $\tfrac{\pi}{4} D^2\,(e_{i+1} - e_i)$ |

## The volume fraction

### What is being approximated

The volume fraction field is the indicator of "inside some bead",
$\Phi(\mathbf r) = \sum_n \Theta(R_b - |\mathbf r - \mathbf r_n|)$. Its
exact bin average over one frame needs every bead position,

$$
\phi_i^{\text{exact}} = \frac{1}{V_i}\int_{\text{bin } i}\Phi\,d^3r
= \frac{1}{V_i}\sum_n s(x_n \to i),
$$

but the probe keeps only the counts $H_j$, so the beads of a bin are treated
alike: the density is piecewise constant on the bins. In the continuum,
$\phi(\mathbf r) = \int \rho(\mathbf r')\,\Theta(R_b - |\mathbf r -
\mathbf r'|)\,d^3r'$, the density convolved with a ball, and with
$\rho$ piecewise constant the value at a point $x$ is
$\sum_j \rho_j\, s(x \to j)$.

### The centred form, `shares='centre'`

The research code evaluated that convolution at the bin centre
$c_i = \tfrac{1}{2}(e_i + e_{i+1})$,

$$
\phi_i^{\text{centre}} = \sum_j \rho_j\, s(c_i \to j),
$$

which is eq. 7 of the thesis appendix "Discrete local distributions in the
cylindrical coordinate system". It returns $\rho v$ exactly for a uniform
fluid, because the row $\sum_j s(c_i \to j) = v$. It does not conserve the
total bead volume along $r$: $\sum_i V_i \phi_i$ differs from $N v$ because
the point value at $c_i$ is not the bin average, by about 2 % of a bead
volume at a bin width $h = a/4$ and 0.5 % at $a/8$, second order in $h$.
Along $z$ with equal slabs the two coincide. This mode is kept so that the
research code's numbers can be reproduced.

### The averaged form, `shares='average'`

The default averages the kernel over the test bead's position in its bin,
with the bin's own measure:

$$
\bar s(i \to j) = \frac{1}{V_i}\int_{e_i}^{e_{i+1}} s(x \to j)\, w(x)\,dx,
\qquad
\phi_i = \sum_j \bar s(i \to j)\, \rho_j .
$$

This is exact under the one assumption already made for $\rho$, that the
beads of a bin are spread uniformly over it. Writing the kernel as an
integral of the indicator, $s(x \to j) = \int_{e_j}^{e_{j+1}} \chi(x, y)\,
w(y)\,dy$ with $\chi(x, y)$ the fraction of the surface at $y$ inside a
bead centred at $x$,

$$
V_i\, \bar s(i \to j) = \int_{e_i}^{e_{i+1}}\!\!\int_{e_j}^{e_{j+1}}
\chi(x, y)\, w(x)\, w(y)\, dy\, dx ,
$$

and since the indicator $\Theta(R_b - |\mathbf r - \mathbf r'|)$ is
symmetric in its two points, so is the integrand, hence

$$
V_i\, \bar s(i \to j) = V_j\, \bar s(j \to i).
$$

Two properties follow from this identity alone, and the tests check both to
round-off in every scheme:

1. **Conservation.** $\sum_i V_i \phi_i = \sum_j \rho_j \sum_i V_i \bar
   s(i \to j) = \sum_j \rho_j V_j \sum_i \bar s(j \to i) = \sum_j (H_j / F)\,
   v$: the beads' total volume comes back.
2. **Uniform fluid.** $\phi_i = \rho \sum_j \bar s(i \to j) = \rho v$.

The averaged form is also the exact bin average of the centred form, so it
is the research code's formula with its one remaining approximation
removed. The plan's scatter form, which moves each bead to its bin centre
and spreads it ($\phi_i = \tfrac{1}{V_i}\sum_j (H_j / F)\, s(c_j \to i)$),
conserves volume but misses the uniform fluid by the same 2 % at $h = a/4$;
it is the other half of the same missing average, and it is not offered.

## The kernels

### Slabs

The volume of a bead below the plane at signed height $u$ from its centre
is the spherical segment $G(u) = \int_{-R_b}^{u} \pi (R_b^2 - t^2)\,dt$
(MathWorld, "Spherical Segment", eq. 14, as
{func}`~polyphys.analyze.measurer.spherical_segment` computes it), and the
share of slab $j$ for a bead at $x$ is $G(e_{j+1} - x) - G(e_j - x)$. The
average over the bead's position is then a second antiderivative,
$G_2(u) = \int_{-\infty}^{u} G$:

$$
G_2(u) = \begin{cases}
0 & u \le -R_b, \\[2pt]
\pi\!\left[\dfrac{R_b^2 (u + R_b)^2}{2} - \dfrac{u^4 - R_b^4}{12}
- \dfrac{R_b^3 (u + R_b)}{3}\right] & |u| \le R_b, \\[8pt]
v\,u & u \ge R_b,
\end{cases}
$$

$$
\bar s(i \to j) = \frac{1}{h_i}\Big[G_2(e_{j+1} - e_i) - G_2(e_{j+1} - e_{i+1})
- G_2(e_j - e_i) + G_2(e_j - e_{i+1})\Big].
$$

The axial direction of a cylinder is periodic: the slabs at one end are
also reached, through their images shifted by the domain length, from the
other end, and the same bracket is summed over the images. For a bead of
diameter 1 in slabs of width 1 the shares are $13\pi/96$ into the bead's
own slab and $\pi/64$ into each neighbour.

### Spherical shells

The volume of a bead centred at distance $x$ from the origin that lies
inside the ball of radius $e$ is the sphere-sphere intersection volume
$W(R_b, e, x)$ (MathWorld, "Sphere-Sphere Intersection", eq. 16, as
{func}`~polyphys.analyze.measurer.sphere_sphere_intersection` computes it),
and the share of shell $j$ is $W(R_b, e_{j+1}, x) - W(R_b, e_j, x)$. In
the lens regime $x\,W$ is a quartic in $x$, so the position average with
weight $4\pi x^2$ is a polynomial on each piece between the breakpoints
$|e - R_b|$ and $e + R_b$ of the two edges, and a Gauss-Legendre rule on
each piece integrates it exactly.

### Cylindrical shells

The kernel is the volume $C(R_b, R, d)$ of a bead whose centre is at axial
distance $d$ inside the coaxial cylinder of radius $R$,
{func}`~polyphys.analyze.density.sphere_cylinder_intersection`. Slicing the
bead by planes normal to the axis, the slice at height $z$ is a disc of
radius $\sqrt{R_b^2 - z^2}$ whose centre is $d$ from the axis, the
cylinder cuts the plane in a disc of radius $R$ about the axis, and the
slice's share is their lens area (MathWorld, "Circle-Circle Intersection",
eq. 15):

$$
C(R_b, R, d) = \int_{-R_b}^{R_b} A\!\left(\sqrt{R_b^2 - z^2},\, R,\, d\right) dz,
$$

$$
A(r, R, d) = r^2 \cos^{-1}\!\frac{d^2 + r^2 - R^2}{2 d r}
+ R^2 \cos^{-1}\!\frac{d^2 + R^2 - r^2}{2 d R}
- \tfrac{1}{2}\sqrt{(-d + r + R)(d + r - R)(d - r + R)(d + r + R)}
$$

for $|R - r| < d < r + R$, with $A = 0$ when the discs are apart and
$\pi \min(r, R)^2$ when one holds the other. The $z$-integral has no closed
form; it is a 32-node Gauss-Legendre rule on each of the three pieces
between the heights where the lens changes regime, and it agrees with an
adaptive quadrature to about $10^{-8}$ of the bead volume. Two limits are
exact: the bead is inside for $R \ge d + R_b$ and outside for
$R \le d - R_b$; on the axis with $R < R_b$ the volume is
$\tfrac{4\pi}{3}[R_b^3 - (R_b^2 - R^2)^{3/2}]$.

The research code used the sphere-sphere volume $W(R_b, R, d)$ in place of
$C$ "for the sake of simplicity". That is `intersection='sphere'`, kept for
reproducing its numbers. Its error, as a fraction of the bead volume and at
the worst axial distance, is about $0.09\,a / R$:

| $R / a$ | 1 | 2 | 5 | 10 | 20 | 40 |
| --- | --- | --- | --- | --- | --- | --- |
| $\max_d \lvert W - C \rvert / v$ | 0.12 | 0.049 | 0.019 | 0.0094 | 0.0047 | 0.0023 |

It is first order in $a / R$ and independent of the bin width, so finer
bins never remove it, and it is largest at the innermost edges, where
`sumrule` takes the `'first'` normalizer of a radial profile. The shares
still telescope to $v$ over all shells, so the error moves volume between
neighbouring shells rather than losing it; but the symmetry
$V_i S_{ij} = V_j S_{ji}$ needs the true indicator, so with the sphere
kernel the conservation of the averaged form holds only to the same order.
(A planar cut, the "strips" of the maintainer's 2021 report, errs by the
same amount with the opposite sign; the mean of the two is second order in
$a / R$, because a cylinder's mean curvature is half a sphere's. The exact
kernel costs well under a second per space, so neither shortcut is
offered.)

## Conventions the kernels leave to the caller

- **The slab diameter.** {func}`~polyphys.analyze.density.bin_volumes` uses
  the diameter it is given. {func}`~polyphys.analyze.density.density_profiles`
  passes $D_\mathrm{eff} = D - a$, the diameter of the region a bead
  centre can reach, which is how the parsers define the bulk values,
  $\rho_c = N_c / (\pi D_{\mathrm{eff},c}^2 L / 4)$ and $\phi_c = \rho_c
  v_c$ (thesis, chapter 3). A uniform crowder fluid then returns $\phi_c$
  exactly, so the scaled profiles of `sumrule` compare with the nominal
  bulk value. The research code used the full $D$ for its slabs, so its
  far-field $z$ values were $((D - a)/D)^2$ times the bulk; its published
  profiles were normalized by their far-field value, which cancels the
  factor.
- **Half the cube.** A spherical shell beyond half the side of a periodic
  cube overstates the volume inside the box, so the density there is
  underestimated. Given `lcube`, `bin_volumes` warns when the outermost
  edge passes half the side by more than one bin, as MAICoS's
  `DensitySphere` does for its `rmax`, and computes the volumes regardless.
- **Frames.** `n_frames` divides the summed counts, so `rho` is a
  per-frame density comparable with $N / V$; the research code did not
  divide.
- **Normalization.** `normalize='sum'` divides by the sum over bins, the
  research code's `normalized=True`; `'integral'` divides by
  $\sum_i p_i V_i$. Profiles that feed {mod}`polyphys.analyze.sumrule` must
  stay unnormalized: its species sum $\sum_\alpha \phi_\alpha / a_\alpha$
  needs each profile's physical scale, and only its own `'first'` and
  `'max'` normalizations are scale-free.
  {func}`~polyphys.analyze.sumrule.space_profiles` warns, naming the file,
  when a profile's `-mean` sums to one, the signature of `normalize='sum'`;
  it cannot tell `'integral'` or a hand normalization from one profile, so
  the rule above is the caller's contract.

## Cost and use

The share matrix depends on the edges and the bead diameter only. Building
it is $O(B k)$ kernel evaluations, with $k \approx a / h + 2$ the band
width (for 300 cylindrical shells with the exact kernel, under a second;
for 6000 periodic slabs, milliseconds), and it is returned as a sparse
`csr_array` of $O(B k)$ entries. Applying it is one sparse product per
profile. So build it once per space and reuse it for every whole, frame, and
species of the same diameter:

```python
import numpy as np
from polyphys.analyze import density

edges = np.arange(0.0, 15.01, 0.1)
volumes = density.bin_volumes(edges, 'cylindrical', 'r', lcyl=300.0)
share = density.share_matrix(edges, 1.0, 'cylindrical', 'r')
for counts in histograms:  # one summed histogram per whole
    rho = density.number_density(counts, volumes, n_frames=14001)
    phi = density.volume_fraction(rho, share)
```

Per-frame counts of shape `(n_frames, n_bins)` go through the same two
calls row by row, which is what a per-bin error estimate at the probe step
will need. `density_profiles` wraps the four calls for one profile.

In the artifact lineage, `rho` and `phi` are whole-level arrays named
`<direction>Rho<Species>` and `<direction>Phi<Species>`; the organizer
averages them over replicates into the `ensAvg` profiles that
{func}`~polyphys.analyze.sumrule.space_profiles` reads. The sum rule is
linear in the profiles, so averaging over wholes and summing over species
commute; its normalization is not, and the published estimator normalizes
the ensemble average.

## References

- Weisstein, Eric W. "Spherical Segment", "Sphere-Sphere Intersection",
  and "Circle-Circle Intersection". From MathWorld--A Wolfram Web
  Resource. <https://mathworld.wolfram.com/>
- A. H. Sadeghi Isfahani, PhD thesis, University of Waterloo (2024),
  chapter 3 and the appendix "Discrete local distributions in the
  cylindrical coordinate system".
- C. Jeon, C. Hyeon, Y. Jung, and B.-Y. Ha, "How are molecular crowding and
  the spatial organization of a biopolymer interrelated", Soft Matter 12,
  9786-9796 (2016). <https://doi.org/10.1039/c6sm01924b>
- MAICoS, `DensitySphere`, for the half-box warning.
  <https://maicos-analysis.org/>
