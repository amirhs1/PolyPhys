# A FENE ring of 8 beads in a periodic cube, in the atom style given by
# -var style (bond, angle, molecular, or full). The ring starts across the
# x = 0 face, so half its beads are remapped with ix = -1.
# Writes ring.<style>.lammpstrj (3 frames: steps 0, 500, 1000) and
# ring.<style>.data (write_data at step 1000, the last frame).
# With style bond it also writes the dump formats in in.formats.

units           lj
dimension       3
boundary        p p p
atom_style      ${style}

variable        randseed equal 7211
variable        topo string "bond/types 1"
variable        molfile string ring8.mol
if "${style} != bond" then &
    "variable topo string 'bond/types 1 angle/types 1 extra/angle/per/atom 3'" &
    "variable molfile string ring8_angles.mol"

region          box block 0 6 0 6 0 6 units box
create_box      1 box ${topo} extra/bond/per/atom 2 extra/special/per/atom 6
molecule        ring ${molfile}
create_atoms    0 single 0.0 3.0 3.0 mol ring 4127 remap yes units box
mass            1 1.0
if "${style} == full" then "set type 1 charge 0.25"

pair_style      lj/cut 1.122462
pair_modify     shift yes
pair_coeff      1 1 1.0 1.0 1.122462
bond_style      fene
bond_coeff      1 30.0 1.5 1.0 1.0
if "${style} != bond" then "angle_style harmonic" "angle_coeff 1 2.0 180.0"
special_bonds   fene
neighbor        0.3 bin
neigh_modify    delay 0 every 1 check yes

velocity        all create 1.0 ${randseed}
fix             F_langevin all langevin 1.0 1.0 1.0 ${randseed}
fix             F_nve all nve
timestep        0.005
thermo          500

dump            D_x all custom 500 ring.${style}.lammpstrj id mol type x y z ix iy iz
dump_modify     D_x format float %.16g
if "${style} == bond" then "include in.formats"
run             1000
write_data      ring.${style}.data
