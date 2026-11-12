# How `ArgumentPrinciple` searches

`ArgumentPrinciple` refines a tree of cells over `(z, p...)` toward the loci of
zeros and poles, reading only the phase of `f` at cell corners. It then links the
surviving cells into branches and labels each one by an argument-principle winding.
This page describes those steps and where they can fail.

## Lattice

The z-box is halved along both axes together, so every cell has the box's aspect, down
to `zres`, and a few levels further for escalation. Each parameter interval between knots
is halved dyadically.

## Scan

The scan refines breadth-first. A cell coarser than `minsep` is split first. Below that, a cell **fires** when one of its z-faces (one per corner of its
parameter extent) has a nonzero winding of the corner quadrants walked around it. 
A step of two quadrants is an ambiguous ±π turn, which also fires. A face enclosing one
simple zero reads +4, one simple pole −4.

A locus moving fast in `p` can cross a cell between its z-faces, so no face fires. The
**grazing guard** catches it: corners spanning at least 3 quadrants mean the sheets
Re f = 0 and Im f = 0 both cross the cell, and they meet only on loci.

Firing cells split in z down to `zres`. A parameter axis splits only where a locus
leaves the cell's z-column over the cell's parameter step and no neighbouring column picks it up. 
That is, the locus moves more than a cell per step, so adjacency could no longer link it.
Parameter resolution therefore follows the loci's speed rather than
`zres`, and samples crowd where `|dz/dp| → ∞` (two loci meeting) or a zero–pole pair
annihilates.

## Linking and labels

Firing cells at `zres` are grouped into connected branches. Two cells are connected when
their z-columns are equal or adjacent and their parameter ranges touch.

A branch is labelled by the winding of `f` around its **outline**: at a parameter slice,
the branch's z-columns grown by one column. The outline holds every locus of the branch,
so a fused opposite partner cancels. It holds no other branch's loci, since the cells of
another branch are at least two columns away. The winding is computed at up to three
slices where one of the branch's faces fires. If every usable vote has one sign, the
median vote is the label.

Each outline edge is integrated by bisection until every half-step turns the phase by
less than π/2 and changes |f| by less than 2×. The phase test alone would miss a
half-step that turns by an extra full turn, which a multiple locus hugging the edge
causes. Such a locus also swings |f| across the half-step, which the modulus test
catches.

Corner faces are not trusted for labels. They only flag a branch that encloses faces of
both signs, which may hold a zero and a pole fused together.

## Escalation

A branch whose faces show both signs, or whose votes disagree in sign or are all
unusable, is refined one z level below `zres`, regrouped, and labelled again. This
repeats up to `zoom` levels. A zero–pole pair
a few sub-cells apart splits into cleanly labelled branches. A pair that stays connected
down to the last level keeps winding 0. An outline that would leave the box withholds
its vote, so a locus within a column of the box edge also escalates, and stays
unresolved if it is still that close at the last level.

## Aliasing

The scan decides from corner quadrants alone whether a cell may hold something. It is
blind in two ways.

- **Pairs.** Two opposite loci in one cell cancel on its faces. A zero–pole pair's phase
  also fades within a few pair-widths, so nothing marks the cell. Refining every cell to
  `minsep` before any is discarded bounds this. A lone simple locus always fires a face
  or the grazing guard, at any `minsep`.
- **Multiplicity.** A k-fold locus turns the phase k times around a face, and four
  corners read that in quadrant steps. A triple zero centred in a face turns 270° per
  edge and reads −4. A 4-fold one turns a whole number of turns per edge, so all corners
  share a quadrant and the face reads 0. `(z − z0)^k` looks the same at every scale, so
  refining does not help. Labels are not affected, because outlines are integrated
  adaptively; but a k-fold locus that no corner reveals is never discovered.
