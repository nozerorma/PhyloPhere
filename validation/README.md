# PhyloPhere validation

Checks PhyloPhere's output against fixtures with a known-truth site list.

The approach: run `ortholog_characterizator` (FUBAR/MEME positive selection)
and PhyloPhere itself on the same fixture, then compare what each recovers
against the fixture's curated truth set -- how much of the expected signal
PhyloPhere sees, and what else it finds beyond it.

## Layout

```
validation/
  tier1/            fixtures + runners (see validation/tier1/README.md)
  truthsets/tier1/  curated known-positive site lists, one per fixture
  papers/           source literature for each fixture (pepc/)
```

## Fixtures

Two fixtures under `tier1/input/`:

- `pepc/` -- PEPC gene, C3/C4 photosynthesis in sedges (Cyperaceae)

Each fixture directory has its own `scripts/build.py` (alignment + tree +
trait table) and `scripts/build_cds.py` (codon-level CDS alignment for
positive-selection analysis).

## Workflow

1. **Build the fixture** -- run the fixture's `scripts/build.py` then
   `scripts/build_cds.py`.
2. **Run ortholog_characterizator** -- `validation/tier1/run_ortholog_characterizator.sh {pepc}`
   runs the `ortholog_characterizator` pipeline (quality -> translation ->
   phylogeny -> positive selection) on the fixture's CDS alignment, giving an
   independent FUBAR/MEME positive-selection call set.
3. **Run PhyloPhere** -- load `gui/templates/tier1_pepc.json` in the GUI (File -> Load template)
   and run the pipeline against the same fixture.
4. **Compare** -- check both outputs against the fixture's truth set in
   `validation/truthsets/tier1/` (`pepc_c4.sites.tsv`):
   which truth sites each method recovers, and what additional sites each
   flags beyond the truth set.

See `validation/tier1/README.md` for the fixture/template/output layout in
detail.
