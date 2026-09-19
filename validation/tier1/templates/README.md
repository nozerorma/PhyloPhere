# Templates

Hard copies of the two Tier 1 GUI project files, kept here for version
tracking:

- [`tier1_pepc.json`](tier1_pepc.json) — PEPC / C4, categorical
- [`tier1_hb.json`](tier1_hb.json) — Hb / high altitude, continuous + categorical

The versions Phylophere's GUI actually loads live in `gui/templates/` (its
"Load template" dialog opens there) — update both copies together. See
`../README.md` for what each one runs.

Plain JSON (`ProjectConfig`, schema_version 1) — hand-editable, diffable, and
loadable either via the GUI's File → Load template or directly by
`gui.project_io.load_project`.
