# Template

`tier1_pepc_c4.json` is a copy of `gui/templates/tier1_pepc_c4.json`, the GUI project that runs the PEPC case (trait `c4`, categorical). The GUI loads its templates from `gui/templates/`, and `render_tier1_scripts.py` renders the run scripts from that file; this copy keeps the case with its fixture and must stay byte-identical to it.

The file is plain JSON (`ProjectConfig`): it can be edited by hand, compared with `diff`, and loaded from the GUI (File, Load template) or with `gui.project_io.load_project`. `../../README.md` states what the template runs.
