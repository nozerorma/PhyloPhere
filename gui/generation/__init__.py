# __init__.py — Script generation for the runner GUI: render context, templates, validation, report regeneration.
# PhyloPhere | gui/generation/
#
# Package of the runner GUI that turns a ProjectConfig into the two shell scripts the GUI
# saves. It has no PySide6 dependency.
#
#   context.py          ProjectConfig → Jinja2 render context
#   render.py           context + templates/*.j2 → batch and single-phenotype scripts
#   validate.py         required-field and path-existence checks before rendering
#   report_registry.py  standalone re-rendering of the HTML reports of a finished run
#   templates/          sbatch_array.sh.j2 (batch runner), run_single.sh.j2 (one phenotype)
#
# Imported by: gui/widgets/main_window.py, gui/widgets/common/regenerate_dialog.py
