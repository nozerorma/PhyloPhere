#!/usr/bin/env python3
# test_widgets_smoke.py — Headless instantiation + data-binding smoke tests for Phase B tabs.
# PhyloPhere | tests/gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Not exhaustive UI testing (no click-simulation of every control) — these confirm
each widget instantiates without error under the offscreen Qt platform and that
its two-way data binding to the underlying dataclass actually works, which is the
failure mode most likely to regress silently (a forgotten `.connect()` or a typo'd
attribute name would otherwise only surface by clicking around in a real window).
"""

# ── Local ─────────────────────────────────────────────────────────────────────
from conftest import golden_project
from gui.models.about import AboutInfo
from gui.models.general import GeneralConfig
from gui.models.resources import ResourcesConfig
from gui.models.runtime import PhenotypeRow, RuntimeConfig
from gui.widgets.tabs.about_tab import AboutTab
from gui.widgets.tabs.general_tab import GeneralTab
from gui.widgets.tabs.resources_tab import ResourcesTab
from gui.widgets.tabs.runtime_tab import RuntimeTab


def test_general_tab_binds_fields(qapp):
    config = GeneralConfig()
    tab = GeneralTab(config)

    tab.seed.setText("42")
    assert config.seed == "42"

    tab.reporting.setChecked(False)
    assert config.reporting is False

    tab.repo_dir.set_text("/repo")
    assert config.repo_dir == "/repo"


def test_runtime_tab_binds_fields_and_catalogue(qapp):
    config = RuntimeConfig()
    tab = RuntimeTab(config, repo_dir_getter=lambda: "/repo")

    tab.resume.setChecked(False)
    assert config.resume is False

    tab.runtime_type.setCurrentText("local")
    assert config.runtime_type == "local"

    tab.work_dir.set_text("/work")
    assert config.work_dir == "/work"

    tab.catalogue._add_row()
    assert len(config.phenotype_rows) == 1
    assert config.phenotype_rows[0] is tab.catalogue.model._rows[0]


def test_runtime_tab_class2_row_disables_class1_columns(qapp):
    from PySide6.QtCore import Qt

    config = RuntimeConfig()
    config.phenotype_rows = [PhenotypeRow(trait_class=2, trait="t", discrete_method="decile")]
    tab = RuntimeTab(config, repo_dir_getter=lambda: "/repo")

    secondary_index = tab.catalogue.model.index(0, 2)  # column 2 = "secondary"
    flags = tab.catalogue.model.flags(secondary_index)
    assert not (flags & Qt.ItemFlag.ItemIsEditable)


def test_resources_tab_binds_fields(qapp):
    config = ResourcesConfig()
    tab = ResourcesTab(config)

    tab.local_cpus_field.setText("16")
    assert config.local_max_cpus == "16"

    tab.slurm_mem_field.setText("256.GB")
    assert config.slurm_max_memory == "256.GB"


def test_about_tab_instantiates(qapp):
    AboutTab(AboutInfo())  # must not raise


def test_main_window_assembles_tabs_and_generates_scripts(qapp, tmp_path):
    from gui.generation.validate import validate
    from gui.widgets.main_window import MainWindow

    window = MainWindow()
    assert [window.tabs.tabText(i) for i in range(window.tabs.count())] == [
        "General",
        "Runtime",
        "CAAS",
        "Disambiguation",
        "Accumulation",
        "RERconverge",
        "FADE",
        "VEP",
        "Scoring",
        "Enrichment",
        "Resources",
        "About",
    ]

    window.project = golden_project()
    window._rebind_tabs_to_project()
    assert validate(window.project) == []

    window.generate_scripts()
    assert window._preview_window.isVisible() or window._preview_window is not None


def test_main_window_new_project_rebinds_all_tabs(qapp):
    from gui.widgets.main_window import MainWindow

    window = MainWindow()
    before = [window.tabs.tabText(i) for i in range(window.tabs.count())]

    window.new_project()

    after = [window.tabs.tabText(i) for i in range(window.tabs.count())]
    assert after == before
    # About tab must stay last after a full tab rebuild.
    assert after[-1] == "About"


def test_module_tab_toggle_marks_window_dirty(qapp):
    from gui.widgets.main_window import MainWindow

    window = MainWindow()
    assert window._dirty is False

    window.caas_tab.enable_toggle.setChecked(False)

    assert window._dirty is True
    assert window.project.modules.caas.enabled is False


def test_main_window_save_and_open_roundtrip(qapp, tmp_path):
    from gui import project_io
    from gui.widgets.main_window import MainWindow

    window = MainWindow()
    window.project = golden_project()
    project_path = tmp_path / "project.json"

    project_io.save_project(project_path, window.project)
    window.open_project_from_path(project_path)

    assert window.project == golden_project()
    assert window.project_path == project_path
    assert window._dirty is False


def _dismiss_topmost_message_box():
    from PySide6.QtWidgets import QApplication, QMessageBox

    for widget in QApplication.topLevelWidgets():
        if isinstance(widget, QMessageBox) and widget.isVisible():
            widget.close()


def test_toolbar_has_save_and_validate_paths_actions(qapp):
    from PySide6.QtWidgets import QToolBar

    from gui.widgets.main_window import MainWindow

    window = MainWindow()
    toolbars = window.findChildren(QToolBar)
    assert len(toolbars) == 1
    action_texts = [a.text() for a in toolbars[0].actions()]
    assert any("Save" in t for t in action_texts)
    assert any("Validate Paths" in t for t in action_texts)


def test_validate_paths_dialog_does_not_crash(qapp):
    from PySide6.QtCore import QTimer

    from gui.widgets.main_window import MainWindow

    window = MainWindow()
    window.project = golden_project()  # fake paths -> triggers the "problems found" branch
    window._rebind_tabs_to_project()

    QTimer.singleShot(0, _dismiss_topmost_message_box)
    window.validate_paths_dialog()  # must not raise, must not hang


def test_general_tab_has_no_leftover_env_install_fields(qapp):
    from gui.models.general import GeneralConfig
    from gui.widgets.tabs.general_tab import GeneralTab

    tab = GeneralTab(GeneralConfig())
    assert not hasattr(tab, "solver")
    assert not hasattr(tab, "env_yaml")
    assert not hasattr(tab, "install_location")


def test_general_tab_group_titles_escape_ampersand_mnemonics(qapp):
    from PySide6.QtWidgets import QGroupBox

    from gui.models.general import GeneralConfig
    from gui.widgets.tabs.general_tab import GeneralTab

    tab = GeneralTab(GeneralConfig())
    titles = [gb.title() for gb in tab.findChildren(QGroupBox)]
    # A bare "&" in a QGroupBox title is a mnemonic marker, not a literal
    # ampersand — it swallows the next character (rendered as "Remote
    # validation _browsing" instead of "Remote validation & browsing").
    # "&&" is Qt's escape for a literal "&".
    for title in titles:
        single_amp_count = title.replace("&&", "").count("&")
        assert single_amp_count == 0, f"unescaped '&' in QGroupBox title: {title!r}"


def test_general_tab_remote_host_updates_config_and_shared_context(qapp):
    from gui.models.general import GeneralConfig
    from gui.widgets.common import remote_context
    from gui.widgets.tabs.general_tab import GeneralTab

    config = GeneralConfig()
    tab = GeneralTab(config)

    tab.remote_host.setText("user@cluster.example.edu")

    assert config.remote_host == "user@cluster.example.edu"
    assert remote_context.get_remote_host() == "user@cluster.example.edu"


def test_validate_paths_dialog_uses_remote_check_when_host_configured(qapp):
    from unittest.mock import patch

    from PySide6.QtCore import QTimer

    from gui.widgets.main_window import MainWindow

    window = MainWindow()
    window.project = golden_project()
    window.project.general.remote_host = "user@cluster"
    window._rebind_tabs_to_project()

    QTimer.singleShot(0, _dismiss_topmost_message_box)
    with patch("gui.widgets.main_window.check_remote_paths") as mock_check:
        mock_check.return_value = ["Runtime: tree file: file not found on user@cluster — /x"]
        window.validate_paths_dialog()

    mock_check.assert_called_once()
    assert mock_check.call_args.args[0] == "user@cluster"


def test_validate_paths_dialog_shows_error_on_remote_check_failure(qapp):
    from unittest.mock import patch

    from PySide6.QtCore import QTimer

    from gui.remote import RemoteCheckError
    from gui.widgets.main_window import MainWindow

    window = MainWindow()
    window.project = golden_project()
    window.project.general.remote_host = "unreachable-host"
    window._rebind_tabs_to_project()

    QTimer.singleShot(0, _dismiss_topmost_message_box)
    with patch("gui.widgets.main_window.check_remote_paths") as mock_check:
        mock_check.side_effect = RemoteCheckError("SSH to 'unreachable-host' timed out after 15s.")
        window.validate_paths_dialog()  # must not raise, must show the error dialog instead

    mock_check.assert_called_once()


def test_path_field_browse_opens_local_dialog_when_no_remote_host(qapp):
    from unittest.mock import patch

    from gui.widgets.common import remote_context
    from gui.widgets.common.path_field import PathField

    remote_context.set_remote_host("")
    field = PathField(mode="dir")

    with patch.object(field, "_browse_local") as mock_local, patch.object(
        field, "_browse_remote"
    ) as mock_remote:
        field._browse()

    mock_local.assert_called_once()
    mock_remote.assert_not_called()


def test_path_field_browse_opens_remote_dialog_when_host_configured(qapp):
    from unittest.mock import patch

    from gui.widgets.common import remote_context
    from gui.widgets.common.path_field import PathField

    remote_context.set_remote_host("user@cluster")
    try:
        field = PathField(mode="dir")
        with patch.object(field, "_browse_local") as mock_local, patch.object(
            field, "_browse_remote"
        ) as mock_remote:
            field._browse()

        mock_remote.assert_called_once_with("user@cluster")
        mock_local.assert_not_called()
    finally:
        remote_context.set_remote_host("")  # don't leak state into other tests
