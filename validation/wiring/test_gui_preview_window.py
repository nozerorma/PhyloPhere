"""The "Generated scripts preview" window closes once the scripts are saved and the "Scripts saved" message is accepted."""
import importlib
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import qt_stub  # noqa: E402


@pytest.fixture(scope="module")
def mw():
    """gui.widgets.main_window imported against the Qt stub; everything imported meanwhile is forgotten afterwards."""
    before = set(sys.modules)
    saved = {n: sys.modules.get(n) for n in ("PySide6", "PySide6.QtCore", "PySide6.QtGui", "PySide6.QtWidgets")}
    qt_stub.install()
    try:
        yield importlib.import_module("gui.widgets.main_window")
    finally:
        for name in set(sys.modules) - before:
            del sys.modules[name]
        for name, module in saved.items():
            if module is not None:
                sys.modules[name] = module


class _Box:
    """QMessageBox whose answers are set per test and whose calls are recorded in `events`."""
    events = []
    answer_yes = True

    class StandardButton:
        Yes, No = "Yes", "No"

    @classmethod
    def question(cls, *args):
        cls.events.append("question")
        return cls.StandardButton.Yes if cls.answer_yes else cls.StandardButton.No

    @classmethod
    def information(cls, *args):
        cls.events.append("information")

    @classmethod
    def warning(cls, *args):
        cls.events.append("warning")

    @classmethod
    def critical(cls, *args):
        cls.events.append("critical")


class _Preview:
    def __init__(self, events):
        self.events = events

    def close(self):
        self.events.append("close")


def _window(mw, monkeypatch, repo_dir, answer_yes=True, with_preview=True):
    _Box.events, _Box.answer_yes = [], answer_yes
    monkeypatch.setattr(mw, "QMessageBox", _Box)
    win = mw.MainWindow.__new__(mw.MainWindow)
    win.project = SimpleNamespace(general=SimpleNamespace(repo_dir=str(repo_dir), remote_host=""))
    win._preview_window = _Preview(_Box.events) if with_preview else None
    return win


def test_the_preview_closes_after_the_scripts_are_saved_and_the_message_is_shown(mw, monkeypatch, tmp_path):
    win = _window(mw, monkeypatch, tmp_path)
    win._save_generated_scripts([("run.sh", "echo hi\n")])
    assert (tmp_path / "run.sh").read_text() == "echo hi\n"
    assert _Box.events == ["question", "information", "close"]
    assert win._preview_window is None


def test_the_preview_stays_open_when_the_user_declines_to_write(mw, monkeypatch, tmp_path):
    win = _window(mw, monkeypatch, tmp_path, answer_yes=False)
    preview = win._preview_window
    win._save_generated_scripts([("run.sh", "x")])
    assert _Box.events == ["question"] and win._preview_window is preview and not (tmp_path / "run.sh").exists()


def test_the_preview_stays_open_when_writing_fails(mw, monkeypatch, tmp_path):
    win = _window(mw, monkeypatch, tmp_path / "missing_dir")
    preview = win._preview_window
    win._save_generated_scripts([("run.sh", "x")])
    assert _Box.events == ["question", "critical"] and win._preview_window is preview


def test_the_preview_stays_open_without_a_repo_dir(mw, monkeypatch, tmp_path):
    win = _window(mw, monkeypatch, "")
    preview = win._preview_window
    win._save_generated_scripts([("run.sh", "x")])
    assert _Box.events == ["warning"] and win._preview_window is preview


def test_saving_without_a_preview_window_is_fine(mw, monkeypatch, tmp_path):
    win = _window(mw, monkeypatch, tmp_path, with_preview=False)
    win._save_generated_scripts([("run.sh", "x")])
    assert _Box.events == ["question", "information"] and win._preview_window is None


def test_a_window_that_never_had_a_preview_attribute_is_fine(mw, monkeypatch, tmp_path):
    class Strict(mw.MainWindow):  # the stub answers any missing attribute; Qt objects raise
        def __getattr__(self, name):
            raise AttributeError(name)

    win = _window(mw, monkeypatch, tmp_path)
    win.__class__ = Strict
    del win._preview_window
    win._save_generated_scripts([("run.sh", "x")])
    assert _Box.events == ["question", "information"]
