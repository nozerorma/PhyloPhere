#!/usr/bin/env python3
# remote_context.py — Shared "current remote host/root dir" state for PathField Browse buttons.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Module-level holder of the current remote host and remote root directory.

PathField instances are created in many places (General, Runtime and
Precomputed tabs, every module tab's path fields, the regenerate dialog). Passing the current
remote host to each constructor would mean touching every call site whenever the
active project changes (New/Open). Instead, PathField reads the host and root
directory from this holder when Browse is clicked. This is safe because the
application has one window and one project at a time.

The holder is a synchronized copy of GeneralConfig.remote_host and
GeneralConfig.remote_root_dir, not a second source of truth: the General tab
(gui/widgets/tabs/general_tab.py) writes it whenever those values change or the
tab is rebuilt.

Imported by: gui/widgets/common/path_field.py, gui/widgets/common/regenerate_dialog.py,
gui/widgets/tabs/general_tab.py, gui/widgets/tabs/runtime_tab.py
"""

_current_remote_host = ""
_current_remote_root_dir = ""


def get_remote_host() -> str:
    return _current_remote_host


def set_remote_host(host: str) -> None:
    global _current_remote_host
    _current_remote_host = host.strip()


def get_remote_root_dir() -> str:
    return _current_remote_root_dir


def set_remote_root_dir(path: str) -> None:
    global _current_remote_root_dir
    _current_remote_root_dir = path.strip()
