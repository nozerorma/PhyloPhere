#!/usr/bin/env python3
# remote_context.py — Shared "current remote host" state for PathField Browse buttons.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
PathField instances are created in ~15 different places across the 12 tabs
(General, Runtime, and each module tab's essential/fallback fields). Threading a
"current remote host" value through every one of those constructors would mean
touching every call site whenever the active project changes (New/Open). Instead,
PathField reads the current remote host from this tiny shared holder at
Browse-click time — cheap, and correct because this is a single-window,
single-project-at-a-time desktop app (see gui/widgets/main_window.py, which
updates this on every project load/switch).

Deliberately not part of GeneralConfig's own object identity — it's read-only
convenience state derived from GeneralConfig.remote_host, not a second source of
truth for it.
"""

_current_remote_host = ""


def get_remote_host() -> str:
    return _current_remote_host


def set_remote_host(host: str) -> None:
    global _current_remote_host
    _current_remote_host = host.strip()
