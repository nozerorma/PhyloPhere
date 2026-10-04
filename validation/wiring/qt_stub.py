"""A stand-in for PySide6 so that GUI modules import and their plain-Python logic can be tested where Qt is not installed.

Every name of PySide6.QtCore/QtGui/QtWidgets is a class that accepts any call and any attribute; nothing is drawn. Tests
that need an observable widget or message box replace that one name with their own fake.
"""
import sys
import types


class _Meta(type):
    def __getattr__(cls, name):
        if name.startswith("__"):
            raise AttributeError(name)
        value = _Meta(name, (_Any,), {})
        type.__setattr__(cls, name, value)
        return value

    def __or__(cls, other):  # Qt flag combinations
        return cls


class _Any(metaclass=_Meta):
    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):  # also lets an instance such as Slot(...) act as a decorator
        return args[0] if len(args) == 1 and callable(args[0]) and not kwargs else _Any()

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        return _Any()

    def __iter__(self):
        return iter(())


def install():
    """Put the stub modules in sys.modules; returns the names added or replaced so the caller can undo it."""
    names = ("PySide6", "PySide6.QtCore", "PySide6.QtGui", "PySide6.QtWidgets")
    for name in names:
        module, cache = types.ModuleType(name), {}

        def getter(attr, cache=cache):
            if attr.startswith("__"):
                raise AttributeError(attr)
            return cache.setdefault(attr, _Meta(attr, (_Any,), {}))

        module.__getattr__ = getter
        sys.modules[name] = module
    return names
