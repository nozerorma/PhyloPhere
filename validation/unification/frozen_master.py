"""Master columns written after the frozen PEPC master (golden/pepc_c4_complete) was made.

The frozen file is the output of an earlier run and is not edited. A test that compares a produced master with it
leaves these columns out; each of them has tests of its own.
"""
NEW_COLUMNS = ("agreement_ambiguous",)


def without_new(df):
    return df.drop(columns=[c for c in NEW_COLUMNS if c in df.columns])


def fields_without_new(fields):
    return [f for f in fields if f not in NEW_COLUMNS]
