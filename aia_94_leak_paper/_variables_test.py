import aastex

import aia_94_leak_paper


def test_variables():
    result = aia_94_leak_paper.variables()
    assert result
    assert all(isinstance(v, aastex.Variable) for v in result)


def test_variables_unique():
    """A name defined twice would silently take whichever value came last."""
    names = [v.name for v in aia_94_leak_paper.variables()]
    assert len(names) == len(set(names))
