import pytest

from checkatlas import atlas
from checkatlas.atlas import list_atlases

given = pytest.mark.parametrize


@given("fn", [atlas(), list_atlases()])
def test_parameterized(fn):
    assert "hello from" in fn()


def test_list_atlases():
    assert list_atlases('.') == "hello from base function"


def test_atlas():
    assert atlas() == "hello from BaseClass"
