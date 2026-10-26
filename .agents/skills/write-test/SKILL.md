---
name: write-test
description: Write or change a test in polyphys/tests/. Use when writing a module's test file, adding a test for a behaviour change or bug fix, or materially changing one; not for doctest Examples in source docstrings.
---

# Write a test

## Before you write

1. For a module's test file, `polyphys/tests/<layer>/test_<module>.py`, start
   once the module is complete, and land the file in the same pull request as
   the module. Start from the test-case list in the approved acceptance
   criteria: for each case the behaviour, the input, and the source of the
   expected value. Write the file in one pass against the finished API. Use
   one fixture set for the module. Do not add cases the list lacks without
   saying so.
2. When the module already has a test file, carry forward every test that
   still states a valid contract, rewritten to the rules below. Delete a test
   only with its reason given in the pull request; the module's coverage must
   not fall.
3. For a bug fix, write the regression test first and watch it fail against
   the unfixed code, then fix the code. For any other behaviour change, the
   claim comes from the change.

## Each test

4. Name it after the behaviour: `test_<unit>_<behaviour>`.
5. Give it a one-line docstring stating the behaviour, invariant, or regression
   it protects, in the present tense: the claim that is false if the test
   fails. No issue numbers; no "Test that…" or restating the name.
6. Test one behaviour, in arrange, act, assert order: set up, call the code
   once, then assert. Several asserts are fine when they check that one
   behaviour.
7. Keep loops and branches out of the test body. Put cases in
   `@pytest.mark.parametrize` with `pytest.param(..., id="...")`, so each case
   has a readable id and fails on its own.
8. Show in the test, or in its `parametrize` cases, every value the result
   depends on; fixtures supply only incidental setup. Put shared data in
   pytest fixtures, never in module-level constants, and use a factory fixture
   with overridable defaults when tests need variants of one object.
9. Give fixtures function scope. Use module or session scope only for setup
   that is expensive and read-only: make a shared array read-only
   (`arr.flags.writeable = False`) and build shared files with
   `tmp_path_factory`. A fixture used by more than one test file goes in a
   `conftest.py`; never import fixtures.
10. Check errors with `pytest.raises(ErrorType, match="...")` and warnings with
    `pytest.warns(WarningType, match="...")`, each around the single call that
    should raise or warn.
11. Take every expected value from an analytic result, a cited reference,
    measured data, or a reproducible computation, and name it in the commit's
    `Ground-truth-source:` trailer. With no such source, test a property
    (symmetry, conservation, invariance) and say so in the docstring.
12. Compare floating-point results with a stated tolerance:
    `numpy.testing.assert_allclose(actual, expected, rtol=..., atol=...)` for
    arrays, which reports the elements that differ, and
    `pytest.approx(expected, rel=..., abs=...)` for scalars. Use
    `numpy.testing.assert_array_equal` for integer, boolean, and string
    arrays. Give `atol` whenever an expected value can be zero, because a
    relative tolerance alone never accepts a non-zero result there. Say where
    each tolerance comes from, such as float64 round-off or a statistical
    error bar, in a comment beside it.

## After you write

13. Show that each new test can fail: change its expected value or break the
    code under test once, run it, watch it fail, and restore it. One case of a
    parametrized test is enough. Say in the pull request that you did.
14. Never weaken an expectation or paste observed output to make a test pass;
    report the failure. A test written after the code checks the criteria,
    not what the code happens to do.
15. Leave other modules' legacy tests as they are.
16. Run the changed test first, then the full suite, as AGENTS.md,
    "Commands", gives them.

## Example

```python
def test_ensemble_average_is_replicate_order_invariant():
    """Averaging replicates gives the same ensemble regardless of input order."""

def test_parse_segment_keeps_last_attribute():
    """A trailing underscore no longer drops a segment name's last attribute."""
```

Not: `"""Test parse() on a segment filename."""` (it restates the code and
states no claim).
