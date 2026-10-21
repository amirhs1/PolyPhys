---
name: write-test
description: Write or change a test in polyphys/tests/. Use when adding a test for a behaviour change or bug fix, or materially changing one; not for doctest Examples in source docstrings.
---

# Write a test

1. Name it after the behaviour: `test_<unit>_<behaviour>`.
2. Give it a one-line docstring stating the behaviour, invariant, or regression
   it protects, in the present tense: the claim that is false if the test
   fails. No issue numbers; no "Test that…" or restating the name.
3. Take every expected value from an analytic result, a cited reference,
   measured data, or a reproducible computation, and name it in the commit's
   `Ground-truth-source:` trailer. With no such source, test a property
   (symmetry, conservation, invariance) and say so in the docstring.
4. Never weaken an expectation or paste observed output to make a test pass;
   report the failure.
5. Leave unrelated legacy tests as they are.
6. Run the changed test first, then the full suite, as AGENTS.md, "Commands",
   gives them.

## Example

```python
def test_ensemble_average_is_replicate_order_invariant():
    """Averaging replicates gives the same ensemble regardless of input order."""

def test_parse_segment_keeps_last_attribute():
    """A trailing underscore no longer drops a segment name's last attribute."""
```

Not: `"""Test parse() on a segment filename."""` (it restates the code and
states no claim).
