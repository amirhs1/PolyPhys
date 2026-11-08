# Study specs

This folder holds the parser and the specs of each of the maintainer's
published studies, for analysis scripts run from a checkout of this
repository. It is part of the repository, not of the installed package
(`pyproject.toml` excludes `projects*`), so PolyPhys holds no table of studies.
The [Studies](../docs/source/studies.md) page shows how to declare a study's
parser and specs; a study elsewhere builds its own in the same way.

| Folder            | Study          |
| ----------------- | -------------- |
| `trans_foci_cyl/` | `TransFociCyl` |
| `trans_foci_cub/` | `TransFociCub` |
| `hns_cub/`        | `HnsCub`       |
| `hns_cyl/`        | `HnsCyl`       |

The parsers and specs of unpublished studies are kept out of this public
repository.

Each folder is a package whose `parser.py` declares the study's `ParserBase`
subclass, whose `spec.py` builds its `ProjectSpec` (`SPEC`), `ProbeSpec`
(`PROBE`), and `AnalysisTable` records (`ANALYSIS`), and whose `__init__.py`
exports the parser class and the three specs; the study's drivers and figure
code will join them there. The `ProjectSpec` holds the facts the study's
analysis needs beyond the parser class, such as the polymer group, the
rounding step of the bulk crowder volume fraction, the attributes added as
columns, the species-to-size map of the density profiles, and the time
series. `ProjectSpec` itself is in `polyphys.spec`.

The recipe below stacks a study's time series into a `space` table. It is not
self-contained: it needs the `ensAvg` files of an existing sweep of the study.
`database` is the `ensAvg` directory of one space, such as
`analysis/<space>-nucleoid-ensAvg`, which holds the files
`<ensemble_long>-<group>-<property>-ensAvg.csv`. The
[Quickstart](../docs/source/quickstart.md) builds a `space` table from a
synthetic sweep, and its examples run as tests.

```python
import sys
sys.path.insert(0, "<path to PolyPhys>/projects")
from hns_cub import SPEC
from polyphys.analyze.pipeline import space_tseries_from_spec

database = "<path to the sweep>/analysis/<space>-nucleoid-ensAvg"
space = space_tseries_from_spec(database, SPEC, "gyrTMon")
```

`test_studies.py` checks every spec against its parser and against the values
of the research code's registry; CI runs it with flake8 and strict mypy.
