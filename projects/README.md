# Study specs

This folder holds the `ProjectSpec` of each of the maintainer's studies: the
facts their analysis needs beyond the parser class, such as the polymer group,
the rounding step of the bulk crowder volume fraction, the attributes added as
columns, the species-to-size map of the density profiles, and the time series.
`ProjectSpec` itself is in `polyphys.manage.projects`.

The folder is part of the repository, not of the installed package
(`pyproject.toml` excludes `projects*`), so PolyPhys holds no table of studies.
A study elsewhere builds its own `ProjectSpec` in the same way.

| Folder            | Study          | Status      |
| ----------------- | -------------- | ----------- |
| `trans_foci_cyl/` | `TransFociCyl` | published   |
| `trans_foci_cub/` | `TransFociCub` | published   |
| `hns_cub/`        | `HnsCub`       | published   |
| `hns_cyl/`        | `HnsCyl`       | provisional |

A provisional spec belongs to a study whose molecular dynamics data is not yet
confirmed. The specs of unpublished studies are kept out of this public
repository. Each folder is a package whose `spec.py` builds the study's
`ProjectSpec` and whose `__init__.py` exports it as `SPEC`; the probe,
analysis, and figure code of a study will join it there. From an analysis
script in a checkout of this repository:

```python
import sys
sys.path.insert(0, "<path to PolyPhys>/projects")
from hns_cub import SPEC

space = space_tseries(database, "gyrTMon", SPEC.parser, SPEC.hierarchy,
                      list(SPEC.physical_attrs), SPEC.group,
                      divisor=SPEC.divisor)
```

`test_studies.py` checks every spec against its parser and against the values
of the research code's registry; CI runs it with flake8 and strict mypy.
