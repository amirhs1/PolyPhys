# API stability

From 1.0.0, PolyPhys follows [Semantic Versioning 2.0.0][semver]. This page
says which names and files that promise covers and how a name is
deprecated and removed.

## What is public

A name is public when it is listed in the `__all__` of a module on the
[API](api.rst) page. The methods and attributes of a public class that its
API page documents, those without a leading underscore, are public with the
class. Some packages, such as `polyphys.probe` and `polyphys.simulate`,
re-export public names of their modules; a re-export is the same object and
has the same status, and the module on the API page is its canonical
location. Every other name is private, even if it has no leading
underscore: a private name may change or disappear in any release, without a
deprecation period.

Each public module is stable or provisional:

- **Stable**: a backward-incompatible change to a stable name comes only in
  a major release.
- **Provisional**: the name may still change in a minor release. Each such
  change gets a `**Breaking:**` entry in the changelog.

| Module                    | Status      |
| ------------------------- | ----------- |
| `polyphys.spec`           | Stable      |
| `polyphys.manage.*`       | Stable      |
| `polyphys.simulate.*`     | Stable      |
| `polyphys.analyze.*`      | Stable      |
| `polyphys.probe.*`        | Provisional |
| `polyphys.viz.*`          | Provisional |

(output_files_target)=

## Output files

The files PolyPhys writes are public as far as the documentation names
them: a file name pattern, a column, or a record key that a page of this
documentation or the docstring of a public function names. A file follows
the status, stable or provisional, of the module whose public function
writes it; a column that a function carries over from a file it reads
keeps the status of the function that first wrote it.

- Adding a file to a directory, a column to a table, or a key to a record
  is backward compatible and may come in a minor release. Code that reads
  PolyPhys's files should select columns and keys by name and ignore the
  files, columns, and keys it does not know.
- Renaming or removing a documented file name pattern, column, or key, or
  changing what a documented column, key, or array means, its unit, or the
  shape of a documented array, is backward incompatible. It gets a
  `**Breaking:**` changelog entry, and for a stable module it comes only in
  a major release.
- The order of a table's columns and the text form of its numbers are not
  part of the promise, even where the documentation lists the columns in
  order.

## How a name is deprecated

A stable public name that is to change or go away is first deprecated, in a
minor release. Semantic Versioning requires the minor version to increase
when public functionality is marked deprecated (clause 7), and its FAQ asks
for at least one minor release that contains the deprecation before the
functionality is removed in a major release [^semver-faq]. So a deprecated
name:

1. keeps working as before;
2. warns with {exc}`DeprecationWarning` when it is used, naming its
   replacement and the release that removes it;
3. carries a `.. deprecated::` directive in its docstring, with the version
   that deprecated it; and
4. gets a changelog entry under "Changed" that starts with
   `**Deprecated:**`, since [Common Changelog][common-changelog] has no
   "Deprecated" category and lists a deprecation under "Changed".

A deprecated name is removed only in the next major release, after at least
one minor release has shipped its warning. The removal gets a
`**Breaking:**` changelog entry with the migration path.

## Before 1.0.0

While the major version was `0`, anything could change in any release
(Semantic Versioning, clause 4). The breaking changes since 0.4.0, and how
to migrate from each, are listed in `CHANGELOG.md`, under "Migrating to
1.0.0".

[semver]: https://semver.org/spec/v2.0.0.html
[common-changelog]: https://common-changelog.org/

[^semver-faq]: Semantic Versioning 2.0.0, "How should I handle deprecating
    functionality?", <https://semver.org/spec/v2.0.0.html>.
