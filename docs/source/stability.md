# API stability

From 1.0.0, PolyPhys follows [Semantic Versioning 2.0.0][semver]. This page
says which names that promise covers and how a name is deprecated and
removed.

## What is public

A name is public when it is listed in the `__all__` of a module on the
[API](api.rst) page. Every other name is private, even if it has no leading
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
