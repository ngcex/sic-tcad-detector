# Releasing ETNA

ETNA is released from GitHub only (no PyPI). A release is a git tag `vX.Y.Z`. Pushing the
tag runs `.github/workflows/release.yml`, which:

- builds the sdist and wheel;
- installs the wheel in a clean venv outside the checkout and smoke-tests it
  (`scripts/wheel_smoke.py --app`);
- publishes a GitHub Release with the matching `CHANGELOG.md` section as notes, the
  artefacts and `SHA256SUMS.txt`.

Only the maintainer creates tags.

## Checklist

1. **Decide the version** ([Semantic Versioning](https://semver.org/)):
   - patch: fixes, with no change to default numbers beyond the stated tolerances;
   - minor: new opt-in features;
   - major: a changed default or a removed or renamed API.

   A changed default that moves pinned numbers (for example a new default damage model)
   counts as major, or as minor with a prominent "Changed" entry. Pick one and say why
   in the changelog.

2. **Bump** `etna/_version.py`. It is the single source of truth; `pyproject.toml` reads
   it.
3. **Update `CITATION.cff`**: `version:` must equal `etna/_version.py`
   (`tests/test_project_metadata.py` checks this).
4. **Update `CHANGELOG.md`**:
   - Move the `## [Unreleased]` entries into a new `## [X.Y.Z] - YYYY-MM-DD` section and
     leave an empty `## [Unreleased]` above it.
   - Update the link references at the bottom.
   - Check the section locally with `uv run python scripts/release_notes.py X.Y.Z`.
5. **Run the gates locally** on the release commit:
   - `uv sync --extra dev`
   - `uv run pytest -n auto -q`: all green.
   - `uv run --no-sync ruff check .`
   - `uvx cffconvert --validate -i CITATION.cff`
   - `uvx --from mkdocs==1.6.1 --with mkdocs-material==9.7.7 mkdocs build --strict`
   - Optional: `uv build` and the wheel smoke test from the `scripts/wheel_smoke.py`
     docstring.
6. **Merge to `main` and push.** Wait until CI (`ci.yml`) and the matrix workflow
   (`matrix.yml`) are green on that commit.
7. **Tag and push only that tag.** Never use `git push --tags` or `--all` from a working
   clone.

   ```bash
   git tag -a vX.Y.Z -m "ETNA X.Y.Z"
   git push origin vX.Y.Z
   ```

8. **Check the workflow run and the release page.** It should have the sdist, the wheel,
   `SHA256SUMS.txt` and the changelog section as notes. If the workflow fails before the
   release is created:
   - fix the cause on `main`;
   - delete the tag locally and on the remote (`git push origin :refs/tags/vX.Y.Z`);
   - tag again.

   Never move a tag after a release has been published or archived.

## DOI (Zenodo)

No DOI exists yet. To get one:

1. The maintainer links the repository in Zenodo's GitHub integration.
2. Zenodo archives every GitHub Release published after that, reads `CITATION.cff`, and
   mints a version DOI plus a concept DOI that covers all versions.
3. Add the DOI to `CITATION.cff` (`identifiers:` / `doi:`) and a DOI badge to
   `README.md`, both in a follow-up commit.

A DOI is permanent, so do not mint test releases.

## What the release does not do

- No PyPI or TestPyPI upload. The `etna` name is taken on PyPI, and the maintainer
  decided on GitHub-only distribution (issue #3, H3).
- No docs deployment. `.github/workflows/docs.yml` is manual until GitHub Pages is
  configured.
