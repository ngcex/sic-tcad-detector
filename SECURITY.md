# Security policy

ETNA is a research simulation library. It runs locally, and the bundled Streamlit app
binds to `localhost` in the documented commands. It is not a network service. Security
issues are still possible, for example:

- unsafe parsing of user-supplied measurement or layout files;
- a vulnerable pinned dependency;
- the Streamlit app being exposed beyond `localhost`.

## Supported versions

Only the latest release, and `main`, receive fixes.

## Reporting a vulnerability

Please **do not open a public issue** for a security problem. Report it privately instead
through GitHub's private vulnerability reporting: open the repository's **Security** tab
and choose **Report a vulnerability**.

> Maintainer note: private vulnerability reporting must be enabled in the repository
> settings (Settings -> Code security -> Private vulnerability reporting). Until it is,
> open a minimal public issue that asks for a private contact channel, and leave out any
> detail of the problem.

Please include:

- the affected version or commit;
- the steps to reproduce;
- the impact as you understand it.

The project is maintained on a best-effort basis by one person. You should get an
acknowledgement within about two weeks, and a fix or mitigation as soon as is practical.
Credit will be given in the changelog unless you ask otherwise.

## Dependencies

- A nightly CI job audits the locked dependency set with `pip-audit` (see
  `.github/workflows/ci.yml`).
- Dependabot proposes updates for `uv.lock` and the GitHub Actions
  (`.github/dependabot.yml`).
