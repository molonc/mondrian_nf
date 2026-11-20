# molonc fork of mondrian_nf

This fork is the Mondrian pipeline that iris runs in production. It diverged
from [mondrian-scwgs/mondrian_nf](https://github.com/mondrian-scwgs/mondrian_nf)
at upstream `v0.1.0` and takes upstream changes only by cherry-pick: upstream
removed the `metadata.yaml` inputs and outputs that iris depends on.

## Releases

- Releases are tagged `vX.Y.Z`, starting at `v1.0.0`, and `manifest.version` in
  `nextflow.config` matches the tag. Test runs use a branch or commit SHA, since
  `nextflow run -r` and iris's `pipeline_source` accept both.
- A PR that changes `manifest.version` is a release. When it merges,
  `.github/workflows/release.yml` refuses a tag that upstream already uses,
  tags the merge commit, and publishes the GitHub Release from the matching
  `### vX.Y.Z` section of `CHANGELOG.md`.
- `params.container_version` in `conf/modules.config` names one
  [molonc/mondrian_utils](https://github.com/molonc/mondrian_utils) release; all
  images come from `ghcr.io/molonc/mondrian_utils/<image>` with that tag.
- A tag deploys nothing. iris runs the tag its `pipeline_source` Variable names,
  set on staging first and on production after the analysis test passes.

Tags from before `v1.0.0` remain for the runs that recorded them. Of those, only
the `v0.1.x` tags were releases; `b0.*`, `32cpu`, and `pythonegg-pwd-fix` were
test builds. Upstream also uses `v0.1.2` to `v0.1.6`, for different commits.
