<!-- Release pull request: beta into main, merged with a merge commit.
     Open it with ?template=release.md on the compare page, or gh pr create --template release.md.
     Title: Release vX.Y.Z — <theme> -->

## Theme
<!-- One or two sentences: the opening paragraph of the release's changelog section. -->

## What the release contains
- Changelog section: `CHANGELOG-v0.md`, [X.Y.Z]
- Answer-changing pull requests: <!-- from the section's "Answer changes" list, or "none" -->
- Upgrading notes for users: <!-- yes, with what users must change, or "none" -->

## Before merging
- [ ] Freeze: since release preparation began, only fixes for this release's milestone have merged into beta, and the checks below were rerun after each
- [ ] The "build(release): vX.Y.Z" pull request has merged into beta (version, changelog date, CITATION.cff, README status)
- [ ] CI passes on beta's tip
- [ ] Full test suite on every compiler available (gfortran, ifx; nvfortran if available) in the dev and debug builds
- [ ] Examples rerun on this commit and their outputs archived
- [ ] The release milestone has no open issues

## After merging
- [ ] Annotated tag vX.Y.Z on the merge commit
- [ ] GitHub release, with the changelog section as its notes (Zenodo archives it and issues the DOI)
- [ ] beta fast-forwarded to main
- [ ] Announcement in Discussions → Announcements
