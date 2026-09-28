# Before the transfer to a NOAA organization

State and plan, not policy. Written 2026-09-23. This is the checklist for everything that has to
happen **before** the repository changes owners. The transfer itself, and what follows it, are
in `PLAN-adoption-and-NOAA-transfer.md` sections 2–4. The release steps are in
`inst/RELEASE-CHECKLIST.md`, and what is queued now is in `SESSION_HANDOFF.md`. This file points
to both rather than repeating them.

**Owner tags:**

- **[Grant]** needs a decision, credentials, or another person.
- **[agent]** can be done in a working session and then reviewed.
- **[either]** could be done by either.

Work the stages in order. Stage A can run in parallel with B and C, because it is mostly
waiting on other people.

## Where things stand (2026-09-27)

- **The transfer is done**: the repository is `afsc-assessments/Rceattle` as of 2026-09-26, and
  Stage D's link sweep shipped at 5.45.0. Stages A, E and F are settings work that outlives it.
- `dev` is at 5.45.0. `main` is at 5.33.0.
- **The newest tag is 5.28.0.** `main` sits 42 commits past it and is untagged.
- `golden` (the `deep-checks` workflow) fails on `main`. Windows `R-CMD-check` fails
  intermittently. Both are recorded in `SESSION_HANDOFF.md` and `TRAPS.md`.
- There is no `LICENSE` file, no `.github/CODEOWNERS`, and no `.mailmap`.
- There are 24 remote branches, including three `depricated-*` branches and several merged
  feature branches.
- 36 lines in 11 tracked files hard-code `grantdadams` (not counting `NEWS.md`).
- `review/pr158-tier0` merged as PR #159 on 2026-09-24; its follow-on corrections are on
  `fix/release-doc-corrections`. Check `git status` before starting.

---

## Stage A. Decisions and permissions [Grant]

Nothing in Stage D can be finalized until A1 is settled.

- [ ] **A1. Choose the destination organization.**
  - Recommended: `afsc-assessments`. Its reasoning and the alternatives are in the PLAN,
    section 0.
  - Email the organization owners. Ask for permission to create repositories there, and
    confirm the repository name `Rceattle` is free.
- [ ] **A2. Find out whether that organization is part of NMFS GitHub Enterprise Cloud (GHEC).**
  - If it is, get the NMFS GHEC user agreement signed by your supervisor, and request
    organization membership through the NMFS form.
  - Confirm that SSO will not lock out outside collaborators (UW, international users).
- [ ] **A3. Check the organization's GitHub Actions policy.** The workflows use these
  third-party actions:
  - `r-lib/actions/*`
  - `JamesIves/github-pages-deploy-action@v4`
  - `codecov/codecov-action`

  If any are blocked, ask for them to be allowed, or plan to swap to first-party actions
  (`actions/deploy-pages` for the site).
- [ ] **A4. License.**
  - Ask the NMFS GitHub Governance Team whether GPL (>= 2) meets NAO 201-118 for this package.
  - Their guidance defaults to Apache 2.0, but Rceattle links TMB (GPL), and Apache 2.0 cannot
    be combined with GPL-2.
  - Keep GPL unless they object. Relicensing would need every copyright holder's consent.
- [ ] **A5. Name a co-maintainer** and get their agreement. Cole is the natural fit.
  - They will be a second owner on the repository and the second reviewer in CODEOWNERS.
  - Tell Melissa, because it answers the continuity ask in the one-pager.
- [ ] **A6. Decide whether `Rceattle-models` moves too** (recommended: yes, same organization).
  If it does, it needs its own short pass through Stages E and F.
- [ ] **A7. Pick the transfer date.**
  - It must fall after Stage B is done.
  - Avoid the two weeks before a Plan Team document deadline, and any period when the GOA
    pollock or arrowtooth teams are running final models.

**Done when:** each of A1–A7 has a written answer in this file.

---

## Stage B. Make `main` releasable and tagged

Tags, releases, and redirects all carry through the transfer. **Tag first** so users have a
clean pin that works on both sides of the move.

- [ ] **B1. Land the in-flight work [either].**
  - Finish, merge, or park `review/pr158-tier0` and any open PRs into `dev`.
  - Do not start the release from a dirty tree.
- [ ] **B2. Ship the 5.34.0–5.45.0 release [either, Grant publishes]**, following
  `inst/RELEASE-CHECKLIST.md` exactly:
  1. Open and merge the `dev` -> `main` PR.
  2. Tag **the merge commit** with the bare version, which is the DESCRIPTION version --
     `5.45.0`, **no `v`** (checklist section 3). It is not 5.41.0: #160 and the review of
     #158 both landed after the release PR was written.
  3. Publish the GitHub Release. The site rebuilds only on `release: published`.
  4. **Confirm pkgdown actually ran.** The event has failed to fire before; if it did not, run
     `gh workflow run pkgdown.yaml --ref main`.
  5. Dispatch `deep-checks`.

  Expect `golden` to be red for the known reason, and do not chase it in this release.
- [ ] **B3. Make `golden` robust [agent, Grant reviews].**
  - The problem: `goa_ss` has a second local minimum, 52.9 negative log-likelihood units
    higher, and a one-ULP gradient change can send `nlminb` there (`TRAPS.md`).
  - The fix: warm-start the four reference fits from the pinned `par`, or fit from two starts
    and keep the lower minimum.
  - This is a harness change only. It must not move any pinned objective.
  - **Done when:** `deep-checks` `golden` is green on `main`.

  Ship it as a patch release (5.45.1) if it lands after B2, so the version being transferred
  has a green release gate.
- [x] **B4. Apply the tag convention everywhere [agent].** Already satisfied, verified
  2026-09-23: `R/0-rceattle_class.R:12` and `man/print.Rceattle.Rd:28` read `@X.Y.Z`, and
  `README.md:43` reads `@5.45.0`. Both are bare, per the convention.
- [x] **B5. Check the other install lines [agent].** Verified with `git grep -n
  "install_github"`: every install command in the README, the vignettes, `examples/` and the Rd
  files is either untagged or uses the bare convention. `README.md:43` names `5.45.0`, which
  becomes valid the moment B2's tag is pushed — so **push the tag, or that line is wrong.**
  It read `5.41.0` until the review of #158; that version is never tagged, so the line named a
  reference `install_github()` could not resolve. Re-check it whenever the version moves
  again: this line is only correct against the version actually tagged.

---

## Stage C. Files the repository should have before it moves [agent]

Make these one PR into `dev`, then a patch release. None of them depend on the destination
organization.

- [ ] **C1. License file.** Run `usethis::use_gpl_license(version = 2, include_future = TRUE)`.
  - This writes `LICENSE.md` and adds it to `.Rbuildignore`. CRAN does not want a copy of the
    GPL inside the tarball, but NOAA guidance wants a license file in the repository; this
    satisfies both.
  - Confirm `DESCRIPTION` still reads `License: GPL (>= 2)`.
- [ ] **C2. `.github/CODEOWNERS`.** Give `*` two owners, Grant and the co-maintainer from A5,
  using their GitHub usernames.
  - Do **not** require code-owner review yet. That can only be switched on once there really
    are two reviewers, or it blocks every PR.
- [ ] **C3. `.mailmap`.** Map Grant's four commit identities to one (`adamsgd@uw.edu`,
  `grantadams60091@gmail.com` as both "Grant Adams" and "Grant.Adams", `grant.adams@noaa.gov`),
  and `grantdadams` to the same person.
  - **Done when:** `git shortlog -sne` shows one line for Grant.
- [ ] **C4. README elements the NMFS guide expects for an R package:**
  - badges, description, install, documentation link, issues link, authors, citation, and the
    disclaimer. Most of these already exist.
  - Add a short "Citation" section that points to `citation("Rceattle")` (there is already an
    `inst/CITATION`).
  - Optional: add a `CITATION.cff` so GitHub shows "Cite this repository".
- [ ] **C5. Obvious README fixes that don't depend on the destination:** change `/tree/master/`
  and `/blob/master/` to `main`.

---

## Stage D. Prepare the link-update branch, but do not merge it [agent]

Build this once A1 is settled. Merge it only after the transfer (PLAN section 3).

- [x] **D1. DONE at 5.45.0** (done in the release branch, not a separate `chore/` branch).
      Changed every `grantdadams`
  reference to the new owner, in these files:
  - `DESCRIPTION`: `URL:` and `BugReports:`
  - `_pkgdown.yml`: `url:` (the new Pages address) and the three navbar and home links
  - `README.md`: badges (R-CMD-check, test-coverage, Codecov), install lines, examples links,
    wiki links
  - `CONTRIBUTING.md`
  - `inst/RELEASE-CHECKLIST.md`, including the section 4 reproducible-install command
  - `R/0-rceattle_class.R` (then regenerate `man/Rceattle-package.Rd` and
    `man/print.Rceattle.Rd` with roxygen; do not hand-edit Rd)
  - `vignettes/introduction.Rmd`
  - `vignettes/articles/developer-guide.Rmd`
  - `examples/Install_Rceattle.R`

  Leave `NEWS.md` alone; it is history.
- [x] **D2. DONE at 5.45.0.** Only this file and `PLAN-adoption-and-NOAA-transfer.md`
      still say `grantdadams`, and both do so to describe the move. All four rewritten
      URLs return 200; `grantdadams.github.io/Rceattle` returns 404, so the old
      website link was already dead when this shipped. Original check: `git grep -n grantdadams -- . ':!NEWS.md'` should return
  nothing, or only lines deliberately kept, each with a comment explaining why.
- [ ] **D3. Record the matching edits in the other repositories**, but do not make them yet:
  - `../Rceattle-models/Rceattle install.R` and its README
  - the READMEs of `../GOA-ATF-ESP` and `../GOA-multispecies-assessment`
  - `inst/dev/SIBLING-REPOS.md`
  - the readiness workbook and decks in `~/Claude/Projects/FIMS overview/`

---

## Stage E. Clean up branches [Grant decides, agent executes]

Branches move with the repository, and the first thing the new organization sees is the branch
list. Nothing is lost if each branch gets an archive tag first.

- [ ] **E1. Decide each remote branch:** keep, archive (tag, then delete), or delete. Current
  list:

  | Group | Branches |
  |---|---|
  | Keep | `main`, `dev`, `gh-pages` |
  | Parked by decision (`SESSION_HANDOFF.md`) | `sel-penalty-form`, `dsem-v5-integration` |
  | Old lines | `depricated-ICES2024-and-CJFAS2025`, `depricated-ceattle_classic`, `depricated-ceattle_classic_reorganized`, `dev-DSEM`, `dev-RTMB`, `dev-cod-bridge`, `Hake_test` |
  | Probably merged | `chore/record-inert-guards`, `chore/remove-qar1`, `docs/consolidate-dev-notes`, `docs/contributor-path`, `docs/correct-dev-notes`, `docs/handoff-release-state`, `feat/osa-cdf-method`, `osa-cdf-method`, `feat/sel-nonparametric-integrable`, `fix/ease-of-use`, `fix/silent-wrong-numbers` |
  | Unmerged remedy | `ci/macos-libomp` |

- [ ] **E2. Archive and delete [agent, after E1].**
  - For each archive: `git tag archive/<branch> origin/<branch>`, push the tag, then delete the
    remote branch.
  - Confirm each "probably merged" branch with `git branch -r --merged origin/dev` before
    deleting it.
  - **Done when:** the remote shows only the branches marked keep or parked.

---

## Stage F. Account settings [Grant]

These are GitHub settings; no code changes.

- [ ] **F1. Set an account successor** (GitHub Settings -> Account -> Successor settings).
  - The NMFS guide asks this of anyone holding NOAA work in a personal account.
  - It protects the repository even if the transfer slips.
- [ ] **F2. Make the NOAA address the primary and notification email** on the GitHub account,
  and set `git config --global user.email grant.adams@noaa.gov` for NOAA work.
- [ ] **F3. Check that two-factor authentication is on,** and that the profile lists the NOAA
  affiliation.
- [ ] **F4. Record the external services tied to the personal account,** so each can be
  re-linked after the move:
  - Codecov (the `CODECOV_TOKEN` secret moves with the repo, but the Codecov app must be
    installed on the organization)
  - GitHub Pages settings
  - any webhooks under Settings -> Webhooks
  - any personal access tokens used in scripts

---

## Stage G. Communication, drafted but not sent [either]

- [ ] **G1. Draft the user notice.** Three sentences:
  - the new URL
  - old install commands and clones still redirect, but point clones at the new URL
  - the tag to pin
- [ ] **G2. Recipient list:**
  - the GOA arrowtooth and GOA pollock author teams
  - Cole and Matt
  - the FIMS team contact
  - the authors of `Rceattle-models` stocks
  - anyone who has opened an issue
  - Melissa (FYI)
- [ ] **G3. Standing note now the transfer is done:** never recreate or fork anything at
  `grantdadams/Rceattle`. Reusing the old name permanently deletes GitHub's redirects, and every
  pinned `install_github("grantdadams/Rceattle@X.Y.Z")` stops resolving. Use branches in the
  organization repo, or rename any fork.

---

## Go / no-go before clicking Transfer

- [ ] A1–A7 are answered, and the organization has granted repository-creation permission.
- [ ] `main` is tagged at the latest release, the GitHub Release is published, and the site
      rebuilt from it.
- [ ] `golden` is green on `main`, or its known failure is explicitly accepted for this transfer.
- [ ] `LICENSE.md`, `CODEOWNERS`, and `.mailmap` are on `main`.
- [ ] Every README install command runs as written.
- [ ] `chore/move-to-<org>` is ready and reviewed, and D2 comes back clean.
- [ ] Branches are cleaned up (E2), or deliberately left.
- [ ] A successor is set (F1), and the external services are listed (F4).
- [ ] The user notice (G1) is drafted, and the date (A7) is on the calendar.

When every box is ticked, follow PLAN section 2.
