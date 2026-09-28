# Plan: adoption fixes and transfer to a NOAA organization

State and plan, not policy. **The transfer happened on 2026-09-26: the repository is now
`afsc-assessments/Rceattle`, and section 1's link sweep shipped at 5.43.0.** Sections 2 and 3 are
kept as the record of what moved and what broke. Written 2026-09-23 from a review of `dev` (5.41.0),
the public site (5.33.0), `../Rceattle-models`, and the AFSC platform readiness workbook. It sits **after** the
release sequence in `SESSION_HANDOFF.md` and does not replace it. Companion documents, for
leadership and for Grant, are in `~/Claude/Projects/FIMS overview/`:
`Rceattle_2027_leadership_onepager.docx` and `Rceattle_working_plan_bridges_and_friction.docx`.

**Why.** Rceattle can run 26 of the 27 age-structured AFSC assessments now. GOA arrowtooth passed
Plan Team and SSC in 2024, and GOA pollock is next. Two things stand in the way of wider use:

- **It is hard to start.** New users cannot see output before installing, cannot import an SS3
  model, and must compile from source.
- **It looks like one person's project.** It lives on a personal account, and about 1,213 of
  1,222 commits since Sept 2024 are Grant's. NMFS guidance says mission-critical code should not
  live in a personal account (see "Sources").

Moving it to a NOAA organization is the cheapest fix for the second problem and makes the first
easier (org-level r-universe, CODEOWNERS, shared admin).

---

## 0. Decisions Grant makes first

These gate everything below. Agents: do not pick these.

1. **Destination organization.** The options:
   - **`afsc-assessments` (recommended).** It already hosts the AFSC assessment code that
     analysts use (`GOApollock`, `ebswp`, `AMAK`, `afscdata`, `sara`, the assessment SOP), so
     the package sits where its users already look.
   - **`nmfs-fish-tools` (Fisheries Integrated Toolbox).** More national visibility, but more
     formal. Rceattle can be listed in the toolbox later without moving there.
   - **A new `Rceattle` organization.** This gives shared admin but no institutional signal,
     so it fixes little.

   Ask the `afsc-assessments` owners before starting. Two things to confirm: whether it is an
   NMFS GitHub Enterprise Cloud (GHEC) organization, and whether it allows the third-party
   Actions the workflows use (`r-lib/actions`, `JamesIves/github-pages-deploy-action`,
   `codecov/codecov-action`). If it is GHEC, members need the signed GHEC user agreement and
   SSO. Non-NOAA contributors (UW, international users) can still be outside collaborators,
   and a public repo stays readable by everyone.

2. **License.** There is **no `LICENSE` file**; only `License: GPL (>= 2)` in `DESCRIPTION`.
   - NMFS guidance recommends Apache 2.0, but Rceattle links TMB (GPL), and Apache 2.0 cannot
     be combined with GPL-2. Relicensing would also need every copyright holder's consent.
   - **Recommendation:** keep GPL (>= 2) and add the standard `LICENSE` file. Confirm with the
     NMFS GitHub Governance Team that GPL satisfies NAO 201-118 for this package.

3. **Co-maintainer.** Name one person with admin rights and review duty.
   - Cole is the natural candidate: he has commits and shares the DSEM recruitment/M/growth
     lane.
   - This also answers leadership's continuity question and is the ask in the one-pager.

4. **`Rceattle-models`.** NMFS guidance says to transfer what a successor would need, and
   several stocks' bridging history lives there. Recommend moving it to the same organization
   and keeping its "not the operational models" note.

5. **Timing.** Transfer after the pending `dev` -> `main` release is tagged and after `golden`
   is robust (section 1). Avoid the two weeks before a Plan Team document deadline. Installs
   keep working through the move (see section 2), but the docs site URL changes.

---

## 1. Before the transfer: finish what `SESSION_HANDOFF.md` already queues

**The working checklist, with owners and a go/no-go list, is `TODO-pre-transfer.md`.** What
follows is the summary. In this order:

1. **Ship the 5.34.0–5.41.0 release.** Tag the merge commit on `main` and publish a GitHub
   Release (`SESSION_HANDOFF.md`, "The release sequence").
2. **Make `golden` robust.** It fails on `main` at 5.33.0 because `goa_ss` lands in a second
   local minimum. Until this is fixed, `deep-checks` cannot gate a release.
   - This matters more for adoption than any docs work: "bomb-proof" is the pitch.
   - The fix is a harness change (warm start or best of two), so it cannot move a fitted number.
3. **Apply the tag convention everywhere.** Done, verified 2026-09-27: `inst/RELEASE-CHECKLIST.md`
   section 3 sets it -- bare `X.Y.Z` with no `v`, where only the first tag, `v4.3.0`, differs --
   and `R/0-rceattle_class.R:12` now reads `@X.Y.Z`. What remains is the tag itself: nothing
   tagged sits above 5.28.0 while `main` is 5.33.0.
4. **Fix the pin in the README.** Done at 5.43.0. It read
   `remotes::install_github("grantdadams/Rceattle@4.3.0")`, which fails -- the tag is `v4.3.0`,
   and it was years out of date. `README.md:43` now pins `afsc-assessments/Rceattle@5.43.0`, valid
   the moment the release tag is pushed.
5. **Add `LICENSE`** (decision 2) and **`.github/CODEOWNERS`** listing Grant and the
   co-maintainer.
   - Do not turn on required reviews on `main`/`dev` until the co-maintainer exists, or you will
     block yourself.
6. **Account hygiene** (GitHub settings, no code):
   - Set a successor on the personal account now; the NMFS guide asks for this and it protects
     the repo even before the move.
   - Commit with the NOAA email (most commits use `adamsgd@uw.edu`).
   - Add a `.mailmap` so Grant's four identities count as one.
7. **Prepare the link-update branch but do not merge it yet.** 36 lines in 11 tracked files
   reference `grantdadams` (NEWS.md excluded):
   - `DESCRIPTION` (URL, BugReports)
   - `_pkgdown.yml` (url, navbar links)
   - `README.md` (badges, install lines, examples links)
   - `CONTRIBUTING.md`
   - `inst/RELEASE-CHECKLIST.md`
   - `R/0-rceattle_class.R`
   - `man/Rceattle-package.Rd` and `man/print.Rceattle.Rd` (regenerate; do not hand-edit)
   - `vignettes/introduction.Rmd`
   - `vignettes/articles/developer-guide.Rmd`
   - `examples/Install_Rceattle.R`

   Check the list with `git grep -n grantdadams -- . ':!NEWS.md'`.

---

## 2. The transfer itself (about 30 minutes)

Settings -> Danger Zone -> Transfer -> the chosen organization. Keep the name `Rceattle`. Grant
needs permission to create repositories in that organization.

**Moves with the repo:**

- issues, PRs, wiki, stars, watchers
- releases and tags
- Actions secrets, including `CODECOV_TOKEN`
- web and git redirects: `remotes::install_github("grantdadams/Rceattle")` and existing clones
  keep working through the redirect

**Breaks:**

- **The GitHub Pages site is not redirected.** `grantdadams.github.io/Rceattle` stops serving (it
  404s today), and the site moves to `afsc-assessments.github.io/Rceattle`.
- **Codecov** needs the Codecov app installed on the organization, and possibly a new token.
- **Organization Actions policy** may block the third-party Actions (decision 1). If the deploy
  action is blocked, `actions/deploy-pages` is the first-party replacement.

**Trap: never create a repository or fork named `grantdadams/Rceattle` now that the move has
happened.** GitHub permanently deletes the redirects when the old name is reused, and every pinned
`install_github("grantdadams/Rceattle@X.Y.Z")` in a consumer script stops resolving. Anyone forking
back to a personal account for PR work must rename the fork (for example `Rceattle-fork`), or
better, use branches in the organization repo.

---

## 3. Immediately after the transfer

1. **Merge the link-update branch** from section 1, step 7, and set `_pkgdown.yml` `url:` to the
   new Pages address. Re-enable Pages on the organization repo and confirm the site builds.
2. **Point every local clone at the new URL:**
   `git remote set-url origin https://github.com/<org>/Rceattle.git`
   Do this for Grant's clones and for the sibling repos in `inst/dev/SIBLING-REPOS.md`. Update
   `../Rceattle-models/Rceattle install.R`.
3. **Update external links:** the Rceattle-models README, the GOA-ATF-ESP and
   GOA-multispecies-assessment READMEs, the readiness workbook, the leadership decks, and
   `grantdadams.wordpress.com`.
4. **Tell users once, briefly:**
   - the GOA arrowtooth and GOA pollock author teams
   - Cole and Matt (DSEM lanes)
   - the FIMS team
   - anyone who has filed an issue

   The message: the new URL, that old install commands still redirect, and the tag to pin.
5. **Add Rceattle to the organization profile README**, and consider listing it in the
   Fisheries Integrated Toolbox.

---

## 4. Package roadmap after the move

Detail and evidence are in the friction-audit doc. `CONTRIBUTOR-EXPERIENCE.md` covers people
**extending** the package; this list covers people **adopting** it. Its item 0 still applies
here: **ask two or three would-be users where they stopped, then re-order this list.**

### Quick wins (about a week, mostly configuration and wording)

1. **Make the site show output.**
   - The problem: `pkgdown.yaml` never sets `RCEATTLE_EVAL_VIGNETTES`, so every article on the
     live site is code with no plots or tables.
   - The fix: set it in the pkgdown job (`vignettes.yaml` has the runtimes), or precompute slow
     articles with the `.Rmd.orig` pattern.
   - **Acceptance:** the live Introduction page shows figures.
2. **Update the package identity.** The `DESCRIPTION` Title and Description are 2016-era text
   ("three target species", ration functions), and the site title says "Multispecies Stock
   Assessment in R". Reuse the README tagline, which leads with single- and multispecies
   operational use.
3. **Add an "AFSC assessments" article.**
   - Lead with the SSC-accepted GOA arrowtooth bridge (`GOAatf2023`), then `Atka2022`,
     `NorthernRockfish2022`, and `GOAcod`, all of which already ship with the package.
   - Link `Rceattle-models` from the README and the site; neither links it today.
4. **Fix the `nages` label.** `stock-synthesis-conversion.Rmd` line 42 maps SS `Nages` to
   `nages` and calls it "Maximum age class". CONTRIBUTING rule 5 says `nages` is a count of
   bins, and the two agree only when `minage = 1`. Add a worked example with an SS3 model that
   starts at age 0.
5. **Use string switch values in user-facing docs** (`estimateMode = "Hindcast"`, not `0`). The
   aliases already exist.
6. **README links:** change `/tree/master/` to `/tree/main/`.
7. **Release notes:** open each release note with three to five "what changes for you" bullets,
   and put the same bullets in the GitHub Release. `NEWS.md` is 7,355 lines, too long for a
   user to scan.

### One to two months

8. **Publish on r-universe** under the organization, as FIMS does, so Windows and macOS users
   get prebuilt binaries instead of compiling TMB. It needs a `packages.json` registry repo in
   the organization. This overlaps `CONTRIBUTOR-EXPERIENCE.md` item H.
9. **`read_ss3()` plus `bridge_compare()`.** This is the largest adoption lever.
   - `read_ss3()`: use `r4ss::SS_read` (in Suggests) to build a `data_list`. Nine of the 27
     AFSC stocks are SS3.
   - `bridge_compare()`: produce the SSB, recruitment, and likelihood comparison against the
     source model that a Plan Team bridging appendix needs.
   - Then write the same pair for the shared ADMB codebases: AMAK (Atka mackerel and AI pollock)
     and the flatfish `fm.tpl` (northern rock sole and yellowfin sole).
   - **This adds API and needs Grant's sign-off on the design before any code.**
10. **Rebuild the SS3 vignette** around an AFSC stock (the GOA cod bridge exists in
    Rceattle-models) and add a bridge-check section built on item 9.
11. **CRAN submission.** `cran-comments.md` was last measured at 5.15.0. Remove `Remotes:`
    before submitting. The NMFS guide names CRAN or rOpenSci review as the standard review for R
    software used for advice.

### Ongoing: operational trust

12. **Pin behaviour, not just API, for accepted models.**
    - `SIBLING-REPOS.md` records two silent changes: `initMode = 1` changed meaning, and
      environmental SRR indices were inert from 4.4.0 through 5.31.0 while NEWS said they worked.
    - Before recommending Rceattle more widely, add a golden fit of the SSC-accepted GOA
      arrowtooth configuration, and of GOA pollock once it is accepted, so behaviour drift
      turns a test red.
    - **Confirm the accepted arrowtooth model never used environmental SRR indices.**
13. **A "Behaviour changes" section in every release,** separate from "Breaking changes".
14. **Show one deprecation message per session** instead of the current default of none
    (`Rceattle.warn_deprecated_args = FALSE`), so old scripts keep running but users learn.
15. **Check that the README's wiki links** (Onboarding, Workflow for updating) still resolve.
    If they duplicate CONTRIBUTING and the developer guide, fold them into the site.

---

## Acceptance for the whole plan

- [ ] Repo, site, and issues live under a NOAA organization with at least two owners.
- [ ] `LICENSE`, `CODEOWNERS`, and a tagged GitHub Release for every version on `main`.
- [ ] `golden` passes on `main` and gates releases.
- [ ] The live site shows figures, and the README install and pin commands work as written.
- [ ] An analyst who has never used Rceattle installs a binary, finds the arrowtooth example,
      and fits it without asking Grant. (Same success test as `CONTRIBUTOR-EXPERIENCE.md`:
      one observable event.)

## Sources

- NMFS Open Science GitHub Guide, sections 7.3–7.5 and 10.1–10.3 (repositories under individual
  accounts, organization settings, licenses): https://nmfs-opensci.github.io/GitHub-Guide/
- GitHub Docs, "Transferring a repository" (what moves, Pages not redirected, redirect deletion):
  https://docs.github.com/en/repositories/creating-and-managing-repositories/transferring-a-repository
- `afsc-assessments` organization: https://github.com/afsc-assessments
