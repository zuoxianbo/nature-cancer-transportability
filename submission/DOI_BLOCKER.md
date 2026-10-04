# Code availability DOI — BLOCKER

**Status**: BLOCKED. DOI cannot be minted from this machine. Credentials required from the author.

## Current placeholders (5 locations)

| # | File | Location | Current text |
|---|---|---|---|
| 1 | `NatureCancer_Final_Submission_Manuscript.docx` | P52 | "The analysis code, manifest and environment specification are prepared for deposition in a versioned public repository; the permanent repository identifier will be inserted before submission." |
| 2 | `NatureCancer_Reporting_Summary_Final_Template.docx` | P18 | (repository identifier placeholder) |
| 3 | `NatureCancer_PreSubmission_Editorial_Reviewer_Revision_Notes.docx` | P31, P56 | (repository identifier placeholders) |
| 4 | `NatureCancer_Cover_Letter_Final_Template.docx` | P9 | (repository identifier placeholder) |
| 5 | `README.txt` | bullet 2 | (repository identifier placeholder) |

## What is missing on this machine

- **GitHub CLI** `gh` — not installed (`command -v gh` → empty).
- **Zenodo CLI** `zenodo` — not installed.
- **`GITHUB_TOKEN` env** — unset.
- **`ZENODO_TOKEN` env** — unset.
- **SSH key / keychain entry** for `github.com:zuoxianbo/...` — none.
- Existing template `zenodo_metadata.json` — -era, needs to be updated to  title/abstract (template below).
- Existing `scripts/github_sync.sh` — -era paths, needs to be updated to point at `_20261004_Nature_Cancer/` (see `github_sync_.sh` in this folder).

## Two acceptable paths — pick one

### Path A — GitHub-only (simplest, ~2 minutes)

The DOI is the repository URL / commit SHA of the public repo `https://github.com/zuoxianbo/nature-cancer-immune-state`. Fill in **P52** with the actual URL (or commit SHA tag ``) once you push.

```
brew install gh
gh auth login                        # interactive browser flow
./__work/_final/github_sync_.sh
```

Then replace placeholder P52 with:
```
The analysis code, frozen manifest and environment specification are available at
https://github.com/zuoxianbo/nature-cancer-immune-state (commit tag ,
SHA <FILL_AFTER_PUSH>).
```

### Path B — Zenodo minted DOI (recommended for Nature)

Nature prefers a DOI to a raw URL. This needs a Zenodo API token.

1. Create a Zenodo API token at `https://zenodo.org/settings/tokens` (scopes: `deposit`).
2. Export `ZENODO_TOKEN=...`.
3. Run the minting script `zenodo_mint_.py` (in this folder).

That script will:
- POST `https://zenodo.org/api/deposit/depositions` to create an empty deposit.
- PUT the metadata from `zenodo_metadata_.json`.
- Upload the package ZIP `_20261004_Nature_Cancer.zip` (or file-by-file) to the deposit bucket.
- POST `https://zenodo.org/api/deposit/depositions/{id}/actions/publish` to publish.
- Print the returned DOI (`https://doi.org/10.5281/zenodo.<id>`).

Then replace placeholder P52 with:
```
The analysis code, frozen manifest and environment specification are available at
https://doi.org/10.5281/zenodo.<ID> (Zenodo record, published <DATE>).
```

## Insert-script for the 5 docx placeholders

`__work/_final/insert_doi_.py` reads `NEW_DOI.txt` from this folder and rewrites all 5 placeholders. Run it after minting.

## Author information required before minting

The `zenodo_metadata_.json` template in this folder has `[AUTHOR TO CONFIRM]` markers for:
- Full name (family, given)
- ORCID
- Institutional affiliation
- Funding acknowledgement (grant number)

These four fields cannot be filled in automatically — they are personal facts.
