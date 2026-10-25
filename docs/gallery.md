# Publish reviewed simulation figures

`gallery/` is the public static website. Its landing page links to separate
result collections, beginning with `gallery/maintenance-01/`. The Pages workflow
uploads only this directory. It checks the assets but does not run Julia,
submit cluster jobs, read model checkpoints, or render simulations in Actions.

**Unity completion does not automatically update this website.** A completed
stage must be collected, transferred, rendered and published. The homepage
records publication status rather than live scheduler status.

`gallery/stabilization-validation-01/` records the three 400-second CPU/GPU
checks from Unity array 66081490. Their collection archive was verified against
SHA256 `54e39721cc04f9307098823d8e0dc72f98348f6277e8480abd924d4e556a879b`.
They are software validation, separate from the longer scientific experiments.
For `purpose = "smoke"`, the renderer verifies the saved validation report,
success marker, configuration/source/dependency/validator hashes, and reported
initial/final comparison tolerances before displaying validation success.
Only an allowlisted summary of that report is published; raw run files remain
in the research archive.

`gallery/rostami-validation-01/` records the three 512², 400-second
CPU/GPU checks from Unity array 66091773. The archive SHA256 is
`23002cacb2628194b81209b726c7618aed93b668f1f6b348bff88c3ebfd1ccd0`.
The maximum final relative L2 field difference is `5.074531993358771e-16`.
The jet-only and reference cyclone cases have different upper flows; this
collection is software validation, not a controlled stabilization comparison
or a reproduction of the paper's mature hexagon. Longer runs are separate.

## From a cluster run to a local preview

Run the simulations with the existing cluster launchers. For a completed Unity
experiment stage, use the corresponding `experiments_01.sh collect` command
described in [the experiment protocol](experiments-01.md). Download that archive
to your laptop and extract it into a local data directory **outside this public
checkout**. Keep the original archive, configs, provenance and raw JLD2 files
there for scientific analysis. GitHub authentication stays on your laptop;
compute nodes need no GitHub token, SSH publishing key, or repository credentials.

The renderer requires Python 3.8+ and the packages listed in
`gallery-requirements.txt`. From the repository root, an optional dedicated
environment outside the checkout keeps these dependencies separate:

```sh
python3 -m venv ../saturn-gallery-venv
source ../saturn-gallery-venv/bin/activate
python3 -m pip install -r gallery-requirements.txt
```

From the local repository root, render the selected `results` directories into
a **new** collection directory. Replace the example input paths with the actual
downloaded run directories:

```sh
python3 scripts/build_gallery.py \
  /path/to/downloaded/run-a/results \
  /path/to/downloaded/run-b/results \
  --output gallery/experiments-01 \
  --title "Lower-layer experiments: initial screen"
python3 scripts/check_gallery.py gallery
python3 -m http.server 8000 --bind 127.0.0.1 --directory gallery
```

Add `--mp4` to the render command if you also want MP4 playback. This requires
a local `ffmpeg` executable; `--ffmpeg /path/to/ffmpeg` selects one explicitly.
The movie displays the actual saved frames at two frames per playback second,
with no temporal interpolation. The interactive frame viewer is available
without ffmpeg. Choose the movie option before rendering into the new output
directory; the renderer refuses to overwrite an existing collection.

Open `http://127.0.0.1:8000/experiments-01/` to review the collection; stop the
preview server with Ctrl-C. Inspect the figures, units, times, captions,
parameter labels, completion status and scientific limitations. A rendered
movie is not evidence of numerical convergence or a stable hexagon. Review
the model and physical-regime caveats alongside each comparison.

Add a relative link from `gallery/index.html` to the new collection, for example
`experiments-01/index.html`, and preview the landing page too. Each collection
is self-contained. Preserve the existing directories and landing-page links;
do not replace the whole gallery when adding a result. For revised rendering,
use a new name such as `experiments-01-v2` and choose which version the landing
page links to after review.

## Commit and publish explicitly

Run these from the repository root on `main`, after configuring the intended
GitHub repository as `origin`. Check the branch and remote before proceeding:

```sh
git status --short --branch
git remote -v
python3 scripts/check_gallery.py gallery
git add -- gallery/experiments-01/ gallery/index.html
git diff --cached --stat -- gallery/experiments-01/ gallery/index.html
git diff --cached -- gallery/experiments-01/ gallery/index.html
```

Review the staged files and the browser preview. The explicit paths stage only
this addition and its landing-page link. Once they are ready, commit only those
paths and push deliberately:

```sh
git commit --only -m "Add reviewed initial-screen gallery" -- \
  gallery/experiments-01/ gallery/index.html
git push origin main
```

The path-limited commit leaves unrelated staged changes out of this commit.
The push publishes all commits currently ahead of `origin/main`, so review
that branch history before pushing. No publishing script commits, pushes,
deletes old collections, or submits a simulation automatically.

## Enable Pages once the repository exists

In the GitHub repository, select **Settings → Pages → Build and deployment →
Source: GitHub Actions**. Then a push to `main` that changes `gallery/**` or the
Pages workflow triggers deployment. **Actions → Publish gallery → Run workflow**
can deploy the current `main` manually. Inspect the workflow result and its
`github-pages` deployment URL; local preparation alone does not establish that
Pages is enabled or that a remote deployment succeeded.

The workflow follows GitHub's [custom Pages workflow documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages):
checkout v6, configure-pages v5, upload-pages-artifact v4 and deploy-pages v4,
with `contents: read`, `pages: write`, `id-token: write`, and the `github-pages`
environment. Only committed `gallery/` assets enter the deployed artifact.

## What the local and Actions check permits

`scripts/check_gallery.py` is read-only and uses Python's standard library.
It accepts ordinary HTML/CSS/JavaScript, JSON, CSV, TOML, text and Markdown, images,
and MP4/WebM files. It rejects unsupported formats including raw JLD2/HDF5
states and archives, symbolic or hard links, hidden files other than
`.nojekyll`, credential-like filenames, and common private-key/token signatures.
Each file must be at most 90 MB; the full gallery must be smaller than 1 GB.
The renderer's downloadable TOML contains only recognized model configuration
fields. An extensionless `COMPLETED` file is accepted only when its full content
is exactly `Static gallery build completed.` followed by one newline; model
completion markers and other extensionless files are rejected.
These are deliberately conservative project limits within the documented
[Pages limits](https://docs.github.com/en/pages/getting-started-with-github-pages/github-pages-limits)
and [artifact requirements](https://github.com/actions/upload-pages-artifact#artifact-validation).

The checker cannot recognize every credential or private annotation. Review
captions, JSON metadata and text before committing, and retain raw scientific
data outside the repository. Compress or shorten a derived movie if it exceeds
the asset limit; use an appropriate data archive for full simulation output.

Run the focused checker tests locally with:

```sh
python3 test/test_gallery_checks.py
```
