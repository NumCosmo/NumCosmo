# Reviewing contributions to NumCosmo

This guide is for whoever reviews a pull request, a maintainer or an LLM agent acting for
one. [CONTRIBUTING.md](CONTRIBUTING.md) and [TESTING.md](TESTING.md) say how a change must
be made. A review checks that, and first decides whether the change is needed and in what
form. Scope decisions belong to the maintainer: the review proposes, it does not cut.

The sections below follow the order of review. Each states a principle and the reason for
it, then what the reviewer checks. Stop at the first step that changes the shape of the
pull request, since the later steps apply to the reshaped change.

## 1. Necessity

Every part of the repository has a cost after it is merged: build and test time, review
of later changes, and maintenance when the code it depends on changes. A contribution is
accepted for what it adds that the repository does not already provide, and that addition
must justify its cost.

For each component of the pull request (a function, a test, a page, a helper module),
state in one sentence what it adds. Search the existing code, tests and documentation for
the same quantity or operation before deciding. A component that only repeats what exists
is removed; a component that repeats it and adds one part is reduced to that part. When a
component adds build or test time, measure the time and weigh it against what the
component adds.

## 2. Proportion

Each part of a change should serve the purpose of the file it is in. Repeated text and
code make a reader read the same thing several times, and the copies drift apart when one
of them is edited. Material whose subject belongs elsewhere hides the subject of the
file.

Check for content stated more than once (across files as well as within one), content
that restates what another part of the repository already establishes, and content about
a subject outside the file's purpose. Shared material is defined once and linked; code
used by several tests, tools or pages becomes one module. A proposal to trim names what
stays, what goes, and the size of each.

## 3. One concern per commit

A commit that does one thing can be reviewed, bisected and reverted on its own. The
history on `master` describes the change as it is, not how it was developed.

Check that each commit has one concern and builds on its own, and that the commits come in
dependency order: library code with its tests and regenerated stubs, then code that uses
it, then documentation. Fixup commits are folded into the commits they fix, and merge
commits from `master` are removed by rebasing. Separate commits in one pull request are
enough when the maintainer merges with a merge commit. A separate pull request is
warranted when a piece is useful and reviewable on its own, or when it needs a decision of
its own, such as a change of results.

## 4. Claims

A reader acts on the numbers and explanations in a contribution without re-deriving them.
Each must therefore be reproducible, and each stated cause must have been tested.

- A number in a commit message, comment, docstring or page is computed by code the reader
  can run, or comes with the command that measured it. Documentation pages compute their
  numbers when they render.
- A stated cause ("X is large because of Y") is accepted only after a measurement that
  varies Y. A cause also predicts every other case with the same property; check those
  cases, not only the one that motivated the explanation.
- A tolerance comes from a measured value and the test states how it was set. A tolerance
  is never widened to make a test pass.
- A change of results is measured against a build of `master` on the affected functions:
  the number of values that differ and the largest relative difference.

## 5. Placement

Code that is not where the repository's conventions put it may not run, may not be found,
and may not be maintained.

- Tests mirror the source path and carry the marker of their lane; a test placed in
  another lane's file runs only when that lane runs.
- A check that uses only the C API is a C test, without a Python copy.
- Tests and pages call the library; a re-implementation of a formula the library evaluates
  can disagree with it and tests nothing about it.
- A private interface of another code is used only together with a test that pins it, so
  that a change in that code fails a test instead of a computation.
- Shared code is a `numcosmo_py` module that takes its inputs as arguments; it does not
  read global variables of its caller, and it is not reached through a `sys.path` insert.

## 6. Writing

Documentation and comments are read by people who hold much of the library in mind at
once, so every sentence should be readable once, at speed. Apply "Writing style for
comments and documentation" in CONTRIBUTING.md and run
`.github/scripts/check_doc_style.sh`. In addition, flag sentences that only announce the
next one, labels invented for things that have a technical name, prose wider than the
file's wrapping width, and text about how the change was developed rather than what it
does.

## 7. Commits and attribution

Check commit messages and the pull request description against "Commit messages" in
CONTRIBUTING.md, in particular bullets that carry several changes, long bullets, prose
paragraphs and development history.

Authors are the people who answer for the content. When commits are rewritten during
review, the original author stays the author, and others who wrote part of the final
content are added with `Co-authored-by:` trailers. A review judges the claims and who
answers for them, never which tools produced the text: whether a text looks
machine-written is not a finding.

## Reporting a review

Rank the findings in the order of the sections above. For each finding give the file and
line, what is wrong, the evidence (a command, a measurement or the code that contradicts
it), and the proposed fix. Keep what was measured separate from what was inferred.
