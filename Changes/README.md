
## Changelog

This directory contains *news fragments* which are short files that contain a
small *Markdown*-formatted text that will be added to the next release page.

Make sure to use full sentences with correct case and punctuation, and please
try to use backticks. The fragments are bulletized, so just write what change
has been introduced, no heading!

Each file should be named like `<PULL REQUEST>.<TYPE>.md`, where
`<PULL REQUEST>` is a pull request number, and `<TYPE>` is one of:

- `feat`: New user facing features, new options or usage
- `build`: New changes/features/bug-fixes for the build-system
- `fix`: A fix for the code-base, could be a bugfix, or behavioral.
- `change`: Changes to API, and other operational changes.
- `doc`: Changes to the documentation
- `highlight`: Adds a highlight bullet point to use as a possibly highlight

It is possible to use the same pull request number for several categories
(e.g. `123.feat.md` and `123.change.md`), each containing a separate entry.
For example a new feature might change the API of other related functions.

Most categories should be formatted as paragraphs with a heading (just like
a `git commit` message).
So for example: ``123.feat.md`` would have the content::

    Enabled *new feature* in the *old solver*

    The *new feature* option is now available for *old solver*.
    To use it, add this to your input fdf file: `OldSolver.NewFeature ...`

`highlight` is usually formatted as bulled points making the fragment
`- This is a highlight`.

If you are unsure what pull request type to use, don't hesitate to ask in your
MR.


### Utilities

If the changes are related to utilities, then please follow the recipe above,
but put them in the `utility` sub-directory. In this way they'll show up as
a separate section in the final release notes.


## Developers

``towncrier`` (a Python package) is required to build the docs; it will be automatically run when
you build the docs locally. You can also run ``towncrier
build --draft --version 1.18`` if you want to get a preview of how your change
will look in the final release notes.

> [!NOTE]
> This README was adapted from the sisl changelog readme, which again was
> adapted from numpy:
> This README was adapted from the numpy changelog readme, which again was
> adapted from pytest:
> This README was adapted from the pytest changelog readme under the terms of
> the MIT licence.
