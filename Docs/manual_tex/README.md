Execute:

```
bash get_manual.sh
```

from this directory to generate a siesta.rst file from siesta user manual written in LaTeX.

There are several things involved to make the conversion. 
`pandoc` is used to convert .tex to .rst.
However, the out-of-the-box conversion is not perfect and we need to do several things to smooth the edges:

- Preprocess tex files to cut tex that is too complex for `pandoc` to parse. This is done in `cut_pandoc_file.py`.
- Redefine some latex commands to make the job simpler for `pandoc`. This is done in `sphinx_commands.tex`.
- Preprocess tex files to substitute some command calls. This is done inside `get_manual.sh` using `sed` calls.
- While converting, filter the conversion to process some of the structures that `pandoc` parses. This is done in `pandoc_filter.py` using the `pandocfilters` python package.
- Post process the output rst to remove artifacts of pandoc. This is done in `get_manual.sh` using `sed` calls.

The `rst_siesta_manual.tex` file is the skeleton that will merge all the sections.
In fact, the conversion is done like `pandoc rst_siesta_manual.tex --filter pandoc_filter.py -o siesta.rst`
