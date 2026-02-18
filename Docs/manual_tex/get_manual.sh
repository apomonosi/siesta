
TEX_ROOT=$(pwd)/$(dirname $0)
SIESTA_DOCS=${TEX_ROOT}/..
echo "Pandoc manual conversion. Using TEX_ROOT=$TEX_ROOT"
echo "Siesta Docs sources in tex in: ${SIESTA_DOCS}"
ls ${SIESTA_DOCS}/tex

# The BSD and GNU sed commands inplace option works differently,
# therefore we need to use a temporary file
SED_TEMP_FILE="${TEX_ROOT}/tmp_file"
sed_in_place() {
    # Get the last argument, which is the name of the file to modify
    for file in "$@"; do true; done

    # Use sed to perform the command and write to a temp file, then move it back
    sed "$@" > ${SED_TEMP_FILE} && mv ${SED_TEMP_FILE} "$file"
}

prepare_tex_file(){
    # Makes sure a tex file is ready for conversion.
    # First argument is the file to clean
    # Second argument is optional. If it is "helper", we will interpret that we are processing a helper file.

    # If second argument is true
    if [ "$2" = helper ]; then
        sed_in_place \
            -e 's/\\begin{fdfentry}{\(.*\)}\[\(.*\)\]<\(.*\)>/\\begin{fdfentry}{\1}{\2}{\3}/g' \
            -e 's/\\begin{fdfentry}{\(.*\)}\[\(.*\)\]/\\begin{fdfentry}{\1}{\2}{}/g' \
            -e 's/\\begin{fdfentry}{\(.*\)}<\(.*\)>/\\begin{fdfentry}{\1}{}{\2}/g' \
            -e 's/\\begin{fdfentry}{\(.*\)}/\\begin{fdfentry}{\1}{}{}/g' \
            -e 's/\\begin{fdfexample}/\\begin{verbatim}/g' \
            -e 's/\\end{fdfexample}/\\end{verbatim}/g' \
            -e 's/\\begin{codeexample}/\\begin{verbatim}/g' \
            -e 's/\\end{codeexample}/\\end{verbatim}/g' \
            -e 's/\\begin{shellexample}/\\begin{verbatim}/g' \
            -e 's/\\end{shellexample}/\\end{verbatim}/g' \
            -e 's/\\begin{output}/\\begin{verbatim}/g' \
            -e 's/\\end{output}/\\end{verbatim}/g' \
            $1
    else
        sed_in_place -e 's/\\begin{fdfentry}{\(.*\)}\[\(.*\)\]<\(.*\)>/\\begin{fdfentry}{\1}{\2}{\3}/g' \
            -e 's/\\begin{fdfentry}{\(.*\)}\[\(.*\)\]/\\begin{fdfentry}{\1}{\2}{}/g' \
            -e 's/\\begin{fdfentry}{\(.*\)}<\(.*\)>/\\begin{fdfentry}{\1}{}{\2}/g' \
            -e 's/\\begin{fdfentry}{\(.*\)}/\\begin{fdfentry}{\1}{}{}/g' \
            -e 's/\\begin{fdfexample}/\\begin{verbatim}/g' \
            -e 's/\\end{fdfexample}/\\end{verbatim}/g' \
            -e 's/\\begin{codeexample}/\\begin{verbatim}/g' \
            -e 's/\\end{codeexample}/\\end{verbatim}/g' \
            -e 's/\\begin{shellexample}/\\begin{verbatim}/g' \
            -e 's/\\end{shellexample}/\\end{verbatim}/g' \
            -e 's/\\begin{output}/\\begin{verbatim}/g' \
            -e 's/\\end{output}/\\end{verbatim}/g' \
            -e 's/\\fdfindex\*/\\fdfindexstar/g' \
            -e 's/\\fdf\*/\\fdfstar/g' \
            $1
    fi

}

# Function that replaces placeholders that have been used to avoid pandoc
# messing with special characters that we want to write literally.
postprocess_rst_specialchars(){
    sed_in_place \
        -e 's/PLACEHOLDERSIESTA_PANDOC_BACKSLASHPLACEHOLDER/\\/g' \
        -e 's/PLACEHOLDERSIESTA_PANDOC_BACKTICKPLACEHOLDER/`/g' \
        "$1"
}

# Function that processes the fdfparam definitions in the rst files adding
# a prefix to the fdf-param variables. This prevents the RTD from having
# duplicated definitions of the same fdfparam variable.
postprocess_fdfparam_rst(){
    sed_in_place -e "s/fdfparam:/fdfparam:${2}/g" "$1"
}


rm -rf ${TEX_ROOT}/tex
# Download
echo "COPYING SIESTA MANUAL'S TEX FILES..."
cp -rp ${SIESTA_DOCS}/tex ${TEX_ROOT}
cp -rp ${SIESTA_DOCS}/tbtrans.tex ${TEX_ROOT}/tex/tbtrans.tex
cp -rp ${SIESTA_DOCS}/siesta.bib ${TEX_ROOT}/siesta.bib

echo "PREPROCESSING TEX FILES..."

# If we don't do this, there is a "siestatextsize" at the beggining of the file
# (because of pandoc's inability to parse complex tex files)
sed_in_place 's/ siestatextsize//' ${TEX_ROOT}/tex/helpers/setup.tex

# Cut the helper tex files if needed to avoid tex that pandoc can't parse
for f in ${TEX_ROOT}/tex/helpers/*; do
    python ${TEX_ROOT}/cut_pandoc_file.py $f
done

# Preprocess the tex files as needed by our redefined latex commands and pandoc's limitations.
for f in ${TEX_ROOT}/tex/helpers/*.tex; do
    prepare_tex_file "$f" helper
done

for f in ${TEX_ROOT}/tex/*.tex ${TEX_ROOT}/tex/sections/*.tex ${TEX_ROOT}/tex/sections/*/*.tex; do
    prepare_tex_file "$f"
done

# Incorporate custom sphinx commands into tbtrans.tex.
# We don't need to do this for the siesta manual because we use a custom tex template.
sed_in_place '15s/^/\\input{.\/sphinx_commands.tex}\n/' ${TEX_ROOT}/tex/tbtrans.tex
# Remove the SIESTA version input from tbtrans.tex to avoid warnings.
sed_in_place 's/\\input{..\/SIESTA.version}//g' ${TEX_ROOT}/tex/tbtrans.tex


# These are for specific variables that bring problems with pandoc
# due to custom latex commands.
sed_in_place -e 's/fdfindexstar{Command line options:-V}/ref{fdfparam:command-line-options:-V}/g' \
             ${TEX_ROOT}/tex/sections/TranSIESTA/Options.tex
sed_in_place -e 's/fdfindex{TS!Elec.<>!Bloch}/ref{fdfparam:ts.elec.<>.bloch}/g' \
             ${TEX_ROOT}/tex/sections/TranSIESTA/Kpoint_sampling.tex

# We need to replace only the second instance here:
sed_in_place -e 's/fdfindex{TS!Elec.<>!Bloch}/ref{fdfparam:ts.elec.<>.bloch}/g' \
             ${TEX_ROOT}/tex/sections/TranSIESTA/Electrodes.tex
sed_in_place -e '0,/ref{fdfparam:ts.elec.<>.bloch}/{s/ref{fdfparam:ts.elec.<>.bloch}/fdfindex{TS!Elec.<>!Bloch}/}' \
             ${TEX_ROOT}/tex/sections/TranSIESTA/Electrodes.tex

# Same for tbtrans
sed_in_place -e 's/fdfindex{TBT.Elec.<>!Bloch}/ref{fdfparam:tbt.elec.<>.bloch}/g' \
             ${TEX_ROOT}/tex/tbtrans.tex
sed_in_place -e '0,/ref{fdfparam:tbt.elec.<>.bloch}/{s/ref{fdfparam:tbt.elec.<>.bloch}/fdfindex{TBT.Elec.<>!Bloch}/}' \
             ${TEX_ROOT}/tex/tbtrans.tex

echo "CONVERTING TEX TO RST WITH PANDOC..."

# Run pandoc
# We need this cd because the \input{} commands in the tex files
# are interpreted relative to the current directory.
cd ${TEX_ROOT}
pandoc ${TEX_ROOT}/rst_siesta_manual.tex --filter ${TEX_ROOT}/pandoc_filter.py --citeproc -o ${TEX_ROOT}/../siesta.rst
pandoc ${TEX_ROOT}/tex/tbtrans.tex --filter ${TEX_ROOT}/pandoc_filter.py --citeproc -o ${TEX_ROOT}/../tbtrans.rst

echo "POSTPROCESSING RST FILES..."
# Replace any placeholders that we have used to write special characters.
postprocess_rst_specialchars ${TEX_ROOT}/../siesta.rst
postprocess_rst_specialchars ${TEX_ROOT}/../tbtrans.rst

# Postprocess rsts to avoid duplicated fdfparam definitions in the RTD.
postprocess_fdfparam_rst ${TEX_ROOT}/../tbtrans.rst "tbt"
sed_in_place -e 's/tbtsave.hs/save.hs/g' ${TEX_ROOT}/../tbtrans.rst
sed_in_place -e 's/tbtts.hs.save/ts.hs.save/g' ${TEX_ROOT}/../tbtrans.rst
sed_in_place -e 's/tbtcdf.save/cdf.save/g' ${TEX_ROOT}/../tbtrans.rst

# Make sure that there is a top level title in the tbtrans manual
sed_in_place '1s/^/TBTrans User Guide\n******************\n/' ${TEX_ROOT}/../tbtrans.rst
