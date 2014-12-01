#!/usr/bin/env bash

# Commands to produce the matrix file for the fake estimate prediction
# Run without any argument to get help.
# Required inputs: ntuples produced by MeasureFakeRate2.
#
# davide.gerbaudo@gmail.com
# 2014-08-08

readonly SCRIPT_NAME=$(basename $0)
# see http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in
readonly PROGDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
readonly SIGNAL_REGION="razor0j" #"emu" #ssinc1j" # depends on what you are working; probably 'emu', 'ssinc1j', ...
# the two options below are needed only when scanning the tight defs for Liz
readonly TIGHTDEF_OPT="" #"--tight-def fakeu.lepIsTight_07"
readonly TIGHTDEF_SUFFIX="" # "_tight_07"

readonly TAG_NTUPLE="Oct_20"
readonly TAG_COMPOSITIONS="razor_mdr00_Oct_20"
readonly TAG_EFFICIENCIES="Jul_26"
readonly TAG_SCALE_FACTORS="Jul_26" # still used for el conv, others are flat for these tests
readonly TAG_WEIGHTED_AVG="razor_mdr00_Nov_26"

function help {
	echo -e "These are the steps to produce the fake matrix file"
    echo -e "export TAG=Jul_26 \t # Setup a tag, used for input and outputs"
    echo -e "${SCRIPT_NAME} efficiencies \t # Compute efficiencies"
    echo -e "${SCRIPT_NAME} compositions \t # Compute compositions"
    echo -e "${SCRIPT_NAME} scalefactors \t # Compute data/mc scale factors (not implemented yet...)"
    echo -e "${SCRIPT_NAME} build        \t # Build weighted matrix"
    echo -e "Note: current tag is ${TAG}"
    echo -e "All these scripts need the input ntuples, which can be generated with"
    echo -e "submitJobs.py --fakerate -o --tag=${TAG} -e '(Sherpa_CT10_lllnu_WZ|noslep.*_WH)' --other-opt=\" --write-tuple \" --submit"
    echo -e "submitJobs.py --fakerate -o --tag=${TAG} -s 'PowhegPythia8_AU2CT10_WZ_W' --alsoplaceholders --other-opt=\" --write-tuple \" --submit"
}

function check_env {
    : ${TAG:?"Need to set a valid TAG."}
}

function compute_efficiencies {
    echo -e "computing efficiencies...this might take some time if we need to fill the histograms"
    local FILL="--fill-histos" # force fill
    local IN_DIR="out/fakerate/"
    local OUT_DIR="out/fake/efficiencies_${TAG_EFFICIENCIES}"
    local COMMON_OPT="--verbose --tag ${TAG_NTUPLE} --input-dir ${IN_DIR} --output-dir ${OUT_DIR}"
    COMMON_OPT="${COMMON_OPT} ${FILL} ${TIGHTDEF_OPT}"
    mkdir -p ${OUT_DIR}
    echo -e "\n--- electron --- `date +%F-%T`\n"
    time python/compute_eff_from_ntuple.py ${COMMON_OPT} --lepton el --mode real  2>&1  | tee ${OUT_DIR}/el_real.txt
    time python/compute_eff_from_ntuple.py ${COMMON_OPT} --lepton el --mode conv  2>&1  | tee ${OUT_DIR}/el_conv.txt
    time python/compute_eff_from_ntuple.py ${COMMON_OPT} --lepton el --mode hflf  2>&1  | tee ${OUT_DIR}/el_hflf.txt
    echo -e "\n--- muon --- `date +%F-%T`\n"
    time python/compute_eff_from_ntuple.py ${COMMON_OPT} --lepton mu --mode real  2>&1  | tee ${OUT_DIR}/mu_real.txt
    time python/compute_eff_from_ntuple.py ${COMMON_OPT} --lepton mu --mode hflf  2>&1  | tee ${OUT_DIR}/mu_hflf.txt
    echo -e "Done: log files in `ls ${OUT_DIR}/*.txt`"
}

function compute_compositions {
    local FILL="--fill-histos" # force fill
    local IN_DIR="out/fakerate/"
    local OUT_DIR="out/fake/compositions_${TAG_COMPOSITIONS}"
    local REGION="${SIGNAL_REGION}"
    local OPT="" # "--syst-fudge" # this and the one below if you are doing fake composition variation
    #local OUT_DIR=${OUT_DIR}_syst_shift
    mkdir -p ${OUT_DIR}
    local COMMON_OPT="${FILL} ${OPT} -v --tag ${TAG_NTUPLE} --region ${REGION} --input-dir ${IN_DIR} --output-dir ${OUT_DIR}"
    echo -e "\n--- electron --- `date +%F-%T`\n"
    time python/compute_fake_compositions.py ${COMMON_OPT} --lepton el 2>&1 | tee ${OUT_DIR}/el.txt
    echo -e "\n--- muon --- `date +%F-%T`\n"
    time python/compute_fake_compositions.py ${COMMON_OPT} --lepton mu 2>&1 | tee ${OUT_DIR}/mu.txt
    echo -e "Done: log files in `ls ${OUT_DIR}/*.txt`"
}

function compute_scalefactors {
    local FILL="--fill-histos" # force fill
    local IN_DIR="out/fakerate/"
    local OUT_DIR="out/fake/scalefactors_${TAG_SCALE_FACTORS}"
    local COMMON_OPT="${FILL} ${OPT} -v --tag ${TAG_NTUPLE} --input-dir ${IN_DIR} --output-dir ${OUT_DIR}"
    COMMON_OPT="${COMMON_OPT} ${TIGHTDEF_OPT}"
    mkdir -p ${OUT_DIR}
    # todo: implement real lepton case
    echo -e "\n--- electron --- `date +%F-%T`\n"
    time python/compute_fake_scale_factor.py ${COMMON_OPT} --lepton el --region hflf 2>&1 | tee ${OUT_DIR}/el_hflf.txt
    time python/compute_fake_scale_factor.py ${COMMON_OPT} --lepton el --region conv 2>&1 | tee ${OUT_DIR}/el_conv.txt
    #time python/compute_fake_scale_factor.py ${COMMON_OPT} --lepton el --region real 2>&1 | tee ${OUT_DIR}/el_real.txt
    echo -e "\n--- muon --- `date +%F-%T`\n"
    time python/compute_fake_scale_factor.py ${COMMON_OPT} --lepton mu --region hflf 2>&1 | tee ${OUT_DIR}/mu_hflf.txt
    #time python/compute_fake_scale_factor.py ${COMMON_OPT} --lepton mu --region real 2>&1 | tee ${OUT_DIR}/mu_real.txt
}

function build_matrix {
    local SYS="" #"_syst_shift" # uncomment this to pick up the shifted compositions
    local IN_COMPOSITION_DIR="out/fake/compositions_${TAG_COMPOSITIONS}/"
    local IN_EFFICICENCY_DIR="out/fake/efficiencies_${TAG_EFFICIENCIES}"
    local IN_SCALEFACTOR_DIR="out/fake/scalefactors_${TAG_SCALE_FACTORS}"
    local OUT_DIR="out/fake/weigtedmatrix_${TAG_WEIGHTED_AVG}"
    local REGION="${SIGNAL_REGION}"
    local COMMON_OPT="${OPT} -v --output-dir ${OUT_DIR} --region ${REGION}"

    mkdir -p ${OUT_DIR}
    python/build_fake_matrices.py \
        ${COMMON_OPT} \
        --lepton el \
        --comp-histos   ${IN_COMPOSITION_DIR}/el/el_composition_histos.root \
        --eff-histos    ${IN_EFFICICENCY_DIR}/mcconv_histos_el_eff.root \
        --eff-histos    ${IN_EFFICICENCY_DIR}/mcqcd_histos_el_eff.root \
        --scale-factors ${IN_SCALEFACTOR_DIR}/conv/el/conv_el_scale_histos.root \
        2>&1 | tee ${OUT_DIR}/build_fake_matrices_el.txt

        # --scale-factors ${IN_SCALEFACTOR_DIR}/hflf/el/hflf_el_scale_histos.root \


    python/build_fake_matrices.py \
        ${COMMON_OPT} \
        --lepton mu \
        --comp-histos   ${IN_COMPOSITION_DIR}/mu/mu_composition_histos.root \
        --eff-histos    ${IN_EFFICICENCY_DIR}/mcqcd_histos_mu_eff.root \
        2>&1 | tee ${OUT_DIR}/build_fake_matrices_mu.txt

        # --scale-factors ${IN_SCALEFACTOR_DIR}/hflf/mu/hflf_mu_scale_histos.root \

    echo -e "Putting all the histograms in the same output file"
    echo -e "TODO"
    # python/pick_objects_
}

#-------------------------
# main
#-------------------------

check_env
if [ -z "$1" ]
then
    help
elif [ "$1" == "efficiencies" ]
then
    echo "Computing efficiencies"
    compute_efficiencies
elif [ "$1" == "compositions" ]
then
    echo "Computing compositions"
    compute_compositions
elif [ "$1" == "scalefactors" ]
then
    echo "Computing data/mc scale factors"
    # echo "not implemented yet...using fixed values"
    compute_scalefactors
elif [ "$1" == "build" ]
then
    echo "Building weighted matrix"
    build_matrix
else
    echo "Unknown option '$*'"
fi

