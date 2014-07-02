#!/usr/bin/env bash

# Commands to produce and verify HFT trees.
# Run without any argument to get help.
#
# davide.gerbaudo@gmail.com
# 2014-03-21

SCRIPT_NAME=$0

function help {
	echo -e "These are the steps to produce the trees for HistFitter"
	echo -e "export TAG=Feb_21 \t # Setup a tag"
	echo -e "${SCRIPT_NAME} filltrees  \t # Submit the jobs to generate the trees"
	echo -e "${SCRIPT_NAME} maketar    \t # Create the tar file for HistFitter"
	echo -e "${SCRIPT_NAME} fillhistos \t # Fill the histograms"
	echo -e "${SCRIPT_NAME} makeplots  \t # Produce plots and summary tables"
    echo -e "${SCRIPT_NAME} test       \t # sanity checks"
}

FAKE_SYSTEMATICS="NOM"
FAKE_SYSTEMATICS+=" EL_FR_DOWN EL_FR_UP"
FAKE_SYSTEMATICS+=" EL_RE_DOWN EL_RE_UP"
FAKE_SYSTEMATICS+=" MU_FR_DOWN MU_FR_UP"
FAKE_SYSTEMATICS+=" MU_RE_DOWN MU_RE_UP"
FAKE_SYSTEMATICS+=" EL_FRAC_DO EL_FRAC_UP"
FAKE_SYSTEMATICS+=" MU_FRAC_DO MU_FRAC_UP"

MC_SYSTEMATICS="NOM"
MC_SYSTEMATICS+=" EER_DN EER_UP"
MC_SYSTEMATICS+=" EES_LOW_DN EES_LOW_UP"
MC_SYSTEMATICS+=" EES_MAT_DN EES_MAT_UP"
MC_SYSTEMATICS+=" EES_PS_DN EES_PS_UP"
MC_SYSTEMATICS+=" EES_Z_DN EES_Z_UP"
MC_SYSTEMATICS+=" ID_DN ID_UP"
MC_SYSTEMATICS+=" JER"
MC_SYSTEMATICS+=" JES_DN JES_UP"
MC_SYSTEMATICS+=" MS_DN MS_UP"
MC_SYSTEMATICS+=" RESOST"
MC_SYSTEMATICS+=" SCALEST_DN SCALEST_UP"

function testvars {
    echo "Printing all systematic variations"
    echo "Fake variations:"
	for FAKE_SYS in ${FAKE_SYSTEMATICS} ; do echo "    ${FAKE_SYS}" ; done
    echo "MC variations:"
	for MC_SYS in ${MC_SYSTEMATICS} ; do echo "   ${MC_SYS}" ; done
}

function filltrees {
    local SUBMIT_OPTION="" && [[ $#>0 ]] && SUBMIT_OPTION="${1}"
# 	./python/submitJobs.py --susyplot -o  -t ${TAG} -s 'Herwigpp_UEEE3_CTEQ6L1_DGnoSL_TB10_05' --other-opt "--with-hft --with-syst" ${SUBMIT_OPTION}
	./python/submitJobs.py --susyplot -o  -t ${TAG} -e 'period' --other-opt "--with-hft --with-syst" ${SUBMIT_OPTION}
	./python/submitJobs.py --susyplot -o  -t ${TAG} -s 'period' --other-opt "--with-hft"             ${SUBMIT_OPTION}
	./python/submitJobs.py --fakepred -o  -t ${TAG} -s 'period' --other-opt "--with-hft"             ${SUBMIT_OPTION}
    if [[ "$SUBMIT_OPTION" == "--submit" ]]
        then
        ./python/email_me_when_youre_done.py \
        --message="hft trees with tag ${TAG} done; now you can go to `pwd` and type '${SCRIPT_NAME} maketar; ${SCRIPT_NAME} fillhistos'" &
        echo "Jobs submitted, now wait for the email"
    else
        echo "This was a dry run; use '${SCRIPT_NAME} filltrees --submit' to actually submit the jobs"
    fi
}

function cleantarlist {
    local HFT_FILE_LIST=$1
	rm    ${HFT_FILE_LIST}
	touch ${HFT_FILE_LIST}
}

function collectfake {
    printf "Collecting fake"
    local HFT_FILE_LIST=$1
	for FAKE_SYS in ${FAKE_SYSTEMATICS}
      do
      ls out/fakepred/${FAKE_SYS}_*.root >> ${HFT_FILE_LIST}
    done
    printf "...now %d files\n" $(cat ${HFT_FILE_LIST} | wc -l)
}

function collectdata {
    printf "Collecting data"
    local HFT_FILE_LIST=$1
    ls out/susyplot/NOM_period*.root >> ${HFT_FILE_LIST}
    printf "...now %d files\n" $(cat ${HFT_FILE_LIST} | wc -l)
}

function collectmc {
    printf "Collecting MC"
    local HFT_FILE_LIST=$1
    for MC_SYS in ${MC_SYSTEMATICS}
      do
      find out/susyplot/ -maxdepth 1 -name "${MC_SYS}_*.root" | sort >> ${HFT_FILE_LIST}
    done
    printf "...now %d files\n" $(cat ${HFT_FILE_LIST} | wc -l)
}

function maketar {
    local HFT_FILE_LIST=$1
    local TAR_FILE=$2
	echo "creating tar file ~/tmp/hft_${TAG}.tgz"
	tar czf ~/tmp/hft_${TAG}.tgz --files-from ${HFT_FILE_LIST}
}

function fillhistos {
	./python/check_hft_trees.py \
 --verbose \
 --input-gen  out/susyplot/merged/ \
 --input-fake out/fakepred/merged/ \
 --output-dir out/hft/ \
 --batch
	./python/email_me_when_youre_done.py \
 --message="histograms filled (tag ${TAG}); now you can go to `pwd` and type 'hft.sh makeplots'" &
}

function makeplots {
	./python/check_hft_trees.py -v --input-dir out/hft/ --output-dir out/hft/plots_all_syst  | tee plot_all_syst_${TAG}.log
}

#-------------------------
# main
#-------------------------

if [ -z "$1" ]
then
    help
elif [ "$1" == "test" ]
then
        testvars
elif [ "$1" == "filltrees" ]
then
    filltrees "${2}"
elif [ "$1" == "maketar" ]
then
    HFT_FILE_LIST="hft_files_${TAG}.txt"
    TAR_FILE="~/tmp/hft_$${TAG}.tgz"
    echo "Files that will go in the tar: ${HFT_FILE_LIST}"
    cleantarlist ${HFT_FILE_LIST}
    collectfake  ${HFT_FILE_LIST}
    collectdata  ${HFT_FILE_LIST}
    collectmc    ${HFT_FILE_LIST}
    maketar      ${HFT_FILE_LIST} ${TAR_FILE}
elif [ "$1" == "fillhistos" ]
then
    echo "Filling histograms"
    fillhistos
elif [ "$1" == "makeplots" ]
then
    echo "Making plots"
    makeplots
else
    echo "Unknown option '$*'"
fi