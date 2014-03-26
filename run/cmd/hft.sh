#!/usr/bin/env bash

# Steps to produce and verify HFT trees.
# Run without any argument to get help.
#
# davide.gerbaudo@gmail.com
# 2014-03-21

function help {
	echo "These are the steps to produce the trees for HistFitter"
	echo "- Setup a tag : 'export TAG=Feb_21'"
	echo "- Submit the jobs to generate the trees : 'hft.sh submit'"
	echo "- Create the tar file for HistFitter :    'hft.sh maketar'"
	echo "- Fill the histograms :                   'hft.sh fillhistos'"
	echo "- Produce plots and summary tables :      'hft.sh makeplots'"
    echo "- sanity checks:                          'hft.sh test'"
}

FAKE_SYSTEMATICS="NOM"
FAKE_SYSTEMATICS+=" EL_FR_DOWN EL_FR_UP"
FAKE_SYSTEMATICS+=" EL_RE_DOWN EL_RE_UP"
FAKE_SYSTEMATICS+=" MU_FR_DOWN MU_FR_UP"
FAKE_SYSTEMATICS+=" MU_RE_DOWN MU_RE_UP"

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
#FORREAL=false
#ifeq (${FORREAL},true)
#	SUBMIT_OPTION = --submit
#	SUBMIT_NOTIFY =  ./python/email_me_when_youre_done.py --message "hft trees with tag $TAG done; now you can go to `pwd` and type 'make tar; make histos'" &
#	SUBMIT_MSG = Jobs submitted, now wait for the email
#else
#	SUBMIT_OPTION =
#	SUBMIT_NOTIFY =
#	SUBMIT_MSG = "This was a dry run; use 'make submit FORREAL=true' if you actually want to submit the jobs"
#endif
#
#submit:
#	@echo ${SUBMIT_OPTION}
#	./python/submitJobs.py --susyplot -o  -t ${TAG} -e 'period' --other-opt "--with-hft --with-syst" $(SUBMIT_OPTION)
#	./python/submitJobs.py --susyplot -o  -t ${TAG} -s 'period' --other-opt "--with-hft"             $(SUBMIT_OPTION)
#	./python/submitJobs.py --fakepred -o  -t ${TAG} -s 'period' --other-opt "--with-hft"             $(SUBMIT_OPTION)
#	${SUBMIT_NOTIFY}
#	@echo $(SUBMIT_MSG)
#
## the list of files that will go in the tar
#HFT_FILE_LIST="hft_files_${TAG}.txt"

function cleantarlist {
    local HFT_FILE_LIST=$1
	rm    ${HFT_FILE_LIST}
	touch ${HFT_FILE_LIST}
}

function mergefake {
    printf "Merge fake"
    local HFT_FILE_LIST=$1
	for FAKE_SYS in ${FAKE_SYSTEMATICS}
      do
      rm   out/fakepred/merged/${FAKE_SYS}_${TAG}.root
      hadd out/fakepred/merged/${FAKE_SYS}_${TAG}.root out/fakepred/${FAKE_SYS}_*.root > /tmp/mergefake_${FAKE_SYS}_${TAG}.log
      echo "out/fakepred/merged/${FAKE_SYS}_${TAG}.root" >> ${HFT_FILE_LIST}
    done
    printf "...now %d files\n" $(cat ${HFT_FILE_LIST} | wc -l)
}

function mergedata {
    printf "Merge data"
    local HFT_FILE_LIST=$1
	rm    out/susyplot/merged/NOM_period_PhysCont_${TAG}.root
	hadd  out/susyplot/merged/NOM_period_PhysCont_${TAG}.root out/susyplot/NOM_period*.root > /tmp/mergedata_${TAG}.log
	echo "out/susyplot/merged/NOM_period_PhysCont_${TAG}.root" >> ${HFT_FILE_LIST}
    printf "...now %d files\n" $(cat ${HFT_FILE_LIST} | wc -l)
}

function collectmc {
    printf "Collect MC"
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

#fillhistos:
#	./python/check_hft_trees.py --verbose --input-gen  out/susyplot/merged/ --input-fake out/fakepred/merged/ --output-dir out/hft/   --batch
#	./python/email_me_when_youre_done.py --message="histograms filled; now you can go to `pwd` and type 'make makeplots'" &
#
#makeplots:
#	./python/check_hft_trees.py -v --input-dir out/hft/ --output-dir out/hft/plots_all_syst  | tee plot_all_syst.log

#-------------------------
# main
#-------------------------

if [ -z "$1" ]
then
    help
elif [ "$1" == "test" ]
then
        testvars
elif [ "$1" == "submit" ]
then
    echo
elif [ "$1" == "maketar" ]
then
    HFT_FILE_LIST="hft_files_${TAG}.txt"
    TAR_FILE="~/tmp/hft_$${TAG}.tgz"
    echo "Files that will go in the tar: ${HFT_FILE_LIST}"
    cleantarlist ${HFT_FILE_LIST}
    mergefake    ${HFT_FILE_LIST}
    mergedata    ${HFT_FILE_LIST}
    collectmc    ${HFT_FILE_LIST}
    maketar      ${HFT_FILE_LIST} ${TAR_FILE}
elif [ "$1" == "fillhistos" ]
then
    echo
elif [ "$1" == "makeplots" ]
then
    echo
else
    echo "Unknown option '$*'"
fi