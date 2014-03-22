#!/usr/bin/make -f

# Steps to produce and verify HFT trees.
# For instructions, type
#   hft.txt help
#
# Note to self:
# remember the differences between make variables and shell variables
# ($ vs $$, see
# http://stackoverflow.com/questions/19475037/function-and-difference-between-and-in-makefile)
#
# davide.gerbaudo@gmail.com
# 2014-03-21

help:
	@echo These are the steps to produce the trees for HistFitter
	@echo - Setup a tag : "export TAG=Feb_21"
	@echo - Submit the jobs to generate the trees : "make submit"
	@echo - Create the tar file for HistFitter : "make tar"
	@echo - Fill the histograms : "make fillhistos"
	@echo - Produce plots and summary tables : "make makeplots"

FAKE_SYSTEMATICS=NOM
FAKE_SYSTEMATICS+=EL_FR_DOWN EL_FR_UP
FAKE_SYSTEMATICS+=EL_RE_DOWN EL_RE_UP
FAKE_SYSTEMATICS+=MU_FR_DOWN MU_FR_UP
FAKE_SYSTEMATICS+=MU_RE_DOWN MU_RE_UP

MC_SYSTEMATICS=NOM
MC_SYSTEMATICS+=EER_DN EER_UP
MC_SYSTEMATICS+=EES_LOW_DN EES_LOW_UP
MC_SYSTEMATICS+=EES_MAT_DN EES_MAT_UP
MC_SYSTEMATICS+=EES_PS_DN EES_PS_UP
MC_SYSTEMATICS+=EES_Z_DN EES_Z_UP
MC_SYSTEMATICS+=ID_DN ID_UP
MC_SYSTEMATICS+=JER
MC_SYSTEMATICS+=JES_DN JES_UP
MC_SYSTEMATICS+=MS_DN MS_UP
MC_SYSTEMATICS+=RESOST
MC_SYSTEMATICS+=SCALEST_DN SCALEST_UP

testvars:
	for FAKE_SYS in ${FAKE_SYSTEMATICS} ; do echo $${FAKE_SYS} ; done
	for MC_SYS in ${MC_SYSTEMATICS} ; do echo $${MC_SYS} ; done

FORREAL=false
ifeq (${FORREAL},true)
	SUBMIT_OPTION = --submit
	SUBMIT_NOTIFY =  ./python/email_me_when_youre_done.py --message "hft trees with tag $TAG done; now you can go to `pwd` and type 'make tar; make histos'" &
	SUBMIT_MSG = Jobs submitted, now wait for the email
else
	SUBMIT_OPTION =
	SUBMIT_NOTIFY =
	SUBMIT_MSG = "This was a dry run; use 'make submit FORREAL=true' if you actually want to submit the jobs"
endif

submit:
	@echo ${SUBMIT_OPTION}
	./python/submitJobs.py --susyplot -o  -t ${TAG} -e 'period' --other-opt "--with-hft --with-syst" $(SUBMIT_OPTION)
	./python/submitJobs.py --susyplot -o  -t ${TAG} -s 'period' --other-opt "--with-hft"             $(SUBMIT_OPTION)
	./python/submitJobs.py --fakepred -o  -t ${TAG} -s 'period' --other-opt "--with-hft"             $(SUBMIT_OPTION)
	${SUBMIT_NOTIFY}
	@echo $(SUBMIT_MSG)

# the list of files that will go in the tar
HFT_FILE_LIST="hft_files_${TAG}.txt"

cleantarlist:
	rm    ${HFT_FILE_LIST}
	touch ${HFT_FILE_LIST}

mergefake:
	for FAKE_SYS in ${FAKE_SYSTEMATICS} ; \
 do \
 rm   out/fakepred/merged/$${FAKE_SYS}_${TAG}.root ; \
 hadd out/fakepred/merged/$${FAKE_SYS}_${TAG}.root out/fakepred/$${FAKE_SYS}_*.root ; \
 echo "out/fakepred/merged/$${FAKE_SYS}_${TAG}.root" >> ${HFT_FILE_LIST} ; \
 done

mergedata:
	rm    out/susyplot/merged/NOM_period_PhysCont_$${TAG}.root
	hadd  out/susyplot/merged/NOM_period_PhysCont_$${TAG}.root out/susyplot/NOM_period*.root
	echo "out/susyplot/merged/NOM_period_PhysCont_$${TAG}.root" >> ${HFT_FILE_LIST}

collectmc:
	for MC_SYS in ${MC_SYSTEMATICS} ; \
 do \
 echo "$${MC_SYS}_*.root" ; \
 find out/susyplot/ -maxdepth 1 -name "$${MC_SYS}_*.root" >> ${HFT_FILE_LIST} ; \
 done
	echo Mc files added to $(HFT_FILE_LIST)

tar: cleantarlist mergefake mergedata collectmc
	@echo "creating tar file ~/tmp/hft_$${TAG}.tgz"
	tar czf ~/tmp/hft_$${TAG}.tgz --files-from ${HFT_FILE_LIST}

fillhistos:
	./python/check_hft_trees.py --verbose --input-gen  out/susyplot/merged/ --input-fake out/fakepred/merged/ --output-dir out/hft/   --batch
	./python/email_me_when_youre_done.py --message="histograms filled; now you can go to `pwd` and type 'make makeplots'" &

makeplots:
	./python/check_hft_trees.py -v --input-dir out/hft/ --output-dir out/hft/plots_all_syst  | tee plot_all_syst.log
