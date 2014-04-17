#!/bin/sh


function run_lep_eff {
    local SCRIPT="python/compute_lepton_efficiency.py"
    local OUT_BASE="out/fakerate/merged"
    local LOG_BASE="log/fakerate"
    time  ${SCRIPT} -v -t ${TAG} --output-dir ${OUT_BASE}/el_eff_conv_${TAG} -l el -m conv -f 2>&1 ${LOG_BASE}/el_eff_conv_${TAG}.log
    time  ${SCRIPT} -v -t ${TAG} --output-dir ${OUT_BASE}/el_eff_hflf_${TAG} -l el -m hflf -f 2>&1 ${LOG_BASE}/el_eff_hflf_${TAG}.log
    time  ${SCRIPT} -v -t ${TAG} --output-dir ${OUT_BASE}/el_eff_hf_${TAG}   -l el -m hf   -f 2>&1 ${LOG_BASE}/el_eff_hf_${TAG}.log
    time  ${SCRIPT} -v -t ${TAG} --output-dir ${OUT_BASE}/el_eff_lf_${TAG}   -l el -m lf   -f 2>&1 ${LOG_BASE}/el_eff_lf_${TAG}.log
    time  ${SCRIPT} -v -t ${TAG} --output-dir ${OUT_BASE}/el_eff_real_${TAG} -l el -m real -f 2>&1 ${LOG_BASE}/el_eff_real_${TAG}.log

    time  ${SCRIPT} -v -t ${TAG} --output-dir ${OUT_BASE}/mu_eff_hflf_${TAG} -l mu -m hflf -f 2>&1 ${LOG_BASE}/mu_eff_hflf_${TAG}.log
    time  ${SCRIPT} -v -t ${TAG} --output-dir ${OUT_BASE}/mu_eff_hf_${TAG}   -l mu -m hf   -f 2>&1 ${LOG_BASE}/mu_eff_hf_${TAG}.log
    time  ${SCRIPT} -v -t ${TAG} --output-dir ${OUT_BASE}/mu_eff_lf_${TAG}   -l mu -m lf   -f 2>&1 ${LOG_BASE}/mu_eff_lf_${TAG}.log
    time  ${SCRIPT} -v -t ${TAG} --output-dir ${OUT_BASE}/mu_eff_real_${TAG} -l mu -m real -f 2>&1 ${LOG_BASE}/mu_eff_real_${TAG}.log
}

run_lep_eff
