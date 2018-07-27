#!/bin/bash

prefix=`sed "/^prefix = /!d;s/.* = //" Makefile`

LD_LIBRARY_PATH="$prefix/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH

valgrind=""
insure=""
for o in "$@"; do
    case "$o" in
	--valgrind)
	    valgrind=1
	    ;;
	--insure)
	    insure=1
	    ;;
	--run-5-p-p-mb-towerstat)
	    macro="Run_towerstat.C"
	    data_set="pp200"
	    ;;
	--run-5-p-p-ert-towerstat)
	    macro="Run_towerstat.C"
	    data_set="pp200ERT"
	    ;;
	--run-5-cu-cu-mb-towerstat)
	    macro="Run_towerstat.C"
	    data_set="CuCu200"
	    ;;
	--run-5-cu-cu-ert-towerstat)
	    macro="Run_towerstat.C"
	    data_set="CuCu200ERT"
	    ;;
	--run-5-p-p-mb-rec)
	    macro="Run_jetextract_run_5_p_p_mb_rec.C"
	    data_set="pp200"
	    ;;
	--run-5-p-p-mb-rec-no-recal)
	    macro="Run_jetextract_run_5_p_p_mb_rec_norecal.C"
	    data_set="pp200"
	    ;;
	--run-5-p-p-ert-rec)
	    macro="Run_jetextract_run_5_p_p_ert_rec.C"
	    data_set="pp200ERT"
	    ;;
	--run-5-p-p-ert-rec-no-recal)
	    macro="Run_jetextract_run_5_p_p_ert_rec_norecal.C"
	    data_set="pp200ERT"
	    ;;
	--run-5-cu-cu-mb-rec)
	    macro="Run_jetextract_run_5_cu_cu_mb_rec.C"
	    data_set="CuCu200"
	    ;;
	--run-5-cu-cu-ert-rec)
	    macro="Run_jetextract_run_5_cu_cu_ert_rec.C"
	    data_set="CuCu200ERT"
	    ;;
	--run-5-cu-cu-mb-fixed)
	    macro="Run_jetextract_run_5_cu_cu_mb_fixed.C"
	    data_set="CuCu200"
	    ;;
    esac
done

if [[ -z "$data_set" ]]; then
    echo "`basename \"$0\"`: error: no data set specified" 1>&2
fi

if [[ -n "$valgrind" ]]; then
	root="valgrind --num-callers=20 --tool=memcheck --leak-check=full --error-limit=no --log-file=valgrind.log --suppressions=$ROOTSYS/root.supp --leak-resolution=high root.exe"
	nevent=10
elif [[ -n "$insure" ]]; then
	root="root_insure.exe"
	nevent=10
else
	root="root.exe"
	nevent=1000
fi
cd `dirname "$0"`/../../pat/macro
exec $root -b 'RunMyMacro.C("'"$macro"'", "test_run.root", '"$nevent"', "'"$data_set"'")'
