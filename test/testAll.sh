#!/bin/bash
#test script for metadamage, 
#


PRG=../metaDMG-cpp
BAM1=./data/f570b1db7c.dedup.filtered.bam

LOG=${0}.log
echo Using Logfile: ${LOG}
rm -f ${LOG} *.bin #remove old logfile and binary tempfile

RVAL=0

echo "Testing Existence of ${PRG}"
if [[ ! -f "${PRG}" ]]; then
    echo "Problem finding program: ${PRG}"
    RVAL=1
fi

echo "Testing Existence of ${BAM1}"
if [[ ! -f "${BAM1}" ]]; then
    echo "Problem finding file: ${PRG}"
    RVAL=$((2+RVAL))
fi

echo "Sorting bamfile"
BAM="$(dirname ${BAM1})/$(basename ${BAM1} .bam).rname.bam"
CMD="samtools sort -n ${BAM1} -o ${BAM}"
${CMD}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((4+RVAL))
fi

mkdir -p output

echo "Running getdamage global"
CMD="${PRG} getdamage --run_mode 0 --min_length 35 --print_length 5 --out_prefix output/test_getdamage_global ${BAM}"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((8+RVAL))
fi

echo "Running getdamage local"
CMD="${PRG} getdamage -r 1 -l 35 -p 5 -o output/test_getdamage_local ${BAM}"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((16+RVAL))
fi

echo "Running lca"
CMD="${PRG} lca --bam ${BAM} --names data/names.dmp.gz --nodes data/nodes.dmp.gz --acc2tax data/acc2taxid.map.gz --sim_score_low 0.9 --sim_score_high 1.0 --edit_dist_min 0 --edit_dist_max 10 --min_mapq 0 --how_many 35 --weight_type 1 --fix_ncbi 0 --out_prefix output/test_lca"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((32+RVAL))
fi

echo "Running aggregate"
CMD="${PRG} aggregate output/test_lca.bdamage.gz --nodes data/nodes.dmp.gz --names data/names.dmp.gz --lcastat output/test_lca.stat.gz --out_prefix output/test_aggregate"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((64+RVAL))
fi

echo "Running dfit local (single threaded)"
CMD="${PRG} dfit output/test_lca.bdamage.gz --names data/names.dmp.gz --nodes data/nodes.dmp.gz --showfits 2 --nopt 2 --nbootstrap 2 --seed 12345 --lib ds --out output/test_dfit_local"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((64+RVAL))
fi
# Remove 'ncall' column and round values, since it fails on GitHub tests
#zcat output/test_dfit_local.dfit.gz | cut -f 1-6,8- | head -n 10 | numfmt -d $'\t' --header --format='%.2f' --field=2- --invalid=ignore > output/test_dfit_local.dfit.fix

echo HEEJEJEJEJEJEJEJEJEJ
zcat output/test_dfit_local.dfit.gz |head  |wc -c
echo HEEJEJEJEJEJEJEJEJEJ
zcat output/test_dfit_local.dfit.gz |tail  |wc -c
exit 1
echo "Running dfit local (10 threaded)"
CMD="${PRG} dfit output/test_lca.bdamage.gz --threads 10 --names data/names.dmp.gz --nodes data/nodes.dmp.gz --showfits 2 --nopt 2 --nbootstrap 2 --seed 12345 --lib ds --out output/test_dfit_local_10threads"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((128+RVAL))
fi
# Remove 'ncall' column and round values, since it fails on GitHub tests
zcat output/test_dfit_local_10threads.dfit.gz | cut -f 1-6,8- | head -n 10 | numfmt -d $'\t' --header --format='%.2f' --field=2- --invalid=ignore | sort -r > output/test_dfit_local_10threads.dfit.fix

echo "Running dfit global"
CMD="${PRG} dfit output/test_getdamage_global.bdamage.gz --showfits 2 --seed 12345 --lib ds --out output/test_dfit_global"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((128+RVAL))
fi
# Remove 'ncall' column since it fail on GitHub tests
zcat output/test_dfit_global.dfit.gz | cut -f 1-6,8- > output/test_dfit_global.dfit.fix

echo "Running printoptions"
CMD="${PRG} print output/test_getdamage_local.bdamage.gz"
${CMD} 1>output/test_getdamage_local.bdamage.tsv 2>>${LOG}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((256+RVAL))
fi

CMD="${PRG} print output/test_getdamage_global.bdamage.gz"
${CMD} 1>output/test_getdamage_global.bdamage.tsv 2>>${LOG}
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((256+RVAL))
fi

CMD="${PRG} print_ugly output/test_lca.bdamage.gz --out_prefix output/test_lca"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((512+RVAL))
fi

CMD="${PRG} print_ugly output/test_lca.bdamage.gz --names data/names.dmp.gz --nodes data/nodes.dmp.gz --lcastat output/test_lca.stat.gz --out_prefix output/test_lca_taxa"
${CMD} >> ${LOG} 2>&1
if [[ $? -ne 0 ]]; then
    echo "Problem running command: ${CMD}"
    RVAL=$((1024+RVAL))
fi


echo "Validating checksums"
echo "========================"
gunzip -f output/*.gz
md5sum -c output.md5 
if [[ $? -ne 0 ]]; then
    echo "Problem with md5sums"
    RVAL=$((2048+${RVAL}))
fi
echo "=====RVAL:${RVAL}======="

if [[ ${RVAL} -ne 0 ]];then
    echo "====StartOfLog==="
    cat ${LOG}
    echo "====EndOfLog==="
    exit 1 #exit codes are capped at 255
fi

exit 0
