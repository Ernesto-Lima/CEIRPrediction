#!/bin/bash
rm -rf *# *~ outputData &> /dev/null
if [ -d "./outputQueso/" ]; then
    echo "Directory already exists."
else
    mkdir outputQueso
fi

input_n=$1
model_n=$2
position=$3
echo "Model........${model_n}"
echo "Input data...${input_n}"

#======================================================================
echo "model_n = ${model_n}" > model.in
if [[ ${input_n} -eq '123' ]]; then
    echo "n_files = 3" >> model.in
    echo "input_n = 1" >> model.in
    echo "read_file = 0" >> model.in
    echo "all_tv = '../data/data_mean_std_control_tv.dat ../data/data_mean_std_hypio_tv.dat ../data/data_mean_std_norio_tv.dat'" >> model.in
    echo "all_hf = '../data/data_mean_std_control_hf.dat ../data/data_mean_std_hypio_hf.dat ../data/data_mean_std_norio_hf.dat'" >> model.in
elif [[ ${input_n} -eq '456' ]]; then
    echo "n_files = 3" >> model.in
    echo "all_tv = '../data/data_mean_std_e_tv.dat ../data/data_mean_std_hypioe_tv.dat ../data/data_mean_std_norioe_tv.dat'" >> model.in
    echo "all_hf = '../data/data_mean_std_e_hf.dat ../data/data_mean_std_hypioe_hf.dat ../data/data_mean_std_norioe_hf.dat'" >> model.in
    echo "input_n = 4" >> model.in
    echo "read_file = 1" >> model.in
    echo "pars_file = './outputQueso/parameters_d123_m6.dat'" >> model.in
elif [[ ${input_n} -eq '45' ]]; then
    echo "n_files = 2" >> model.in
    echo "all_tv = '../data/data_mean_std_e_tv.dat ../data/data_mean_std_hypioe_tv.dat'" >> model.in
    echo "all_hf = '../data/data_mean_std_e_hf.dat ../data/data_mean_std_hypioe_hf.dat'" >> model.in
    echo "input_n = 4" >> model.in
    echo "read_file = 1" >> model.in
    echo "pars_file = './outputQueso/parameters_d123_m6.dat'" >> model.in
elif [[ ${input_n} -eq '46' ]]; then
    echo "n_files = 2" >> model.in
    echo "all_tv = '../data/data_mean_std_e_tv.dat ../data/data_mean_std_norioe_tv.dat'" >> model.in
    echo "all_hf = '../data/data_mean_std_e_hf.dat ../data/data_mean_std_norioe_hf.dat'" >> model.in
    echo "input_n = 4" >> model.in
    echo "read_file = 1" >> model.in
    echo "pars_file = './outputQueso/parameters_d123_m6.dat'" >> model.in
elif [[ ${input_n} -eq '6' ]]; then
    echo "n_files = 1" >> model.in
    echo "all_tv = '../data/data_mean_std_norioe_tv.dat'" >> model.in
    echo "all_hf = '../data/data_mean_std_norioe_hf.dat'" >> model.in
    echo "input_n = 6" >> model.in
    echo "read_file = 2" >> model.in
    echo "position = ${position}" >> model.in
    echo "pars_file = './outputQueso/parameters_d123_m6.dat ./outputQueso/parameters_d45_m10.dat'" >> model.in
else
    echo "needs to define"
fi

make run
mv ./outputData/rawChain_ml.m ./outputQueso/rawChain_d${input_n}_m${model_n}.m
mv ./outputData/display_sub0.txt ./outputQueso/display_d${input_n}_m${model_n}.txt

#rm -rf *# *~ outputData *.o bc_model job_output.txt model.in &> /dev/null
cd ./outputQueso/
    for chain in raw*.m; do
        suffix=$(echo ${chain} | cut -d_ -f2- | cut -d. -f1)
        echo ${suffix}
        lines=$(grep zeros ${chain} | head -n1 | cut -d"(" -f2 | cut -d, -f1 | awk '{printf "%d\n",$1+1}')
        sed -n -e "2,${lines}p" ${chain} | cut -d"[" -f2 > lastchain_${suffix}.dat
        let iline=lines+3
        let eline=lines+lines+1
        sed -n -e "${iline},${eline}p" ${chain} | cut -d"[" -f2 > loglike_${suffix}.dat
        grep evidence display_${suffix}.txt | tail -n1 | cut -d= -f2 | cut -d, -f1 | awk '{printf "%f\n",$1+0.0}' > evidence_${suffix}.txt
        grep evidence display_${suffix}.txt | tail -n1 | cut -d= -f2 | cut -d, -f1 | awk '{printf "%f\n",$1+0.0}'
        grep factor display_${suffix}.txt | tail -n1 | cut -d, -f3
    done
cd ..
