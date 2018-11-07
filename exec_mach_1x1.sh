#!/bin/bash

processes=(1 2 16 64)
nodes=(1 2 8)
exec_modes=(0 1)
matrix_sizes=(2050 4098)
max_execs=5

for exec_mode in "${exec_modes[@]}"
do
    echo "-------------------------------"
    echo "|  Execution mode $exec_mode  |"
    echo "-------------------------------"

    for num_nodes in "${nodes[@]}"
    do
        echo "-------------------------------"
        echo "| Number of nodes $num_nodes  |"
        echo "-------------------------------"

        for num_processes in "${processes[@]}"
        do
            echo "--------------------------------------"
            echo "| Number of processes $num_processes |"
            echo "--------------------------------------"

            for size in "${matrix_sizes[@]}"
            do
                echo "-------------------------------"
                echo "|       Matrix size $size     |"
                echo "-------------------------------"

                num_exec=0
                while [[ $num_exec -lt $max_execs ]]
                do
                    echo "-------------------------------"
                    echo "|  Execution number $num_exec |"
                    echo "-------------------------------"

                    ppn=$num_processes/$num_nodes
                    mpiexec -f machines_1x1.txt -np $num_nodes -ppn $ppn ./gs_mpi $size $exec_mode
                    num_exec=$(( num_exec+1 ))
                done
            done
        done
    done
done
