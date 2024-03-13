# Once you are happy with the program you may drop this step and move the compiled executable to another directory.
./compiler.sh


# This moves the old simulation files in a subfolder such that Mathematica notebook ReadOut.nb is simpler to handle.
mv *.tsv oldfiles/
mv *.log oldfiles/

./dislocation-solver --r0 0.5 --K0 25 --lambda0 8.1 --lambdax 10.0 --g 4.0  --passes 10000 --printConfiguration True &

 wait

# You may add several of the following lines (ending with '&')
# to run simulations in parallel. However do not add more lines
# than nodes on your computing machine.
# If you want to have many simulations following one another
# you may use the format separated by wait
#
# ./dislocation-solver ... &
#            .
#            .
# ./dislocation-solver ... &
#
# wait
#
# ./dislocation-solver ... &
#            .
#            .
# ./dislocation-solver ... &
