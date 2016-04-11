diff expected.txt <(for i in $(seq 0 100)
do
    ../sr --random -n 100 --np $i -a --seed $i --timeout 1000000 --maxiter 100 --gen-type 9
done)
