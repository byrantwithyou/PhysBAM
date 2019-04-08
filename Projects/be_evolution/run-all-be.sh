for dt in 0.01 0.005 0.002 0.001; do
#    echo ./be_evolution 720 -o out-720-$dt -last_frame 2400 -no_self -rd_stiffness 1e6 -kry_it 20 -kry_tol 1e-4 -no_collisions_in_solve -newton_it 100 -dt $dt -resolution 5 -grad_ls -stiffen 10 -dampen 1 -rr -rd -dd -rd_ccd -rr_ccd
#    echo ./be_evolution 730 -o out-730-$dt -no_self -rd_stiffness 1000 -rd_friction 0.3 -kry_it 20 -kry_tol 1e-4 -no_collisions_in_solve -newton_it 10000 -dt $dt -rr -rd -grad_ls -stiffen 1 -dampen 0 -newton_tol .1 -angle_tol .0001 -rd_ccd -rr_ccd
    echo ./be_evolution 740 -o out-740d-$dt -no_self -rd_stiffness 1e2 -rd_friction 0.1 -kry_it 20 -kry_tol 1e-4 -no_collisions_in_solve -newton_it 20 -dt $dt -resolution 1 -rd -grad_ls -rd_ccd -rr_ccd -newton_tol .1
    echo ./be_evolution 740 -o out-740s-$dt -no_self -rd_stiffness 1e2 -rd_friction 0.12 -kry_it 20 -kry_tol 1e-4 -no_collisions_in_solve -newton_it 20 -dt $dt -resolution 1 -rd -grad_ls -rd_ccd -rr_ccd -newton_tol .1
#    echo ./be_evolution 750 -o out-750s-$dt -no_self -rd_stiffness 1000 -rd_friction 0.05 -kry_it 20 -kry_tol 1e-4 -no_collisions_in_solve -newton_it 10000 -dt $dt -rr -rd -grad_ls -stiffen 1 -dampen 0 -newton_tol .1 -angle_tol .0001
#    echo ./be_evolution 750 -o out-750f-$dt -no_self -rd_stiffness 1000 -rd_friction 0.002 -kry_it 20 -kry_tol 1e-4 -no_collisions_in_solve -newton_it 10000 -dt $dt -rr -rd -grad_ls -stiffen 1 -dampen 0 -newton_tol .1 -angle_tol .0001
#    echo ./be_evolution 751 -o out-751-$dt -no_self -rd_stiffness 1000 -rd_friction 0.03 -kry_it 20 -kry_tol 1e-4 -no_collisions_in_solve -newton_it 10000 -dt $dt -rr -rd -grad_ls -stiffen 1 -dampen 0 -newton_tol .1 -angle_tol .0001
#    echo ./be_evolution 130 -o out-130-$dt -no_self -rd_stiffness 100 -rd_friction 0.6 -kry_it 20 -kry_tol 1e-4 -no_collisions_in_solve -newton_it 10000 -dt $dt -rr -rd -dd -grad_ls -stiffen 1 -dampen 0 -newton_tol .01 -angle_tol .0001 -last_frame 200 -rd_ccd -rr_ccd
done | sed 's/$/ | grep Frame/' | xargs -d '\n' -P 16 -n 1 bash -c

