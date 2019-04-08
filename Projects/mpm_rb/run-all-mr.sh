for dt in 0.01 0.005 0.002 0.001; do
#echo ./mpm_rb -3d 21 -o out-21h-$dt -double -scale_E .1 -rd_stiffness 1e2 -max_dt $dt -regular_seeding -rd_friction 0.3 -rd -rr -grad_ls -last_frame 400 -rr_ccd -rd_ccd -solver_iterations 20 -angle_tol .0001
#echo ./mpm_rb -3d 21 -o out-21l-$dt -double -scale_E .1 -rd_stiffness 1e2 -max_dt $dt -regular_seeding -rd_friction 0.05 -rd -rr -grad_ls -rr_ccd -rd_ccd -solver_iterations 20 -angle_tol .0001
#echo ./mpm_rb -3d 13 -o out-13l-$dt -double -rd_stiffness 1e3 -max_dt $dt -T .1 -T .1 -T .2 -rd_friction .1 -grad_ls -rr -rr_ccd -rd_ccd -solver_iterations 20 -angle_tol .0001 -newton_tolerance .01
#echo ./mpm_rb -3d 13 -o out-13h-$dt -double -rd_stiffness 1e3 -max_dt $dt -T .1 -T .1 -T .2 -rd_friction .125 -grad_ls -rr -rr_ccd -rd_ccd -solver_iterations 20 -angle_tol .0001 -newton_tolerance .01
#echo ./mpm_rb -3d 8 -o out-8l-$dt -double -scale_E .1 -rd_stiffness 10 -max_dt $dt -T .1 -T .1 -T .2 -regular_seeding -rd_friction 0.1 -di -grad_ls -newton_tolerance .01 -rr_ccd -rd_ccd -solver_iterations 20 -angle_tol .0001
#echo ./mpm_rb -3d 8 -o out-8h-$dt -double -scale_E .1 -rd_stiffness 10 -max_dt $dt -T .1 -T .1 -T .2 -regular_seeding -rd_friction 0.125 -di -grad_ls -newton_tolerance .01 -rr_ccd -rd_ccd -solver_iterations 20 -angle_tol .0001
#echo ./mpm_rb -3d 9 -o out-9l-$dt -double -scale_E .1 -rd_stiffness 1000 -max_dt $dt -T .1 -T .1 -T .2 -regular_seeding -rd_friction 0.1 -rd -grad_ls -newton_tolerance .1 -rr_ccd -rd_ccd -solver_iterations 20 -angle_tol .0001
#echo ./mpm_rb -3d 9 -o out-9h-$dt -double -scale_E .1 -rd_stiffness 1000 -max_dt $dt -T .1 -T .1 -T .2 -regular_seeding -rd_friction 0.125 -rd -grad_ls -newton_tolerance .1 -rr_ccd -rd_ccd -solver_iterations 20 -angle_tol .0001
echo ./mpm_rb -3d 71 -o out-71h-$dt -double -scale_E .1 -rd_stiffness 1e2 -max_dt $dt -regular_seeding -rd_friction 0.25 -grad_ls -rr -rd -di -solver_iterations 20 -angle_tol .0001
echo ./mpm_rb -3d 71 -o out-71l-$dt -double -scale_E .1 -rd_stiffness 1e2 -max_dt $dt -regular_seeding -rd_friction 0.1 -grad_ls -rr -rd -di -solver_iterations 20 -angle_tol .0001
#echo ./mpm_rb -3d 72 -o out-72h-$dt -double -rd_stiffness 1e1 -max_dt $dt -rd_friction 0.08 -rr -grad_ls -solver_iterations 20 -angle_tol .0001
#echo ./mpm_rb -3d 72 -o out-72l-$dt -double -rd_stiffness 1e1 -max_dt $dt -rd_friction 0.03 -rr -grad_ls -solver_iterations 20 -angle_tol .0001
#echo ./mpm_rb -3d 12 -o out-12h-$dt -double -rd_stiffness 1e1 -max_dt $dt -T .2 -rd_friction 0.3 -rr -grad_ls -rr_ccd -rd_ccd -solver_iterations 20 -angle_tol .0001
#echo ./mpm_rb -3d 12 -o out-12l-$dt -double -rd_stiffness 1e1 -max_dt $dt -T .2 -rd_friction 0.01 -rr -grad_ls -rr_ccd -rd_ccd -solver_iterations 20 -angle_tol .0001
done | sed 's/$/ | grep -i Frame/' | xargs -d '\n' -P 16 -n 1 bash -c
