[ -e "test.txt" ] && rm -f test.txt
for i in {16,24,32,40,48,56,64,80,96,112,128,144,160,176,192,208,224,240,256}; do
echo -n "$i " >> test.txt;
../multigrid_embedded_poisson_3d_v2_nocona --dimension=3 --scalar=double --n-cell="$i^D" --level-set="scale[0.25](trefoil_knot(0.75))" --min-dist-to-vertex=0.02 --example="interface(({2},{0}), ({1},{1}), aggregate)" --solver=petsc-cg --rel-tol=1e-14 --infty-norm-error-filename=test.txt;
done
# u^-    = x + 2*y + 3*z
# beta^- = 1 + 3*pow<2>(x) + 2*pow<2>(y) + pow<2>(z)
# u^+    = pow<2>(x) + 2*pow<2>(y) + 3*pow<2>(z)
# beta^+ = 1
