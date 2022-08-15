@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1u "path_in.evolve" "path_out.evolve"
evolve_batch s 0u "path_out.evolve" "path_out.txt"
pause