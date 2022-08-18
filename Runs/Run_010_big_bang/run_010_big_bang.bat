@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1u "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_010_big_bang\big_bang_10000.evolve" "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_010_big_bang\big_bang_10001.evolve"
evolve_batch s 0u "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_010_big_bang\big_bang_10001.evolve" "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_010_big_bang\big_bang_10001.txt"
pause