@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 10000u "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_012_big_bang\big_bang_1.evolve" "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_012_big_bang\big_bang_10001.evolve"
evolve_batch s 0u "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_012_big_bang\big_bang_10001.evolve" "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_012_big_bang\big_bang_10001.txt"
pause