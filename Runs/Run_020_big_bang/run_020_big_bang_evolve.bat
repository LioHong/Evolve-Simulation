@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1000u "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_020_big_bang\big_bang_0.evolve" "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_020_big_bang\big_bang_1.evolve"
evolve_batch s 0u "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_020_big_bang\big_bang_1.evolve" "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_020_big_bang\big_bang_1.txt"
pause