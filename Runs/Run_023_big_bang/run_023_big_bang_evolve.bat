@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1000u "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_023_big_bang\big_bang_8001.evolve" "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_023_big_bang\big_bang_9001.evolve"
evolve_batch s 0u "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_023_big_bang\big_bang_9001.evolve" "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_023_big_bang\big_bang_9001.txt"
pause