@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 100u "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_019_big_bang\big_bang_801.evolve" "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_019_big_bang\big_bang_901.evolve"
evolve_batch s 0u "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_019_big_bang\big_bang_901.evolve" "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_019_big_bang\big_bang_901.txt"
pause