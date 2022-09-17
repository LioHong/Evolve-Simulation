@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 2000u "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_024_big_bang\big_bang_6001.evolve" "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_024_big_bang\big_bang_8001.evolve"
evolve_batch s 0u "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_024_big_bang\big_bang_8001.evolve" "C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_024_big_bang\big_bang_8001.txt"
pause