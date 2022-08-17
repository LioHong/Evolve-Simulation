@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1u "C:\Users\Lio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_006_big_bang\big_bang_0.evolve" "C:\Users\Lio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_006_big_bang\big_bang_1.evolve"
evolve_batch s 0u "C:\Users\Lio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_006_big_bang\big_bang_1.evolve" "C:\Users\Lio Hong\Documents\LioHong\Evolve-Simulation\\Runs\Run_006_big_bang\big_bang_1.txt"
pause