@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1000000u "%~dp0smol.evolve" "%~dp0smol_1e06.evolve"
evolve_batch s 0u "%~dp0smol_1e06.evolve" "%~dp0smol_1e06.txt"
pause