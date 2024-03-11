@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1000000u "%~dp0smol02.evolve" "%~dp0smol02_1e06.evolve"
evolve_batch s 0u "%~dp0smol02_1e06.evolve" "%~dp0smol02_1e06.txt"
pause