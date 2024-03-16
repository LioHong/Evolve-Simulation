@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1u "%~dp0smol02_199999.evolve" "%~dp0smol02_200000.evolve"
evolve_batch s 0u "%~dp0smol02_200000.evolve" "%~dp0smol02_200000.txt"
pause