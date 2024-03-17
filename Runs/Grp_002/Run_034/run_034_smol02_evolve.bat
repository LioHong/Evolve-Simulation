@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1u "%~dp0smol02_299999.evolve" "%~dp0smol02_300000.evolve"
evolve_batch s 0u "%~dp0smol02_300000.evolve" "%~dp0smol02_300000.txt"
pause