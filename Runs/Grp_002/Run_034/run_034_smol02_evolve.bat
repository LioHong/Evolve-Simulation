@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1u "%~dp0smol02_200999.evolve" "%~dp0smol02_201000.evolve"
evolve_batch s 0u "%~dp0smol02_201000.evolve" "%~dp0smol02_201000.txt"
pause