@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1000u "%~dp0smol02_20000.evolve" "%~dp0smol02_21000.evolve"
evolve_batch s 0u "%~dp0smol02_21000.evolve" "%~dp0smol02_21000.txt"
pause