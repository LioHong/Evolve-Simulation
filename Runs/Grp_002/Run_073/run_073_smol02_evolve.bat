@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1u "%~dp0smol02_1300459.evolve" "%~dp0smol02_1300460.evolve"
evolve_batch s 0u "%~dp0smol02_1300460.evolve" "%~dp0smol02_1300460.txt"
pause