@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 0u "%~dp0p_in.evolve" "%~dp0p_out.txt"
pause