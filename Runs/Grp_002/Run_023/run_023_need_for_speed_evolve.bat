@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1u "%~dp0need_for_speed_49999.evolve" "%~dp0need_for_speed_50000.evolve"
evolve_batch s 0u "%~dp0need_for_speed_50000.evolve" "%~dp0need_for_speed_50000.txt"
pause