@echo off
cd "C:\Program Files (x86)\Evolve"
evolve_batch s 1u "%~dp0need_for_speed_37499.evolve" "%~dp0need_for_speed_37500.evolve"
evolve_batch s 0u "%~dp0need_for_speed_37500.evolve" "%~dp0need_for_speed_37500.txt"
pause