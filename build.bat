setlocal enableDelayedExpansion

for /f "eol=: tokens=3 delims= " %%a in ('find "VERSION = " Parse_Device_Log-GUI-CPIW.py') do (
   set "VERSION=%%a"
   echo !VERSION!
)
set "VERSION=%VERSION:~1%"
set "VERSION=%VERSION:~0,-1%"

::set VERSION=3.0.16.0

PyInstaller --onefile --hidden-import=libscrc --distpath .\build\CPIW_Parser_%VERSION%\ -n "CPIW_File_Parser_"%VERSION% Parse_Device_Log-GUI-CPIW.py 

chdir .\build\CPIW_Parser_%VERSION%\

"C:\Program Files\7-Zip\7z.exe" a CPIW_File_Parser-%VERSION%.zip .

chdir ..\..