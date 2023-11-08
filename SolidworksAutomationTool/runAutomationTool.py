# this is a prototype of running the automation tool with python
import subprocess
print("running the automation tool executable")
# store the path as a raw string so that the backslashes won't be misinterpreted as escape characters
automationToolPath = r"C:\SolidworksAutomationTool\bin\Debug\net7.0\SolidworksAutomationTool.exe"
# ask the OS to run this executable in a new console window
subprocess.Popen(automationToolPath, creationflags= subprocess.CREATE_NEW_CONSOLE)