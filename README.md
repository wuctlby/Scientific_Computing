# Scientific_Computing
This is a repository for the course of Scientific Computing, specified with tools for CERN ALICE

## WSL
### [What is WSL](https://learn.microsoft.com/en-us/windows/wsl/about)

Windows Subsystem for Linux (WSL) is a feature of Windows that allows you to run a Linux environment on your Windows machine, without needing a separate virtual machine or dual booting.

### Required `Virtualization Technology` for WSL
Step 1:

  - Restart your computer and enter BIOS/UEFI settings(typically by pressing F2, Del, or Esc during startup).
  
Step 2:

  - Locate the virtualization option (e.g., Intel Virtualization Technology, VT-x, AMD-V, or SVM Mode) and enable it.
  
Step 3:

  - Save changes and exit BIOS/UEFI.


### [How to install WSL](https://learn.microsoft.com/en-us/windows/wsl/install)
Step 1:
- #### Enable the WSL feature
  ```powershell
  dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
  ```
Step 2:
- #### Install the default Linux distribution(`Ubuntu`)
  ```powershell   
  wsl --install
  ```
or
- #### Install a specific distribution
  
  See the list of available Linux distributions of WSL:
  ```
  wsl --list --online
  ```
  Install a specific Linux distribution of WSL:
  ```
  wsl --install -d <Distribution Name>
  ```
or
- #### Install other distributions, like `Almalinux 9`
  Open the Microsoft Store, search `Almalinux 9`, and download.

## Docker
### [What is Docker](https://docs.docker.com/get-started/docker-overview/)
Docker is an open platform for developing, shipping, and running applications. Docker enables you to separate your applications from your infrastructure so you can deliver software quickly.
#Namely#, Docker would package the code, environment dependencies, and the operation function, and one could easily run the software on different platforms.

### [How to install Docker in Windows](https://docs.docker.com/desktop/setup/install/windows-install/)
Referring to this [link](https://docs.docker.com/desktop/setup/install/windows-install/), one can find the download link of Docker for Windows.

Download ==> open the `.exe` ==> click on the next step ==> log out the Windows to complete the installation.
