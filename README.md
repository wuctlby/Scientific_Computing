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

**(1/3) WSL 💯**

## Docker
### [What is Docker](https://docs.docker.com/get-started/docker-overview/)
Docker is an open platform for developing, shipping, and running applications. Docker enables you to separate your applications from your infrastructure so you can deliver software quickly.
**Namely**, Docker would package the code, environment dependencies, and the operation function, and one could easily run the software on different platforms.

### [How to install Docker in Windows](https://docs.docker.com/desktop/setup/install/windows-install/)
Referring to this [link](https://docs.docker.com/desktop/setup/install/windows-install/), one can find the download link of Docker for Windows.

Download ==> click on the `.exe` ==> click on the next step ==> log out of Windows to complete the installation.

Currently, the Docker for Windows supports enabling WSL automatically. Only the distribution needs to be chosen.

1. Enable `Settings` > `Use WSL 2 based engine...` (enabled by default).
2. Enable `Settings` > `Resources` > `WSL Integration`.

### [How to use Docker](https://docs.docker.com/get-started/introduction/)
Referring to this introduction [link](https://docs.docker.com/get-started/introduction/), one can find the basic concept and utilization of Docker.

**But this readme is more about how to install and combine these tools, after all the tools (based on WSL, Docker, and VS Code) were installed, some basic methods of using will be introduced.**

**(2/3) Docker 💯**

## VS Code
### [What is VS Code](https://www.youtube.com/watch?v=B-s71n0dHUk)
VS Code is a free and lightweight code editor. It is powerful, supports most available hardware and platforms, and contains a lot of extensions. **Namely**, one could easily use it to develop almost all programming languages.

### [How to install VS Code](https://code.visualstudio.com/download)
Referring to this [link](https://code.visualstudio.com/download), one can download the VS Code for Windows.

Download ==> click on the `.exe` ==> click on the next step (**be sure to choose the `Add to PATH` in `Select Additional Tasks`**) until on install ==> reboot to make available the writing of the `PATH` of VS Code into the system path.

**Must be sure to choose `Add to PATH` and reboot, then the folder in WSL can be easily accessed**.

### [How to use the VS code]
Due to this readme is more about how to install and combine these tools, here only how to install the necessary extensions will be introduced.

#### Open VS code
There are two ways one can open the vs code:
1. Click on the icon of VS code.
2. Enter the command in WSL.
    ```bash
    code . # it doesn't matter without the `.` as usual
    ```
    **Attention: if one has opened the VS Code, one now can access the folder in WSL, read and edit all files. But, actually, what this uses is `Windows-based tools`, one can't run or debug the code with Linux environment and also wouldn't have access to features such as autocomplete, debugging, or lining. Make sure that you have already installed the extensions below, then start working, finally no more issues, even [the cross-OS challenges](https://code.visualstudio.com/docs/remote/wsl)**

#### Extensions
1. WSL (re-open the folder)
2. Dev Containers
3. Docker
4. Programming langues package, like python


