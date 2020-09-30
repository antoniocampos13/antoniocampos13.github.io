---
Title: Setting Up Your Unix Computer for Bioinformatics Analysis
Date: 2020-09-30 18:00
Author: Antonio Victor Campos Coelho
Tags: Bioinformatics
---

## Introduction

In this post I will show how I set up my Unix machine to use Bioinformatics programs and tools. I am currently using Ubuntu 20.04 LTS (Focal Fossa) on a [Windows Subsystem for Linux (WSL2)](https://www.digitalocean.com/community/posts/trying-the-new-wsl-2-its-fast-windows-subsystem-for-linux) on Windows 10, so no GUI today!

## Preparing the system

First, it is recommended that we upgrade the system. Open the command line terminal in your machine and copy and paste or type the following commands, pressing Enter after each one (make sure you type your password correctly whenever asked):

```bash
sudo apt-get update
sudo apt-get upgrade
```

Then we must install some useful libraries, especially to be sure that all future libraries we need will be installed and work properly. Some of these (e.g. default-jdk, the Java libraries), may already be installed in your system, but just to ensure:

```bash
sudo apt-get install -y curl unzip build-essential ncurses-dev
sudo apt-get install -y byacc zlib1g-dev python-dev git cmake
sudo apt-get install -y default-jdk ant
```

## Installing (mini)conda

Now we will install [miniconda](https://conda.io/miniconda.html). What is miniconda? Miniconda is a simplified version of Conda, an environment management system. Every program we install on your computers depend on other programs to work. So if a program X needs a program Y to work, it may stop working if Y gets an update that for some reason is incompatible with the original X program.

Thus, environments were developed to solve this kind of problem, because they serve to isolate groups of programs, ensuring only compatible versions of software are working together. Therefore, miniconda serves to create and manage environments. The best practice is that we should create one environment for one specific use. In my case, I installed miniconda to create a environment and populate it with tools used for several Bioinformatics analysis. Other people can create environments for other uses with specific programs needed and so on. Other advantage of miniconda is that the configuration files for environments can be shared with others, ensuring **backup** and **reproducibility**.

Without further ado, let's finally install miniconda. Since I am using a Unix with Python 3.7.7 pre-installed, the version of the installer [is this one](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh). Check the [installation page](https://docs.conda.io/en/latest/miniconda.html#linux-installers) if you have a different Python version.

You can download the installer from your browser or via command line:

```bash
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Then, go to the folder where the installer was downloaded and run the script:

```bash
bash Miniconda3-latest-Linux-x86_64.sh

./Miniconda3-latest-Linux-x86_64.sh # same effect
```

When the installation finishes, we must initialize conda:

```bash
miniconda3/condabin/conda init
```

Close the terminal and open it again. Now miniconda must be ready to use. Check by typing and pressing Enter:

```bash
conda
```

Then, I added two **channels**. Channels are ["the locations where packages are stored"](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html). Miniconda has the `defaults` channel pre-configured. The two channels in question are dedicated to Bioinformatics and Data analysis programs, which may not be present in the default channels, so we must add them.

## Configuring miniconda channels

Once again in the terminal enter the following commands:

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
```

Miniconda sets up priorities in the list of channels it receives. When we need to install some program, miniconda will search in the higher-priority channels first, then in the channels with lower-priority. "Different channels can have the same package" and you can "safely put channels at the bottom of your channel list to provide additional packages that are not in the default channels" as stated in the [official website](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html). The flag `--add` adds the respective channels (`bioconda` and `conda-forge`) to the **top** of the priorities list. If you want to give lower priority, putting them in the **bottom** of the list, use the `--append` command instead. Thus, according to the command above, the order of channel priorities in our new miniconda installation will be: `conda-forge`, `bioconda` and lastly, `defaults`.

## Create an environment for Bioinformatics programs

Now that miniconda is configured, we will create the environment that will receive them. I will name it `bioenv`. You can choose whatever name you like! 

```bash
conda create -y --name bioenv python=3.8
```

## Activating and deactivating an environment

With the `bioenv` created, we must **activate** it:

```bash
conda activate bioenv
```

We need to perform this step every time we want to use the programs that we will install in this environment. If you do not need to use the environment for the moment, simply **deactivate** it:

```bash
conda deactivate
```

Simply **activate** it again when needed.

## Installing programs
<!-- COLOCAR LINK DO REPOSITÓRIO AQUI E NO CODE BLOCK-->
Now we can finally install our programs. Activate the environment again (only if you have deactivated it). Download the `bioenv.txt` file in my GitHub repository. This file contains a selection of most used Bioinformatics programs (hat tip to [Dr. István Albert](https://www.biostarhandbook.com/index.html))

```bash
cat bioenv.txt | xargs conda install -y
```

## Backing up and restoring your environment configuration

Miniconda has a special command to backup your environment configuration. **Activate** (if needed) the environment you want to backup and enter the command:

```bash
conda env export > bioenv.yml
```

It will result in a `YAML` file in the current working folder containing all configurations in your environment. Again, I named the file `bioenv.yml` but you can choose whatever you like. Note that if you already have a `bioenv.yml` in your directory, it will be overwritten, so be careful.

To restore this environment in your computer, or on other computer, first install miniconda again, and then use the command:

```bash
conda env create -f bioenv.yml
```

The `-f` flag means you are creating an environment using the configurations in the `bioenv.yml` file. The first line of the yml file sets the new environment's name, so you can change it in the file if you like. It will also restore the channels configured in the previous installation of miniconda.

## Conclusion

This is how I configured my system so I could use the major Bioinformatics tools out there. In summary, I:

* Prepared an Unix (Ubuntu) system;
* Installed miniconda, a environment manager;
* Configured channels so I could retrieve desired software;
* Created an environment, showed how to activate and deactivate it, and finally installed software in it;
* Showed how to backup your environment for safekeeping or share with others.

In future posts I will demo some uses of the programs I installed in the new environment.

## References

[Trying the New WSL 2. It's Fast! (Windows Subsystem for Linux) | DigitalOcean](https://www.digitalocean.com/community/posts/trying-the-new-wsl-2-its-fast-windows-subsystem-for-linux)

[Miniconda](https://conda.io/miniconda.html)

[Miniconda installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh)

[Miniconda & Conda documentation](https://docs.conda.io/en/latest/miniconda.html#linux-installers)

[Managing channels; conda 4.8.4.post65+1a0ab046 documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html)

[The Biostar Handbook: 2nd Edition](https://www.biostarhandbook.com/index.html)