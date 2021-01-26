# StageM2
Accès aux codes du stage de M2

## Relate Install :

[Follow the link](https://myersgroup.github.io/relate/)

Follow installation guide from the website
<strong>Install relate in the home directory and rename the directory as relate </strong>

## Conda env with tsinfer:

This part requier **conda** or **miniconda**, as well as **python3.5** or higher.
Please make sure to install before continue

Download the **environment.yml** file 
```shell
conda env create -f environment.yml
```
This will install a full operational environment for conda, as well as **tsinfer**
You might want to give a specific name for the environment, to do so change the first line starting with **name :**
Then you will need to activate the environment everytime using :
```shell
conda activate name
```

## Bash :

This software require parallel function, if not installed, please do :
```shell
sudo apt install parallel
```

To silence parallel, you can write once :
```shell
parallel --citation
```
## Slim installation :
Now everything should be setup
