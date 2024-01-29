# Frog: FROm molecular dynamics to second harmonic Generation 

FROG is a software originally designed to compute the second harmonic generation response of liquids from Molecular Dynamics (MD) at the quantum level. To achieve this, it uses the PE scheme of the DALTON software which provides optical response of molecules within an electrostatic embedding.

Here are the main FROG functionalities:

  + Perform usual structural analysis such as density, molecular orientation, Radial Distribution Function or H-bond.
  + Compute automatically the optical response (polarizability, first and second hyperpolarizability) of molecules at the quantum level in an embedded environment.

Its main advantages are:

  + Can open many types of MD, works for pure liquid or mixture.
  + Can deal with 2D liquid interfaces, but also work in bulk.
  + The user can quite easily define the parameters for a molecule description or analysis.
  + Has a whole part designed to help you to treat thousands of QM calculations in a cluster.

And for honesty's sake, its drawback:

  + Not that easy to use at first.
  + Can be slow.

# Installation

The code is available in this repository in the *Frog* directory. The packaged version to get is in the main directory and is named ```Frog-version.tar.gz```. 
You can also get the code on Zenodo: https://zenodo.org/records/5998193

Frog is available as a Python3 package. To install it, first make sure your ```pip``` and ```numpy``` version are up-to-date:

```
python -m pip install --upgrade pip
python -m pip install --upgrade numpy
```

Assuming that python refers to your Python3 version. 
Then, install Frog using pip thanks to the local package:

```
python -m pip install localization_of_the_file/Frog-version.tar.gz 
```

Note that Frog requires ```MDAnalysis``` and ```pytim```. These packages should be installed automatically using the previous command. If not, or if you have problems with them, you can use the following command to make sure the latest version is installed:

```
python -m pip install --upgrade --force-reinstall MDAnalysis
python -m pip install --upgrade --force-reinstall pytim 
```

If you are using Windows, you can follow this procedure by using Anaconda with the QtConsole (do not use ```python -m```, use directly ```pip``` instead).

# Check installation

Download this repository, for instance using git:

```
git clone https://github.com/glb96/Frog.git
```

Or download this project and untar it. 
To check your installation, run the regtests available in the ```Regtests``` directory. To do so, go to this directory and run in the shell:

```
bash RUN_REGTESTS.sh
```

It should take about a few minutes depending on your machine. 
If everything is okay, you should have as final print: *If you are reading this, it means that FROG regtests have been "successful* . 
Note that you can also have a look at these inputs to see how to use Frog options. 

If you struggle to run the tests using the bash script, you can open a jupyter notebook from the tutorials and check that you can import Frog. 

# Wiki and tutorials

The wiki is available at https://glb96.github.io/Frog . 
You can also download the ```docs``` directory for offline usages.
In this case, open ```index.html``` with your preferred browser. For instance, with a Linux shell:

```
firefox index.html
```

The tutorials are contained in the Tutorials directory as well as some useful scripts. 
The wiki contains an extensive presentation of Frog objects and how to get started.
A more detailed presentation of the tutorials is also available in the wiki. 

