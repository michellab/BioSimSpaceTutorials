Author: Lester Hedges<br>
Email:&nbsp;&nbsp; lester.hedges@bristol.ac.uk

# Nodes: _Interoperable workflow components_

The companion notebook for this section can be found [here](https://github.com/michellab/BioSimSpaceTutorials/blob/4844562e7d2cd0b269cead56562ec16a3dfaef7c/01_introduction/04_writing_nodes.ipynb)

So far we have been working with BioSimSpace in a rather ad hoc fashion. While this intereactive exploration is a great way of learning and prototyping ideas, it is not a good way of producing a reproducible and interoperable script that can be shared with others. For example, we created processes that used specific packages such as AMBER and GROMACS. If a user didn't have these available on their system, then the script simply wouldn't work. We also used hard-coded paths to input files. This means the user would have to edit the paths each time they ran the script with different input, which would quickly become tedious.

In order to solve this problem, a core concept of BioSimSpace is the interoperable workflow component, or _node_. These are robust and portable Python scripts that typically do a small, well-defined piece of work. All inputs and outputs from the node are validated and the node is written in a such a way that it is _independent_ of the underlying software packages, i.e. the same script can work with a range of different packages. In addition, nodes are aware of the environment in which they are run, so can be used interactively, from the command-line, or within a workflow engine.

While it is possible to write a node directly as a Python script, we suggest that the best way of writing one is inside of a [Jupyter](http://jupyter.org) notebook. As you've already seen, the interactive notebook environment provides a fantastic way of prototyping and documenting your node and will allow a user to interact with it directly on a remote cloud server, such as [notebook.biosimspace.org](https://notebook.biosimspace.org). The notebook can provide a complete record of your work, inlcuding documentation, visualisation, and graphs. When you are happy with the node, you can download it as a regular Python script (by clicking on `File/Download As/Python` in JupyterHub or `File/Export Notebook As/Export Notebook to Executable Script` in JupyterLab) and run it directly from the command-line on your workstation, laptop, or on a high-performance computing cluster. Any interactive BioSimSpace elements, such as molecular visualisations, will simply be ignored when run this way.


## An example: Minimisation

In the rest of the notebook you'll learn how to use BioSimSpace to write a robust and interoperable workflow node to perform energy minimisation on a molecular system.

As always, we'll first need to import BioSimSpace:


```python
import BioSimSpace as BSS
```

We begin by creating a `Node` object. This is the core of our molecular workflow component. It defines what it does, what input is needed, and the output that is produced.


```python
node = BSS.Gateway.Node("A node to perform energy minimisation and save the minimised molecular configuration to file.")
```

We'll now set the author and license of the node. When nodes are run the the authorship can be queried so that people can get credit for their work. Eventually, BioSimSpace nodes also will also contain built in tracking information to determine how many times they are run.


```python
node.addAuthor(name="Lester Hedges", email="lester.hedges@bristol.ac.uk", affiliation="University of Bristol")
node.setLicense("GPLv3")
```

Nodes require inputs. To specify inputs we use the `BSS.Gateway` package, which is used as a bridge between BioSimSpace and the outside world. This will allow us to document the inputs, define their type, and specify any constraints on their allowed values. Here we will need a set of files that define the molecular system, and an integer that indicates the number of minimisation steps to perform.


```python
node.addInput("files", BSS.Gateway.FileSet(
    help="A set of molecular input files.")
)

node.addInput("steps", BSS.Gateway.Integer(
    help="The number of minimisation steps.",
    minimum=0,
    maximum=1000000,
    default=10000)
)

node.addInput("engine", BSS.Gateway.String(
    help="The molecular dynamics engine",
    allowed=BSS.MD.engines(),
    default="auto")
)
```

Note that the input requirements `steps` and `engine` have default values, so are optional.

We now need to define the output of the node. In this case we will return a set of files representing the minimised molecular system.


```python
node.addOutput("minimised", BSS.Gateway.FileSet(help="The minimised molecular system."))
```

When working interactively within a Jupyter notebook we need a way to allow users to set the input requirements. The `node.showControls` method will display a graphical user interface (GUI), from which inputs can be set. All of the elements for this GUI are automatically generated by the `addInput` and `addOutput` functions above. As you'll see in the next section, if we were to run the same node from the command-line, we would instead get an automatically generated [argparse](https://docs.python.org/3/library/argparse.html) parser.

Note that the GUI requires active user input. All input requirements that don't have a default value _must_ be set before the node can proceed. If you try to query the node for one of the user values then an error will be raised. For bounded integer inputs you can use a slider to set the value, or type in the input box and press enter.

When working interactively you will typically be running on a remote server where you won't have access to the local filesystem. In this case you'll need to upload files for any of the `File` or `FileSet` input requirements. The GUI below will provide buttons that allow you to browse your own filesystem and select files. Since Jupyter has a limit of 5MB for file transfers, we provide support for compressed formats, such as `.zip` or `.tar.gz`. (A single archive can contain a set of files, allowing you to set a single value for a `FileSet` requirement.) We've provided some example input files that can be used in the training notebooks, which are available to download from the links below. These can then be re-uploaded using the GUI.

AMBER: [ala.crd](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/amber/ala/ala.crd), [ala.top](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/amber/ala/ala.top)

GROMACS: [kigaki.gro](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/gromacs/kigaki/kigaki.gro), [kigaki.top](https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/gromacs/kigaki/kigaki.top)


When uploading files the name of the current file(s) will replace the `Upload` button. If you need to change the file, simply click on the button again and choose a new file.


```python
node.showControls()
```

![Notebook GUI](https://github.com/michellab/BioSimSpaceTutorials/blob/c06201e9464732df5fd64fa560779ef333f59651/01_introduction/assets/04_gui.png)


Once all requirements are set then we can acces the values using the `node.getInput` method. The first time this is called the `node` will automatically validate all of the input and report the user if any errors were found.

We'll now create a molecular system using the input files uploaded by the user. As in the previous section, we don't need to specify the format of the files, since this is automatically determined by BioSimSpace. (BioSimSpace has support for a wide range of formats and can convert between many formats too.)


```python
system = BSS.IO.readMolecules(node.getInput("files"))
```

As learned in the previus notebook, in order to run a minimisation we need to define a protocol. This can be done using the `BSS.Protocol` package. Here we will create a "best practice" minimisation protocol, overriding the number of steps with the input from the user.


```python
protocol = BSS.Protocol.Minimisation(steps=node.getInput("steps"))
```

We now have everything that is required to run a minimisation. To do so, we use the `BSS.MD` package to find an appropriate molecular dynamics package on our current environment. What package is found will depend upon both the system and protocol, as well as the hardware that is available to the user. (For example, the user can choose to find packages with GPU support.)

Note that this is different to the previous section, where we specifically launched AMBER and GROMACS processes ourselves. This is what makes the node interoperable, i.e. it will work regardles of what MD packages are installed. (As long as we find a package that supports minimisation and supports a molecular file format to which we can convert the input system.) By adding the optional `engine` requirement we have also allowed the user to override the `auto` setting if they prefer to use a specific engine.

(By default, the `run` function automatically starts the process so it will be running as once you execute the cell below.)


```python
process = BSS.MD.run(system, protocol, engine=node.getInput("engine"))
```

We now wait for the process to finish, then check whether there were any errors before continuing. If errors were raised, then we raise an exception and print the last 10 lines of stdout and stderr to the user.


```python
process.wait()

if process.isError():
    print(process.stdout(10))
    print(process.stdout(10))
    raise RuntimeError("The process exited with an error!")
```

When the process has finished running we can get the minimised molecular configuration. We will save this to file using the same format as the original system, and set the `minimised` output requirement to the list of file names that were written.


```python
node.setOutput("minimised",
    BSS.IO.saveMolecules("minimised", process.getSystem(), system.fileFormat()))
```

Finally, we validate that the node completed succesfully. This will check that all output requirements are satisfied and that no errors were raised by the user. Any file outputs will be available for the user to download as a compressed archive.

Note that the validation will fail until the cell above finishes running.


```python
node.validate()
```




<a href='https://github.com/michellab/BioSimSpaceTutorials/blob/0800c84f845d7f7863bb916aaeb0d3bba3bf4137/01_introduction/assets/output.zip' target='_blank' download='output.zip'>output.zip</a><br>



Once we are satisfied with our node we can choosed to download it as a regular Python script that can be run from the command-line.

In JupyterHub, click on: `File/Download As/Python`\
In JupyterLab, click on: `File/Export Notebook As/Export Notebook to Executable Script`

That's it, you've now succesfully executed your first BioSimSpace node!