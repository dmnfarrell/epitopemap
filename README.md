epitopemap
==========

### A web application for visualizing mhc binding predictions

**Note:** It is now recommended that you use the epitopepredict tool instead. http://dmnfarrell.github.io/epitopepredict

#### Background

This web2py application is designed to allow binding predictions to be run for multiple proteins, such as a whole bacterial genome.

#### Installation

* Install web2py
* Unzip web2py and place in the folder where you wish to run it
* Download eptiopemap by cloning this repo
* place under web2py/applications
* Start the server using ```python web2py.py -i localhost -a password -p 8000 -K epitopemap -X```
* Go to http://localhost:8000/epitopemap in your browser

#### Other software

The ncbi-blast+ tools and muscle are needed for conservation analysis. On Ubuntu type:

`sudo apt-get install ncbi-blast+ muscle`

#### Python dependencies

Setup currently requires you to install the Python libraries yourself. This is a cleaner approach than providing them with the application and is now simple on linux. In addition we recommend you use easy_install or pip to install the packages rather than the OS package manager (e.g. apt-get) but both methods should work.
On Ubuntu type the following on the command line to install the Python modules:

```
sudo apt-get install python-pip
sudo pip install numpy pandas matplotlib biopython bokeh mpld3
```

#### Usage

See the help pages at http://dmnfarrell.github.io/epitopemap/help.html or inside the web application for details.

