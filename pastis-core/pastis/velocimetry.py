import numpy as np


def initialize_RV(datadict):
    """
    Function to initialize PASTIS_RV
    """
    global stepccf, RV_CCF_ALL
    # global RVdatadict
    print("Initializing RV module...")

    for key in datadict.keys():
        if datadict[key]["type"] == "RV":
            if "vspan" in datadict[key]["data"]:
                print("... Vspan detected ...")
                print("... Oversampling the CCF ...")
                stepccf = 0.001
                RV_CCF_ALL = np.arange(-250.0, 250.0, stepccf)
            else:
                stepccf = 0.1
                RV_CCF_ALL = np.arange(-250.0, 250.0, stepccf)

    print("... DONE! \n")
    return  # RVdatadict
