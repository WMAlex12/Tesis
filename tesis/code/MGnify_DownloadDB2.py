import requests
import pandas as pd
import torch
import torch.nn as nn
from numba import jit, cuda
import lib.variable_utils



# Device configuration
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(torch.cuda.is_available())

from jsonapi_client import Session
import pandas as pd
from package import *
import importlib

from lib.variable_utils import get_variable_from_link_or_input

# You can also just directly set the accession variable in code, like this:
# accession = "MGYS00005292"
accession = get_variable_from_link_or_input('MGYS', 'Study Accession', 'MGYS00005384')

with Session("https://www.ebi.ac.uk/metagenomics/api/v1") as mgnify:
    resources = map(lambda r: r.json, mgnify.iterate(api_endpoint))
    resources = pd.json_normalize(resources)
    resources.to_csv(f"{api_endpoint}.csv")
resources
