"""
This file contains the classes needed to store the data
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class Atom:
    name: str
    id: int
    x: float
    y: float
    z: float


@dataclass
class Molecule:
    name: str
    id: int
    contains: tuple
    rdf = np.ndarray
