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
    x_traj = np.ndarray
    y_traj = np.ndarray
    z_traj = np.ndarray


@dataclass
class Molecule:
    name: str
    id: int
    contains: tuple
