"""
JSON-based I/O for Simcoon solver configuration.

This module provides functions to load and save simulation configurations
using JSON format for material properties and loading paths.

The JSON format provides:
- Clear, self-documenting structure
- Easy programmatic generation and modification
- Better validation and error messages
- Compatibility with web APIs and modern tooling

Example JSON material file:
```json
{
    "name": "ELISO",
    "props": {"E": 70000, "nu": 0.3, "alpha": 1e-5},
    "nstatev": 1,
    "orientation": {"psi": 0, "theta": 0, "phi": 0}
}
```

Example JSON path file:
```json
{
    "initial_temperature": 293.15,
    "blocks": [
        {
            "type": "mechanical",
            "control_type": "small_strain",
            "corate_type": "jaumann",
            "ncycle": 1,
            "steps": [
                {
                    "time": 1.0,
                    "Dn_init": 10,
                    "Dn_mini": 1,
                    "Dn_inc": 100,
                    "DEtot": [0.01, 0, 0, 0, 0, 0],
                    "Dsigma": [0, 0, 0, 0, 0, 0],
                    "control": ["strain", "stress", "stress", "stress", "stress", "stress"],
                    "DT": 0
                }
            ]
        }
    ]
}
```
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Union, Any, TYPE_CHECKING

import numpy as np

# Re-export micromechanics classes and functions for backward compatibility
from .micromechanics import (
    MaterialOrientation,
    GeometryOrientation,
    Phase,
    Layer,
    Ellipsoid,
    Cylinder,
    Section,
    load_phases_json,
    save_phases_json,
    load_layers_json,
    save_layers_json,
    load_ellipsoids_json,
    save_ellipsoids_json,
    load_cylinders_json,
    save_cylinders_json,
    load_sections_json,
    save_sections_json,
)

# Lazy imports for solver classes - only needed for path loading
# This allows material I/O to work without simcoon._core
if TYPE_CHECKING:
    from .solver import Block, Step, StepMeca, StepThermomeca


def _get_solver_classes():
    """Lazy import of solver classes."""
    from .solver import (
        Block, Step, StepMeca, StepThermomeca,
        StateVariables, StateVariablesM, StateVariablesT,
        CONTROL_TYPES, CORATE_TYPES
    )
    return {
        'Block': Block,
        'Step': Step,
        'StepMeca': StepMeca,
        'StepThermomeca': StepThermomeca,
        'StateVariables': StateVariables,
        'StateVariablesM': StateVariablesM,
        'StateVariablesT': StateVariablesT,
        'CONTROL_TYPES': CONTROL_TYPES,
        'CORATE_TYPES': CORATE_TYPES
    }


# =============================================================================
# Material Loading
# =============================================================================

def load_material_json(filepath: Union[str, Path]) -> Dict[str, Any]:
    """
    Load material properties from a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to the JSON material file

    Returns
    -------
    dict
        Dictionary with keys:
        - 'name': UMAT name (str)
        - 'props': material properties as numpy array
        - 'nstatev': number of internal state variables (int)
        - 'orientation': dict with 'psi', 'theta', 'phi' in degrees

    Examples
    --------
    >>> mat = load_material_json('material.json')
    >>> print(mat['name'], mat['props'])
    ELISO [70000.0, 0.3, 1e-05]
    """
    with open(filepath, 'r') as f:
        data = json.load(f)

    # Convert props dict to array if needed
    if isinstance(data.get('props'), dict):
        props = np.array(list(data['props'].values()), dtype=float)
    else:
        props = np.array(data.get('props', []), dtype=float)

    # Get orientation (default to zero)
    orientation = data.get('orientation', {})
    psi = orientation.get('psi', 0.0)
    theta = orientation.get('theta', 0.0)
    phi = orientation.get('phi', 0.0)

    return {
        'name': data.get('name', 'ELISO'),
        'props': props,
        'nstatev': data.get('nstatev', 1),
        'orientation': {
            'psi': psi,
            'theta': theta,
            'phi': phi
        }
    }


def save_material_json(filepath: Union[str, Path], name: str, props: Union[np.ndarray, Dict[str, float]],
                       nstatev: int = 1, psi: float = 0.0, theta: float = 0.0, phi: float = 0.0,
                       prop_names: List[str] = None):
    """
    Save material properties to a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to save the JSON file
    name : str
        UMAT name (e.g., 'ELISO', 'EPICP')
    props : np.ndarray or dict
        Material properties. If dict, keys are property names.
    nstatev : int
        Number of internal state variables
    psi, theta, phi : float
        Euler angles in degrees
    prop_names : list of str, optional
        Property names if props is an array

    Examples
    --------
    >>> save_material_json('material.json', 'ELISO',
    ...                    {'E': 70000, 'nu': 0.3, 'alpha': 1e-5}, nstatev=1)
    """
    if isinstance(props, np.ndarray):
        if prop_names:
            props_dict = {name: float(val) for name, val in zip(prop_names, props)}
        else:
            props_dict = {f'prop_{i}': float(val) for i, val in enumerate(props)}
    else:
        props_dict = {k: float(v) for k, v in props.items()}

    data = {
        'name': name,
        'props': props_dict,
        'nstatev': nstatev,
        'orientation': {
            'psi': psi,
            'theta': theta,
            'phi': phi
        }
    }

    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)


# =============================================================================
# Path Loading
# =============================================================================

def load_path_json(filepath: Union[str, Path]) -> Dict[str, Any]:
    """
    Load simulation path from a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to the JSON path file

    Returns
    -------
    dict
        Dictionary with keys:
        - 'initial_temperature': float
        - 'blocks': list of Block objects ready for Solver

    Examples
    --------
    >>> path = load_path_json('path.json')
    >>> solver = Solver(blocks=path['blocks'])
    >>> history = solver.solve()
    """
    with open(filepath, 'r') as f:
        data = json.load(f)

    initial_temp = data.get('initial_temperature', 293.15)
    blocks = []

    for block_data in data.get('blocks', []):
        block = _parse_block(block_data)
        blocks.append(block)

    return {
        'initial_temperature': initial_temp,
        'blocks': blocks
    }


def _parse_block(block_data: Dict[str, Any]):
    """Parse a block from JSON data."""
    classes = _get_solver_classes()
    Block = classes['Block']

    # Get block type
    block_type = block_data.get('type', 'mechanical')

    # Get control types
    control_type = block_data.get('control_type', 'small_strain')
    corate_type = block_data.get('corate_type', 'jaumann')

    # Parse steps
    steps = []
    for step_data in block_data.get('steps', []):
        step = _parse_step(step_data, block_type)
        steps.append(step)

    return Block(
        steps=steps,
        nstatev=block_data.get('nstatev', 0),
        umat_name=block_data.get('umat_name', 'ELISO'),
        umat_type=block_type,
        props=np.array(block_data.get('props', []), dtype=float) if 'props' in block_data else None,
        control_type=control_type,
        corate_type=corate_type,
        ncycle=block_data.get('ncycle', 1)
    )


def _parse_step(step_data: Dict[str, Any], block_type: str):
    """Parse a step from JSON data."""
    classes = _get_solver_classes()
    StepMeca = classes['StepMeca']
    StepThermomeca = classes['StepThermomeca']

    # Common parameters
    Dn_init = step_data.get('Dn_init', 1)
    Dn_mini = step_data.get('Dn_mini', 1)
    Dn_inc = step_data.get('Dn_inc', 100)
    time = step_data.get('time', 1.0)

    # Control mode
    control = step_data.get('control', ['strain'] * 6)

    # Strain/stress targets
    DEtot = np.array(step_data.get('DEtot', [0] * 6), dtype=float)
    Dsigma = np.array(step_data.get('Dsigma', [0] * 6), dtype=float)

    if block_type == 'thermomechanical':
        DT = step_data.get('DT', 0.0)
        Q = step_data.get('Q', 0.0)
        thermal_control = step_data.get('thermal_control', 'temperature')

        return StepThermomeca(
            Dn_init=Dn_init,
            Dn_mini=Dn_mini,
            Dn_inc=Dn_inc,
            control=control,
            time=time,
            DEtot_end=DEtot,
            Dsigma_end=Dsigma,
            DT_end=DT,
            Q_end=Q,
            thermal_control=thermal_control
        )
    else:
        return StepMeca(
            Dn_init=Dn_init,
            Dn_mini=Dn_mini,
            Dn_inc=Dn_inc,
            control=control,
            time=time,
            DEtot_end=DEtot,
            Dsigma_end=Dsigma
        )


def save_path_json(filepath: Union[str, Path], blocks: List,
                   initial_temperature: float = 293.15):
    """
    Save simulation path to a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to save the JSON file
    blocks : list of Block
        List of Block objects
    initial_temperature : float
        Initial temperature in Kelvin

    Examples
    --------
    >>> step = StepMeca(DEtot_end=np.array([0.01, 0, 0, 0, 0, 0]),
    ...                 control=['strain', 'stress', 'stress', 'stress', 'stress', 'stress'])
    >>> block = Block(steps=[step], umat_name='ELISO')
    >>> save_path_json('path.json', [block])
    """
    classes = _get_solver_classes()
    StepMeca = classes['StepMeca']
    StepThermomeca = classes['StepThermomeca']

    blocks_data = []

    for block in blocks:
        block_data = {
            'type': block.umat_type,
            'control_type': block.control_type,
            'corate_type': block.corate_type,
            'ncycle': block.ncycle,
            'umat_name': block.umat_name,
            'nstatev': block.nstatev,
            'steps': []
        }

        if block.props is not None and len(block.props) > 0:
            block_data['props'] = block.props.tolist()

        for step in block.steps:
            step_data = {
                'time': step.time,
                'Dn_init': step.Dn_init,
                'Dn_mini': step.Dn_mini,
                'Dn_inc': step.Dn_inc,
                'control': step.control
            }

            if isinstance(step, StepMeca):
                step_data['DEtot'] = step.DEtot_end.tolist()
                step_data['Dsigma'] = step.Dsigma_end.tolist()

            if isinstance(step, StepThermomeca):
                step_data['DT'] = step.DT_end
                step_data['Q'] = step.Q_end
                step_data['thermal_control'] = step.thermal_control

            block_data['steps'].append(step_data)

        blocks_data.append(block_data)

    data = {
        'initial_temperature': initial_temperature,
        'blocks': blocks_data
    }

    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)


# =============================================================================
# Combined Loading (convenience function)
# =============================================================================

def load_simulation_json(material_file: Union[str, Path],
                         path_file: Union[str, Path]) -> Dict[str, Any]:
    """
    Load both material and path files and create configured blocks.

    This is a convenience function that loads material properties and path
    configuration, then assigns the material to all blocks.

    Parameters
    ----------
    material_file : str or Path
        Path to the JSON material file
    path_file : str or Path
        Path to the JSON path file

    Returns
    -------
    dict
        Dictionary with keys:
        - 'blocks': list of Block objects with material assigned
        - 'initial_temperature': float
        - 'material': dict with material properties

    Examples
    --------
    >>> sim = load_simulation_json('material.json', 'path.json')
    >>> solver = Solver(blocks=sim['blocks'])
    >>> sv = StateVariablesM(nstatev=sim['material']['nstatev'])
    >>> sv.T = sim['initial_temperature']
    >>> history = solver.solve(sv)
    """
    material = load_material_json(material_file)
    path = load_path_json(path_file)

    # Assign material to blocks that don't have explicit props
    for block in path['blocks']:
        if block.props is None or len(block.props) == 0:
            block.props = material['props']
        if block.nstatev == 0:
            block.nstatev = material['nstatev']
        if block.umat_name == 'ELISO' and material['name'] != 'ELISO':
            block.umat_name = material['name']

    return {
        'blocks': path['blocks'],
        'initial_temperature': path['initial_temperature'],
        'material': material
    }


