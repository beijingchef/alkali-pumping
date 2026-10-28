"""Bounded caches for transient densities and compact display observables."""
import numpy as np
import streamlit as st
from ..physics import time_evolution, optical_readout, probe_scan
from ..physics.constants import ATOMS, MU_B_OVER_H_HZ_PER_T, ground_hyperfine_lande_g
from ..version import PHYSICS_MODEL_VERSION

LARGE_DENSITY_BYTES = 8 * 1024**2
SMALL_EVOLUTION_ENTRIES = 2
LARGE_EVOLUTION_ENTRIES = 1
DISPLAY_ENTRIES = 4


def estimated_density_bytes(config_A, config_B, cell, total_time):
    configs = [config_A]
    if config_B is not None and config_B['atom_name'] != config_A['atom_name']:
        configs.append(config_B)
    fastest, sizes = 0., []
    for config in configs:
        atom = ATOMS[config['atom_name']]
        upper = atom['I'] + .5
        g = max(abs(ground_hyperfine_lande_g(config['atom_name'], atom, F)) for F in (upper-1, upper))
        fastest = max(fastest, 2*upper*g*MU_B_OVER_H_HZ_PER_T*1e-9*abs(cell['static_field_nT']))
        sizes.append(int(2*(2*atom['I']+1)))
    samples = min(20001, max(801, int(np.ceil(20*fastest*total_time))+1))
    return samples * sum(n*n for n in sizes) * np.dtype(complex).itemsize


@st.cache_data(max_entries=SMALL_EVOLUTION_ENTRIES, show_spinner=False)
def _small_evolution(config_A, config_B, beams, cell, t1, t2, model_version, solver_version, physics_version):
    return time_evolution.compute_density_evolution(config_A, config_B, beams, cell, t1, t2)


@st.cache_data(max_entries=LARGE_EVOLUTION_ENTRIES, show_spinner=False)
def _large_evolution(config_A, config_B, beams, cell, t1, t2, model_version, solver_version, physics_version):
    return time_evolution.compute_density_evolution(config_A, config_B, beams, cell, t1, t2)


def evolution_cached(config_A, config_B, beams, cell, t1, t2):
    cache = (_large_evolution if estimated_density_bytes(config_A, config_B, cell, t2) >= LARGE_DENSITY_BYTES
             else _small_evolution)
    return cache(config_A, config_B, beams, cell, t1, t2,
                 time_evolution.TIME_EVOLUTION_MODEL_VERSION,
                 time_evolution.TIME_EVOLUTION_SOLVER_VERSION, PHYSICS_MODEL_VERSION)


@st.cache_data(max_entries=DISPLAY_ENTRIES, show_spinner=False)
def display_cached(species, zeeman_axis, moment_axis):
    view = time_evolution.density_display_view(species, zeeman_axis, moment_axis)
    # Cache populations and moments, never another full complex trajectory.
    view.pop('density_matrix')
    return view


@st.cache_data(max_entries=4, show_spinner=False)
def _steady_state(config_A, config_B, beams, cell, solver_version, physics_version):
    return probe_scan.pump_on_steady_state(config_A, config_B, beams, cell)


def steady_cached(config_A, config_B, beams, cell):
    # Streamlit hashes supplied arguments, not omitted default arguments.
    return _steady_state(config_A, config_B, beams, cell,
                         time_evolution.TIME_EVOLUTION_SOLVER_VERSION, PHYSICS_MODEL_VERSION)


@st.cache_data(max_entries=8, show_spinner=False)
def spectrum_cached(species, probe, cell, coefficients, lower, upper, physics_version=PHYSICS_MODEL_VERSION):
    line = probe['line']
    return optical_readout.density_probe_spectrum(
        atom=ATOMS[species['atom_name']], ground_states=species['ground_states'],
        density_matrix=species['density_matrix'][0], q_axis='z', line=line,
        detunings_MHz=np.linspace(lower, upper, 2001),
        **{name: probe[name] for name in ('k_axis', 'azimuth_deg', 'ellipticity_deg', 'path_length_cm')},
        density_cm3=species['density_cm3'], n2_pressure_torr=cell['n2_pressure_torr'],
        temperature_C=cell['temperature_C'], n2_width_MHz_per_torr=coefficients[line]['width'],
        n2_shift_MHz_per_torr=coefficients[line]['shift'],
    )


@st.cache_data(max_entries=4, show_spinner=False)
def probe_cached(species, probe, cell, coefficients,
                 model_version=time_evolution.TIME_EVOLUTION_MODEL_VERSION, physics_version=PHYSICS_MODEL_VERSION):
    line = probe['line']
    return optical_readout.density_probe_signals(
        atom=ATOMS[species['atom_name']], ground_states=species['ground_states'],
        density_matrix=species['density_matrix'], q_axis=species['q_axis'], **probe,
        density_cm3=species['density_cm3'], n2_pressure_torr=cell['n2_pressure_torr'],
        temperature_C=cell['temperature_C'], n2_width_MHz_per_torr=coefficients[line]['width'],
        n2_shift_MHz_per_torr=coefficients[line]['shift'],
    )


@st.cache_data(max_entries=4, show_spinner=False)
def scan_display_cached(species, axis, direction):
    """Cache only a single stationary matrix per display basis, never a trajectory."""
    from ..physics.scan_basis import scan_population_view
    return scan_population_view(species, axis, direction)
