"""Bounded ACC caches separated by physical state, dynamics, and section readout."""
import streamlit as st

from ..physics.system_model import (
    build_species_linear_model, solve_alkali_equilibrium,
    prepare_alkali_response_contexts, compute_alkali_response,
)
from ..version import PHYSICS_MODEL_VERSION


@st.cache_data(max_entries=24, show_spinner=False)
def models_cached(species_A, species_B, beams, common, physics_model_version=PHYSICS_MODEL_VERSION):
    return {
        "A": build_species_linear_model(species_A, beams, common),
        "B": build_species_linear_model(species_B, beams, common) if species_B is not None else None,
    }


@st.cache_data(max_entries=24, show_spinner=False)
def equilibrium_cached(models, common, physics_model_version=PHYSICS_MODEL_VERSION):
    return solve_alkali_equilibrium(models["A"], models["B"], common)


@st.cache_data(max_entries=24, show_spinner=False)
def contexts_cached(equilibrium, common, physics_model_version=PHYSICS_MODEL_VERSION):
    return prepare_alkali_response_contexts(equilibrium, common)


@st.cache_data(max_entries=48, show_spinner=False)
def atomic_cached(equilibrium, label, rf, common, physics_model_version=PHYSICS_MODEL_VERSION):
    return compute_alkali_response(equilibrium, label, rf, common, contexts_cached(equilibrium, common))


@st.cache_data(max_entries=48, show_spinner=False)
def optical_cached(equilibrium, label, rf, probe, common, physics_model_version=PHYSICS_MODEL_VERSION):
    # Atomic observable and plot choices never enter this cache key.
    return compute_alkali_response(
        equilibrium, label, dict(rf, probe=probe), common,
        contexts_cached(equilibrium, common), optical=True,
    )
