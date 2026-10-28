"""Shared spherical-coordinate field controls and condition validation."""
from math import isfinite
import streamlit as st


def validate_field_conditions(values):
    for field, upper in (('static_field_theta_deg', 180), ('static_field_phi_deg', 360)):
        value = values[field]
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not isfinite(value) or not 0 <= value <= upper or int(value) != value):
            raise ValueError(f'{field} must be a whole number of degrees between 0 and {upper}.')
    strength = values['static_field_nT']
    if isinstance(strength, bool) or not isinstance(strength, (int, float)) or not isfinite(strength):
        raise ValueError('static_field_nT must be finite.')


def render_field_controls(prefix=''):
    columns = st.columns(3, gap='xsmall')
    with columns[0]:
        theta = st.number_input(r'$\theta_{\mathbf{B}}$ (°)', min_value=0, max_value=180, step=1, format='%d',
                                key=prefix+'static_field_theta_deg', width='stretch', help='Polar angle from laboratory +z.')
    with columns[1]:
        phi = st.number_input(r'$\phi_{\mathbf{B}}$ (°)', min_value=0, max_value=360, step=1, format='%d',
                              key=prefix+'static_field_phi_deg', width='stretch', help='Azimuth from +x toward +y.')
    with columns[2]:
        strength = st.number_input(r'$|\mathbf{B}|$ (nT)', step=1., format='%g', key=prefix+'static_field_nT', width='stretch',
                                   help='A negative value reverses the direction specified by the angles.')
    return theta, phi, strength
