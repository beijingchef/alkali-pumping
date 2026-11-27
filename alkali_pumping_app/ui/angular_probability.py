"""Interactive coherent-state angular probability surfaces and CSV samples."""

from math import pi

import numpy as np
import pandas as pd
import streamlit as st

from ..physics.spin_coherent import angular_probability_surface


def _surface_component():
    # AppTest starts fresh sessions while keeping imported modules cached. Keep
    # one registration per session so timed redraws reuse the same definition.
    component_key = "_acc_angular_probability_component"
    if component_key in st.session_state:
        return st.session_state[component_key]
    component = st.components.v2.component(
    "acc_angular_probability_surface",
    html="""
<div class="surface-view">
  <canvas width="600" height="480" aria-label="Angular momentum probability surface; drag to rotate"></canvas>
  <div class="surface-hint">Drag to rotate · Scroll to zoom</div>
</div>
""",
    css="""
.surface-view { width: 100%; }
.surface-view canvas {
  display: block; width: 100%; height: auto; aspect-ratio: 5 / 4;
  background: #fff; border: 1px solid #d8dce3; border-radius: 5px;
  cursor: grab; touch-action: none;
}
.surface-view canvas:active { cursor: grabbing; }
.surface-hint { color: #555; font: 12px sans-serif; margin-top: 3px; text-align: center; }
""",
    js="""
export default function ({ data, key, parentElement }) {
  const canvas = parentElement.querySelector("canvas");
  if (!canvas || !data) return;
  const ctx = canvas.getContext("2d");
  const width = canvas.width, height = canvas.height;
  const cameraKey = `acc-angular-camera-${key}`;
  let camera = { yaw: 0.72, pitch: 0.40, zoom: 1 };
  try { camera = { ...camera, ...JSON.parse(sessionStorage.getItem(cameraKey) || "{}") }; }
  catch (_) { /* Use the default view. */ }
  const nTheta = data.n_theta, nPhi = data.n_phi;
  const vertices = [];
  for (let i = 0; i < nTheta; i++) {
    const theta = Math.PI * i / (nTheta - 1);
    for (let j = 0; j < nPhi; j++) {
      const phi = 2 * Math.PI * j / (nPhi - 1);
      const r = data.radii[i * nPhi + j];
      vertices.push([r * Math.sin(theta) * Math.cos(phi),
                     r * Math.sin(theta) * Math.sin(phi), r * Math.cos(theta), r]);
    }
  }
  // Light blue/purple to peach/yellow, with soft directional light.
  const palette = [
    [0.00, [111, 130, 188]], [0.32, [158, 145, 197]],
    [0.66, [223, 170, 155]], [1.00, [251, 230, 181]],
  ];
  function color(value, light) {
    const span = data.radius_max - data.radius_min;
    const t = span < 1e-7 ? 0.55 : Math.max(0, Math.min(1, (value - data.radius_min) / span));
    let a = palette[0], b = palette[palette.length - 1];
    for (let i = 1; i < palette.length; i++) {
      if (t <= palette[i][0]) { a = palette[i - 1]; b = palette[i]; break; }
    }
    const fraction = (t - a[0]) / (b[0] - a[0]);
    const rgb = a[1].map((v, i) => Math.round((v + fraction * (b[1][i] - v)) * light));
    return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
  }
  function draw() {
    ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, width, height);
    const cy = Math.cos(camera.yaw), sy = Math.sin(camera.yaw);
    const cp = Math.cos(camera.pitch), sp = Math.sin(camera.pitch);
    const scale = 0.43 * Math.min(width, height) * camera.zoom / data.limit;
    const projected = vertices.map(([x, y, z, r]) => {
      const xr = cy * x - sy * y, yr = sy * x + cy * y;
      const yd = cp * yr - sp * z, depth = sp * yr + cp * z;
      return [width / 2 + scale * xr, height / 2 - scale * yd, depth, r, xr, yd];
    });
    const triangles = [];
    for (let i = 0; i < nTheta - 1; i++) {
      for (let j = 0; j < nPhi - 1; j++) {
        const a = i * nPhi + j, b = a + 1, c = a + nPhi, d = c + 1;
        triangles.push([a, b, c], [b, d, c]);
      }
    }
    triangles.sort((first, second) =>
      (projected[first[0]][2] + projected[first[1]][2] + projected[first[2]][2]) -
      (projected[second[0]][2] + projected[second[1]][2] + projected[second[2]][2]));
    for (const ids of triangles) {
      const p = ids.map(i => projected[i]);
      const radius = (p[0][3] + p[1][3] + p[2][3]) / 3;
      const cx = (p[0][4] + p[1][4] + p[2][4]) / 3;
      const cy = (p[0][5] + p[1][5] + p[2][5]) / 3;
      const cz = (p[0][2] + p[1][2] + p[2][2]) / 3;
      const norm = Math.hypot(cx, cy, cz) || 1;
      const light = 0.88 + 0.12 * Math.max(0, (0.35 * cx - 0.30 * cy + 0.89 * cz) / norm);
      ctx.fillStyle = color(radius, light);
      ctx.beginPath(); ctx.moveTo(p[0][0], p[0][1]);
      ctx.lineTo(p[1][0], p[1][1]); ctx.lineTo(p[2][0], p[2][1]);
      ctx.closePath(); ctx.fill();
    }
    // Front-facing latitude and longitude curves give the surface a quiet
    // Mathematica-like mesh without showing hidden lines through the body.
    ctx.strokeStyle = "rgba(39, 45, 71, 0.22)"; ctx.lineWidth = 0.7;
    for (let i = 5; i < nTheta - 1; i += 5) {
      ctx.beginPath(); let visible = false;
      for (let j = 0; j < nPhi; j++) {
        const p = projected[i * nPhi + j];
        if (p[2] > 0) {
          if (visible) ctx.lineTo(p[0], p[1]); else ctx.moveTo(p[0], p[1]);
          visible = true;
        } else visible = false;
      }
      ctx.stroke();
    }
    for (let j = 0; j < nPhi; j += 10) {
      ctx.beginPath(); let visible = false;
      for (let i = 0; i < nTheta; i++) {
        const p = projected[i * nPhi + j];
        if (p[2] > 0) {
          if (visible) ctx.lineTo(p[0], p[1]); else ctx.moveTo(p[0], p[1]);
          visible = true;
        } else visible = false;
      }
      ctx.stroke();
    }
    ctx.fillStyle = "#1f2937"; ctx.font = "15px sans-serif";
    ctx.fillText("x′", 28, height - 25); ctx.fillText("y′", 61, height - 25);
    ctx.fillText("z′", 94, height - 25);
    const axes = [[1, 0, 0, "#b35e5e"], [0, 1, 0, "#659657"], [0, 0, 1, "#526eae"]];
    for (const [x, y, z, stroke] of axes) {
      const xr = cy * x - sy * y, yr = sy * x + cy * y;
      ctx.strokeStyle = stroke; ctx.lineWidth = 3;
      ctx.beginPath(); ctx.moveTo(38, height - 60);
      ctx.lineTo(38 + 33 * xr, height - 60 - 33 * (cp * yr - sp * z)); ctx.stroke();
    }
  }
  let dragging = false, lastX = 0, lastY = 0;
  canvas.onpointerdown = event => {
    dragging = true; lastX = event.clientX; lastY = event.clientY;
    canvas.setPointerCapture(event.pointerId);
  };
  canvas.onpointermove = event => {
    if (!dragging) return;
    camera.yaw += (event.clientX - lastX) * 0.009;
    camera.pitch = Math.max(-1.45, Math.min(1.45, camera.pitch + (event.clientY - lastY) * 0.009));
    lastX = event.clientX; lastY = event.clientY;
    sessionStorage.setItem(cameraKey, JSON.stringify(camera)); draw();
  };
  canvas.onpointerup = canvas.onpointercancel = () => { dragging = false; };
  canvas.onwheel = event => {
    event.preventDefault();
    camera.zoom = Math.max(0.5, Math.min(2.5, camera.zoom * Math.exp(-event.deltaY * 0.001)));
    sessionStorage.setItem(cameraKey, JSON.stringify(camera)); draw();
  };
  draw();
}
""",
)
    st.session_state[component_key] = component
    return component


def render_angular_probability_surfaces(result, density_matrix, selected_F, *, key_prefix="acc_surface"):
    """Mount the selected manifold and return its sampled CSV values."""
    states = result["ground_states"]
    theta, phi = np.meshgrid(
        np.linspace(0.0, pi, 31), np.linspace(0.0, 2.0 * pi, 61),
        indexing="ij",
    )
    surface_component = _surface_component()
    manifold_population = sum(float(result["population"][i]) for i, state in enumerate(states)
                              if np.isclose(float(state["F"]), selected_F))
    probability = angular_probability_surface(density_matrix, states, selected_F, theta, phi)
    radius = 4.0 * pi * probability / manifold_population
    st.markdown(f"**F={selected_F:g}**  ($P_F$={manifold_population:.3g})")
    surface_component(
        key=f"{key_prefix}_{selected_F:g}",
        data={
            "n_theta": theta.shape[0], "n_phi": theta.shape[1],
            "radii": radius.ravel().tolist(),
            "radius_min": float(radius.min()), "radius_max": float(radius.max()),
            "limit": max(1.1, float(radius.max()) * 1.05),
        },
    )
    return pd.DataFrame({
        "F": selected_F,
        "theta_deg": np.rad2deg(theta).ravel(),
        "phi_deg": np.rad2deg(phi).ravel(),
        "Q_per_sr": probability.ravel(),
        "conditional_radius": radius.ravel(),
        "P_F": manifold_population,
    })
