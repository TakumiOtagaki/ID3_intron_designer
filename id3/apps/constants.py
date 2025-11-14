"""
Common constants used by the ID3 demo apps.
"""

MODE_CONFIG = {
    'det.soft': {'alpha': 0, 'beta': 0},   # Deterministic, Soft output
    'det.hard': {'alpha': 0, 'beta': 1},   # Deterministic, Hard output
    'sto.soft': {'alpha': 0.1, 'beta': 0},  # Stochastic, Soft output
    'sto.hard': {'alpha': 0.1, 'beta': 1},  # Stochastic, Hard output
}
