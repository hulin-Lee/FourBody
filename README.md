# FourBody
<img src="FourBody.png" alt="Schematic configuration of the four-body system" height="300">

This figure illustrates the four-body system, consisting of a stellar binary and a massive black hole binary (MBHB).

In this simulation, we model a scenario where an intermediate-mass black hole (IMBH) migrates toward a supermassive black hole (SMBH). During the migration, the IMBH may capture a stellar binary into a resonant state, forcing the stellar binary to migrate toward the SMBH. As the process unfolds, the stellar binary may collide, undergo tidal disruption, or survive. This program predicts the fate of the stellar binary based on various input parameters. For more details, see our forthcoming paper (in preparation).

### Customizable Parameters:
- **`--seed`**: Initial seed for the simulation.
- **`--q_MBHB`**: Mass ratio between the IMBH and SMBH, ranging from 0 to 1.
- **`--Mp`**: Mass of the primary star (in solar masses), ranging from 1 to 100.
- **`--Ms`**: Mass of the secondary star (in solar masses), ranging from 0.1 to 100.
- **`--a_CM`**: Initial semi-major axis of the stellar binary's center of mass, ranging from 0.1 to 0.7 times the initial semi-major axis of the MBHB ($a_i$).
- **`--I_CM`**: Initial inclination of the stellar binary’s center of mass (in degrees), ranging from 0 to 90.
- **`--a_b`**: Initial semi-major axis of the relative motion between the two stars, ranging from 100 to 2000 solar radii.
- **`--tau`**: Migration timescale of the IMBH, where the evolution of the MBHB’s semi-major axis follows the equation: $a_{\text{MBHB}} = a_i + (a_i - a_f)e^{-t/\tau}$, with $\tau$ in units of the MBHB’s initial period.

### Output Files:
- **`data.npy`**: Contains the trajectories of the stellar binary’s center of mass and the IMBH.
- **`ending.txt`**: Provides the parameter settings and a brief summary of the stellar binary’s final state.



