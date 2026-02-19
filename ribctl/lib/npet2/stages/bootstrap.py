# ribctl/lib/npet2/stages/bootstrap.py (updated)
from __future__ import annotations

from typing import Any, Dict

import numpy as np

from npet2.core.pipeline import Stage
from npet2.core.ribosome_types import RibosomeProfile
from npet2.core.types import StageContext


class Stage00Inputs(Stage):
    key = "00_inputs"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        return {}

    def run(self, ctx: StageContext) -> None:
        structure_provider = ctx.require("structure_provider")
        data = structure_provider.load_atoms(ctx.rcsb_id)

        atom_xyz = np.asarray(data["atom_xyz"], dtype=np.float32)
        ctx.inputs["atom_xyz"] = atom_xyz
        ctx.inputs["atom_element"] = data.get("atom_element", None)
        ctx.inputs["mmcif_path"] = data["mmcif_path"]

        # Profile: validate it's the right type
        profile = data["profile"]
        if not isinstance(profile, RibosomeProfile):
            profile = RibosomeProfile.model_validate(profile)
        ctx.inputs["profile"] = profile

        # If the provider gave us a biopython structure or RibosomeOps, stash it.
        # Standalone providers won't -- Stage30 will parse mmcif on demand.
        if "ro" in data:
            ro = data["ro"]
            ctx.inputs["ro"] = ro
            ctx.inputs["biopython_structure"] = ro.assets.biopython_structure()
        elif "biopython_structure" in data:
            ctx.inputs["biopython_structure"] = data["biopython_structure"]

        ctx.artifacts["atom_xyz"] = ctx.store.put_numpy(
            name="atom_xyz",
            stage=self.key,
            arr=atom_xyz,
            meta={"shape": list(atom_xyz.shape), "dtype": str(atom_xyz.dtype)},
        )

        mins = atom_xyz.min(axis=0)
        maxs = atom_xyz.max(axis=0)
        ctx.stats["atom_bounds"] = {"min": mins.tolist(), "max": maxs.tolist()}
        ctx.stats["n_atoms"] = int(atom_xyz.shape[0])


class Stage10Landmarks(Stage):
    key = "10_landmarks"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        return {}

    def run(self, ctx: StageContext) -> None:
        landmark_provider = ctx.require("landmark_provider")
        lm = landmark_provider.get_landmarks(ctx.rcsb_id)

        ptc = np.asarray(lm["ptc_xyz"], dtype=np.float32)
        constr = np.asarray(lm["constriction_xyz"], dtype=np.float32)

        if ptc.shape != (3,):
            raise ValueError(f"PTC must be shape (3,), got {ptc.shape}")
        if constr.shape != (3,):
            raise ValueError(f"Constriction must be shape (3,), got {constr.shape}")

        ctx.inputs["ptc_xyz"] = ptc
        ctx.inputs["constriction_xyz"] = constr

        D = float(np.linalg.norm(constr - ptc))
        if D < 1.0:
            raise ValueError(
                f"PTC and constriction are too close ({D:.1f}A) -- check coordinates"
            )

        z_min = -float(ctx.config.cylinder_ptc_extension_A)
        z_max = float(ctx.config.cylinder_height_A)

        ctx.inputs["cylinder_z_min"] = z_min
        ctx.inputs["cylinder_z_max"] = z_max
        ctx.inputs["landmark_distance"] = D

        print(
            f"  [10_landmarks] PTC-Constriction distance={D:.1f}A, "
            f"cylinder z=[{z_min:.1f}, {z_max:.1f}]A"
        )

        ctx.artifacts["ptc"] = ctx.store.put_json(
            name="ptc",
            stage=self.key,
            obj={"location": ptc.tolist()},
            meta={"units": "A"},
        )
        ctx.artifacts["constriction_site"] = ctx.store.put_json(
            name="constriction_site",
            stage=self.key,
            obj={"location": constr.tolist()},
            meta={
                "units": "A",
                "landmark_distance_A": D,
                "cylinder_z_min": z_min,
                "cylinder_z_max": z_max,
            },
        )
